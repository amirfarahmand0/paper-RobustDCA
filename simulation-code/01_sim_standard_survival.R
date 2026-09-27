# Simulation study: standard survival setting.
# Evaluates the KM-based and CIPCW-based net benefit estimators under non-informative,
# linear informative, and nonlinear informative censoring.

library(data.table)
library(survival)
library(survivalROC)
library(mvtnorm)
library(pbapply)
library(future.apply)
library(parallel)
library(dplyr)
library(splines)
library(mgcv)
library(progressr)

# Initialize parallel processing
n_workers <- detectCores() - 1
plan(multisession, workers = n_workers)
cat("Workers used:", n_workers, "\n")
handlers("txtprogressbar")

# 1. Data Generating Mechanisms (DGM) ----

generate_survival_tuned <- function(n = 800, target_event_pct = 0.35, target_t = 1,
                                    delta_beta = 0.45, HRc_Z1, HRc_Z2, target_cens_pct = 0.3) {

  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Z <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]

  W1 <- rexp(n, rate = 0.1)
  W2 <- rexp(n, rate = 0.2)

  beta_Z1 <- 1 + delta_beta
  beta_Z2 <- 1 - delta_beta
  lp_event <- beta_Z1 * Z1 + beta_Z2 * Z2

  gamma1 <- log(HRc_Z1)
  gamma2 <- log(HRc_Z2)
  exp_gamma <- exp(gamma1 * Z1 + gamma2 * Z2)

  optim_event <- function(log_lambda) {
    lambda <- exp(log_lambda)
    prob <- 1 - exp(-lambda * target_t * exp(lp_event))
    abs(mean(prob) - target_event_pct)
  }
  lambda_event <- exp(optim(par = log(0.4), fn = optim_event, method = "Brent", lower = -10, upper = 10)$par)

  optim_cens <- function(log_lambda_c) {
    lambda_c <- exp(log_lambda_c)
    rate_c <- lambda_c * exp_gamma
    rate_e <- lambda_event * exp(lp_event)
    prob_obs_cens <- (rate_c / (rate_c + rate_e)) * (1 - exp(-target_t * (rate_c + rate_e)))
    abs(mean(prob_obs_cens) - target_cens_pct)
  }
  lambda_cens <- exp(optim(par = log(0.1), fn = optim_cens, method = "Brent", lower = -10, upper = 10)$par)

  true_T <- rexp(n, rate = lambda_event * exp(lp_event))
  C <- rexp(n, rate = lambda_cens * exp_gamma)

  obs_time <- pmin(true_T, C)
  delta <- as.integer(true_T <= C)

  data.frame(obs_time = obs_time, delta = delta, true_T = true_T, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

generate_survival_tuned_cens_nonlinear <- function(n = 800, target_event_pct = 0.35, target_t = 1,
                                                   delta_beta = 0.45, target_cens_pct = 0.30) {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Z <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]

  W1 <- rexp(n, rate = 0.1)
  W2 <- rexp(n, rate = 0.2)

  beta_Z1 <- 1 + delta_beta
  beta_Z2 <- 1 - delta_beta
  lp_event <- beta_Z1 * Z1 + beta_Z2 * Z2

  optim_event <- function(log_lambda) {
    lambda <- exp(log_lambda)
    prob <- 1 - exp(-lambda * target_t * exp(lp_event))
    abs(mean(prob) - target_event_pct)
  }
  lambda_event <- exp(optim(par = log(0.4), fn = optim_event, method = "Brent", lower = -10, upper = 10)$par)
  rate_e <- lambda_event * exp(lp_event)

  lp_cens <- 1.5 * tanh(Z1) - 2.0 * exp(-(Z2^2))
  exp_gamma <- exp(lp_cens)

  optim_cens <- function(log_lambda_c) {
    lambda_c <- exp(log_lambda_c)
    rate_c <- lambda_c * exp_gamma
    term1 <- rate_c / (rate_c + rate_e)
    term2 <- 1 - exp(-target_t * (rate_c + rate_e))
    prob_obs_cens <- term1 * term2
    abs(mean(prob_obs_cens) - target_cens_pct)
  }
  lambda_cens <- exp(optim(par = log(0.1), fn = optim_cens, method = "Brent", lower = -10, upper = 10)$par)

  true_T <- rexp(n, rate = rate_e)
  C <- rexp(n, rate = lambda_cens * exp_gamma)

  obs_time <- pmin(true_T, C)
  delta <- as.integer(true_T <= C)

  data.frame(obs_time = obs_time, delta = delta, true_T = true_T, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

generate_data_by_scenario <- function(scenario, n, target_event_pct, target_t, delta_beta, target_cens_pct) {
  if (scenario == "informative_simple_cox") {
    return(generate_survival_tuned(n = n, target_event_pct = target_event_pct, target_t = target_t,
                                   delta_beta = delta_beta, HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = target_cens_pct))
  }
  if (scenario == "non_informative_cox") {
    return(generate_survival_tuned(n = n, target_event_pct = target_event_pct, target_t = target_t,
                                   delta_beta = delta_beta, HRc_Z1 = 1, HRc_Z2 = 1, target_cens_pct = target_cens_pct))
  }
  if (scenario == "informative_nonlinear") {
    return(generate_survival_tuned_cens_nonlinear(n = n, target_event_pct = target_event_pct, target_t = target_t,
                                                  delta_beta = delta_beta, target_cens_pct = target_cens_pct))
  }
  stop("Unknown scenario")
}

# 2. Model Training & Risk Prediction ----

fit_black_box_model <- function(train_data) {
  coxph(Surv(obs_time, delta) ~ Z1 + Z2, data = train_data)
}

predict_risk <- function(model, data, target_t = 1) {
  lp <- predict(model, newdata = data, type = "lp")
  bh <- basehaz(model, centered = TRUE)

  H0_t <- if (target_t <= max(bh$time)) {
    approx(bh$time, bh$hazard, xout = target_t, method = "linear", rule = 2)$y
  } else {
    tail(bh$hazard, 1)
  }

  1 - exp(-H0_t * exp(lp))
}

# 3. Net Benefit Components ----

true_components <- function(data, t, z, risk_score) {
  dt <- as.data.table(data)
  D <- dt$true_T <= t
  pos <- risk_score > z
  c(prevalence = mean(D), sensitivity = mean(pos[D]), specificity = 1 - mean(pos[!D]))
}

km_components <- function(data, t, z, risk_score) {
  dt <- as.data.table(data)
  fit <- survfit(Surv(obs_time, delta) ~ 1, data = dt)

  S_t <- if (t %in% fit$time) {
    fit$surv[fit$time == t]
  } else {
    approx(fit$time, fit$surv, xout = t, method = "linear", rule = 2)$y
  }
  prev <- 1 - S_t

  roc <- survivalROC(Stime = dt$obs_time, status = dt$delta, marker = risk_score,
                     predict.time = t, cut.values = z, method = "KM")
  c(prevalence = prev, sensitivity = roc$TP[2], specificity = 1 - roc$FP[2])
}

fit_all_censoring_weights <- function(data, t) {
  dt <- as.data.table(data)
  out <- list()

  # IPCW True Cox
  fit_cox <- coxph(Surv(obs_time, 1 - delta) ~ Z1 + Z2, data = dt)
  bh_cox <- basehaz(fit_cox, centered = TRUE)
  setorder(bh_cox, time)
  get_L0_cox <- function(u) approx(bh_cox$time, bh_cox$hazard, xout = u, method = "linear", rule = 2)$y
  lp_cox <- predict(fit_cox, newdata = dt, type = "lp")
  out$cox <- list(obs = pmax(exp(-get_L0_cox(dt$obs_time) * exp(lp_cox)), 0.02),
                  t = pmax(exp(-get_L0_cox(t) * exp(lp_cox)), 0.02))

  # IPCW Misspecified Cox
  fit_miscox <- coxph(Surv(obs_time, 1 - delta) ~ 1 + Z1, data = dt)
  bh_miscox <- basehaz(fit_miscox, centered = TRUE)
  setorder(bh_miscox, time)
  get_L0_miscox <- function(u) approx(bh_miscox$time, bh_miscox$hazard, xout = u, method = "linear", rule = 2)$y
  lp_miscox <- predict(fit_miscox, newdata = dt, type = "lp")
  out$miscox <- list(obs = pmax(exp(-get_L0_miscox(dt$obs_time) * exp(lp_miscox)), 0.02),
                     t = pmax(exp(-get_L0_miscox(t) * exp(lp_miscox)), 0.02))

  # IPCW GAM
  fit_gam <- gam(obs_time ~ s(Z1, k = 4, bs = "cr") + s(Z2, k = 4, bs = "cr"), family = cox.ph(), data = dt, weights = (1 - delta))
  out$gam <- list(obs = pmax(predict(fit_gam, newdata = dt, type = "response"), 0.02))
  dt_t <- copy(dt)
  dt_t$obs_time <- t
  out$gam$t <- pmax(predict(fit_gam, newdata = dt_t, type = "response"), 0.02)

  return(out)
}

fast_ipcw_components <- function(data, t, z, risk_score, G_hat_obs, G_hat_t) {
  dt <- as.data.table(data)
  weight_event <- dt$delta / G_hat_obs
  indicator_event <- (dt$obs_time <= t)
  prev <- mean(weight_event * indicator_event)

  indicator_pos <- (risk_score > z)
  num_se <- sum(weight_event * indicator_event * indicator_pos)
  den_se <- sum(weight_event * indicator_event)
  Se <- if (den_se == 0) 0 else num_se / den_se

  indicator_surv <- (dt$obs_time > t)
  indicator_neg <- (risk_score <= z)
  num_sp <- sum(indicator_surv * indicator_neg / G_hat_t)
  den_sp <- sum(indicator_surv / G_hat_t)
  Sp <- if (den_sp == 0) 1 else num_sp / den_sp

  c(prevalence = prev, sensitivity = Se, specificity = Sp)
}

# 4. Simulation Execution Logic ----

compute_nb <- function(components, z) {
  prev <- components["prevalence"]
  sens <- components["sensitivity"]
  spec <- components["specificity"]
  sens * prev - (1 - spec) * (1 - prev) * (z / (1 - z))
}

compute_true_nb <- function(model, large_data, target_t = 1, pt = 0.2) {
  risk_score <- predict_risk(model, large_data, target_t)
  comps <- true_components(large_data, t = target_t, z = pt, risk_score = risk_score)
  compute_nb(comps, z = pt)
}

run_single_replication <- function(model, valid_data, target_t, thresholds, B, NB_true_list) {
  n <- nrow(valid_data)
  risk_score <- predict_risk(model, valid_data, target_t)
  weights_main <- fit_all_censoring_weights(valid_data, target_t)
  res_main_list <- list()

  for (pt in thresholds) {
    nb_t <- NB_true_list[[as.character(pt)]]
    res_main_list[[as.character(pt)]] <- data.frame(
      threshold = pt,
      NB_true = nb_t,
      NB_km = compute_nb(km_components(valid_data, target_t, pt, risk_score), pt),
      NB_ipcw_cox = compute_nb(fast_ipcw_components(valid_data, target_t, pt, risk_score, weights_main$cox$obs, weights_main$cox$t), pt),
      NB_ipcw_miscox = compute_nb(fast_ipcw_components(valid_data, target_t, pt, risk_score, weights_main$miscox$obs, weights_main$miscox$t), pt),
      NB_ipcw_gam = compute_nb(fast_ipcw_components(valid_data, target_t, pt, risk_score, weights_main$gam$obs, weights_main$gam$t), pt)
    )
  }
  res_main_df <- do.call(rbind, res_main_list)

  boot_res <- vector("list", B)
  for (b in 1:B) {
    success <- FALSE
    attempts <- 0

    while (!success && attempts < 10) {
      boot_idx <- sample(1:n, n, replace = TRUE)
      boot_data <- valid_data[boot_idx, ]
      boot_risk <- risk_score[boot_idx]

      fail_reason <- NULL
      weights_boot <- tryCatch({
        fit_all_censoring_weights(boot_data, target_t)
      }, error = function(e) {
        fail_reason <<- paste("error:", conditionMessage(e))
        NULL
      }, warning = function(w) {
        fail_reason <<- paste("warning:", conditionMessage(w))
        NULL
      })

      if (!is.null(weights_boot)) {
        boot_pts <- list()
        for (pt in thresholds) {
          boot_pts[[as.character(pt)]] <- c(
            compute_nb(km_components(boot_data, target_t, pt, boot_risk), pt),
            compute_nb(fast_ipcw_components(boot_data, target_t, pt, boot_risk, weights_boot$cox$obs, weights_boot$cox$t), pt),
            compute_nb(fast_ipcw_components(boot_data, target_t, pt, boot_risk, weights_boot$miscox$obs, weights_boot$miscox$t), pt),
            compute_nb(fast_ipcw_components(boot_data, target_t, pt, boot_risk, weights_boot$gam$obs, weights_boot$gam$t), pt)
          )
        }
        boot_res[[b]] <- boot_pts
        success <- TRUE
      } else {
        message(sprintf("[boot retry %d/10] %s", attempts + 1, fail_reason))
      }
      attempts <- attempts + 1
    }

    if (!success) {
      boot_pts <- list()
      for (pt in thresholds) {
        boot_pts[[as.character(pt)]] <- rep(NA, 4)
      }
      boot_res[[b]] <- boot_pts
    }
  }

  ci_list <- list()
  for (pt in thresholds) {
    pt_char <- as.character(pt)
    mat <- do.call(rbind, lapply(boot_res, function(x) x[[pt_char]]))
    cis <- apply(mat, 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE)

    ci_list[[pt_char]] <- data.frame(
      threshold = pt,
      ci_km_lower = cis[1, 1], ci_km_upper = cis[2, 1],
      ci_ipcw_cox_lower = cis[1, 2], ci_ipcw_cox_upper = cis[2, 2],
      ci_ipcw_miscox_lower = cis[1, 3], ci_ipcw_miscox_upper = cis[2, 3],
      ci_ipcw_gam_lower = cis[1, 4], ci_ipcw_gam_upper = cis[2, 4]
    )
  }

  final_df <- merge(res_main_df, do.call(rbind, ci_list), by = "threshold")

  final_df$bias_km <- final_df$NB_km - final_df$NB_true
  final_df$bias_ipcw_cox <- final_df$NB_ipcw_cox - final_df$NB_true
  final_df$bias_ipcw_miscox <- final_df$NB_ipcw_miscox - final_df$NB_true
  final_df$bias_ipcw_gam <- final_df$NB_ipcw_gam - final_df$NB_true

  final_df$sqerr_km <- (final_df$NB_km - final_df$NB_true)^2
  final_df$sqerr_ipcw_cox <- (final_df$NB_ipcw_cox - final_df$NB_true)^2
  final_df$sqerr_ipcw_miscox <- (final_df$NB_ipcw_miscox - final_df$NB_true)^2
  final_df$sqerr_ipcw_gam <- (final_df$NB_ipcw_gam - final_df$NB_true)^2

  return(final_df)
}

# 5. Main Execution Loop ----

common_params <- list(target_event_pct = 0.3, target_t = 5, delta_beta = 0)
scenario_names <- c("non_informative_cox", "informative_simple_cox", "informative_nonlinear")
thresholds <- c(0.01, 0.1, 0.25, 0.5, 0.75)
censoring_rates <- c(0.10, 0.20, 0.30)

n_reps <- 1000
B_boot <- 500
all_results <- list()

for (scen_name in scenario_names) {
  for (cens_rate in censoring_rates) {
    cat(sprintf("\n[Running scenario: %s | Censoring: %.1f%%]\n", scen_name, cens_rate * 100))

    # 1. Train Model
    set.seed(123)
    train_data <- generate_survival_tuned(n = 1000, target_event_pct = common_params$target_event_pct,
                                          target_t = common_params$target_t, delta_beta = common_params$delta_beta,
                                          HRc_Z1 = 1, HRc_Z2 = 1, target_cens_pct = cens_rate)
    model <- fit_black_box_model(train_data)

    # 2. Generate Full Population
    set.seed(456)
    large_data <- generate_data_by_scenario(scenario = scen_name, n = 1000000,
                                            target_event_pct = common_params$target_event_pct,
                                            target_t = common_params$target_t, delta_beta = common_params$delta_beta,
                                            target_cens_pct = cens_rate)

    pop_event_rate <- mean(large_data$obs_time <= common_params$target_t & large_data$delta == 1)
    n_test <- ceiling(100 / pop_event_rate)

    # 3. Calculate True NB on Full 1,000,000 Cohort
    NB_true_list <- list()
    for (pt in thresholds) {
      NB_true_list[[as.character(pt)]] <- compute_true_nb(model, large_data, target_t = common_params$target_t, pt = pt)
    }

    # 4. Parallel Replications
    with_progress({
      p <- progressor(steps = n_reps)
      results <- future_lapply(1:n_reps, function(r) {

        success_main <- FALSE
        res <- NULL
        attempts_main <- 0

        # The while loop protects the main estimate (up to 10 tries)
        while (!success_main && attempts_main < 10) {
          valid_idx <- sample(seq_len(nrow(large_data)), n_test, replace = FALSE)

          fail_reason_main <- NULL
          # Try to run the main estimate and bootstraps
          res <- tryCatch({
            run_single_replication(model = model, valid_data = large_data[valid_idx, ],
                                   target_t = common_params$target_t, thresholds = thresholds,
                                   B = B_boot, NB_true_list = NB_true_list)
          }, error = function(e) {
            fail_reason_main <<- paste("error:", conditionMessage(e))
            NULL
          }, warning = function(w) {
            fail_reason_main <<- paste("warning:", conditionMessage(w))
            NULL
          })

          # Stop retrying once the replication succeeds
          if (!is.null(res)) {
            success_main <- TRUE
          } else {
            message(sprintf("[main retry %d/10] %s", attempts_main + 1, fail_reason_main))
          }
          attempts_main <- attempts_main + 1
        }

        # If it failed all 10 times, generate an NA dataframe for this replication
        if (!success_main) {
          res <- data.frame(
            threshold = thresholds,
            NB_true = unlist(NB_true_list),
            NB_km = NA, NB_ipcw_cox = NA, NB_ipcw_miscox = NA, NB_ipcw_gam = NA,
            ci_km_lower = NA, ci_km_upper = NA, ci_ipcw_cox_lower = NA, ci_ipcw_cox_upper = NA,
            ci_ipcw_miscox_lower = NA, ci_ipcw_miscox_upper = NA, ci_ipcw_gam_lower = NA, ci_ipcw_gam_upper = NA,
            bias_km = NA, bias_ipcw_cox = NA, bias_ipcw_miscox = NA, bias_ipcw_gam = NA,
            sqerr_km = NA, sqerr_ipcw_cox = NA, sqerr_ipcw_miscox = NA, sqerr_ipcw_gam = NA
          )
        }

        p()
        return(res)
      }, future.seed = 123L, future.scheduling = 1)
    })

    # 5. Format output
    results_df <- do.call(rbind, results)
    results_df$scenario <- scen_name
    results_df$target_cens_pct <- cens_rate
    results_df$target_t <- common_params$target_t

    results_df <- results_df %>% mutate(
      cover_km = (NB_true >= ci_km_lower & NB_true <= ci_km_upper),
      cover_ipcw_cox = (NB_true >= ci_ipcw_cox_lower & NB_true <= ci_ipcw_cox_upper),
      cover_ipcw_miscox = (NB_true >= ci_ipcw_miscox_lower & NB_true <= ci_ipcw_miscox_upper),
      cover_ipcw_gam = (NB_true >= ci_ipcw_gam_lower & NB_true <= ci_ipcw_gam_upper)
    )

    all_results[[paste0(scen_name, "_cens_", gsub("\\.", "", as.character(cens_rate)))]] <- results_df
  }
}

combined_results <- bind_rows(all_results)
cat("\nStandard survival simulation script complete.\n")
