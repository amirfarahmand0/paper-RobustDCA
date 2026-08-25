# ==============================================================================
# Script Name: 02_sim_competing_risks.R
# Purpose:     Monte Carlo simulation evaluating robust Net Benefit (NB)
#              estimators under informative censoring for competing risk outcomes.
# ==============================================================================

library(survival)
library(data.table)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(mgcv)
library(parallel)
library(future.apply)
library(progressr)

# Initialize parallel processing
n_workers <- detectCores() - 1
plan(multisession, workers = n_workers)
cat("Workers used:", n_workers, "\n")
handlers("txtprogressbar")

# ==============================================================================
# 1. Helper & Calibration Functions
# ==============================================================================

validate_cr_inputs <- function(target_any_event_pct, competing_fraction, target_cif1, target_cif2, target_cens_pct) { invisible(TRUE) }

tune_event_lambdas <- function(lp1, lp2, target_cif1, target_cif2, target_t) {
  objective <- function(log_lambdas) {
    lambda1 <- exp(log_lambdas[1])
    lambda2 <- exp(log_lambdas[2])
    r1 <- lambda1 * exp(lp1)
    r2 <- lambda2 * exp(lp2)
    rsum <- r1 + r2

    one_minus_surv <- 1 - exp(-rsum * target_t)
    F1 <- (r1 / rsum) * one_minus_surv
    F2 <- (r2 / rsum) * one_minus_surv

    rel_err1 <- (mean(F1) - target_cif1) / target_cif1
    rel_err2 <- (mean(F2) - target_cif2) / target_cif2
    rel_err1^2 + rel_err2^2
  }

  opt <- optim(par = c(log(0.10), log(0.10)), fn = objective, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-14))
  lambda1 <- exp(opt$par[1])
  lambda2 <- exp(opt$par[2])

  r1 <- lambda1 * exp(lp1)
  r2 <- lambda2 * exp(lp2)
  rsum <- r1 + r2
  one_minus_surv <- 1 - exp(-rsum * target_t)

  list(lambda1 = lambda1, lambda2 = lambda2, achieved_cif1 = mean((r1 / rsum) * one_minus_surv), achieved_cif2 = mean((r2 / rsum) * one_minus_surv), opt = opt)
}

tune_censoring_lambda <- function(rate1, rate2, lp_cens, target_cens_pct, target_t) {
  objective <- function(log_lambda_c) {
    lambda_c <- exp(log_lambda_c)
    rc <- lambda_c * exp(lp_cens)
    rsum <- rate1 + rate2 + rc
    prob_cens <- (rc / rsum) * (1 - exp(-rsum * target_t))
    (mean(prob_cens) - target_cens_pct)^2
  }

  opt <- optim(par = log(0.10), fn = objective, method = "Brent", lower = -10, upper = 10)
  lambda_cens <- exp(opt$par)
  rc <- lambda_cens * exp(lp_cens)
  rsum <- rate1 + rate2 + rc

  list(lambda_cens = lambda_cens, achieved_cens_pct = mean((rc / rsum) * (1 - exp(-rsum * target_t))), opt = opt)
}

# ==============================================================================
# 2. Data Generating Mechanisms (DGM)
# ==============================================================================

generate_competing_risk_tuned <- function(n, target_any_event_pct = 0.30, competing_fraction = 1/3, target_t = 4,
                                          delta_beta = 0, beta2_Z1 = -0.50, beta2_Z2 = 0.75,
                                          censoring_type = c("noninformative", "linear", "nonlinear"),
                                          HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = 0.30) {

  censoring_type <- match.arg(censoring_type)
  target_cif2 <- target_any_event_pct * competing_fraction
  target_cif1 <- target_any_event_pct * (1 - competing_fraction)

  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Z <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]; Z2 <- Z[, 2]

  W1 <- rexp(n, rate = 0.1)
  W2 <- rexp(n, rate = 0.2)

  beta1_Z1 <- 1 + delta_beta
  beta1_Z2 <- 1 - delta_beta
  lp1 <- beta1_Z1 * Z1 + beta1_Z2 * Z2
  lp2 <- beta2_Z1 * Z1 + beta2_Z2 * Z2

  event_tuning <- tune_event_lambdas(lp1, lp2, target_cif1, target_cif2, target_t)
  rate1 <- event_tuning$lambda1 * exp(lp1)
  rate2 <- event_tuning$lambda2 * exp(lp2)

  lp_cens <- switch(censoring_type,
                    noninformative = rep(0, n),
                    linear = log(HRc_Z1) * Z1 + log(HRc_Z2) * Z2,
                    nonlinear = 1.5 * tanh(Z1) - 2.0 * exp(-(Z2^2)))

  cens_tuning <- tune_censoring_lambda(rate1, rate2, lp_cens, target_cens_pct, target_t)
  rate_cens <- cens_tuning$lambda_cens * exp(lp_cens)

  T1 <- rexp(n, rate = rate1)
  T2 <- rexp(n, rate = rate2)
  C <- rexp(n, rate = rate_cens)

  true_time <- pmin(T1, T2)
  true_status <- ifelse(T1 <= T2, 1L, 2L)
  obs_time <- pmin(T1, T2, C)
  status <- ifelse(C < pmin(T1, T2), 0L, ifelse(T1 <= T2, 1L, 2L))

  data.frame(obs_time = obs_time, status = status, delta1 = as.integer(status == 1L),
             delta2 = as.integer(status == 2L), delta_any = as.integer(status != 0L),
             T1 = T1, T2 = T2, C = C, true_time = true_time, true_status = true_status,
             Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

generate_data_by_scenario_cr <- function(scenario, n, target_any_event_pct = 0.30, competing_fraction = 1/3, target_t = 4,
                                         delta_beta = 0, beta2_Z1 = -0.50, beta2_Z2 = 0.75, target_cens_pct = 0.30) {
  scenario <- match.arg(scenario, choices = c("non_informative_cox", "informative_simple_cox", "informative_nonlinear"))

  if (scenario == "non_informative_cox") {
    return(generate_competing_risk_tuned(n = n, target_any_event_pct = target_any_event_pct, competing_fraction = competing_fraction, target_t = target_t, delta_beta = delta_beta, beta2_Z1 = beta2_Z1, beta2_Z2 = beta2_Z2, censoring_type = "noninformative", target_cens_pct = target_cens_pct))
  }
  if (scenario == "informative_simple_cox") {
    return(generate_competing_risk_tuned(n = n, target_any_event_pct = target_any_event_pct, competing_fraction = competing_fraction, target_t = target_t, delta_beta = delta_beta, beta2_Z1 = beta2_Z1, beta2_Z2 = beta2_Z2, censoring_type = "linear", HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = target_cens_pct))
  }
  if (scenario == "informative_nonlinear") {
    return(generate_competing_risk_tuned(n = n, target_any_event_pct = target_any_event_pct, competing_fraction = competing_fraction, target_t = target_t, delta_beta = delta_beta, beta2_Z1 = beta2_Z1, beta2_Z2 = beta2_Z2, censoring_type = "nonlinear", target_cens_pct = target_cens_pct))
  }
}

# ==============================================================================
# 3. Model Training & Risk Prediction
# ==============================================================================

fit_black_box_model_cr <- function(train_data) {
  list(cox1 = coxph(Surv(obs_time, status == 1) ~ Z1 + Z2, data = train_data),
       cox2 = coxph(Surv(obs_time, status == 2) ~ Z1 + Z2, data = train_data))
}

step_hazard_eval <- function(bh, query_times, left_continuous = TRUE) {
  if (nrow(bh) == 0) return(rep(0, length(query_times)))
  sapply(query_times, function(u) {
    idx <- if (left_continuous) which(bh$time < u) else which(bh$time <= u)
    if (length(idx) == 0) 0 else bh$hazard[max(idx)]
  })
}

predict_risk_cr <- function(model, data, target_t = 1) {
  bh1 <- basehaz(model$cox1, centered = TRUE); bh1 <- bh1[order(bh1$time), ]
  bh2 <- basehaz(model$cox2, centered = TRUE); bh2 <- bh2[order(bh2$time), ]

  keep <- bh1$time <= target_t
  event_times1 <- bh1$time[keep]
  if (length(event_times1) == 0) return(rep(0, nrow(data)))

  Lambda01_at_events <- bh1$hazard[keep]
  Lambda01_left <- c(0, head(Lambda01_at_events, -1))
  dLambda01 <- Lambda01_at_events - Lambda01_left
  Lambda02_left <- step_hazard_eval(bh2, event_times1, left_continuous = TRUE)

  lp1 <- predict(model$cox1, newdata = data, type = "lp")
  lp2 <- predict(model$cox2, newdata = data, type = "lp")

  S_left_mat <- exp(-(outer(exp(lp1), Lambda01_left, "*") + outer(exp(lp2), Lambda02_left, "*")))
  CIF1 <- rowSums(S_left_mat * outer(exp(lp1), dLambda01, "*"))

  pmin(pmax(CIF1, 0), 1)
}

# ==============================================================================
# 4. Net Benefit Components
# ==============================================================================

true_components_cr <- function(data, t, z, risk_score) {
  D1 <- (data$true_time <= t) & (data$true_status == 1L)
  pos <- risk_score > z
  c(prevalence = mean(D1), sensitivity = mean(pos[D1]), specificity = 1 - mean(pos[!D1]))
}

fit_all_censoring_weights_cr <- function(data, t) {
  dt <- as.data.table(data)
  out <- list()

  # IPCW True Cox
  fit_cox <- coxph(Surv(obs_time, 1 - delta_any) ~ Z1 + Z2, data = dt)
  bh_cox <- basehaz(fit_cox, centered = TRUE)
  setorder(bh_cox, time)
  get_L0_cox <- function(u) approx(bh_cox$time, bh_cox$hazard, xout = u, method = "linear", rule = 2)$y
  lp_cox <- predict(fit_cox, newdata = dt, type = "lp")
  out$cox <- list(obs = pmax(exp(-get_L0_cox(dt$obs_time) * exp(lp_cox)), 0.02), t = pmax(exp(-get_L0_cox(t) * exp(lp_cox)), 0.02))

  # IPCW Misspecified Cox
  fit_miscox <- coxph(Surv(obs_time, 1 - delta_any) ~ 1 + Z1, data = dt)
  bh_miscox <- basehaz(fit_miscox, centered = TRUE)
  setorder(bh_miscox, time)
  get_L0_miscox <- function(u) approx(bh_miscox$time, bh_miscox$hazard, xout = u, method = "linear", rule = 2)$y
  lp_miscox <- predict(fit_miscox, newdata = dt, type = "lp")
  out$miscox <- list(obs = pmax(exp(-get_L0_miscox(dt$obs_time) * exp(lp_miscox)), 0.02), t = pmax(exp(-get_L0_miscox(t) * exp(lp_miscox)), 0.02))

  # IPCW GAM
  fit_gam <- gam(obs_time ~ s(Z1, k = 4, bs = "cr") + s(Z2, k = 4, bs = "cr"), family = cox.ph(), data = dt, weights = (1 - delta_any))
  out$gam <- list(obs = pmax(predict(fit_gam, newdata = dt, type = "response"), 0.02))
  dt_t <- copy(dt); dt_t$obs_time <- t
  out$gam$t <- pmax(predict(fit_gam, newdata = dt_t, type = "response"), 0.02)

  return(out)
}

fast_ipcw_components_cr <- function(data, t, z, risk_score, G_hat_obs, G_hat_t) {
  dt <- as.data.table(data)

  weight_event <- dt$delta1 / G_hat_obs
  indicator_event <- (dt$obs_time <= t)
  prev <- mean(weight_event * indicator_event)

  num_se <- sum(weight_event * indicator_event * (risk_score > z))
  Se <- if (sum(weight_event * indicator_event) == 0) 0 else num_se / sum(weight_event * indicator_event)

  weight_neg <- ifelse(dt$obs_time > t, 1 / G_hat_t, ifelse((dt$obs_time <= t) & (dt$status == 2L), 1 / G_hat_obs, 0))
  num_sp <- sum(weight_neg * (risk_score <= z))
  Sp <- if (sum(weight_neg) == 0) 1 else num_sp / sum(weight_neg)

  c(prevalence = prev, sensitivity = Se, specificity = Sp)
}

# ==============================================================================
# 5. Simulation Execution Logic
# ==============================================================================

compute_nb <- function(components, z) {
  components["sensitivity"] * components["prevalence"] - (1 - components["specificity"]) * (1 - components["prevalence"]) * (z / (1 - z))
}

compute_true_nb_cr <- function(model, large_data, target_t = 1, pt = 0.2) {
  compute_nb(true_components_cr(large_data, t = target_t, z = pt, risk_score = predict_risk_cr(model, large_data, target_t)), z = pt)
}

compute_aj_nb_cr_fast <- function(valid_data, target_t = 1, pt = 0.2, risk_score) {
  positive <- risk_score >= pt
  q <- mean(positive)

  if (!any(positive)) { return(0) }

  d <- valid_data[positive, , drop = FALSE]
  status_cr <- factor(d$status, levels = c(0, 1, 2), labels = c("censor", "cause1", "cause2"))

  fit <- survival::survfit(survival::Surv(d$obs_time, status_cr) ~ 1)

  cause1 <- which(fit$states == "cause1")
  k <- which(fit$time <= target_t)
  F1 <- if (length(k) == 0) 0 else fit$pstate[max(k), cause1]

  tp <- q * F1
  fp <- q * (1 - F1)
  as.numeric(tp - fp * pt / (1 - pt))
}

run_single_replication_cr <- function(model, valid_data, target_t, thresholds, B, NB_true_list) {
  n <- nrow(valid_data)
  risk_score <- predict_risk_cr(model, valid_data, target_t)
  weights_main <- fit_all_censoring_weights_cr(valid_data, target_t)
  res_main_list <- list()

  for (pt in thresholds) {
    nb_t <- NB_true_list[[as.character(pt)]]
    res_main_list[[as.character(pt)]] <- data.frame(
      threshold = pt,
      NB_true = nb_t,
      NB_aj = compute_aj_nb_cr_fast(valid_data, target_t, pt, risk_score),
      NB_ipcw_cox = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$cox$obs, weights_main$cox$t), pt),
      NB_ipcw_miscox = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$miscox$obs, weights_main$miscox$t), pt),
      NB_ipcw_gam = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$gam$obs, weights_main$gam$t), pt)
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
        fit_all_censoring_weights_cr(boot_data, target_t)
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
            compute_aj_nb_cr_fast(boot_data, target_t, pt, boot_risk),
            compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$cox$obs, weights_boot$cox$t), pt),
            compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$miscox$obs, weights_boot$miscox$t), pt),
            compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$gam$obs, weights_boot$gam$t), pt)
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
      ci_aj_lower = cis[1, 1], ci_aj_upper = cis[2, 1],
      ci_ipcw_cox_lower = cis[1, 2], ci_ipcw_cox_upper = cis[2, 2],
      ci_ipcw_miscox_lower = cis[1, 3], ci_ipcw_miscox_upper = cis[2, 3],
      ci_ipcw_gam_lower = cis[1, 4], ci_ipcw_gam_upper = cis[2, 4]
    )
  }

  final_df <- merge(res_main_df, do.call(rbind, ci_list), by = "threshold")

  final_df$bias_aj <- final_df$NB_aj - final_df$NB_true
  final_df$bias_ipcw_cox <- final_df$NB_ipcw_cox - final_df$NB_true
  final_df$bias_ipcw_miscox <- final_df$NB_ipcw_miscox - final_df$NB_true
  final_df$bias_ipcw_gam <- final_df$NB_ipcw_gam - final_df$NB_true

  final_df$sqerr_aj <- (final_df$NB_aj - final_df$NB_true)^2
  final_df$sqerr_ipcw_cox <- (final_df$NB_ipcw_cox - final_df$NB_true)^2
  final_df$sqerr_ipcw_miscox <- (final_df$NB_ipcw_miscox - final_df$NB_true)^2
  final_df$sqerr_ipcw_gam <- (final_df$NB_ipcw_gam - final_df$NB_true)^2

  return(final_df)
}

# ==============================================================================
# 6. Main Execution Loop
# ==============================================================================

common_params_cr <- list(target_any_event_pct = 0.40, competing_fraction = 1 / 4, target_t = 5, delta_beta = 0)
scenario_names_cr <- c("non_informative_cox", "informative_simple_cox", "informative_nonlinear")
thresholds_cr <- c(0.01, 0.1, 0.25, 0.5, 0.75)
censoring_rates_cr <- c(0.10, 0.20, 0.30)

n_reps <- 1000
B_boot <- 500
all_results_cr <- list()

for (scen_name in scenario_names_cr) {
  for (cens_rate in censoring_rates_cr) {
    cat(sprintf("\n[Running scenario: %s | Censoring: %.1f%%]\n", scen_name, cens_rate * 100))

    # 1. Train Model
    set.seed(123)
    train_data <- generate_competing_risk_tuned(n = 1000, target_any_event_pct = common_params_cr$target_any_event_pct,
                                                competing_fraction = common_params_cr$competing_fraction,
                                                target_t = common_params_cr$target_t, delta_beta = common_params_cr$delta_beta,
                                                censoring_type = "noninformative", target_cens_pct = cens_rate)
    model <- fit_black_box_model_cr(train_data)

    # 2. Generate Full Population
    set.seed(456)
    large_data <- generate_data_by_scenario_cr(scenario = scen_name, n = 1000000,
                                               target_any_event_pct = common_params_cr$target_any_event_pct,
                                               competing_fraction = common_params_cr$competing_fraction,
                                               target_t = common_params_cr$target_t, delta_beta = common_params_cr$delta_beta,
                                               target_cens_pct = cens_rate)

    pop_event_rate <- mean(large_data$obs_time <= common_params_cr$target_t & large_data$status == 1L)
    n_test <- ceiling(100 / pop_event_rate)

    # 3. Calculate True NB directly on Full 1,000,000 Cohort
    NB_true_list <- list()
    for (pt in thresholds_cr) {
      NB_true_list[[as.character(pt)]] <- compute_true_nb_cr(model, large_data, target_t = common_params_cr$target_t, pt = pt)
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
            run_single_replication_cr(model = model, valid_data = large_data[valid_idx, ],
                                      target_t = common_params_cr$target_t, thresholds = thresholds_cr,
                                      B = B_boot, NB_true_list = NB_true_list)
          }, error = function(e) {
            fail_reason_main <<- paste("error:", conditionMessage(e))
            NULL
          }, warning = function(w) {
            fail_reason_main <<- paste("warning:", conditionMessage(w))
            NULL
          })

          # If it succeeded, break the while loop and save the results
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
            threshold = thresholds_cr,
            NB_true = unlist(NB_true_list),
            NB_aj = NA, NB_ipcw_cox = NA, NB_ipcw_miscox = NA, NB_ipcw_gam = NA,
            ci_aj_lower = NA, ci_aj_upper = NA, ci_ipcw_cox_lower = NA, ci_ipcw_cox_upper = NA,
            ci_ipcw_miscox_lower = NA, ci_ipcw_miscox_upper = NA, ci_ipcw_gam_lower = NA, ci_ipcw_gam_upper = NA,
            bias_aj = NA, bias_ipcw_cox = NA, bias_ipcw_miscox = NA, bias_ipcw_gam = NA,
            sqerr_aj = NA, sqerr_ipcw_cox = NA, sqerr_ipcw_miscox = NA, sqerr_ipcw_gam = NA
          )
        }

        p()
        return(res)
      }, future.seed = 123L, future.scheduling = 1)
    })

    # 5. Format and Save Output
    results_df <- do.call(rbind, results)
    results_df$scenario <- scen_name
    results_df$target_cens_pct <- cens_rate
    results_df$target_t <- common_params_cr$target_t

    results_df <- results_df %>% mutate(
      cover_aj = (NB_true >= ci_aj_lower & NB_true <= ci_aj_upper),
      cover_ipcw_cox = (NB_true >= ci_ipcw_cox_lower & NB_true <= ci_ipcw_cox_upper),
      cover_ipcw_miscox = (NB_true >= ci_ipcw_miscox_lower & NB_true <= ci_ipcw_miscox_upper),
      cover_ipcw_gam = (NB_true >= ci_ipcw_gam_lower & NB_true <= ci_ipcw_gam_upper)
    )

    all_results_cr[[paste0(scen_name, "_cens_", gsub("\\.", "", as.character(cens_rate)))]] <- results_df
  }
}

combined_results_cr <- bind_rows(all_results_cr)
results_dir <- "../results"
if (!dir.exists(results_dir)) dir.create(results_dir)

file_results_cr <- file.path(results_dir, "combined_results_competing.csv")
write.csv(combined_results_cr, file_results_cr, row.names = FALSE)
cat("\nCompeting risks simulation script complete.\n")
