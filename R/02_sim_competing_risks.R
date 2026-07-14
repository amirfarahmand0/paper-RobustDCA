# ==============================================================================
# Script Name: 02_sim_competing_risks.R
# Purpose:     Monte Carlo simulation evaluating robust Net Benefit (NB) 
#              estimators under informative censoring for competing risk outcomes.
# Application: Robust Decision Curve Analysis (DCA)
# ==============================================================================

# ---- 1. Setup and Packages ----
library(survival)
library(data.table)
library(mvtnorm)
library(dplyr)
library(tidyr)
library(mgcv)
library(randomForestSRC)
library(parallel)
library(future.apply)

n_workers <- detectCores() - 2
plan(multisession, workers = n_workers)
cat("Workers used:", n_workers, "\n")

# ---- 2. Input Validation Helper ----
validate_cr_inputs <- function(target_any_event_pct,
                               competing_fraction,
                               target_cif1,
                               target_cif2,
                               target_cens_pct) {
  
  if (!(target_any_event_pct > 0 && target_any_event_pct < 1)) stop("target_any_event_pct must satisfy 0 < target_any_event_pct < 1.")
  if (!(competing_fraction >= 0 && competing_fraction <= 1)) stop("competing_fraction must satisfy 0 <= competing_fraction <= 1.")
  if (!(target_cif1 > 0)) stop("target_cif1 must be > 0.")
  if (!(target_cif2 > 0)) stop("target_cif2 must be > 0.")
  if (!(target_cif1 + target_cif2 < 1)) stop("target_cif1 + target_cif2 must be < 1.")
  if (!(target_cens_pct > 0 && target_cens_pct < 1)) stop("target_cens_pct must satisfy 0 < target_cens_pct < 1.")
  
  invisible(TRUE)
}

# ---- 3. Baseline Hazard Calibration Functions ----
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
    mean_F1 <- mean(F1)
    mean_F2 <- mean(F2)
    rel_err1 <- (mean_F1 - target_cif1) / target_cif1
    rel_err2 <- (mean_F2 - target_cif2) / target_cif2
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

# ---- 4. Competing Risks Data Generating Mechanism (DGM) ----
generate_competing_risk_tuned <- function(n, target_any_event_pct = 0.30, competing_fraction = 1/3, target_t = 4, 
                                          delta_beta = 0, beta2_Z1 = -0.50, beta2_Z2 = 0.75, 
                                          censoring_type = c("noninformative", "linear", "nonlinear"), 
                                          HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = 0.30) {
  censoring_type <- match.arg(censoring_type)
  target_cif2 <- target_any_event_pct * competing_fraction
  target_cif1 <- target_any_event_pct * (1 - competing_fraction)
  validate_cr_inputs(target_any_event_pct, competing_fraction, target_cif1, target_cif2, target_cens_pct)
  
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Z <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]; Z2 <- Z[, 2]
  W1 <- rexp(n, rate = 0.1); W2 <- rexp(n, rate = 0.2)
  
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
                    nonlinear = 1.5 * tanh(Z1) - 2.0 * exp(-(Z2^2))
  )
  
  cens_tuning <- tune_censoring_lambda(rate1, rate2, lp_cens, target_cens_pct, target_t)
  rate_cens <- cens_tuning$lambda_cens * exp(lp_cens)
  
  T1 <- rexp(n, rate = rate1)
  T2 <- rexp(n, rate = rate2)
  C  <- rexp(n, rate = rate_cens)
  
  true_time <- pmin(T1, T2)
  true_status <- ifelse(T1 <= T2, 1L, 2L)
  obs_time <- pmin(T1, T2, C)
  status <- ifelse(C < pmin(T1, T2), 0L, ifelse(T1 <= T2, 1L, 2L))
  
  data.frame(
    obs_time = obs_time, status = status, delta1 = as.integer(status == 1L),
    delta2 = as.integer(status == 2L), delta_any = as.integer(status != 0L),
    T1 = T1, T2 = T2, C = C, true_time = true_time, true_status = true_status,
    Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2
  )
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

# ---- 5. Black-box Prediction & Risk Calculation ----
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

# ---- 6. Net Benefit Component Estimation ----
true_components_cr <- function(data, t, z, risk_score) {
  D1 <- (data$true_time <= t) & (data$true_status == 1L)
  pos <- risk_score > z
  c(prevalence = mean(D1), sensitivity = mean(pos[D1]), specificity = 1 - mean(pos[!D1]))
}

marginal_censoring_survival <- function(dt) { survfit(Surv(obs_time, as.integer(status == 0L)) ~ 1, data = dt) }

km_step_eval <- function(fit, query_times, left_continuous = FALSE) {
  sapply(query_times, function(u) {
    idx <- if (left_continuous) which(fit$time < u) else which(fit$time <= u)
    if (length(idx) == 0) 1 else fit$surv[max(idx)]
  })
}

aj_components_cr <- function(data, t, z, risk_score) {
  dt <- data
  aj_fit <- survfit(Surv(obs_time, factor(status, levels = c(0, 1, 2))) ~ 1, data = dt)
  state_idx <- which(aj_fit$states == "1")
  
  idx <- which(aj_fit$time <= t)
  prev <- if (length(idx) == 0) 0 else aj_fit$pstate[max(idx), state_idx]
  
  g_fit <- marginal_censoring_survival(dt)
  G_at_obs_minus <- pmax(km_step_eval(g_fit, dt$obs_time, left_continuous = TRUE), 0.02)
  G_at_t <- pmax(km_step_eval(g_fit, rep(t, nrow(dt)), left_continuous = FALSE), 0.02)
  
  weight_cause1 <- ifelse(dt$status == 1L & dt$obs_time <= t, 1 / G_at_obs_minus, 0)
  num_se <- sum(weight_cause1 * (risk_score > z))
  Se <- if (sum(weight_cause1) == 0) 0 else num_se / sum(weight_cause1)
  
  weight_neg <- ifelse(dt$obs_time > t, 1 / G_at_t, ifelse((dt$obs_time <= t) & (dt$status == 2L), 1 / G_at_obs_minus, 0))
  num_sp <- sum(weight_neg * (risk_score <= z))
  Sp <- if (sum(weight_neg) == 0) 1 else num_sp / sum(weight_neg)
  
  c(prevalence = prev, sensitivity = Se, specificity = Sp)
}

ipcw_components_cr <- function(data, t, z, risk_score, cens_method = c("cox", "misspecified_cox", "spline_cox", "rsf", "aft")) {
  cens_method <- match.arg(cens_method)
  dt <- as.data.table(data)
  
  get_G_from_cox <- function(model_fit) {
    bh <- basehaz(model_fit, centered = TRUE); setorder(bh, time)
    get_Lambda0 <- function(u) approx(bh$time, bh$hazard, xout = u, method = "linear", rule = 2)$y
    lp <- predict(model_fit, newdata = dt, type = "lp")
    list(G_hat_obs = pmax(exp(-get_Lambda0(dt$obs_time) * exp(lp)), 0.02), G_hat_t = pmax(exp(-get_Lambda0(t) * exp(lp)), 0.02))
  }
  
  if (cens_method == "cox") res <- get_G_from_cox(coxph(Surv(obs_time, 1 - delta_any) ~ Z1 + Z2, data = dt))
  if (cens_method == "misspecified_cox") res <- get_G_from_cox(coxph(Surv(obs_time, 1 - delta_any) ~ 1 + Z1, data = dt))
  
  if (cens_method %in% c("cox", "misspecified_cox")) { G_hat_obs <- res$G_hat_obs; G_hat_t <- res$G_hat_t }
  
  if (cens_method == "spline_cox") {
    cens_fit <- gam(obs_time ~ s(Z1, k = 10) + s(Z2, k = 10), family = cox.ph(), data = dt, weights = (1 - delta_any))
    G_hat_obs <- pmax(predict(cens_fit, newdata = dt, type = "response"), 0.02)
    dt_t <- copy(dt); dt_t$obs_time <- t
    G_hat_t <- pmax(predict(cens_fit, newdata = dt_t, type = "response"), 0.02)
  }
  
  if (cens_method == "aft") {
    aft_fit <- survreg(Surv(obs_time, 1 - delta_any) ~ Z1 + Z2, data = dt, dist = "weibull")
    lp_aft <- predict(aft_fit, newdata = dt, type = "lp")
    get_G_aft <- function(u) exp(-exp((log(pmax(u, 1e-10)) - lp_aft) / aft_fit$scale))
    G_hat_obs <- pmax(get_G_aft(dt$obs_time), 0.02)
    G_hat_t <- pmax(get_G_aft(t), 0.02)
  }
  
  if (cens_method == "rsf") {
    cens_fit <- rfsrc(Surv(obs_time, 1 - delta_any) ~ Z1 + Z2, data = as.data.frame(dt), ntree = 1, nodesize = 3, mtry = 2, nsplit = 10, forest = TRUE)
    times_padded <- c(0, cens_fit$time.interest)
    G_hat_obs <- pmax(sapply(seq_len(nrow(dt)), function(i) approx(times_padded, c(1, cens_fit$survival.oob[i, ]), xout = dt$obs_time[i], method = "constant", f = 0, rule = 2)$y), 0.02)
    G_hat_t <- pmax(apply(cens_fit$survival.oob, 1, function(s) approx(times_padded, c(1, s), xout = t, method = "constant", f = 0, rule = 2)$y), 0.02)
  }
  
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

# ---- 7. Net Benefit & Inference ----
compute_nb <- function(components, z) {
  components["sensitivity"] * components["prevalence"] - (1 - components["specificity"]) * (1 - components["prevalence"]) * (z / (1 - z))
}

compute_true_nb_cr <- function(model, large_data, target_t = 1, pt = 0.2) {
  compute_nb(true_components_cr(large_data, t = target_t, z = pt, risk_score = predict_risk_cr(model, large_data, target_t)), z = pt)
}

compute_aj_nb_cr <- function(model, valid_data, target_t = 1, pt = 0.2) {
  valid_data$risk_score <- predict_risk_cr(model, valid_data, target_t)
  valid_data$status_cr <- factor(valid_data$status, levels = c(0, 1, 2), labels = c("censor", "cause1", "cause2"))
  dca_res <- dcurves::dca(Surv(obs_time, status_cr) ~ risk_score, data = valid_data, time = target_t, thresholds = pt)
  as.numeric(dca_res$dca$net_benefit[dca_res$dca$variable == "risk_score"])
}

compute_ipcw_nb_cr <- function(model, valid_data, target_t = 1, pt = 0.2, cens_method = c("cox", "misspecified_cox", "spline_cox", "rsf", "aft")) {
  compute_nb(ipcw_components_cr(valid_data, t = target_t, z = pt, risk_score = predict_risk_cr(model, valid_data, target_t), cens_method = match.arg(cens_method)), z = pt)
}

bootstrap_ci <- function(valid_data, model, estimator_func, B = 200, target_t = 1, pt = 0.2) {
  boot_nbs <- numeric(B)
  for (b in 1:B) {
    boot_idx <- sample(1:nrow(valid_data), nrow(valid_data), replace = TRUE)
    boot_nbs[b] <- estimator_func(model, valid_data[boot_idx, ], target_t, pt)
  }
  quantile(boot_nbs, probs = c(0.05, 0.95), na.rm = TRUE)
}

run_single_replication_cr <- function(model, NB_true, valid_data, target_t = 4, pt = 0.2, B = 200) {
  NB_aj <- compute_aj_nb_cr(model, valid_data, target_t, pt)
  NB_ipcw_cox <- compute_ipcw_nb_cr(model, valid_data, target_t, pt, cens_method = "cox")
  NB_ipcw_miscox <- compute_ipcw_nb_cr(model, valid_data, target_t, pt, cens_method = "misspecified_cox")
  NB_ipcw_rsf <- compute_ipcw_nb_cr(model, valid_data, target_t, pt, cens_method = "rsf")
  NB_ipcw_aft <- compute_ipcw_nb_cr(model, valid_data, target_t, pt, cens_method = "aft")
  NB_ipcw_spline_cox <- compute_ipcw_nb_cr(model, valid_data, target_t, pt, cens_method = "spline_cox")
  
  ci_aj <- bootstrap_ci(valid_data, model, compute_aj_nb_cr, B, target_t, pt)
  ci_ipcw_cox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb_cr(m, d, tt, ptt, cens_method = "cox"), B, target_t, pt)
  ci_ipcw_miscox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb_cr(m, d, tt, ptt, cens_method = "misspecified_cox"), B, target_t, pt)
  ci_ipcw_rsf <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb_cr(m, d, tt, ptt, cens_method = "rsf"), B, target_t, pt)
  ci_ipcw_aft <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb_cr(m, d, tt, ptt, cens_method = "aft"), B, target_t, pt)
  ci_ipcw_spline_cox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb_cr(m, d, tt, ptt, cens_method = "spline_cox"), B, target_t, pt)
  
  data.frame(
    NB_aj = NB_aj, NB_ipcw_cox = NB_ipcw_cox, NB_ipcw_miscox = NB_ipcw_miscox, NB_ipcw_rsf = NB_ipcw_rsf, NB_ipcw_aft = NB_ipcw_aft, NB_ipcw_spline_cox = NB_ipcw_spline_cox,
    bias_aj = NB_true - NB_aj, bias_ipcw_cox = NB_true - NB_ipcw_cox, bias_ipcw_miscox = NB_true - NB_ipcw_miscox, bias_ipcw_rsf = NB_true - NB_ipcw_rsf, bias_ipcw_aft = NB_true - NB_ipcw_aft, bias_ipcw_spline_cox = NB_true - NB_ipcw_spline_cox,
    sqerr_aj = (NB_aj - NB_true)^2, sqerr_ipcw_cox = (NB_ipcw_cox - NB_true)^2, sqerr_ipcw_miscox = (NB_ipcw_miscox - NB_true)^2, sqerr_ipcw_rsf = (NB_ipcw_rsf - NB_true)^2, sqerr_ipcw_aft = (NB_ipcw_aft - NB_true)^2, sqerr_ipcw_spline_cox = (NB_ipcw_spline_cox - NB_true)^2,
    ci_aj_lower = ci_aj[1], ci_aj_upper = ci_aj[2], ci_ipcw_cox_lower = ci_ipcw_cox[1], ci_ipcw_cox_upper = ci_ipcw_cox[2], ci_ipcw_miscox_lower = ci_ipcw_miscox[1], ci_ipcw_miscox_upper = ci_ipcw_miscox[2], ci_ipcw_rsf_lower = ci_ipcw_rsf[1], ci_ipcw_rsf_upper = ci_ipcw_rsf[2], ci_ipcw_aft_lower = ci_ipcw_aft[1], ci_ipcw_aft_upper = ci_ipcw_aft[2], ci_ipcw_spline_cox_lower = ci_ipcw_spline_cox[1], ci_ipcw_spline_cox_upper = ci_ipcw_spline_cox[2]
  )
}

# ---- 8. Main Simulation Loop ----
common_params_cr <- list(target_any_event_pct = 0.3, competing_fraction = 1 / 3, target_t = 4, delta_beta = 0)
scenario_names_cr_new <- c("non_informative_cox", "informative_simple_cox", "informative_nonlinear")
thresholds_cr <- c(0.01, 0.1, 0.25, 0.5, 0.75)
censoring_rates_cr <- c(0.2, 0.4)
all_results_cr <- list()

for (scen_name in scenario_names_cr_new) {
  for (cens_rate in censoring_rates_cr) {
    cat("Running scenario:", scen_name, "with censoring rate:", cens_rate, "\n")
    
    set.seed(123)
    train_data <- generate_competing_risk_tuned(n = 1000, target_any_event_pct = common_params_cr$target_any_event_pct, competing_fraction = common_params_cr$competing_fraction, target_t = common_params_cr$target_t, delta_beta = common_params_cr$delta_beta, censoring_type = "noninformative", target_cens_pct = cens_rate)
    model <- fit_black_box_model_cr(train_data)
    
    set.seed(123)
    large_data <- generate_data_by_scenario_cr(scenario = scen_name, n = 1000000, target_any_event_pct = common_params_cr$target_any_event_pct, competing_fraction = common_params_cr$competing_fraction, target_t = common_params_cr$target_t, delta_beta = common_params_cr$delta_beta, target_cens_pct = cens_rate)
    
    pop_event_rate <- mean(large_data$obs_time <= common_params_cr$target_t & large_data$status == 1L)
    n_test <- ceiling(100 / pop_event_rate)
    
    cat(sprintf("  -> Population cause-1 event rate by t=%.1f: %.2f%%\n", common_params_cr$target_t, pop_event_rate * 100))
    cat(sprintf("  -> Dynamic test sample size set to: %d\n", n_test))
    
    scenario_rate_results <- list()
    
    for (pt in thresholds_cr) {
      cat("  Threshold:", pt, "\n")
      NB_true <- compute_true_nb_cr(model, large_data, target_t = common_params_cr$target_t, pt = pt)
      
      results <- future_lapply(1:100, function(r) {
        valid_idx <- sample(seq_len(nrow(large_data)), n_test, replace = FALSE)
        run_single_replication_cr(model = model, NB_true = NB_true, valid_data = large_data[valid_idx, ], target_t = common_params_cr$target_t, pt = pt, B = 100)
      }, future.seed = TRUE)
      
      results_df <- do.call(rbind, results)
      results_df$NB_true <- NB_true
      results_df$scenario <- scen_name
      results_df$target_cens_pct <- cens_rate
      results_df$threshold <- pt
      results_df$target_t <- common_params_cr$target_t
      
      results_df <- results_df %>% mutate(
        cover_aj = (NB_true >= ci_aj_lower & NB_true <= ci_aj_upper),
        cover_ipcw_cox = (NB_true >= ci_ipcw_cox_lower & NB_true <= ci_ipcw_cox_upper),
        cover_ipcw_miscox = (NB_true >= ci_ipcw_miscox_lower & NB_true <= ci_ipcw_miscox_upper),
        cover_ipcw_rsf = (NB_true >= ci_ipcw_rsf_lower & NB_true <= ci_ipcw_rsf_upper),
        cover_ipcw_aft = (NB_true >= ci_ipcw_aft_lower & NB_true <= ci_ipcw_aft_upper),
        cover_ipcw_spline_cox = (NB_true >= ci_ipcw_spline_cox_lower & NB_true <= ci_ipcw_spline_cox_upper)
      )
      scenario_rate_results[[paste0("thr_", pt)]] <- results_df
    }
    all_results_cr[[paste0(scen_name, "_cens_", gsub("\\.", "", as.character(cens_rate)))]] <- bind_rows(scenario_rate_results)
  }
}

# ---- 9. Data Aggregation & Output ----
combined_results_cr <- bind_rows(all_results_cr)

# Ensure the results directory exists
if (!dir.exists("results")) { dir.create("results") }

# File paths
file_results_cr <- "results/combined_results_competing.csv"

# Check and save raw results
if (!file.exists(file_results_cr)) {
  write.csv(combined_results_cr, file_results_cr, row.names = FALSE)
  cat("\nSaved:", file_results_cr)
} else {
  cat("\nFile already exists, skipping save:", file_results_cr)
}

cat("\nCompeting risks simulation script complete.\n")