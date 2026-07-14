# ==============================================================================
# Script Name: 01_sim_standard_survival.R
# Purpose:     Monte Carlo simulation evaluating robust Net Benefit (NB) 
#              estimators under informative censoring for standard survival outcomes.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Package Loading
# ------------------------------------------------------------------------------
library(data.table)
library(survival)
library(survivalROC)
library(mvtnorm)
library(pbapply)
library(future.apply)
library(parallel)
library(dplyr)
library(splines)
library(randomForestSRC)
library(mgcv)

# Initialize parallel processing to handle heavy Monte Carlo iterations
n_workers <- detectCores() - 2
plan(multisession, workers = n_workers)
cat("Workers used:", n_workers, "\n")

# ------------------------------------------------------------------------------
# 2. Data Generating Mechanisms (DGM)
# ------------------------------------------------------------------------------

# DGM 1: Generates standard survival data with calibrated baseline hazards
# Allows for non-informative or linear informative censoring mechanisms.
generate_survival_tuned <- function(n = 800,
                                    target_event_pct = 0.35,
                                    target_t = 1,
                                    delta_beta = 0.45,
                                    HRc_Z1,
                                    HRc_Z2,
                                    target_cens_pct = 0.3) {
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  Z <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  # Noise covariates to evaluate model misspecification
  W1 <- rexp(n, rate = 0.1)
  W2 <- rexp(n, rate = 0.2)
  
  # Event mechanism
  beta_Z1 <- 1 + delta_beta
  beta_Z2 <- 1 - delta_beta
  lp_event <- beta_Z1 * Z1 + beta_Z2 * Z2
  
  # Censoring mechanism
  gamma1 <- log(HRc_Z1)
  gamma2 <- log(HRc_Z2)
  exp_gamma <- exp(gamma1 * Z1 + gamma2 * Z2)
  
  # Calibration for target event rate
  optim_event <- function(log_lambda) {
    lambda <- exp(log_lambda)
    prob <- 1 - exp(-lambda * target_t * exp(lp_event))
    abs(mean(prob) - target_event_pct)
  }
  
  lambda_event <- exp(
    optim(par = log(0.4), fn = optim_event,
          method = "Brent", lower = -10, upper = 10)$par
  )
  
  # Calibration for target censoring rate
  optim_cens <- function(log_lambda_c) {
    lambda_c <- exp(log_lambda_c)
    rate_c <- lambda_c * exp_gamma
    rate_e <- lambda_event * exp(lp_event)
    term1 <- rate_c / (rate_c + rate_e)
    term2 <- 1 - exp(-target_t * (rate_c + rate_e))
    prob_obs_cens <- term1 * term2
    abs(mean(prob_obs_cens) - target_cens_pct)
  }
  
  lambda_cens <- exp(
    optim(par = log(0.1), fn = optim_cens,
          method = "Brent", lower = -10, upper = 10)$par
  )
  
  # Draw latent times
  true_T <- rexp(n, rate = lambda_event * exp(lp_event))
  C <- rexp(n, rate = lambda_cens * exp_gamma)
  
  obs_time <- pmin(true_T, C)
  delta <- as.integer(true_T <= C)
  
  latent_event <- mean(true_T <= target_t)
  obs_event <- mean(obs_time <= target_t & delta == 1)
  censored <- mean(obs_time <= target_t & delta == 0)
  at_risk <- mean(obs_time > target_t)
  
  cat(sprintf("Target latent event (true_T <= %.1f): %.1f%%\n", target_t, 100 * target_event_pct))
  cat(sprintf("Actual latent event: %.1f%%\n", 100 * latent_event))
  cat(sprintf("Observed events by t: %.1f%%\n", 100 * obs_event))
  cat(sprintf("Censored by t: %.1f%%\n", 100 * censored))
  cat(sprintf("Still at risk at t: %.1f%%\n\n", 100 * at_risk))
  
  data.frame(
    obs_time = obs_time,
    delta = delta,
    true_T = true_T,
    Z1 = Z1,
    Z2 = Z2,
    W1 = W1,
    W2 = W2
  )
}

# DGM 2: Generates survival data with non-linear informative censoring
generate_survival_tuned_cens_nonlinear <- function(
    n = 800,
    target_event_pct = 0.35,
    target_t = 1,
    delta_beta = 0.45,
    target_cens_pct = 0.30
) {
  Sigma <- matrix(c(1, 0.5,
                    0.5, 1), 2, 2)
  Z <- rmvnorm(n, mean = c(0, 0), sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  # Noise covariates for harsh misspecification
  W1 <- rexp(n, rate = 0.1)
  W2 <- rexp(n, rate = 0.2)
  
  # Event model (Linear Cox)
  beta_Z1 <- 1 + delta_beta
  beta_Z2 <- 1 - delta_beta
  lp_event <- beta_Z1 * Z1 + beta_Z2 * Z2
  
  optim_event <- function(log_lambda) {
    lambda <- exp(log_lambda)
    prob <- 1 - exp(-lambda * target_t * exp(lp_event))
    abs(mean(prob) - target_event_pct)
  }
  
  lambda_event <- exp(
    optim(par = log(0.4),
          fn = optim_event,
          method = "Brent",
          lower = -10,
          upper = 10)$par
  )
  
  rate_e <- lambda_event * exp(lp_event)
  
  # Nonlinear censoring truth:Informative and Nonlinear
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
  
  lambda_cens <- exp(
    optim(par = log(0.1),
          fn = optim_cens,
          method = "Brent",
          lower = -10,
          upper = 10)$par
  )
  
  true_T <- rexp(n, rate = rate_e)
  C <- rexp(n, rate = lambda_cens * exp_gamma)
  
  obs_time <- pmin(true_T, C)
  delta <- as.integer(true_T <= C)
  
  latent_event <- mean(true_T <= target_t)
  obs_event <- mean(obs_time <= target_t & delta == 1)
  censored <- mean(obs_time <= target_t & delta == 0)
  at_risk <- mean(obs_time > target_t)
  
  cat("=== NONLINEAR CENSORING DGM ===\n")
  cat(sprintf("Target latent event (true_T <= %.1f): %.1f%%\n",
              target_t, 100 * target_event_pct))
  cat(sprintf("Actual latent event (true_T <= %.1f): %.1f%%\n",
              target_t, 100 * latent_event))
  cat(sprintf("Observed event by t: %.1f%%\n", 100 * obs_event))
  cat(sprintf("Censored by t: %.1f%%\n", 100 * censored))
  cat(sprintf("Still at risk at t: %.1f%%\n", 100 * at_risk))
  
  data.frame(
    obs_time = obs_time,
    delta = delta,
    true_T = true_T,
    Z1 = Z1,
    Z2 = Z2,
    W1 = W1,
    W2 = W2
  )
}

# Wrapper to call the appropriate DGM based on scenario string
generate_data_by_scenario <- function(
    scenario,
    n,
    target_event_pct,
    target_t,
    delta_beta,
    target_cens_pct
) {
  if (scenario == "informative_simple_cox") {
    return(generate_survival_tuned(
      n = n, target_event_pct = target_event_pct, target_t = target_t,
      delta_beta = delta_beta, HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = target_cens_pct
    ))
  }
  
  if (scenario == "informative_cox_mixture") {
    return(generate_survival_tuned_cens_mix3(
      n = n, target_event_pct = target_event_pct, target_t = target_t,
      delta_beta = delta_beta, target_cens_pct = target_cens_pct,
      HRc_Z1_vec = c(0.70, 1.40, 2.20), HRc_Z2_vec = c(1.80, 0.85, 1.25)
    ))
  }
  
  if (scenario == "non_informative_cox") {
    return(generate_survival_tuned(
      n = n, target_event_pct = target_event_pct, target_t = target_t,
      delta_beta = delta_beta, HRc_Z1 = 1, HRc_Z2 = 1, target_cens_pct = target_cens_pct
    ))
  }
  
  if (scenario == "informative_nonlinear") {
    return(generate_survival_tuned_cens_nonlinear(
      n = n, target_event_pct = target_event_pct, target_t = target_t,
      delta_beta = delta_beta, target_cens_pct = target_cens_pct
    ))
  }
  
  if (scenario == "tree_event") {
    return(generate_survival_tree_event(
      n = n, target_event_pct = target_event_pct, target_t = target_t, target_cens_pct = target_cens_pct
    ))
  }
  
  stop("Unknown scenario")
}

# ------------------------------------------------------------------------------
# 3. Model Training & Risk Prediction 
# ------------------------------------------------------------------------------

# Train the baseline Cox model on derivation data
fit_black_box_model <- function(train_data) {
  coxph(Surv(obs_time, delta) ~ Z1 + Z2, data = train_data)
}

# Predict dynamic absolute risk at time target_t
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

# ------------------------------------------------------------------------------
# 4. Net Benefit Component Estimation
# ------------------------------------------------------------------------------

# True components (Oracle, uses uncensored latent true_T)
true_components <- function(data, t, z, risk_score) {
  dt <- as.data.table(data)
  D <- dt$true_T <= t
  pos <- risk_score > z
  
  c(
    prevalence = mean(D),
    sensitivity = mean(pos[D]),
    specificity = 1 - mean(pos[!D])
  )
}

# Standard Kaplan-Meier components (assumes non-informative censoring)
km_components <- function(data, t, z, risk_score) {
  dt <- as.data.table(data)
  
  fit <- survfit(Surv(obs_time, delta) ~ 1, data = dt)
  S_t <- if (t %in% fit$time) {
    fit$surv[fit$time == t]
  } else {
    approx(fit$time, fit$surv, xout = t, method = "linear", rule = 2)$y
  }
  prev <- 1 - S_t
  
  roc <- survivalROC(
    Stime = dt$obs_time, status = dt$delta, marker = risk_score,
    predict.time = t, cut.values = z, method = "KM"
  )
  
  Se <- roc$TP[2]
  Sp <- 1 - roc$FP[2]
  
  c(prevalence = prev, sensitivity = Se, specificity = Sp)
}

# Robust CIPCW components (corrects for informative censoring)
ipcw_components <- function(data, t, z, risk_score,
                            cens_method = c("cox", "misspecified_cox", "spline_cox", "rsf", "aft")) {
  cens_method <- match.arg(cens_method)
  dt <- as.data.table(data)
  
  get_G_from_cox <- function(model_fit) {
    bh <- basehaz(model_fit, centered = TRUE)
    setorder(bh, time)
    get_Lambda0 <- function(u) approx(bh$time, bh$hazard, xout = u, method = "linear", rule = 2)$y
    lp <- predict(model_fit, newdata = dt, type = "lp")
    
    Lambda0_obs <- get_Lambda0(dt$obs_time)
    G_hat_obs <- pmax(exp(-Lambda0_obs * exp(lp)), 0.02)
    Lambda0_t <- get_Lambda0(t)
    G_hat_t <- pmax(exp(-Lambda0_t * exp(lp)), 0.02)
    list(G_hat_obs = G_hat_obs, G_hat_t = G_hat_t)
  }
  
  if (cens_method == "cox") {
    cens_fit <- coxph(Surv(obs_time, 1 - delta) ~ Z1 + Z2, data = dt)
    res <- get_G_from_cox(cens_fit)
    G_hat_obs <- res$G_hat_obs; G_hat_t <- res$G_hat_t
  }
  
  if (cens_method == "misspecified_cox") {
    # Misspecified: omits critical confounder Z2
    cens_fit <- coxph(Surv(obs_time, 1 - delta) ~ 1 + Z1, data = dt)
    res <- get_G_from_cox(cens_fit)
    G_hat_obs <- res$G_hat_obs; G_hat_t <- res$G_hat_t
  }
  
  if (cens_method == "spline_cox") {
    cens_fit <- gam(obs_time ~ s(Z1, k = 10) + s(Z2, k = 10), 
                    family = cox.ph(), data = dt, weights = (1 - delta))
    
    G_hat_obs <- predict(cens_fit, newdata = dt, type = "response")
    G_hat_obs <- pmax(G_hat_obs, 0.02)
    
    dt_t <- copy(dt)
    dt_t$obs_time <- t
    G_hat_t <- predict(cens_fit, newdata = dt_t, type = "response")
    G_hat_t <- pmax(G_hat_t, 0.02)
  }
  
  if (cens_method == "aft") {
    aft_fit <- survreg(Surv(obs_time, 1 - delta) ~ Z1 + Z2, data = dt, dist = "weibull")
    lp_aft <- predict(aft_fit, newdata = dt, type = "lp")
    sigma_aft <- aft_fit$scale
    get_G_aft <- function(u, lp, sigma) {
      u <- pmax(u, 1e-10)
      exp(-exp((log(u) - lp) / sigma))
    }
    G_hat_obs <- pmax(get_G_aft(dt$obs_time, lp_aft, sigma_aft), 0.02)
    G_hat_t <- pmax(get_G_aft(t, lp_aft, sigma_aft),  0.02)
  }
  
  if (cens_method == "rsf") {
    cens_fit <- rfsrc(Surv(obs_time, 1 - delta) ~ Z1 + Z2,
                      data = as.data.frame(dt), ntree = 1, nodesize = 3,
                      mtry = 2, nsplit = 10, forest = TRUE)
    
    # Extract Out-Of-Bag (OOB) survival estimates
    rsf_times <- cens_fit$time.interest
    rsf_surv <- cens_fit$survival.oob 
    
    times_padded <- c(0, rsf_times)
    
    G_hat_obs <- sapply(seq_len(nrow(dt)), function(i) {
      surv_padded <- c(1, rsf_surv[i, ])
      approx(times_padded, surv_padded, xout = dt$obs_time[i], 
             method = "constant", f = 0, rule = 2)$y
    })
    G_hat_obs <- pmax(G_hat_obs, 0.02)
    
    G_hat_t <- apply(rsf_surv, 1, function(s) {
      surv_padded <- c(1, s)
      approx(times_padded, surv_padded, xout = t, 
             method = "constant", f = 0, rule = 2)$y
    })
    G_hat_t <- pmax(G_hat_t, 0.02)
  }
  
  # Calculate Weighted Components
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

# ------------------------------------------------------------------------------
# 5. Net Benefit Integration & Inference
# ------------------------------------------------------------------------------

# Master NB Equation
compute_nb <- function(components, z) {
  prev <- components["prevalence"]
  sens <- components["sensitivity"]
  spec <- components["specificity"]
  sens * prev - (1 - spec) * (1 - prev) * (z / (1 - z))
}

# Method-specific wrappers
compute_true_nb <- function(model, large_data, target_t = 1, pt = 0.2) {
  risk_score <- predict_risk(model, large_data, target_t)
  comps <- true_components(large_data, t = target_t, z = pt, risk_score = risk_score)
  compute_nb(comps, z = pt)
}

compute_km_nb <- function(model, valid_data, target_t = 1, pt = 0.2) {
  risk_score <- predict_risk(model, valid_data, target_t)
  comps <- km_components(valid_data, t = target_t, z = pt, risk_score = risk_score)
  compute_nb(comps, z = pt)
}

compute_ipcw_nb <- function(model, valid_data, target_t = 1, pt = 0.2,
                            cens_method = c("cox", "misspecified_cox", "spline_cox", "rsf", "aft")) {
  cens_method <- match.arg(cens_method)
  risk_score <- predict_risk(model, valid_data, target_t)
  comps <- ipcw_components(valid_data, t = target_t, z = pt, risk_score = risk_score, cens_method = cens_method)
  compute_nb(comps, z = pt)
}

# Calculate 90% Bootstrap CI for empirical coverage assessment
bootstrap_ci <- function(valid_data, model, estimator_func, B = 200, target_t = 1, pt = 0.2) {
  n <- nrow(valid_data)
  boot_nbs <- numeric(B)
  for (b in 1:B) {
    boot_idx <- sample(1:n, n, replace = TRUE)
    boot_data <- valid_data[boot_idx, ]
    boot_nbs[b] <- estimator_func(model, boot_data, target_t, pt)
  }
  quantile(boot_nbs, probs = c(0.05, 0.95), na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# 6. Primary Execution Loops
# ------------------------------------------------------------------------------

run_single_replication <- function(model, NB_true, valid_data, target_t = 1, pt = 0.2, B = 200) {
  
  # Point Estimates
  NB_km <- compute_km_nb(model, valid_data, target_t, pt)
  NB_ipcw_cox <- compute_ipcw_nb(model, valid_data, target_t, pt, cens_method = "cox")
  NB_ipcw_miscox <- compute_ipcw_nb(model, valid_data, target_t, pt, cens_method = "misspecified_cox")
  NB_ipcw_rsf <- compute_ipcw_nb(model, valid_data, target_t, pt, cens_method = "rsf")
  NB_ipcw_aft <- compute_ipcw_nb(model, valid_data, target_t, pt, cens_method = "aft")
  NB_ipcw_spline_cox <- compute_ipcw_nb(model, valid_data, target_t, pt, cens_method = "spline_cox")
  
  # Bias 
  bias_km <- (NB_true - NB_km)
  bias_ipcw_cox <- (NB_true - NB_ipcw_cox)
  bias_ipcw_miscox <- (NB_true - NB_ipcw_miscox) 
  bias_ipcw_rsf <- (NB_true - NB_ipcw_rsf)
  bias_ipcw_aft <- (NB_true - NB_ipcw_aft) 
  bias_ipcw_spline_cox <- (NB_true - NB_ipcw_spline_cox)
  
  # Mean Squared Error
  sqerr_km <- (NB_km - NB_true)^2
  sqerr_ipcw_cox <- (NB_ipcw_cox - NB_true)^2
  sqerr_ipcw_miscox <- (NB_ipcw_miscox - NB_true)^2
  sqerr_ipcw_rsf <- (NB_ipcw_rsf - NB_true)^2
  sqerr_ipcw_aft <- (NB_ipcw_aft - NB_true)^2
  sqerr_ipcw_spline_cox <- (NB_ipcw_spline_cox - NB_true)^2
  
  # Confidence Intervals
  ci_km <- bootstrap_ci(valid_data, model, compute_km_nb, B, target_t, pt)
  ci_ipcw_cox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb(m, d, tt, ptt, cens_method = "cox"), B, target_t, pt)
  ci_ipcw_miscox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb(m, d, tt, ptt, cens_method = "misspecified_cox"), B, target_t, pt)
  ci_ipcw_rsf <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb(m, d, tt, ptt, cens_method = "rsf"), B, target_t, pt)
  ci_ipcw_aft <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb(m, d, tt, ptt, cens_method = "aft"), B, target_t, pt)
  ci_ipcw_spline_cox <- bootstrap_ci(valid_data, model, function(m, d, tt, ptt) compute_ipcw_nb(m, d, tt, ptt, cens_method = "spline_cox"), B, target_t, pt)
  
  data.frame(
    NB_km = NB_km, NB_ipcw_cox = NB_ipcw_cox, NB_ipcw_miscox = NB_ipcw_miscox,
    NB_ipcw_rsf = NB_ipcw_rsf, NB_ipcw_aft = NB_ipcw_aft, NB_ipcw_spline_cox = NB_ipcw_spline_cox,
    
    bias_km = bias_km, bias_ipcw_cox = bias_ipcw_cox, bias_ipcw_miscox = bias_ipcw_miscox,
    bias_ipcw_rsf = bias_ipcw_rsf, bias_ipcw_aft = bias_ipcw_aft, bias_ipcw_spline_cox = bias_ipcw_spline_cox,
    
    sqerr_km = sqerr_km, sqerr_ipcw_cox = sqerr_ipcw_cox, sqerr_ipcw_miscox = sqerr_ipcw_miscox,
    sqerr_ipcw_rsf = sqerr_ipcw_rsf, sqerr_ipcw_aft = sqerr_ipcw_aft, sqerr_ipcw_spline_cox = sqerr_ipcw_spline_cox,
    
    ci_km_lower = ci_km[1], ci_km_upper = ci_km[2],
    ci_ipcw_cox_lower = ci_ipcw_cox[1], ci_ipcw_cox_upper = ci_ipcw_cox[2],
    ci_ipcw_miscox_lower = ci_ipcw_miscox[1], ci_ipcw_miscox_upper = ci_ipcw_miscox[2],
    ci_ipcw_rsf_lower = ci_ipcw_rsf[1], ci_ipcw_rsf_upper = ci_ipcw_rsf[2],
    ci_ipcw_aft_lower = ci_ipcw_aft[1], ci_ipcw_aft_upper = ci_ipcw_aft[2], 
    ci_ipcw_spline_cox_lower = ci_ipcw_spline_cox[1], ci_ipcw_spline_cox_upper = ci_ipcw_spline_cox[2]
  )
}

# --- Main Simulation Loop Execution ---

common_params <- list(target_event_pct = 0.3, target_t = 4, delta_beta = 0)
scenario_names <- c("non_informative_cox", "informative_simple_cox", "informative_nonlinear")
thresholds <- c(0.01, 0.1, 0.25, 0.5, 0.75) 
censoring_rates <- c(0.20, 0.40)
all_results <- list()

for (scen_name in scenario_names) {
  for (cens_rate in censoring_rates) {
    cat("Running scenario:", scen_name, "with censoring rate:", cens_rate, "\n")
    
    # Train model on non-informative dataset
    set.seed(123)
    train_data <- generate_survival_tuned(
      n = 1000, target_event_pct = common_params$target_event_pct, target_t = common_params$target_t,
      delta_beta = common_params$delta_beta, HRc_Z1 = 1, HRc_Z2 = 1, target_cens_pct = cens_rate
    )
    model <- fit_black_box_model(train_data)
    
    # Generate population reference dataset
    set.seed(123)
    large_data <- generate_data_by_scenario(
      scenario = scen_name, n = 1000000, target_event_pct = common_params$target_event_pct,
      target_t = common_params$target_t, delta_beta = common_params$delta_beta, target_cens_pct = cens_rate
    )
    
    # Dynamic Sample Size calculation (targeting 100 expected events)
    pop_event_rate <- mean(large_data$obs_time <= common_params$target_t & large_data$delta == 1)
    n_test <- ceiling(100 / pop_event_rate)
    
    cat(sprintf("  -> Population event rate by t=%.1f: %.2f%%\n", common_params$target_t, pop_event_rate * 100))
    cat(sprintf("  -> Dynamic test sample size set to: %d\n", n_test))
    
    scenario_rate_results <- list()
    
    for (pt in thresholds) {
      cat("  Threshold:", pt, "\n")
      
      NB_true <- compute_true_nb(model, large_data, target_t = common_params$target_t, pt = pt)
      
      # Parallel execution over 100 Monte Carlo replications
      results <- future_lapply(1:100 ,function(r) {
        valid_idx <- sample(seq_len(nrow(large_data)), n_test, replace = FALSE)
        valid_data <- large_data[valid_idx, ]
        run_single_replication(model = model, NB_true = NB_true, valid_data = valid_data, target_t = common_params$target_t, pt = pt, B = 100)
      }, future.seed = TRUE)
      
      results_df <- do.call(rbind, results)
      results_df$NB_true <- NB_true
      results_df$scenario <- scen_name
      results_df$target_cens_pct <- cens_rate
      results_df$threshold <- pt
      results_df$target_t <- common_params$target_t
      
      # Determine empirical coverage mapping
      results_df <- results_df %>%
        mutate(
          cover_km = (NB_true >= ci_km_lower & NB_true <= ci_km_upper),
          cover_ipcw_cox = (NB_true >= ci_ipcw_cox_lower & NB_true <= ci_ipcw_cox_upper),
          cover_ipcw_miscox = (NB_true >= ci_ipcw_miscox_lower & NB_true <= ci_ipcw_miscox_upper),
          cover_ipcw_rsf = (NB_true >= ci_ipcw_rsf_lower & NB_true <= ci_ipcw_rsf_upper),
          cover_ipcw_aft = (NB_true >= ci_ipcw_aft_lower & NB_true <= ci_ipcw_aft_upper),
          cover_ipcw_spline_cox = (NB_true >= ci_ipcw_spline_cox_lower & NB_true <= ci_ipcw_spline_cox_upper)
        )
      
      scenario_rate_results[[paste0("thr_", pt)]] <- results_df
    }
    all_results[[paste0(scen_name, "_cens_", gsub("\\.", "", as.character(cens_rate)))]] <- bind_rows(scenario_rate_results)
  }
}

# ------------------------------------------------------------------------------
# 7. Data Aggregation & Output
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 7. Data Aggregation & Output
# ------------------------------------------------------------------------------

combined_results <- bind_rows(all_results)

# Ensure the results directory exists
if (!dir.exists("results")) { dir.create("results") }

# File paths
file_results <- "results/combined_results_standard.csv"

# Check and save raw results
if (!file.exists(file_results)) {
  write.csv(combined_results, file_results, row.names = FALSE)
  cat("\nSaved:", file_results)
} else {
  cat("\nFile already exists, skipping save:", file_results)
}

cat("\nStandard survival simulation script complete.\n")