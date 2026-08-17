# ==============================================================================
# Script Name: 04_sim_competing_risks_benchmark.R
# Purpose:     Benchmark Competing Risks WITHOUT Oracle
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

n_workers <- detectCores() - 1
plan(multisession, workers = n_workers)
handlers("txtprogressbar") 

tune_event_lambdas <- function(lp1, lp2, target_cif1, target_cif2, target_t) { objective <- function(log_lambdas) { lambda1 <- exp(log_lambdas[1]); lambda2 <- exp(log_lambdas[2]); r1 <- lambda1 * exp(lp1); r2 <- lambda2 * exp(lp2); rsum <- r1 + r2; one_minus_surv <- 1 - exp(-rsum * target_t); F1 <- (r1 / rsum) * one_minus_surv; F2 <- (r2 / rsum) * one_minus_surv; rel_err1 <- (mean(F1) - target_cif1) / target_cif1; rel_err2 <- (mean(F2) - target_cif2) / target_cif2; rel_err1^2 + rel_err2^2 }; opt <- optim(par = c(log(0.10), log(0.10)), fn = objective, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-14)); lambda1 <- exp(opt$par[1]); lambda2 <- exp(opt$par[2]); list(lambda1 = lambda1, lambda2 = lambda2) }
tune_censoring_lambda <- function(rate1, rate2, lp_cens, target_cens_pct, target_t) { objective <- function(log_lambda_c) { rc <- exp(log_lambda_c) * exp(lp_cens); rsum <- rate1 + rate2 + rc; prob_cens <- (rc / rsum) * (1 - exp(-rsum * target_t)); (mean(prob_cens) - target_cens_pct)^2 }; opt <- optim(par = log(0.10), fn = objective, method = "Brent", lower = -10, upper = 10); list(lambda_cens = exp(opt$par)) }
generate_competing_risk_tuned <- function(n, target_any_event_pct = 0.30, competing_fraction = 1/3, target_t = 4, delta_beta = 0, beta2_Z1 = -0.50, beta2_Z2 = 0.75, censoring_type = "linear", HRc_Z1 = 2.65, HRc_Z2 = 2.65, target_cens_pct = 0.30) {
  target_cif2 <- target_any_event_pct * competing_fraction; target_cif1 <- target_any_event_pct * (1 - competing_fraction)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2); Z <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma); Z1 <- Z[, 1]; Z2 <- Z[, 2]
  lp1 <- (1 + delta_beta) * Z1 + (1 - delta_beta) * Z2; lp2 <- beta2_Z1 * Z1 + beta2_Z2 * Z2
  event_tuning <- tune_event_lambdas(lp1, lp2, target_cif1, target_cif2, target_t)
  rate1 <- event_tuning$lambda1 * exp(lp1); rate2 <- event_tuning$lambda2 * exp(lp2)
  lp_cens <- if(censoring_type=="linear") log(HRc_Z1) * Z1 + log(HRc_Z2) * Z2 else rep(0,n)
  cens_tuning <- tune_censoring_lambda(rate1, rate2, lp_cens, target_cens_pct, target_t)
  T1 <- rexp(n, rate = rate1); T2 <- rexp(n, rate = rate2); C <- rexp(n, rate = cens_tuning$lambda_cens * exp(lp_cens))
  obs_time <- pmin(T1, T2, C); status <- ifelse(C < pmin(T1, T2), 0L, ifelse(T1 <= T2, 1L, 2L))
  data.frame(obs_time = obs_time, status = status, delta1 = as.integer(status == 1L), delta_any = as.integer(status != 0L), true_time = pmin(T1, T2), true_status = ifelse(T1 <= T2, 1L, 2L), Z1 = Z1, Z2 = Z2)
}
fit_black_box_model_cr <- function(train_data) { list(cox1 = coxph(Surv(obs_time, status == 1) ~ Z1 + Z2, data = train_data), cox2 = coxph(Surv(obs_time, status == 2) ~ Z1 + Z2, data = train_data)) }
step_hazard_eval <- function(bh, query_times) { sapply(query_times, function(u) { idx <- which(bh$time < u); if (length(idx) == 0) 0 else bh$hazard[max(idx)] }) }
predict_risk_cr <- function(model, data, target_t = 1) {
  bh1 <- basehaz(model$cox1, centered = TRUE); bh1 <- bh1[order(bh1$time), ]; bh2 <- basehaz(model$cox2, centered = TRUE); bh2 <- bh2[order(bh2$time), ]
  keep <- bh1$time <= target_t; event_times1 <- bh1$time[keep]; if (length(event_times1) == 0) return(rep(0, nrow(data)))
  Lambda01_at_events <- bh1$hazard[keep]; Lambda01_left <- c(0, head(Lambda01_at_events, -1)); dLambda01 <- Lambda01_at_events - Lambda01_left; Lambda02_left <- step_hazard_eval(bh2, event_times1)
  lp1 <- predict(model$cox1, newdata = data, type = "lp"); lp2 <- predict(model$cox2, newdata = data, type = "lp")
  S_left_mat <- exp(-(outer(exp(lp1), Lambda01_left, "*") + outer(exp(lp2), Lambda02_left, "*"))); CIF1 <- rowSums(S_left_mat * outer(exp(lp1), dLambda01, "*")); pmin(pmax(CIF1, 0), 1)
}
true_components_cr <- function(data, t, z, risk_score) { D1 <- (data$true_time <= t) & (data$true_status == 1L); pos <- risk_score > z; c(prevalence = mean(D1), sensitivity = mean(pos[D1]), specificity = 1 - mean(pos[!D1])) }
compute_nb <- function(components, z) { components["sensitivity"] * components["prevalence"] - (1 - components["specificity"]) * (1 - components["prevalence"]) * (z / (1 - z)) }
compute_true_nb_cr <- function(model, large_data, target_t = 1, pt = 0.2) { compute_nb(true_components_cr(large_data, t = target_t, z = pt, risk_score = predict_risk_cr(model, large_data, target_t)), z = pt) }
compute_aj_nb_cr_fast <- function(valid_data, target_t, pt, risk_score) { positive <- risk_score >= pt; q <- mean(positive); if (!any(positive)) return(0); d <- valid_data[positive, , drop = FALSE]; status_cr <- factor(d$status, levels = c(0, 1, 2), labels = c("censor", "cause1", "cause2")); fit <- survival::survfit(survival::Surv(d$obs_time, status_cr) ~ 1); cause1 <- which(fit$states == "cause1"); k <- which(fit$time <= target_t); F1 <- if (length(k) == 0) 0 else fit$pstate[max(k), cause1]; as.numeric(q * F1 - q * (1 - F1) * pt / (1 - pt)) }
fit_all_censoring_weights_cr <- function(data, t) { dt <- as.data.table(data); out <- list(); fit_cox <- coxph(Surv(obs_time, 1 - delta_any) ~ Z1 + Z2, data = dt); bh_cox <- basehaz(fit_cox, centered = TRUE); setorder(bh_cox, time); get_L0_cox <- function(u) approx(bh_cox$time, bh_cox$hazard, xout = u, method = "linear", rule = 2)$y; lp_cox <- predict(fit_cox, newdata = dt, type = "lp"); out$cox <- list(obs = pmax(exp(-get_L0_cox(dt$obs_time) * exp(lp_cox)), 0.02), t = pmax(exp(-get_L0_cox(t) * exp(lp_cox)), 0.02)); fit_miscox <- coxph(Surv(obs_time, 1 - delta_any) ~ 1 + Z1, data = dt); bh_miscox <- basehaz(fit_miscox, centered = TRUE); setorder(bh_miscox, time); get_L0_miscox <- function(u) approx(bh_miscox$time, bh_miscox$hazard, xout = u, method = "linear", rule = 2)$y; lp_miscox <- predict(fit_miscox, newdata = dt, type = "lp"); out$miscox <- list(obs = pmax(exp(-get_L0_miscox(dt$obs_time) * exp(lp_miscox)), 0.02), t = pmax(exp(-get_L0_miscox(t) * exp(lp_miscox)), 0.02)); fit_gam <- gam(obs_time ~ s(Z1, k = 10) + s(Z2, k = 10), family = cox.ph(), data = dt, weights = (1 - delta_any)); out$gam <- list(obs = pmax(predict(fit_gam, newdata = dt, type = "response"), 0.02)); dt_t <- copy(dt); dt_t$obs_time <- t; out$gam$t <- pmax(predict(fit_gam, newdata = dt_t, type = "response"), 0.02); return(out) }
fast_ipcw_components_cr <- function(data, t, z, risk_score, G_hat_obs, G_hat_t) { dt <- as.data.table(data); weight_event <- dt$delta1 / G_hat_obs; indicator_event <- (dt$obs_time <= t); prev <- mean(weight_event * indicator_event); Se <- if (sum(weight_event * indicator_event) == 0) 0 else sum(weight_event * indicator_event * (risk_score > z)) / sum(weight_event * indicator_event); weight_neg <- ifelse(dt$obs_time > t, 1 / G_hat_t, ifelse((dt$obs_time <= t) & (dt$status == 2L), 1 / G_hat_obs, 0)); Sp <- if (sum(weight_neg) == 0) 1 else sum(weight_neg * (risk_score <= z)) / sum(weight_neg); c(prevalence = prev, sensitivity = Se, specificity = Sp) }

run_single_replication_cr <- function(model, valid_data, target_t, thresholds, B, NB_true_list) {
  n <- nrow(valid_data); risk_score <- predict_risk_cr(model, valid_data, target_t); weights_main <- fit_all_censoring_weights_cr(valid_data, target_t); res_main_list <- list()
  for (pt in thresholds) { res_main_list[[as.character(pt)]] <- data.frame(threshold = pt, NB_true = NB_true_list[[as.character(pt)]], NB_aj = compute_aj_nb_cr_fast(valid_data, target_t, pt, risk_score), NB_ipcw_cox = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$cox$obs, weights_main$cox$t), pt), NB_ipcw_miscox = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$miscox$obs, weights_main$miscox$t), pt), NB_ipcw_gam = compute_nb(fast_ipcw_components_cr(valid_data, target_t, pt, risk_score, weights_main$gam$obs, weights_main$gam$t), pt)) }
  res_main_df <- do.call(rbind, res_main_list); boot_res <- vector("list", B)
  for (b in 1:B) {
    boot_idx <- sample(1:n, n, replace = TRUE); boot_data <- valid_data[boot_idx, ]; boot_risk <- risk_score[boot_idx]; weights_boot <- fit_all_censoring_weights_cr(boot_data, target_t); boot_pts <- list()
    for (pt in thresholds) { boot_pts[[as.character(pt)]] <- c(compute_aj_nb_cr_fast(boot_data, target_t, pt, boot_risk), compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$cox$obs, weights_boot$cox$t), pt), compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$miscox$obs, weights_boot$miscox$t), pt), compute_nb(fast_ipcw_components_cr(boot_data, target_t, pt, boot_risk, weights_boot$gam$obs, weights_boot$gam$t), pt)) }
    boot_res[[b]] <- boot_pts
  }
  ci_list <- list()
  for (pt in thresholds) { pt_char <- as.character(pt); mat <- do.call(rbind, lapply(boot_res, function(x) x[[pt_char]])); cis <- apply(mat, 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE); ci_list[[pt_char]] <- data.frame(threshold = pt, ci_aj_lower = cis[1, 1], ci_aj_upper = cis[2, 1], ci_ipcw_cox_lower = cis[1, 2], ci_ipcw_cox_upper = cis[2, 2], ci_ipcw_miscox_lower = cis[1, 3], ci_ipcw_miscox_upper = cis[2, 3], ci_ipcw_gam_lower = cis[1, 4], ci_ipcw_gam_upper = cis[2, 4]) }
  final_df <- merge(res_main_df, do.call(rbind, ci_list), by = "threshold")
  return(final_df)
}

cat("\n======================================================\n")
cat("Starting Competing Risks Benchmark (NO ORACLE)...\n")
cat("======================================================\n")

t0 <- Sys.time()
cens_rate <- 0.30
thresholds <- c(0.25)
n_reps <- 1000
B_boot <- 500

set.seed(123)
train_data <- generate_competing_risk_tuned(n = 1000, censoring_type = "noninformative", target_cens_pct = cens_rate)
model <- fit_black_box_model_cr(train_data)

set.seed(456)
large_data <- generate_competing_risk_tuned(n = 1000000, censoring_type = "linear", target_cens_pct = cens_rate)
n_test <- ceiling(100 / mean(large_data$obs_time <= 4 & large_data$status == 1L))

NB_true_list <- list()
for (pt in thresholds) {
  NB_true_list[[as.character(pt)]] <- compute_true_nb_cr(model, large_data, target_t = 4, pt = pt)
}

with_progress({
  p <- progressor(steps = n_reps)
  results <- future_lapply(1:n_reps, function(r) {
    res <- run_single_replication_cr(model, large_data[sample(seq_len(nrow(large_data)), n_test, replace = FALSE), ], 4, thresholds, B_boot, NB_true_list)
    p(); return(res)
  }, future.seed = 123L, future.scheduling = 1) 
})

t1 <- Sys.time()
cat(sprintf("\n[Competing Risks Benchmark Completed] Runtime: %.2f minutes\n", as.numeric(difftime(t1, t0, units="mins"))))
