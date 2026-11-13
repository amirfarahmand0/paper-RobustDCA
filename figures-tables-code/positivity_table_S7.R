# Table S7: positivity diagnostic for the standard survival setting with linear informative censoring.
# For each censoring rate and decision threshold, reports the share of the reference population
# classified as high risk and their true probability of remaining uncensored at 5 years.
# Uses the data-generating functions in simulation-code/ (no simulation replications are run).
# Run from the repository root.
suppressPackageStartupMessages({library(data.table); library(survival); library(mvtnorm)})
root <- if (dir.exists("simulation-code")) "." else ".."
exprs <- parse(file.path(root, "simulation-code", "01_sim_standard_survival.R"))
for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                     identical(e[[3]][[1]], as.name("function"))) eval(e, globalenv())
thresholds <- c(0.01, 0.10, 0.25, 0.50, 0.75)
log_hr_cens <- log(2.65)
res <- NULL
for (cr in c(0.10, 0.20, 0.30)) {
  set.seed(123)   # same training data and prediction model as the simulation
  train <- generate_survival_tuned(n = 1000, target_event_pct = 0.3, target_t = 5, delta_beta = 0,
                                   HRc_Z1 = 1, HRc_Z2 = 1, target_cens_pct = cr)
  model <- fit_black_box_model(train)
  set.seed(456)   # same reference population as the simulation
  pop <- generate_data_by_scenario(scenario = "informative_simple_cox", n = 1e6, target_event_pct = 0.3,
                                   target_t = 5, delta_beta = 0, target_cens_pct = cr)
  # recover the calibrated event and censoring rates from the generator's own calibration
  lp <- pop$Z1 + pop$Z2
  eg <- exp(log_hr_cens * (pop$Z1 + pop$Z2))
  lam_e <- exp(optim(log(0.4), function(a) abs(mean(1 - exp(-exp(a) * 5 * exp(lp))) - 0.3),
                     method = "Brent", lower = -10, upper = 10)$par)
  lam_c <- exp(optim(log(0.1), function(a) {
    rc <- exp(a) * eg; re <- lam_e * exp(lp)
    abs(mean(rc / (rc + re) * (1 - exp(-5 * (rc + re)))) - cr)
  }, method = "Brent", lower = -10, upper = 10)$par)
  G5 <- exp(-lam_c * eg * 5)              # true probability of remaining uncensored at 5 years
  risk <- predict_risk(model, pop, 5)
  for (z in thresholds) {
    hi <- risk > z
    res <- rbind(res, data.frame(censoring_rate = cr, threshold = z, pct_above = 100 * mean(hi),
                                 median_G5_above = median(G5[hi]), pct_below_truncation_above = 100 * mean(G5[hi] < 0.02),
                                 median_G5_below = median(G5[!hi])))
  }
}
dir.create(file.path(root, "tables"), showWarnings = FALSE)
write.csv(res, file.path(root, "tables", "Table_S7_positivity.csv"), row.names = FALSE)
print(format(res, digits = 3))
