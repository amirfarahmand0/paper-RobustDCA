# GBSG case study: informative censoring, censoring model, and decision curve analysis
# with the KM-based and CIPCW-based approaches (Figures 3 and 4, Table S6).
# Data: `rotterdam` and `gbsg` from the R package `survival`.

library(survival)
library(survivalROC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

data(cancer, package = "survival")

# 1. Data preparation ----

# GBSG validation cohort: fractional polynomial terms and prognostic index
gbsg_aligned <- gbsg %>%
  mutate(
    age_trans = age / 100,
    age_fp1   = age_trans^3,
    age_fp2   = (age_trans^3) * log(age_trans),
    sized1    = ifelse(size > 20 & size <= 50, 1, 0),
    sized2    = ifelse(size > 50, 1, 0),
    nodes_fp  = nodes^(-0.5),
    er_trans  = er / 1000
  ) %>%
  mutate(
    PI = 1.07 * age_fp1 + 9.13 * age_fp2 + 0.46 * meno +
      0.23 * sized1 + 0.31 * sized2 - 1.74 * nodes_fp -
      0.34 * er_trans - 0.35 * hormon
  )

# Rotterdam development cohort: node-positive patients.
# RFS is the earlier of recurrence or death. Deaths occurring after recurrence
# follow-up ended (no recurrence, rtime < dtime) are censored at rtime.
rotterdam_aligned <- rotterdam %>%
  filter(nodes > 0) %>%
  mutate(
    late_death = recur == 0 & death == 1 & rtime < dtime,
    rfstime    = ifelse(recur == 1, rtime, pmin(rtime, dtime)),
    status     = ifelse(recur == 1, 1, ifelse(late_death, 0, death)),
    age_trans  = age / 100,
    age_fp1    = age_trans^3,
    age_fp2    = (age_trans^3) * log(age_trans),
    sized1     = ifelse(size == "20-50", 1, 0),
    sized2     = ifelse(size == ">50", 1, 0),
    nodes_fp   = nodes^(-0.5),
    er_trans   = er / 1000
  ) %>%
  mutate(
    PI = 1.07 * age_fp1 + 9.13 * age_fp2 + 0.46 * meno +
      0.23 * sized1 + 0.31 * sized2 - 1.74 * nodes_fp -
      0.34 * er_trans - 0.35 * hormon
  ) %>%
  filter(is.finite(PI))

# 2. Predicted 5-year risk ----

target_t <- 5 * 365.25    # prediction horizon (days)

# Baseline survival at the horizon, estimated in the development cohort
mean_PI_rotterdam <- mean(rotterdam_aligned$PI, na.rm = TRUE)
fit_rotterdam     <- coxph(Surv(rfstime, status) ~ offset(PI - mean_PI_rotterdam),
                           data = rotterdam_aligned)
surv_rotterdam    <- survfit(fit_rotterdam)
S0_t              <- surv_rotterdam$surv[max(which(surv_rotterdam$time <= target_t))]

# Predicted 5-year risk in the validation cohort
gbsg_final <- gbsg_aligned %>%
  mutate(
    abs_risk    = 1 - (S0_t ^ exp(PI - mean_PI_rotterdam)),
    cens_status = ifelse(status == 0, 1, 0)
  )

# 3. Evidence of informative censoring ----

logrank_p <- function(formula, data) {
  test <- survdiff(formula, data = data)
  1 - pchisq(test$chisq, df = length(test$n) - 1)
}

groups <- list(
  All            = gbsg_final,
  Premenopausal  = filter(gbsg_final, meno == 0),
  Postmenopausal = filter(gbsg_final, meno == 1)
)

# Log-rank tests of censoring and RFS by hormonal treatment
logrank_tests <- data.frame(
  group       = names(groups),
  p_censoring = sapply(groups, function(d) logrank_p(Surv(rfstime, cens_status) ~ hormon, d)),
  p_rfs       = sapply(groups, function(d) logrank_p(Surv(rfstime, status) ~ hormon, d)),
  row.names   = NULL
)

# 4. Censoring model and CIPCW weights ----

fit_censoring <- coxph(Surv(rfstime, status == 0) ~ age + meno + size + nodes + er +
                         hormon + hormon:meno, data = gbsg)

# Hazard ratio of hormonal treatment on censoring, by menopausal status
contrast_hr <- function(fit, terms) {
  L   <- as.numeric(names(coef(fit)) %in% terms)
  est <- sum(L * coef(fit))
  se  <- sqrt(drop(t(L) %*% vcov(fit) %*% L))
  exp(c(HR = est, lower = est - 1.96 * se, upper = est + 1.96 * se))
}

censoring_hr <- rbind(
  Premenopausal  = contrast_hr(fit_censoring, "hormon"),
  Postmenopausal = contrast_hr(fit_censoring, c("hormon", "meno:hormon"))
)

# Proportional hazards assumption (Schoenfeld residuals)
ph_test <- cox.zph(fit_censoring)

# G(T_i | X_i): probability of remaining uncensored up to the observed time
G_hat_obs <- pmax(exp(-predict(fit_censoring, type = "expected")), 0.01)

# G(tau | X_i): probability of remaining uncensored up to the horizon
gbsg_t         <- gbsg_final
gbsg_t$rfstime <- target_t
G_hat_t        <- pmax(exp(-predict(fit_censoring, newdata = gbsg_t, type = "expected")), 0.01)

indicator_event <- (gbsg_final$rfstime <= target_t) & (gbsg_final$status == 1)
indicator_surv  <- (gbsg_final$rfstime > target_t)

weight_event <- 1 / G_hat_obs
weight_surv  <- 1 / G_hat_t

weights_used <- c(weight_event[indicator_event], weight_surv[indicator_surv])

# 5. Decision curve analysis ----

thresholds <- seq(0.14, 0.23, by = 0.001)
N          <- nrow(gbsg_final)

# Prevalence at the horizon
prev_km   <- 1 - summary(survfit(Surv(rfstime, status) ~ 1, data = gbsg_final),
                         times = target_t)$surv
prev_ipcw <- sum(weight_event * indicator_event) / N

nb_results <- lapply(thresholds, function(z) {

  # KM-based approach
  roc_km <- survivalROC(Stime = gbsg_final$rfstime, status = gbsg_final$status,
                        marker = gbsg_final$abs_risk, predict.time = target_t,
                        cut.values = z, method = "KM")
  Se_km <- roc_km$TP[2]
  Sp_km <- 1 - roc_km$FP[2]

  NB_km     <- Se_km * prev_km - (1 - Sp_km) * (1 - prev_km) * (z / (1 - z))
  NB_all_km <- prev_km - (1 - prev_km) * (z / (1 - z))

  # CIPCW-based approach
  indicator_pos <- gbsg_final$abs_risk > z
  TP_ipcw <- sum(weight_event * indicator_event * indicator_pos) / N
  FP_ipcw <- sum(weight_surv  * indicator_surv  * indicator_pos) / N

  NB_ipcw     <- TP_ipcw - FP_ipcw * (z / (1 - z))
  NB_all_ipcw <- prev_ipcw - (1 - prev_ipcw) * (z / (1 - z))

  data.frame(threshold = z,
             NB_km = NB_km, NB_all_km = NB_all_km,
             NB_ipcw = NB_ipcw, NB_all_ipcw = NB_all_ipcw)
})

dca_df <- do.call(rbind, nb_results) %>%
  mutate(
    diff_km        = NB_km   - NB_all_km,      # Model - Treat All, KM
    diff_ipcw      = NB_ipcw - NB_all_ipcw,    # Model - Treat All, CIPCW
    diff_model     = NB_ipcw - NB_km,          # CIPCW - KM, Model
    diff_treat_all = NB_all_ipcw - NB_all_km   # CIPCW - KM, Treat All
  )

# Preferred strategy at each threshold
decision_check <- dca_df %>%
  mutate(
    Decision_KM = case_when(
      NB_km > NB_all_km & NB_km > 0      ~ "Model",
      NB_all_km >= NB_km & NB_all_km > 0 ~ "Treat All",
      TRUE                               ~ "Treat None"
    ),
    Decision_CIPCW = case_when(
      NB_ipcw > NB_all_ipcw & NB_ipcw > 0      ~ "Model",
      NB_all_ipcw >= NB_ipcw & NB_all_ipcw > 0 ~ "Treat All",
      TRUE                                     ~ "Treat None"
    ),
    Decision_Changed = Decision_KM != Decision_CIPCW
  )

# Threshold ranges (in %) over which a condition holds
threshold_ranges <- function(condition, thr = decision_check$threshold) {
  if (!any(condition)) return("none")
  runs  <- rle(condition)
  ends  <- cumsum(runs$lengths)
  start <- ends - runs$lengths + 1
  paste(sprintf("%.1f-%.1f", 100 * thr[start[runs$values]], 100 * thr[ends[runs$values]]),
        collapse = ", ")
}

strategies <- c("Model", "Treat All", "Treat None")

decision_summary <- data.frame(
  strategy  = strategies,
  KM        = sapply(strategies, function(s) threshold_ranges(decision_check$Decision_KM == s)),
  CIPCW     = sapply(strategies, function(s) threshold_ranges(decision_check$Decision_CIPCW == s)),
  row.names = NULL
)

# 6. Results ----

# Cohort
print(data.frame(
  n                  = N,
  rfs_events         = sum(gbsg_final$status == 1),
  censored           = sum(gbsg_final$status == 0),
  max_follow_up_yrs  = round(max(gbsg_final$rfstime) / 365.25, 1),
  events_by_5y       = sum(indicator_event),
  event_free_at_5y   = sum(indicator_surv)
))

print(table(menopausal_status  = factor(gbsg_final$meno,   0:1, c("Pre", "Post")),
            hormonal_treatment = factor(gbsg_final$hormon, 0:1, c("No", "Yes"))))

# Log-rank tests by hormonal treatment
print(logrank_tests, digits = 3)

# Censoring model
print(summary(fit_censoring))
print(round(censoring_hr, 2))
print(ph_test)

# CIPCW weights
print(round(c(median = median(weights_used),
              p99    = unname(quantile(weights_used, 0.99)),
              max    = max(weights_used)), 2))

# Prevalence at 5 years
print(round(c(KM = prev_km, CIPCW = prev_ipcw), 3))

# Preferred strategy by threshold range (%)
print(decision_summary)

# Thresholds (%) at which the preferred strategy differs between KM and CIPCW
print(threshold_ranges(decision_check$Decision_Changed))

print(decision_check %>%
        filter(Decision_Changed) %>%
        select(threshold, Decision_KM, Decision_CIPCW, diff_km, diff_ipcw))

# Difference in NB, CIPCW minus KM (range across thresholds)
nb_difference <- rbind(Model       = range(dca_df$diff_model),
                       `Treat All` = range(dca_df$diff_treat_all))
colnames(nb_difference) <- c("min", "max")
print(signif(nb_difference, 3))

# 7. Figures ----

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(1.6, "cm"),
        strip.background = element_rect(fill = "grey95"))

threshold_axis <- scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                     breaks = seq(0.14, 0.23, 0.01))

# Figure 3: RFS and cumulative probability of censoring by hormonal treatment ----

hormon_colours  <- c("Hormonal treatment" = "#1B9E77", "No hormonal treatment" = "#D95F02")
hormon_linetype <- c("Hormonal treatment" = "solid",   "No hormonal treatment" = "dashed")

km_by_group <- function(event_var, cumulative) {
  do.call(rbind, lapply(c(0, 1), function(m) {

    d <- gbsg_final %>%
      mutate(event = .data[[event_var]]) %>%
      filter(meno == m) %>%
      mutate(group = factor(hormon, c(1, 0), names(hormon_colours)))

    fit <- summary(survfit(Surv(rfstime / 365.25, event) ~ group, data = d), censored = TRUE)
    km  <- data.frame(time = fit$time, surv = fit$surv, group = sub("group=", "", fit$strata))
    km  <- rbind(data.frame(time = 0, surv = 1, group = names(hormon_colours)), km)

    km %>% mutate(
      y         = if (cumulative) 1 - surv else surv,
      menopause = c("Premenopausal", "Postmenopausal")[m + 1],
      p         = logrank_p(Surv(rfstime, event) ~ hormon, d)
    )
  })) %>%
    mutate(menopause = factor(menopause, c("Premenopausal", "Postmenopausal")),
           group     = factor(group, names(hormon_colours)))
}

km_panel <- function(km, y_label, label_y) {

  labels <- distinct(km, menopause, p) %>%
    mutate(text = ifelse(p < 0.001, "Log-rank p < 0.001", sprintf("Log-rank p = %.3f", p)))

  ggplot(km, aes(time, y, colour = group, linetype = group)) +
    geom_step(linewidth = 0.7) +
    geom_text(data = labels, aes(x = 0.2, y = label_y, label = text),
              inherit.aes = FALSE, hjust = 0, size = 3.4, colour = "grey20") +
    facet_wrap(~ menopause) +
    scale_colour_manual(values = hormon_colours, name = NULL) +
    scale_linetype_manual(values = hormon_linetype, name = NULL) +
    scale_x_continuous(breaks = 0:7, limits = c(0, 7.3)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Years since surgery", y = y_label) +
    theme_paper
}

figure_censoring <-
  (km_panel(km_by_group("status", cumulative = FALSE),
            "Recurrence-free survival", 0.05) /
     km_panel(km_by_group("cens_status", cumulative = TRUE),
              "Cumulative probability of censoring", 0.95)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

print(figure_censoring)

# Figure 4: impact of accounting for informative censoring on the decision ----

threshold_step <- median(diff(thresholds))

flip_regions <- decision_check %>%
  mutate(flip_group = cumsum(Decision_Changed != lag(Decision_Changed, default = FALSE))) %>%
  filter(Decision_Changed) %>%
  group_by(flip_group) %>%
  summarise(xmin = max(0.14, min(threshold) - threshold_step / 2),
            xmax = min(0.23, max(threshold) + threshold_step / 2),
            .groups = "drop")

panel_a_data <- dca_df %>%
  select(threshold, KM = diff_km, CIPCW = diff_ipcw) %>%
  pivot_longer(-threshold, names_to = "approach", values_to = "difference") %>%
  mutate(approach = factor(approach, c("CIPCW", "KM")))

panel_b_data <- dca_df %>%
  select(threshold, Model = diff_model, `Treat All` = diff_treat_all) %>%
  pivot_longer(-threshold, names_to = "strategy", values_to = "difference")

panel_a <- ggplot() +
  geom_rect(data = flip_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = "Decision changed"),
            alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = panel_a_data,
            aes(threshold, difference, colour = approach, linetype = approach),
            linewidth = 0.9) +
  scale_colour_manual(values = c(CIPCW = "#1B9E77", KM = "#7570B3"), name = NULL) +
  scale_linetype_manual(values = c(CIPCW = "solid", KM = "dashed"), name = NULL) +
  scale_fill_manual(values = c("Decision changed" = "grey60"), name = NULL) +
  threshold_axis +
  labs(x = "Decision threshold",
       y = expression(NB[Model] - NB["Treat All"])) +
  theme_paper

panel_b <- ggplot(panel_b_data, aes(threshold, difference, colour = strategy)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(Model = "#1B9E77", `Treat All` = "#D95F02"), name = NULL) +
  threshold_axis +
  labs(x = "Decision threshold",
       y = expression(NB[CIPCW] - NB[KM])) +
  theme_paper

figure_decision <- (panel_a | panel_b) + plot_annotation(tag_levels = "A")

print(figure_decision)

# 8. Table S6: censoring model (Supporting Information, Section S4) ----

# ER reported per 100 fmol/L for readability (the model is unchanged)
term_info <- data.frame(
  term  = c("age", "meno", "size", "nodes", "er", "hormon", "meno:hormon"),
  label = c("Age", "Postmenopausal (vs premenopausal)", "Tumour size",
            "Number of positive lymph nodes", "Oestrogen receptor",
            "Hormonal treatment (premenopausal women)",
            "Hormonal treatment x postmenopausal (interaction)"),
  unit  = c("per year", "", "per mm", "per node", "per 100 fmol/L", "yes vs no", ""),
  scale = c(1, 1, 1, 1, 100, 1, 1)
)

beta <- coef(fit_censoring)
se   <- sqrt(diag(vcov(fit_censoring)))
ph   <- ph_test$table

model_rows <- term_info |>
  mutate(
    estimate     = beta[term] * scale,
    std_error    = se[term] * scale,
    HR           = exp(estimate),
    lower        = exp(estimate - 1.96 * std_error),
    upper        = exp(estimate + 1.96 * std_error),
    p_value      = 2 * pnorm(-abs(estimate / std_error)),
    schoenfeld_p = ph[match(term, rownames(ph)), "p"]
  )

# Hazard ratio of hormonal treatment on censoring by menopausal status (linear contrasts)
contrast_row <- function(terms, label) {
  L   <- as.numeric(names(beta) %in% terms)
  est <- sum(L * beta)
  s   <- sqrt(drop(t(L) %*% vcov(fit_censoring) %*% L))
  data.frame(term = paste(terms, collapse = " + "), label = label, unit = "yes vs no",
             HR = exp(est), lower = exp(est - 1.96 * s), upper = exp(est + 1.96 * s),
             p_value = 2 * pnorm(-abs(est / s)), schoenfeld_p = NA)
}

table_s6 <- bind_rows(
  model_rows |> select(term, label, unit, HR, lower, upper, p_value, schoenfeld_p),
  contrast_row("hormon", "Hormonal treatment, premenopausal women"),
  contrast_row(c("hormon", "meno:hormon"), "Hormonal treatment, postmenopausal women"),
  data.frame(term = "GLOBAL", label = "Global proportional hazards test", unit = "",
             HR = NA, lower = NA, upper = NA, p_value = NA,
             schoenfeld_p = ph["GLOBAL", "p"])
) |>
  mutate(across(c(HR, lower, upper), ~ round(.x, 3)),
         across(c(p_value, schoenfeld_p), ~ signif(.x, 3)))

print(table_s6)

# Session information ----

print(sessionInfo())

