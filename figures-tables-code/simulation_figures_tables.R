# Simulation figures (Figures 1, 2, S1-S4) and tables (Tables S2-S5).
# Reads simulation-results/ and prints each figure and table. Run from the repository root.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# 0. Paths and settings ----
# Works whether the script is run from the repository root or from figures-tables-code/
root <- if (dir.exists("simulation-results")) "." else ".."
path_results <- file.path(root, "simulation-results")

nominal       <- 0.90     # nominal level of the bootstrap confidence intervals
z_nominal     <- qnorm(1 - (1 - nominal) / 2)
main_rates    <- c(0.10, 0.30)          # censoring rates shown in Figures 1-2
all_rates     <- c(0.10, 0.20, 0.30)

# Figure option
show_mean <- FALSE        # TRUE adds a marker for the mean error (the bias) to each boxplot

# 1. Read results and reshape to one row per replication x method ----
res_std <- read.csv(file.path(path_results, "combined_results_standard.csv"))
res_cr  <- read.csv(file.path(path_results, "combined_results_competing.csv"))

# Method keys used in the simulation output
method_keys <- c("ref", "ipcw_cox", "ipcw_miscox", "ipcw_gam")
method_labels <- function(setting) {
  c(ref         = if (setting == "Standard") "KM" else "AJ",
    ipcw_cox    = "CIPCW (Cox)",
    ipcw_miscox = "CIPCW (Misspecified Cox)",
    ipcw_gam    = "CIPCW (GAM)")
}

mechanism_labels <- c(
  non_informative_cox    = "Non-informative",
  informative_simple_cox = "Linear informative",
  informative_nonlinear  = "Nonlinear informative"
)

to_long <- function(d, setting) {
  ref_key <- if (setting == "Standard") "km" else "aj"
  labs    <- method_labels(setting)
  bind_rows(lapply(method_keys, function(k) {
    col <- if (k == "ref") ref_key else k
    tibble(
      setting   = setting,
      mechanism = unname(mechanism_labels[d$scenario]),
      cens_rate = d$target_cens_pct,
      threshold = d$threshold,
      method    = unname(labs[k]),
      method_key = k,
      nb_true   = d$NB_true,
      nb_hat    = d[[paste0("NB_", col)]],
      ci_lower  = d[[paste0("ci_", col, "_lower")]],
      ci_upper  = d[[paste0("ci_", col, "_upper")]]
    )
  }))
}

long <- bind_rows(to_long(res_std, "Standard"), to_long(res_cr, "Competing")) |>
  mutate(err = nb_hat - nb_true,
         # Bootstrap SE approximated from the 90% interval width (bootstrap SD not stored)
         se_boot = (ci_upper - ci_lower) / (2 * z_nominal))

# 2. Performance measures with Monte Carlo standard errors (Morris et al. 2019) ----
perf <- long |>
  group_by(setting, mechanism, cens_rate, threshold, method_key, method) |>
  summarise(
    n_sim    = sum(!is.na(nb_hat)),
    n_failed = sum(is.na(nb_hat)),
    nb_true  = first(nb_true),
    bias     = mean(err, na.rm = TRUE),
    emp_se   = sd(nb_hat, na.rm = TRUE),
    mse      = mean(err^2, na.rm = TRUE),
    mcse_mse = sqrt(sum((err^2 - mse)^2, na.rm = TRUE) / (n_sim * (n_sim - 1))),
    coverage = mean(ci_lower <= nb_true & nb_true <= ci_upper, na.rm = TRUE),
    mod_se   = sqrt(mean(se_boot^2, na.rm = TRUE)),
    var_se2  = var(se_boot^2, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    mcse_bias      = emp_se / sqrt(n_sim),
    rel_bias       = 100 * bias / nb_true,
    mcse_rel_bias  = 100 * mcse_bias / abs(nb_true),
    mcse_emp_se    = emp_se / sqrt(2 * (n_sim - 1)),
    rmse           = sqrt(mse),
    mcse_rmse      = mcse_mse / (2 * rmse),
    mcse_coverage  = sqrt(coverage * (1 - coverage) / n_sim),
    rel_err_mod_se = 100 * (mod_se / emp_se - 1),
    mcse_rel_err_mod_se = 100 * (mod_se / emp_se) *
      sqrt(var_se2 / (4 * n_sim * mod_se^4) + 1 / (2 * (n_sim - 1)))
  ) |>
  select(-var_se2, -mse, -mcse_mse) |>
  arrange(setting, mechanism, cens_rate, threshold, method_key)

# Report failed replications
failed <- perf |> filter(n_failed > 0)
if (nrow(failed) > 0) {
  message("Replications with missing estimates:")
  print(failed |> select(setting, mechanism, cens_rate, threshold, method, n_failed))
} else {
  message("No failed replications: all estimates available.")
}

# 3. Figures ----
cols   <- c("KM" = "#7570b3", "AJ" = "#7570b3", "CIPCW (Cox)" = "#d95f02",
            "CIPCW (Misspecified Cox)" = "#e7298a", "CIPCW (GAM)" = "#1b9e77")
ltys   <- c("KM" = "solid", "AJ" = "solid", "CIPCW (Cox)" = "dashed",
            "CIPCW (Misspecified Cox)" = "dotted", "CIPCW (GAM)" = "dotdash")
# Shapes add a second visual cue for colour-blind readers
shapes <- c("KM" = 16, "AJ" = 16, "CIPCW (Cox)" = 17,
            "CIPCW (Misspecified Cox)" = 15, "CIPCW (GAM)" = 18)

theme_fig <- theme_bw(base_size = 9) +
  theme(
    strip.background   = element_rect(fill = "grey88", colour = "grey30", linewidth = 0.4),
    strip.text         = element_text(face = "bold", size = 7.5),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    axis.text          = element_text(size = 7),
    axis.title         = element_text(size = 8.5),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 8),
    legend.key.width   = unit(1.8, "lines"),
    plot.tag           = element_text(face = "bold", size = 11)
  )

# Select and label the rows for one figure
prep_fig <- function(d, setting_, mechs, rates, keys, mech_display, rate_suffix = " Censoring") {
  lab <- method_labels(setting_)[keys]
  d |>
    filter(setting == setting_, mechanism %in% mechs, cens_rate %in% rates,
           method_key %in% keys) |>
    mutate(
      method    = factor(method, levels = unname(lab)),
      mechanism = factor(unname(mech_display[mechanism]), levels = unname(mech_display[mechs])),
      rate      = factor(paste0(100 * cens_rate, "%", rate_suffix),
                         levels = paste0(100 * rates, "%", rate_suffix)),
      thr       = factor(paste0(round(100 * threshold), "%"),
                         levels = paste0(round(100 * sort(unique(threshold))), "%"))
    )
}

make_figure <- function(setting_, mechs, rates, keys, mech_display, layout, file,
                        width, height) {
  # Shorter row labels in the stacked layout
  rate_suffix <- if (layout == "stacked") " cens." else " Censoring"
  r <- prep_fig(long, setting_, mechs, rates, keys, mech_display, rate_suffix)
  p <- prep_fig(perf, setting_, mechs, rates, keys, mech_display, rate_suffix)
  lv <- levels(p$method)
  sc_colour <- scale_colour_manual(values = cols[lv], breaks = lv)
  sc_fill   <- scale_fill_manual(values = cols[lv], breaks = lv)
  sc_lines  <- list(scale_linetype_manual(values = ltys[lv], breaks = lv),
                    scale_shape_manual(values = shapes[lv], breaks = lv))
  fg    <- facet_grid(rate ~ mechanism)
  dodge <- position_dodge(width = 0.35)
  # Legend taken from Panel A
  lines_points <- list(
    geom_line(aes(linetype = method), position = dodge, linewidth = 0.5, show.legend = FALSE),
    geom_point(aes(shape = method), position = dodge, size = 1.7, show.legend = FALSE)
  )

  # Panel A: estimation error
  pA <- ggplot(r, aes(thr, err, fill = method, colour = method)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.7,
                 position = position_dodge(width = 0.8), alpha = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    fg + sc_colour + sc_fill + theme_fig +
    labs(x = "Decision threshold", y = "Estimated NB - True NB")
  if (show_mean) {   # optional marker for the mean error (the bias) in each box
    pA <- pA + stat_summary(aes(group = method), fun = mean, geom = "point", shape = 23,
                            size = 1.1, fill = "white", stroke = 0.6,
                            position = position_dodge(width = 0.8), show.legend = FALSE)
  }

  # Panel B: RMSE
  pB <- ggplot(p, aes(thr, rmse, colour = method, group = method)) +
    lines_points + fg + sc_colour + sc_lines +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    theme_fig + labs(x = "Decision threshold", y = "RMSE")

  # Panel C: coverage, with the Monte Carlo interval around 90%
  n_rep   <- median(p$n_sim)
  mc_half <- 1.96 * sqrt(nominal * (1 - nominal) / n_rep)
  pC <- ggplot(p, aes(thr, coverage, colour = method, group = method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = nominal - mc_half,
             ymax = nominal + mc_half, fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = nominal, linetype = "dashed", colour = "grey45") +
    lines_points + fg + sc_colour + sc_lines +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = function(x) paste0(100 * x, "%")) +
    theme_fig + labs(x = "Decision threshold", y = "Empirical coverage")

  fig <- if (layout == "two_row") {
    pA / (pB | pC) + plot_layout(heights = c(1, 1), guides = "collect")
  } else {
    pA / pB / pC + plot_layout(heights = c(1.15, 1, 1), guides = "collect")
  }
  fig <- fig + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")

  # width and height: intended size in inches if the figure is saved
  message(file)
  print(fig)
  invisible(fig)
}

main_mechs <- c("Non-informative", "Linear informative")
sens_mechs <- c("Non-informative", "Linear informative", "Nonlinear informative")
main_keys  <- c("ref", "ipcw_cox", "ipcw_miscox")
sens_keys  <- c("ref", "ipcw_cox", "ipcw_gam")
display_main <- c("Non-informative" = "Non-informative", "Linear informative" = "Informative")
display_sens <- c("Non-informative" = "Non-informative", "Linear informative" = "Linear informative",
                  "Nonlinear informative" = "Nonlinear informative")

# Main text (10% and 30% censoring)
make_figure("Standard",  main_mechs, main_rates, main_keys, display_main, "two_row",
            "Figure_1_standard",  width = 7, height = 8)
make_figure("Competing", main_mechs, main_rates, main_keys, display_main, "two_row",
            "Figure_2_competing", width = 7, height = 8)
# Supporting Information, Section S2 (all censoring rates)
make_figure("Standard",  main_mechs, all_rates, main_keys, display_main, "two_row",
            "Figure_S1_standard_all_rates",  width = 7, height = 10)
make_figure("Competing", main_mechs, all_rates, main_keys, display_main, "two_row",
            "Figure_S2_competing_all_rates", width = 7, height = 10)
# Supporting Information, Section S3 (sensitivity analysis)
make_figure("Standard",  sens_mechs, all_rates, sens_keys, display_sens, "stacked",
            "Figure_S3_standard_sensitivity",  width = 7, height = 9.5)
make_figure("Competing", sens_mechs, all_rates, sens_keys, display_sens, "stacked",
            "Figure_S4_competing_sensitivity", width = 7, height = 9.5)

# 4. Tables S2-S5 ----
print_perf_table <- function(setting_, blocks, name) {
  rows <- list()
  for (b in blocks) for (cr in all_rates) {
    d <- perf |>
      filter(setting == setting_, mechanism == b$mechanism, cens_rate == cr,
             method_key %in% b$keys) |>
      mutate(method_key = factor(method_key, levels = b$keys)) |>
      arrange(threshold, method_key)
    rows[[length(rows) + 1]] <- d |> mutate(method_key = as.character(method_key))
  }
  message(name)
  print(as.data.frame(bind_rows(rows)))
}

main_blocks <- list(
  list(mechanism = "Non-informative",    keys = main_keys),
  list(mechanism = "Linear informative", keys = main_keys))
# Sensitivity tables: nonlinear censoring (all approaches) and the GAM under the other mechanisms
sens_blocks <- list(
  list(mechanism = "Nonlinear informative", keys = sens_keys),
  list(mechanism = "Non-informative",       keys = "ipcw_gam"),
  list(mechanism = "Linear informative",    keys = "ipcw_gam"))

print_perf_table("Standard",  main_blocks, "Table_S2_standard_main")
print_perf_table("Competing", main_blocks, "Table_S3_competing_main")
print_perf_table("Standard",  sens_blocks, "Table_S4_standard_sensitivity")
print_perf_table("Competing", sens_blocks, "Table_S5_competing_sensitivity")

message("Done.")
