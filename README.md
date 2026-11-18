# paper-RobustDCA

Code to reproduce the simulation results, figures, tables, and case study results from the paper
"Robust Decision Curve Analysis for Survival Risk Prediction Models in the Presence of Informative Censoring"
by Amirhosein Farahmand, Mauricio Lopez-Mendez, Mohsen Sadatsafavi, Abdollah Safari, and Tae Yoon Lee.

Preprint: [arXiv:2611.21953](http://arxiv.org/abs/2611.21953)

## Organization

All scripts are run from the root of this repository.

### `simulation-code/`

Code to run the Monte Carlo simulation study. Each script runs 1,000 replications with 500 bootstrap
resamples per replication, for three censoring mechanisms (non-informative, linear informative,
nonlinear informative) and three censoring rates (10%, 20%, 30%). The simulation results used in
the paper are provided in `simulation-results/`.

- `01_sim_standard_survival.R`: standard survival setting (Figures 1 and S1; Tables S2 and S4).
- `02_sim_competing_risks.R`: competing-risks setting (Figures 2 and S2; Tables S3 and S5).

Random seeds are set inside the scripts (`set.seed(123)` for the training data, `set.seed(456)` for the
reference population, and `future.seed = 123L` for the replications).

**Computing time.** The simulations were run in parallel on a machine with 80 CPU cores
(79 parallel workers); each simulation script took approximately 2.4 hours.

### `simulation-results/`

Output of the simulation scripts, with one row per replication, scenario, censoring rate, and decision
threshold (`combined_results_standard.csv` and `combined_results_competing.csv`).

### `figures-tables-code/`

- `simulation_figures_tables.R`: reads the files in `simulation-results/` and produces the simulation
  figures (Figures 1, 2, S1, and S2) and tables (Tables S2–S5). Figures 1 and 2 show all three
  censoring rates, and the main performance measures are relative bias and relative RMSE (as a
  percentage of the true NB). Tables S2–S5 report both relative and absolute measures, saved as CSV
  and as LaTeX table rows. Runs in under a minute.
  Performance measures and Monte Carlo standard errors follow Morris, White and Crowther
  (*Statistics in Medicine*, 2019).
- `positivity_table_S7.R`: positivity diagnostic (Table S7). Uses the data-generating functions of the
  standard survival simulation (no replications are run). Runs in about a minute.

### `case-study-code/`

- `gbsg_case_study.R`: GBSG case study. Computes predicted 5-year risks from the Royston–Altman model,
  assesses informative censoring, fits the censoring model, estimates net benefit with the KM-based and
  CIPCW-based approaches, and produces Figures 3 and 4 (saved to `figures/`) and Table S6.

The GBSG and Rotterdam data are included in the R package `survival` (`survival::gbsg`,
`survival::rotterdam`), so no download is needed.

### `figures/` and `tables/`

Figures in the paper (PDF and PNG) and tables in the Supporting Information (CSV). Table S1 describes the
simulation design and is written by hand; the scripts above print all other figures and tables.

## Software

The proposed CIPCW-based estimator of net benefit is implemented in our fork of the `dcurves` R package,
including a vignette illustrating its use:
<https://github.com/amirfarahmand0/dcurves/tree/feature/cipcw>.

Analyses were run in R 4.6.0. Required packages:

- Simulation: `data.table`, `survival`, `survivalROC`, `mvtnorm`, `mgcv`, `future.apply`,
  `progressr`, `dplyr`, `tidyr`
- Figures and tables: `dplyr`, `tidyr`, `ggplot2`, `patchwork` (optional: `ragg`)
- Case study: `survival`, `survivalROC`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`
