# Robust Decision Curve Analysis for Survival Risk Prediction Models in the Presence of Informative Censoring

This repository contains the R code, simulation datasets, and case study files necessary to reproduce the findings in the paper **"Robust Decision Curve Analysis for Survival Risk Prediction Models in the Presence of Informative Censoring."** 

This project introduces a robust Net Benefit (NB) estimator based on Conditional Inverse Probability of Censoring Weighting (CIPCW) to correct for informative censoring bias in Decision Curve Analysis (DCA) for both standard time-to-event and competing risks outcomes.

---

## Repository Structure

The repository is organized to separate computationally heavy simulation scripts from result aggregation, visualization, and clinical application.

```text
papercode-RobustDCA/
│
├── README.md                              
├── 03_results_summary_and_plots.Rmd       <-- Aggregates CSV results and generates all manuscript tables and figures.
│
├── R/                                     
│   ├── 01_sim_standard_survival.R         <-- Core Monte Carlo simulation for standard survival outcomes.
│   └── 02_sim_competing_risks.R           <-- Core Monte Carlo simulation for competing risks outcomes.
│
├── results/                               
│   ├── combined_results_standard.csv      <-- Pre-computed results from the standard survival simulation.
│   └── combined_results_competing.csv     <-- Pre-computed results from the competing risks simulation.
│
├── figures/                               <-- Destination folder for high-resolution output plots.
│
└── case_study/                            
    └── case_study_GBSG.Rmd                <-- Step-by-step application of CIPCW-DCA on the GBSG breast cancer dataset.