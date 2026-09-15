# Repository Guidelines

## Current Mission

The manuscript draft is in its final stages. Our current mission is to tidy the code and organize files for release in a public repository. Prioritize readability, portable paths, documented dependencies, and reproducible run instructions. Preserve the methods, parameters, and results supporting the manuscript, and flag changes that could affect them. Identify manuscript-related scripts and outputs before removing or consolidating files; ask for clarification when their role is unclear.

## Project Structure & Module Organization

This repository analyzes sleep-spindle coupling in adult and child cohorts, primarily with R and XGBoost.

- Root R scripts handle tuning (`classification_model.R`), channel comparisons (`CvsF.R`), statistics, and visualization.
- `script/` contains pooled leave-one-subject-out (LOSO) evaluation, cohort/channel comparisons, and SHAP analysis.
- `ebm/` contains Python Explainable Boosting Machine experiments, plotting helpers, and related R/HPC scripts.
- `results/` stores generated tables, model summaries, and figures; `results/final/` holds selected final outputs.
- `old/` contains archived scripts, data, and results. Place new work alongside the active scripts.

## Build, Test, and Development Commands

There is no package build or dependency lockfile. Install packages declared by the relevant script's `library()` calls or Python imports before running it. Run these commands from the repository root:

- `Rscript -e 'parse(file="classification_model.R"); parse(file="script/LOSO.R")'` checks R syntax without executing analyses.
- `Rscript -e 'library(readr); source("classification_model.R")'` runs the 50-trial adult XGBoost tuning search. Loading `readr` supplies the script's unqualified `read_csv()` call.
- `Rscript script/LOSO.R` runs pooled subject-level evaluation and exports metrics and a figure; it requires both cohort CSVs and `tuning_results_adults.csv`.

Scripts reference `~/spindles/` and often write into the current working directory. Check input/output paths before execution. Some statistical scripts require data objects already loaded in R.

## Coding Style & Naming Conventions

Use two-space indentation in R and four spaces in Python. Prefer `<-` for R assignments, `snake_case` for new functions and variables, and uppercase constants such as `ADULT_SAMP_FREQ`. Preserve dataset column names such as `peakLoc` and `coupling_label`. Follow existing section comments and descriptive cohort-based output names. No formatter or linter is configured.

## Testing Guidelines

No automated test framework, naming convention, or coverage target exists. Parse changed R scripts, then exercise affected logic on a small reproducible subset. For LOSO changes, verify patient separation between training and validation, label encoding, and feature exclusions. Compare ROC-AUC, PR-AUC, and output schemas with prior results; record seeds and package versions.

## Commit & Pull Request Guidelines

History uses short descriptive subjects without a formal prefix convention. Write concise action-oriented messages, such as `Fix pooled LOSO feature selection`. PRs should explain the analysis change, affected cohorts/scripts, validation commands, and metric differences. Link relevant issues and include before/after figures for visualization changes. Respect `.gitignore` exclusions for CSV/RDS data, session files, and cluster logs; avoid unrelated generated artifacts.
