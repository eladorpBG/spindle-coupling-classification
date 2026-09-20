# Spindle Coupling Classification

R analyses of sleep-spindle coupling in adult and child cohorts. The project uses spindle features to predict coupling status with XGBoost, compares performance across cohorts and electrode locations, and summarizes feature contributions with SHAP.

## Repository structure

```text
script/                  Tuning, manuscript analysis, and visualization scripts
results/final/           Selected manuscript tables, metrics, and figures
```

### Script guide

| Script                                                           | Purpose                                                                           | Main outputs                                                       |
| ---------------------------------------------------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| [tune_xgboost_params.R](script/tune_xgboost_params.R)             | Adult, child, or both-cohort XGBoost hyperparameter search                        | `tuning_results_adults.csv`, `tuning_results_children.csv`     |
| [LOSO.R](script/LOSO.R)                                           | Pooled leave-one-subject-out evaluation                                           | `pooled_loso_metrics.csv`, summary statistics, pooled AUC figure |
| [ChildvsAdult.R](script/ChildvsAdult.R)                           | Compare pooled, child, and adult models using subject-grouped ten-fold validation | Subject metrics, SHAP values, cohort-comparison heatmaps           |
| [CvsF.R](script/CvsF.R)                                           | Central/frontal model comparisons using subject-grouped ten-fold validation       | Subject metrics and central/frontal heatmaps                       |
| [3vs4.R](script/3vs4.R)                                           | Left/right electrode model comparisons using subject-grouped ten-fold validation  | Subject metrics and lateral heatmaps                               |
| [shap.R](script/shap.R)                                           | Plot cohort-specific SHAP summaries from saved values                             | `combined_shap_summary.png`                                      |
| [statistical_analysis.R](script/statistical_analysis.R)           | Descriptive statistics, standardized mean differences, and spindle density        | Summary table, SMD and density figures                             |
| [statistical_visualization.R](script/statistical_visualization.R) | Plot coupled/uncoupled feature means for both cohorts                             | `combined_mean_values_plot.png`                                  |
| [corr_plots.R](script/corr_plots.R)                               | Feature-correlation matrices for both cohorts                                     | Individual and combined correlation figures                        |

## Requirements

The scripts are intended to be run in **RStudio Server**. Open a script and click **Source** after setting the working directory and any settings at the top of that script. The reference analysis environment used **R 4.4.3 (2025-02-28)** on **x86_64-pc-linux-gnu**. There is no dependency lockfile; the attached package versions from that R session are recorded below.

Install missing packages in the RStudio Console. The following covers packages referenced by the current manuscript and tuning scripts:

```r
install.packages(c(
  "caret", "data.table", "dplyr", "fastDummies", "ggcorrplot",
  "ggplot2", "ggpubr", "igraph", "lmtest", "patchwork", "pROC",
  "PRROC", "randomForest", "readr", "reshape2", "ROCR", "rsample",
  "rstatix", "SHAPforxgboost", "splitstackshape", "stringr",
  "survey", "tableone", "tibble", "tidymodels", "tidyr", "tidyverse",
  "xgboost", "zoo"
))
```

### Reference package versions

These versions were reported as *attached packages* in the reference R session:

| Package        | Version | Package      | Version |
| -------------- | ------- | ------------ | ------- |
| yardstick      | 1.3.2   | workflowsets | 1.1.1   |
| workflows      | 1.3.0   | tune         | 2.0.0   |
| tailor         | 0.1.0   | rsample      | 1.3.1   |
| recipes        | 1.3.1   | parsnip      | 1.3.3   |
| modeldata      | 1.5.1   | infer        | 1.0.9   |
| dials          | 1.4.2   | scales       | 1.4.0   |
| broom          | 1.0.10  | tidymodels   | 1.4.1   |
| lmtest         | 0.9-40  | fastDummies  | 1.7.4   |
| randomForest   | 4.7-1.2 | caret        | 6.0-94  |
| lattice        | 0.22-5  | reshape2     | 1.4.4   |
| survey         | 4.4-2   | survival     | 3.8-3   |
| Matrix         | 1.7-2   | tableone     | 0.13.2  |
| rstatix        | 0.7.2   | ggpubr       | 0.6.0   |
| ggcorrplot     | 0.1.4.1 | igraph       | 2.2.1   |
| patchwork      | 1.3.2   | zoo          | 1.8-12  |
| SHAPforxgboost | 0.1.3   | PRROC        | 1.3.1   |
| pROC           | 1.18.5  | data.table   | 1.16.2  |
| xgboost        | 1.7.8.1 | lubridate    | 1.9.3   |
| forcats        | 1.0.0   | stringr      | 1.5.1   |
| dplyr          | 1.1.4   | purrr        | 1.1.0   |
| readr          | 2.1.5   | tidyr        | 1.3.1   |
| tibble         | 3.2.1   | ggplot2      | 4.0.0   |
| tidyverse      | 2.0.0   | ROCR         | 1.0-11  |

## Required inputs

Provide the following local files:

| File                       | Used by                   |
| -------------------------- | ------------------------- |
| `adults_data_2025.csv`   | Adult and pooled analyses |
| `children_data_2025.csv` | Child and pooled analyses |

Each cohort CSV contains one row per spindle, with these columns:

```text
patient_id, spindle_idx, detSample, startSample, endSample,
peakAmp, peakLoc, peakFreq, sigmaPower, energy, duration, numCycles,
numBumps, symmetry, freqGradient, fano, raisingSlope, droppingSlope,
refrPeriod, corrCoef, channel, coupling_label
```

`coupling_label` encodes uncoupled (`0`) and coupled (`1`) spindles. Model scripts use `patient_id` to group subjects and exclude identifying/event-index columns from predictors. The main model analyses remove occipital channels, impute missing `refrPeriod` values by patient, drop `numCycles`, and convert `peakLoc` from samples to seconds using 200 Hz for adults and 256 Hz for children.

## Running the analyses

Use RStudio Server with the repository root as the working directory. In the Console, check `getwd()` and, if needed, run `setwd("~/spindles")` with your actual checkout path. Open each script in the editor and click **Source**.

### Model evaluation and SHAP

With both cohort CSVs and tuning results available, source `script/LOSO.R` for pooled evaluation. To generate the cohort comparison and SHAP summary, then source `script/ChildvsAdult.R` followed by `script/shap.R` in the same session. `shap.R` reads `ChildvsAdults_shap_values.rds` produced by the cohort comparison.

For electrode comparisons, open `script/CvsF.R` or `script/3vs4.R`, set `COHORT` to `"adults"` or `"children"`, and click Source. Results go to `results/`; set `TUNING_DIR <- "results/tuning"` if using new tuning outputs. Run each cohort separately; both scripts default to adults and leave overlap filtering off.

### Descriptive statistics and feature means

Load the cohort data in the Console. Then source `script/statistical_analysis.R`, followed by `script/statistical_visualization.R` in the same RStudio session. The visualization script consumes `final_table` and `vars_continuous` created by the statistics script. It replaces that script's initial mean-values figure with the manuscript-style layout.

### Correlations

Before clicking Source on `script/corr_plots.R`, define `fixed_variable_order` as a character vector of the manuscript's displayed feature names in the intended order. The script uses this object but currently does not define it. The original ordering must be supplied to reproduce the saved figure.

### Hyperparameter tuning

Open `script/tune_xgboost_params.R`, set `COHORT` to `"adults"`, `"children"`, or `"both"`, check `DATA_DIR` and `OUTPUT_DIR`, and click Source.
