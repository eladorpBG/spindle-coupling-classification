# Random-search tuning for the manuscript's XGBoost models.

library(xgboost)
library(splitstackshape)
library(zoo)

## RStudio settings ------------------------------------------------------------

COHORT <- "adults"               # "adults", "children", or "both"
DATA_DIR <- "."                  # Or an absolute path, e.g. "~/spindles"
OUTPUT_DIR <- file.path(DATA_DIR, "results", "tuning")
N_TRIALS <- 50
SEED <- 123

## Preprocessing ---------------------------------------------------------------

prepare_tuning_data <- function(data, cohort = c("adults", "children")) {
  cohort <- match.arg(cohort)
  required <- c("patient_id", "channel", "coupling_label", "refrPeriod",
                "peakLoc", "numCycles")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0) {
    stop("Missing input columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyNA(data[c("patient_id", "channel", "coupling_label")])) {
    stop("patient_id, channel, and coupling_label must not contain missing values.")
  }
  if (!all(data$coupling_label %in% c(0, 1))) {
    stop("coupling_label must contain only 0 (uncoupled) and 1 (coupled).")
  }
  if (!is.numeric(data$refrPeriod) || !is.numeric(data$peakLoc)) {
    stop("refrPeriod and peakLoc must be numeric.")
  }

  # Remove occipital channels and normalize references.
  data$coupling_label <- factor(data$coupling_label, levels = c(0, 1))
  data <- data[!grepl("O", data$channel), , drop = FALSE]
  if (nrow(data) == 0) stop("No spindles remain after channel filtering.")
  data$channel <- sub("-M[12]$", "", data$channel)

  data$refrPeriod <- zoo::na.aggregate(data$refrPeriod, by = data$patient_id)
  data <- subset(data, select = -numCycles)
  sampling_frequency <- c(adults = 200, children = 256)[[cohort]]
  data$peakLoc <- data$peakLoc / sampling_frequency
  data
}

make_tuning_matrix <- function(data) {
  # Select predictors by name, independent of input column order.
  metadata <- c("ID", "row_id", "patient_id", "spindle_idx", "detSample",
                "startSample", "endSample", "Cohort", "channel", "coupling_label")
  predictors <- setdiff(names(data)[vapply(data, is.numeric, logical(1))], metadata)
  if (length(predictors) == 0) stop("No numeric predictors found.")
  as.matrix(data[, predictors, drop = FALSE])
}

## Search and output -----------------------------------------------------------

sample_tuning_params <- function() {
  list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = runif(1, 0.01, 0.3),
    max_depth = sample(3:10, 1),
    min_child_weight = sample(1:5, 1),
    subsample = runif(1, 0.6, 1),
    colsample_bytree = runif(1, 0.6, 1),
    nrounds = sample(c(500, 1000, 2500), 1)
  )
}

write_tuning_results <- function(results, output_file) {
  # Manuscript scripts expect these two names as the final columns.
  export <- results
  names(export)[names(export) == "cv_auc"] <- "as.numeric(score[1])"
  names(export)[names(export) == "best_iteration"] <- "as.numeric(score[2])"
  write.csv(export, output_file, row.names = FALSE)
}

tune_xgboost_params <- function(cohort = c("adults", "children", "both"),
                               data_dir = ".", output_dir = ".",
                               n_trials = 50, seed = 123) {
  cohort <- match.arg(cohort)
  if (length(n_trials) != 1 || !is.finite(n_trials) ||
      n_trials < 1 || n_trials != floor(n_trials)) {
    stop("n_trials must be a positive integer.")
  }

  if (cohort == "both") {
    cohorts <- c("adults", "children")
    input_files <- file.path(path.expand(data_dir),
                             paste0(cohorts, "_data_2025.csv"))
    missing_files <- input_files[!file.exists(input_files)]
    if (length(missing_files) > 0) {
      stop("Input files not found: ", paste(missing_files, collapse = ", "),
           ". Check DATA_DIR in the RStudio settings.", call. = FALSE)
    }

    results <- lapply(cohorts, function(current_cohort) {
      tune_xgboost_params(current_cohort, data_dir, output_dir, n_trials, seed)
    })
    names(results) <- cohorts
    return(invisible(results))
  }

  input_file <- file.path(path.expand(data_dir), paste0(cohort, "_data_2025.csv"))
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file,
         ". Check data_dir (DATA_DIR in the RStudio settings).", call. = FALSE)
  }

  # Validate the output destination before starting the tuning search.
  output_dir <- path.expand(output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }
  if (file.access(output_dir, mode = 2) != 0) {
    stop("Output directory is not writable: ", output_dir, call. = FALSE)
  }
  output_file <- file.path(output_dir, paste0("tuning_results_", cohort, ".csv"))
  if (dir.exists(output_file) ||
      (file.exists(output_file) && file.access(output_file, mode = 2) != 0)) {
    stop("Output file is not writable: ", output_file, call. = FALSE)
  }

  message("Tuning ", cohort, " using ", normalizePath(input_file))
  data <- prepare_tuning_data(read.csv(input_file), cohort)
  set.seed(seed)

  # A 66% sample within each label/patient stratum, followed by row-level five-fold CV.
  split <- splitstackshape::stratified(
    data, c("coupling_label", "patient_id"), size = 0.66, bothSets = TRUE
  )
  train <- as.data.frame(split[[1]])
  labels <- as.numeric(train$coupling_label) - 1
  n_folds <- 5
  if (any(table(factor(labels, levels = c(0, 1))) < n_folds)) {
    stop("The tuning sample needs at least five observations of each class.")
  }
  dtrain <- xgboost::xgb.DMatrix(make_tuning_matrix(train), label = labels)
  imbalance_weight <- sum(labels == 0) / sum(labels == 1)

  trials <- vector("list", n_trials)
  for (i in seq_len(n_trials)) {
    params <- sample_tuning_params()
    # nrounds controls CV; booster and class weight are model parameters.
    cv_params <- params[names(params) != "nrounds"]
    cv_params$booster <- "dart"
    cv_params$scale_pos_weight <- imbalance_weight
    cv <- xgboost::xgb.cv(
      params = cv_params,
      data = dtrain,
      nrounds = params$nrounds,
      nfold = n_folds,
      stratified = TRUE,
      early_stopping_rounds = 10,
      maximize = TRUE,
      verbose = FALSE
    )

    auc <- cv$evaluation_log$test_auc_mean
    if (length(auc) == 0 || any(!is.finite(auc))) {
      stop("Trial ", i, " returned missing or non-finite validation AUC.")
    }
    best_iteration <- which.max(auc)
    trials[[i]] <- cbind(as.data.frame(params),
                         cv_auc = auc[best_iteration],
                         best_iteration = best_iteration)
    message(sprintf("Trial %d/%d: AUC = %.4f, rounds = %d",
                    i, n_trials, auc[best_iteration], best_iteration))
  }

  results <- do.call(rbind, trials)
  write_tuning_results(results, output_file)
  message("Saved tuning results to ", output_file)
  message("Best trial:")
  print(results[which.max(results$cv_auc), , drop = FALSE])
  invisible(results)
}

## Run -------------------------------------------------------------------------

tuning_results <- tune_xgboost_params(
  cohort = COHORT,
  data_dir = DATA_DIR,
  output_dir = OUTPUT_DIR,
  n_trials = N_TRIALS,
  seed = SEED
)
