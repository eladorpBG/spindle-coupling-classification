# Compare central and frontal channel XGBoost models within one cohort.

library(ggpubr)
library(readr)
library(rstatix)
library(tableone)
library(survey)
library(reshape2)
library(caret)
library(randomForest)
library(zoo)
library(pROC)
library(PRROC)
library(ROCR)
library(xgboost)
library(fastDummies)
library(lmtest)
library(tidymodels)
library(dplyr)
library(data.table)
library(igraph)

## RStudio settings ----

# Select "adults" or "children".
COHORT <- "adults"
DATA_DIR <- "."
TUNING_DIR <- "."
OUTPUT_DIR <- "results"
APPLY_FILTER <- FALSE
ADULT_SAMP_FREQ <- 200
CHILD_SAMP_FREQ <- 256
set.seed(123)

COHORT <- match.arg(COHORT, c("adults", "children"))
INPUT_FILE <- file.path(path.expand(DATA_DIR), paste0(COHORT, "_data_2025.csv"))
TUNING_FILE <- file.path(path.expand(TUNING_DIR), paste0("tuning_results_", COHORT, ".csv"))
METRICS_FILE <- file.path(path.expand(OUTPUT_DIR), paste0(COHORT, "_cvsf_subject_metrics.rds"))
FIGURE_FILE <- file.path(path.expand(OUTPUT_DIR), paste0("combined_", COHORT, "_cvsf_heatmaps.png"))

if (!file.exists(INPUT_FILE)) stop("Input file not found: ", INPUT_FILE)
if (!file.exists(TUNING_FILE)) stop("Tuning file not found: ", TUNING_FILE)
if (!dir.exists(path.expand(OUTPUT_DIR))) stop("Output directory not found: ", OUTPUT_DIR)

## Evaluation helpers ----

eval_model_roc <- function(model, test_mat, test_labels){
  xgb_preds <- predict(model, test_mat, reshape = TRUE)
  print_confusion_mat(xgb_preds, test_labels)
  model_auc <- print_AUC_plot_ROC(xgb_preds, test_labels)
  return(model_auc)
}

eval_model_pr <- function(model, test_mat, test_labels){
  xgb_preds <- predict(model, test_mat, reshape = TRUE)

  # Factor labels are converted to 0/1 weights for PR AUC.
  pr_result <- pr.curve(
    scores.class0 = xgb_preds, 
    weights.class0 = as.numeric(test_labels) - 1, 
    curve = TRUE
  )

  return(pr_result$auc.integral)
}
print_confusion_mat <- function(xgb_preds, test_labels){
  df_xgb_preds <- as.data.frame(xgb_preds)
  colnames(df_xgb_preds) <- "1"
  df_xgb_preds$"0" <- 1 - df_xgb_preds$"1"
  df_xgb_preds$PredictedClass <- apply(df_xgb_preds, 1, function(y) colnames(df_xgb_preds)[which.max(y)])
  confusion_matrix <- confusionMatrix(as.factor(df_xgb_preds$PredictedClass), test_labels)
  return(confusion_matrix)
}

print_AUC_plot_ROC <- function(xgb_preds, test_labels){
  xgb_prediction <- prediction(xgb_preds, test_labels);
  xgb_auc <- performance(xgb_prediction, measure = "auc")@y.values[[1]]
  print(xgb_auc)
  roc_curve <- performance(xgb_prediction, measure="tpr", x.measure="fpr")
  png()
  plot(roc_curve)
  dev.off()
  return(xgb_auc)
}

get_ci_from_ttest <- function(x) {
  if (length(na.omit(x)) < 2) {
    return(c(mean = mean(x), lower_ci = NA, upper_ci = NA))
  }

  test <- t.test(x)

  return(c(
    mean = unname(test$estimate),
    lower_ci = test$conf.int[1],
    upper_ci = test$conf.int[2]
  ))
}

## Load and preprocess data ----

data4elec <- read.csv(INPUT_FILE)
samp_freq <- c(adults = ADULT_SAMP_FREQ, children = CHILD_SAMP_FREQ)[[COHORT]]

data4elec <- data4elec %>% 
  arrange(patient_id, startSample) %>% 
  mutate(row_id = row_number())

dt <- as.data.table(data4elec)
setkey(dt, patient_id, startSample, endSample)

overlaps <- foverlaps(dt, dt, type = "any", which = FALSE)

valid_links <- overlaps[row_id < i.row_id, .(row_id, i.row_id)]

# Connected overlaps share a cluster for optional event filtering.
g <- graph_from_data_frame(valid_links, directed = FALSE, vertices = data4elec$row_id)
comps <- components(g)

if (APPLY_FILTER) {
  # When enabled, retain one event per patient and overlap cluster.
  df_filtered <- data4elec %>%
    mutate(cluster_id = comps$membership) %>%
    group_by(patient_id, cluster_id) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-cluster_id, -row_id)
} else {
  df_filtered <- data4elec
}
cat("Original rows:", nrow(data4elec), "\n")
cat("Filtered rows:", nrow(df_filtered), "\n")

data4elec <- df_filtered
data4elec$coupling_label <- as.factor(data4elec$coupling_label)
data4elec <- subset(data4elec, !grepl("O", channel))
data4elec$channel <- ifelse(grepl("F", data4elec$channel), "F", data4elec$channel)
data4elec$channel <- ifelse(grepl("C", data4elec$channel), "C", data4elec$channel)
if (!all(c("F", "C") %in% data4elec$channel) ||
    any(!(data4elec$channel %in% c("F", "C")))) {
  stop("Expected frontal and central channels after preprocessing.")
}
data4elec <- tibble::rowid_to_column(data4elec, "ID")
data4elec$refrPeriod <- na.aggregate(data4elec$refrPeriod, by = data4elec$patient_id)
data4elec <- subset(data4elec, select= -numCycles)
data4elec$peakLoc <- data4elec$peakLoc / samp_freq

## Grouped cross-validation ----

k <- 10
folds <- rsample::group_vfold_cv(data4elec, group = patient_id, v = k)

final_auc_mat <- matrix(0, nrow = 3, ncol = 3)

tune_results <- read_csv(TUNING_FILE)

max_score_idx = which.max(tune_results$`as.numeric(score[1])`)
best_params <- tune_results[max_score_idx,-(ncol(tune_results)-1)]
best_params_list <- as.list(best_params[,-ncol(best_params)])

all_indices <- seq_len(nrow(data4elec))
all_oof_preds <- data.frame()

for (i in 1:k) {
  cat("Processing fold", i, "\n")

  current_split <- folds$splits[[i]]

  train_indices <- current_split$in_id
  test_indices <- setdiff(all_indices, train_indices)

  data4model <- data4elec[ , !names(data4elec) %in% c("patient_id", "spindle_idx", "ID", "detSample", "startSample", "endSample")]
  train4elec <- data4model[train_indices, ]
  test4elec <- data4model[test_indices, ]

  F_train <- train4elec[train4elec$channel == "F", ]
  C_train <-train4elec[train4elec$channel == "C", ]
  min_samples <- min(nrow(train4elec), nrow(F_train), nrow(C_train))
  all_inds <- createDataPartition(train4elec$coupling_label, p = min_samples /  nrow(train4elec), list = FALSE, times = 1)
  all_train <- train4elec[all_inds,]
  F_inds <- createDataPartition(F_train$coupling_label, p = min_samples /  nrow(F_train), list = FALSE, times = 1)
  F_train <- F_train[F_inds,]
  C_inds <- createDataPartition(C_train$coupling_label, p = min_samples /  nrow(C_train), list = FALSE, times = 1)
  C_train <- C_train[C_inds,]
  F_test <- test4elec[test4elec$channel == "F", ]
  C_test <- test4elec[test4elec$channel == "C", ]

  test_patient_ids_all <- data4elec$patient_id[test_indices]
  test_patient_ids_F <- test_patient_ids_all[test4elec$channel == "F"]
  test_patient_ids_C <- test_patient_ids_all[test4elec$channel == "C"]

  all_train_mat <- as.matrix(all_train[,sapply(all_train, is.numeric)])
  F_train_mat <- as.matrix(F_train[,sapply(F_train, is.numeric)])
  C_train_mat <- as.matrix(C_train[,sapply(C_train, is.numeric)])

  all_test_mat <- as.matrix(test4elec[,sapply(test4elec, is.numeric)])
  F_test_mat <- as.matrix(F_test[,sapply(F_test, is.numeric)])
  C_test_mat <- as.matrix(C_test[,sapply(C_test, is.numeric)])

  xgb_all_train <-  xgb.DMatrix(data = all_train_mat, label = as.numeric(all_train$coupling_label) -1)
  xgb_F_train <- xgb.DMatrix(data = F_train_mat, label = as.numeric(F_train$coupling_label) -1)
  xgb_C_train <- xgb.DMatrix(data = C_train_mat, label = as.numeric(C_train$coupling_label) -1)
  xgb_all_test <-  xgb.DMatrix(data = all_test_mat, label = as.numeric(test4elec$coupling_label) -1)
  xgb_F_test <- xgb.DMatrix(data = F_test_mat, label = as.numeric(F_test$coupling_label) -1)
  xgb_C_test <- xgb.DMatrix(data = C_test_mat, label = as.numeric(C_test$coupling_label) -1)

  imbalance_weight_all <-  sum(all_train$coupling_label == 0) / sum(all_train$coupling_label == 1)
  imbalance_weight_F <- sum(F_train$coupling_label == 0) / sum(F_train$coupling_label == 1)
  imbalance_weight_C <- sum(C_train$coupling_label == 0) / sum(C_train$coupling_label == 1)

  xgb_model_all <- xgb.train(
    params = best_params_list,
    data = xgb_all_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_all
  )

  xgb_model_F <- xgb.train(
    params = best_params_list,
    data = xgb_F_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_F
  )
  xgb_model_C <- xgb.train(
    params = best_params_list,
    data = xgb_C_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_C
  )

  model_names <- c("All", "F", "C")
  models <- list(xgb_model_all, xgb_model_F, xgb_model_C)

  datasets <- list(
    list(name = "All", data = xgb_all_test, labels = test4elec$coupling_label, pids = test_patient_ids_all),
    list(name = "F", data = xgb_F_test, labels = F_test$coupling_label, pids = test_patient_ids_F),
    list(name = "C", data = xgb_C_test, labels = C_test$coupling_label, pids = test_patient_ids_C)
  )

  for (l in seq_along(models)) {
    for (j in seq_along(datasets)) {
      preds <- predict(models[[l]], datasets[[j]]$data, reshape = TRUE)

      temp_df <- data.frame(
        patient_id = datasets[[j]]$pids,
        label = as.numeric(datasets[[j]]$labels) - 1,
        pred = preds,
        train_model = model_names[l],
        test_set = datasets[[j]]$name
      )

      all_oof_preds <- bind_rows(all_oof_preds, temp_df)
    }
  }
}

## Subject-level metrics ----

# AUC is undefined for patients with only one outcome class.
subject_metrics <- all_oof_preds %>%
  group_by(train_model, test_set, patient_id) %>%
  summarize(
    n_pos = sum(label == 1),
    n_neg = sum(label == 0),
    ROCAUC = if (n_pos > 0 && n_neg > 0) {
      as.numeric(pROC::auc(label, pred, quiet = TRUE))
    } else { NA_real_ },
    PRAUC = if (n_pos > 0 && n_neg > 0) {
      PRROC::pr.curve(scores.class0 = pred, weights.class0 = label)$auc.integral
    } else { NA_real_ },
    .groups = "drop"
  )

skipped_patients <- subject_metrics %>% filter(is.na(ROCAUC))

if (nrow(skipped_patients) > 0) {
  cat("\n======================================================\n")
  cat(" WARNING: SKIPPED PATIENTS (Missing one or both classes)\n")
  cat("======================================================\n")

  for (i in 1:nrow(skipped_patients)) {
    row <- skipped_patients[i, ]
    cat(sprintf("Skipped Patient ID: %-10s | Train: %-5s | Test: %-5s | Positives: %-3d | Negatives: %-3d\n", 
                row$patient_id, row$train_model, row$test_set, row$n_pos, row$n_neg))
  }
  cat("======================================================\n\n")
} else {
  cat("\nAll patients successfully evaluated (no missing classes).\n\n")
}

saveRDS(subject_metrics, file = METRICS_FILE)

## Summarize model comparisons ----

process_auc_data_cvsf <- function(metrics_df, metric_col) {
  summary_df <- metrics_df %>%
    filter(!is.na(.data[[metric_col]])) %>%
    group_by(train_model, test_set) %>%
    summarize(
      mean_val = get_ci_from_ttest(.data[[metric_col]])["mean"],
      Lower_CI = get_ci_from_ttest(.data[[metric_col]])["lower_ci"],
      Upper_CI = get_ci_from_ttest(.data[[metric_col]])["upper_ci"],
      .groups = "drop"
    ) %>%
    rename(Model = train_model, Test = test_set, !!metric_col := mean_val)

  desired_order <- c("Central", "Frontal", "Pooled")

  summary_df <- summary_df %>%
    mutate(
      Model = case_when(
        Model == "C"   ~ "Central",
        Model == "F"   ~ "Frontal",
        Model == "All" ~ "Pooled",
        TRUE ~ as.character(Model)
      ),
      Test = case_when(
        Test == "C"   ~ "Central",
        Test == "F"   ~ "Frontal",
        Test == "All" ~ "Pooled",
        TRUE ~ as.character(Test)
      )
    ) %>%
    mutate(
      Model = factor(Model, levels = desired_order), 
      Test  = factor(Test, levels = desired_order)
    )

  return(summary_df)
}

subject_metrics <- readRDS(METRICS_FILE)

auc_data_roc <- process_auc_data_cvsf(subject_metrics, "ROCAUC")
auc_data_pr  <- process_auc_data_cvsf(subject_metrics, "PRAUC")

## Plot model comparisons ----

roc_plot <- ggplot(auc_data_roc, aes(x = Test, y = Model, fill = ROCAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 

  geom_text(aes(label = sprintf("%.3f", ROCAUC)), 
            color = "white", size = 7, fontface = "bold", nudge_y = 0.15) + 
  geom_text(aes(label = sprintf("[%.3f - %.3f]", Lower_CI, Upper_CI)), 
            color = "white", size = 4.5, nudge_y = -0.15) + 

  scale_fill_gradient(low = "#0072B2", high = "#D55E00", name = "ROC AUC", limits = c(0.5, 0.8)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "Test Dataset", y = "Model Trained On") +
  theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 26),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10) 
  )

pr_plot <- ggplot(auc_data_pr, aes(x = Test, y = Model, fill = PRAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 

  geom_text(aes(label = sprintf("%.3f", PRAUC)), 
            color = "white", size = 7, fontface = "bold", nudge_y = 0.15) + 
  geom_text(aes(label = sprintf("[%.3f - %.3f]", Lower_CI, Upper_CI)), 
            color = "white", size = 4.5, nudge_y = -0.15) + 

  scale_fill_gradient(low = "#0072B2", high = "#D55E00", name = "PR AUC", limits = c(0.1, 0.5)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "Test Dataset", y = "") + 
  theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 26),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10) 
  )

## Export plots ----

combined_plot <- roc_plot + pr_plot + 
  plot_annotation(tag_levels = list(c('C', 'D'))) &
  theme(plot.tag = element_text(size = 28, face = "bold"))

ggsave(FIGURE_FILE, plot = combined_plot, width = 18, height = 8, dpi = 300)
