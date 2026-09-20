# Compare child, adult, and pooled XGBoost models across cohorts.

library(ggpubr)
library(rstatix)
library(tableone)
library(survey)
library(reshape2)
library(caret)
library(zoo)
library(pROC)
library(PRROC)
library(ROCR)
library(dplyr)
library(xgboost)
library(fastDummies)
library(lmtest)
library(readr)
library(SHAPforxgboost)

set.seed(123)
ADULT_SAMP_FREQ <- 200
CHILD_SAMP_FREQ <- 256

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

## Load and preprocess cohorts ----

children_data <- read.csv("~/spindles/children_data_2025.csv")
children_data$coupling_label <- as.factor(children_data$coupling_label)
children_data <- subset(children_data, !grepl("O", channel))
children_data$channel <- sub("-M[12]$", "", children_data$channel)
children_data <- tibble::rowid_to_column(children_data, "ID")
children_data$refrPeriod <- na.aggregate(children_data$refrPeriod, by = children_data$patient_id)
children_data <- subset(children_data, select= -numCycles)
children_data$peakLoc <- children_data$peakLoc / CHILD_SAMP_FREQ

adults_data <- read.csv("~/spindles/adults_data_2025.csv")
adults_data$coupling_label <- as.factor(adults_data$coupling_label)
adults_data <- subset(adults_data, !grepl("O", channel))
adults_data$channel <- sub("-M[12]$", "", adults_data$channel)
adults_data <- tibble::rowid_to_column(adults_data, "ID")
adults_data$refrPeriod <- na.aggregate(adults_data$refrPeriod, by = adults_data$patient_id)
adults_data <- subset(adults_data, select= -numCycles)
adults_data$peakLoc <- adults_data$peakLoc / ADULT_SAMP_FREQ

## Grouped cross-validation ----

k <- 10
child_folds <- rsample::group_vfold_cv(children_data, group = patient_id, v = k)
adult_folds <- rsample::group_vfold_cv(adults_data, group = patient_id, v = k)

child_tuning_results <- read_csv("~/spindles/tuning_results_children.csv")
child_max_score_idx = which.max(child_tuning_results$`as.numeric(score[1])`)
child_best_params <- child_tuning_results[child_max_score_idx,-(ncol(child_tuning_results)-1)]
child_best_params_list <- as.list(child_best_params[,-ncol(child_best_params)])

adult_tuning_results <- read_csv("~/spindles/tuning_results_adults.csv")
adult_max_score_idx = which.max(adult_tuning_results$`as.numeric(score[1])`)
adult_best_params <- adult_tuning_results[adult_max_score_idx,-(ncol(adult_tuning_results)-1)]
adult_best_params_list <- as.list(adult_best_params[,-ncol(adult_best_params)])

final_auc_mat <- matrix(0, nrow = 3, ncol = 3)

child_all_indices <- seq_len(nrow(children_data))
adult_all_indices <- seq_len(nrow(adults_data))

all_oof_preds <- data.frame()
all_shap_values <- list()

for (i in 1:k) {
  cat("Processing fold", i, "\n")

  child_current_split <- child_folds$splits[[i]]
  child_train_indices <- child_current_split$in_id
  child_test_indices <- setdiff(child_all_indices, child_train_indices)

  adult_current_split <- adult_folds$splits[[i]]
  adult_train_indices <- adult_current_split$in_id
  adult_test_indices <- setdiff(adult_all_indices, adult_train_indices)

  child_data4model <- children_data[ , !names(children_data) %in% c("patient_id", "spindle_idx", "ID", "detSample", "startSample", "endSample")]
  child_train <- child_data4model[child_train_indices, ]
  child_test <- child_data4model[child_test_indices, ]

  adult_data4model <- adults_data[ , !names(adults_data) %in% c("patient_id", "spindle_idx", "ID", "detSample", "startSample", "endSample")]
  adult_train <- adult_data4model[adult_train_indices, ]
  adult_test <- adult_data4model[adult_test_indices, ]

  all_train <- rbind(child_train, adult_train)
  all_test <- rbind(child_test, adult_test)

  child_test_pids <- children_data$patient_id[child_test_indices]
  adult_test_pids <- adults_data$patient_id[adult_test_indices]

  all_test_pids <- c(child_test_pids, adult_test_pids)

  all_train_mat <- as.matrix(all_train[,sapply(all_train, is.numeric)])
  child_train_mat <- as.matrix(child_train[,sapply(child_train, is.numeric)])
  adult_train_mat <- as.matrix(adult_train[,sapply(adult_train, is.numeric)])

  all_test_mat <- as.matrix(all_test[,sapply(all_test, is.numeric)])
  child_test_mat <- as.matrix(child_test[,sapply(child_test, is.numeric)])
  adult_test_mat <- as.matrix(adult_test[,sapply(adult_test, is.numeric)])

  xgb_all_train <-  xgb.DMatrix(data = all_train_mat, label = as.numeric(all_train$coupling_label) -1)
  xgb_child_train <- xgb.DMatrix(data = child_train_mat, label = as.numeric(child_train$coupling_label) -1)
  xgb_adult_train <- xgb.DMatrix(data = adult_train_mat, label = as.numeric(adult_train$coupling_label) -1)
  xgb_all_test <-  xgb.DMatrix(data = all_test_mat, label = as.numeric(all_test$coupling_label) -1)
  xgb_child_test <- xgb.DMatrix(data = child_test_mat, label = as.numeric(child_test$coupling_label) -1)
  xgb_adult_test <- xgb.DMatrix(data = adult_test_mat, label = as.numeric(adult_test$coupling_label) -1)

  datasets <- list(
    list(name = "All",      data = xgb_all_test,   mat = all_test_mat,   labels = all_test$coupling_label,   pids = all_test_pids),
    list(name = "Children", data = xgb_child_test, mat = child_test_mat, labels = child_test$coupling_label, pids = child_test_pids),
    list(name = "Adults",   data = xgb_adult_test, mat = adult_test_mat, labels = adult_test$coupling_label, pids = adult_test_pids)
  )

  imbalance_weight_all <-  sum(all_train$coupling_label == 0) / sum(all_train$coupling_label == 1)
  imbalance_weight_child <- sum(child_train$coupling_label == 0) / sum(child_train$coupling_label == 1)
  imbalance_weight_adult <- sum(adult_train$coupling_label == 0) / sum(adult_train$coupling_label == 1)

  xgb_model_all <- xgb.train(
    params = adult_best_params_list,
    data = xgb_all_train,
    nrounds = adult_best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_all
  )

  xgb_model_child <- xgb.train(
    params = child_best_params_list,
    data = xgb_child_train,
    nrounds = child_best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_child
  )
  xgb_model_adult <- xgb.train(
    params = adult_best_params_list,
    data = xgb_adult_train,
    nrounds = adult_best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_adult
  )

  model_names <- c("All", "Children", "Adults")
  models <- list(xgb_model_all, xgb_model_child, xgb_model_adult)

  datasets <- list(
    list(name = "All",      data = xgb_all_test,   mat = all_test_mat,   labels = all_test$coupling_label,   pids = all_test_pids),
    list(name = "Children", data = xgb_child_test, mat = child_test_mat, labels = child_test$coupling_label, pids = child_test_pids),
    list(name = "Adults",   data = xgb_adult_test, mat = adult_test_mat, labels = adult_test$coupling_label, pids = adult_test_pids)
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

      if (model_names[l] == datasets[[j]]$name) {
        shap_fold <- shap.prep(xgb_model = models[[l]], X_train = datasets[[j]]$mat)
        all_shap_values[[paste0("fold_", i)]][[model_names[l]]] <- shap_fold
      }
    }
  }

}

saveRDS(all_shap_values, "ChildvsAdults_shap_values.rds")

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
    cat(sprintf("Skipped Patient ID: %-10s | Train: %-8s | Test: %-8s | Positives: %-3d | Negatives: %-3d\n", 
                row$patient_id, row$train_model, row$test_set, row$n_pos, row$n_neg))
  }
  cat("======================================================\n\n")
} else {
  cat("\nAll patients successfully evaluated (no missing classes).\n\n")
}

saveRDS(subject_metrics, file = "ChildvsAdults_subject_metrics.rds")

## Summarize model comparisons ----

process_auc_data <- function(metrics_df, metric_col) {
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

  desired_order <- c("Pooled", "Children", "Adults")

  summary_df <- summary_df %>%
    mutate(
      Model = case_when(Model == "All" ~ "Pooled", TRUE ~ as.character(Model)),
      Test = case_when(Test == "All" ~ "Pooled", TRUE ~ as.character(Test))
    ) %>%
    mutate(
      Model = factor(Model, levels = rev(desired_order)),
      Test  = factor(Test, levels = rev(desired_order))
    )

  summary_df$label <- sprintf("%.3f\n(%.3f-%.3f)", summary_df[[metric_col]], summary_df$Lower_CI, summary_df$Upper_CI)

  return(summary_df)
}

subject_metrics <- readRDS("ChildvsAdults_subject_metrics.rds")

auc_data_roc <- process_auc_data(subject_metrics, "ROCAUC")
auc_data_pr  <- process_auc_data(subject_metrics, "PRAUC")

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
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10) 
  )

## Export plots ----

combined_plot <- roc_plot + pr_plot + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 24, face = "bold"))

ggsave("combined_child_vs_adult_heatmaps.png", plot = combined_plot, width = 18, height = 8, dpi = 300)
