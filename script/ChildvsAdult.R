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

eval_model_roc <- function(model, test_mat, test_labels){
  xgb_preds <- predict(model, test_mat, reshape = TRUE)
  print_confusion_mat(xgb_preds, test_labels)
  model_auc <- print_AUC_plot_ROC(xgb_preds, test_labels)
  return(model_auc)
}

eval_model_pr <- function(model, test_mat, test_labels){
  # 1. Get predictions
  xgb_preds <- predict(model, test_mat, reshape = TRUE)
  
  
  # 3. Calculate PR AUC
  # Note: If test_labels is already a 0/1 numeric vector (standard for xgboost), 
  # you don't need the `as.numeric() - 1`. If it's a factor, keep it.
  # I've used test_labels directly here assuming it's numeric.
  pr_result <- pr.curve(
    scores.class0 = xgb_preds, 
    weights.class0 = as.numeric(test_labels) - 1, 
    curve = TRUE
  )
  
  # 6. Return just the AUC metric
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
  # Ensure there are at least 2 data points to perform a test
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

## Load data
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

# Number of folds for cross-validation
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

# Create an empty matrix for AUC values
final_auc_mat <- matrix(0, nrow = 3, ncol = 3)

child_all_indices <- seq_len(nrow(children_data))
adult_all_indices <- seq_len(nrow(adults_data))

# Create an empty dataframe to hold all predictions
all_oof_preds <- data.frame()
all_shap_values <- list()

for (i in 1:k) {
  cat("Processing fold", i, "\n")
  
  # Get the specific split object for the i-th fold
  child_current_split <- child_folds$splits[[i]]
  # Get train and test indices from the split object
  child_train_indices <- child_current_split$in_id
  child_test_indices <- setdiff(child_all_indices, child_train_indices)
  
  # Get the specific split object for the i-th fold
  adult_current_split <- adult_folds$splits[[i]]
  # Get train and test indices from the split object
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
  
  # Keep track of patient IDs for the test sets
  child_test_pids <- children_data$patient_id[child_test_indices]
  adult_test_pids <- adults_data$patient_id[adult_test_indices]
  
  # all_test is an rbind of child and adult, so we concatenate the IDs in the same order
  all_test_pids <- c(child_test_pids, adult_test_pids)
  
  # prop.table(table(train4elec$channel))
  # prop.table(table(all_test$channel))
  # 
  # prop.table(table(child_train$coupling_label))
  # prop.table(table(adult_train$coupling_label))
  # 
  # prop.table(table(child_test$coupling_label))
  # prop.table(table(adult_test$coupling_label))
  # 
  # 
  ## adultlassification model - XGBoost
  
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
  
  ## Train model 
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
  
  
  # Define models and datasets in lists for looping
  model_names <- c("All", "Children", "Adults")
  models <- list(xgb_model_all, xgb_model_child, xgb_model_adult)
  
  datasets <- list(
    list(name = "All",      data = xgb_all_test,   mat = all_test_mat,   labels = all_test$coupling_label,   pids = all_test_pids),
    list(name = "Children", data = xgb_child_test, mat = child_test_mat, labels = child_test$coupling_label, pids = child_test_pids),
    list(name = "Adults",   data = xgb_adult_test, mat = adult_test_mat, labels = adult_test$coupling_label, pids = adult_test_pids)
  )
  
  # Fill the out-of-fold prediction dataframe
  for (l in seq_along(models)) {
    for (j in seq_along(datasets)) {
      # Generate predictions
      preds <- predict(models[[l]], datasets[[j]]$data, reshape = TRUE)
      
      # Create tracking dataframe for this specific combination
      temp_df <- data.frame(
        patient_id = datasets[[j]]$pids,
        label = as.numeric(datasets[[j]]$labels) - 1,
        pred = preds,
        train_model = model_names[l],
        test_set = datasets[[j]]$name
      )
      
      # Append to master dataframe
      all_oof_preds <- bind_rows(all_oof_preds, temp_df)
      
      if (model_names[l] == datasets[[j]]$name) {
        shap_fold <- shap.prep(xgb_model = models[[l]], X_train = datasets[[j]]$mat)
        all_shap_values[[paste0("fold_", i)]][[model_names[l]]] <- shap_fold
      }
    }
  }
  
}

saveRDS(all_shap_values, "ChildvsAdults_shap_values.rds")

# ==============================================================================
# CALCULATE PER-SUBJECT METRICS & PRINT SKIPPED
# ==============================================================================

subject_metrics <- all_oof_preds %>%
  group_by(train_model, test_set, patient_id) %>%
  summarize(
    n_pos = sum(label == 1),
    n_neg = sum(label == 0),
    # Only calculate if the subject has at least one example of both classes
    ROCAUC = if (n_pos > 0 && n_neg > 0) {
      as.numeric(pROC::auc(label, pred, quiet = TRUE))
    } else { NA_real_ },
    PRAUC = if (n_pos > 0 && n_neg > 0) {
      PRROC::pr.curve(scores.class0 = pred, weights.class0 = label)$auc.integral
    } else { NA_real_ },
    .groups = "drop"
  )

# Identify and print skipped patients
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

# ==============================================================================
# 1. HELPER FUNCTION: Clean data and create labels
# ==============================================================================
process_auc_data <- function(metrics_df, metric_col) {
  # Group by model and test set, calculate CIs while ignoring NAs
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
  
  # Rename to "Pooled" and enforce factor order
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
  
  # Format string label for heatmaps (Optional depending on how you use geom_text)
  summary_df$label <- sprintf("%.3f\n(%.3f-%.3f)", summary_df[[metric_col]], summary_df$Lower_CI, summary_df$Upper_CI)
  
  return(summary_df)
}

# ==============================================================================
# 2. GENERATE PLOT DATA
# ==============================================================================

subject_metrics <- readRDS("ChildvsAdults_subject_metrics.rds")

auc_data_roc <- process_auc_data(subject_metrics, "ROCAUC")
auc_data_pr  <- process_auc_data(subject_metrics, "PRAUC")

# ==============================================================================
# 3. BUILD INDIVIDUAL PLOTS
# ==============================================================================
# Plot A: ROC AUC
roc_plot <- ggplot(auc_data_roc, aes(x = Test, y = Model, fill = ROCAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 
  # 1. The Mean AUC: Large, bold, and shifted slightly UP
  geom_text(aes(label = sprintf("%.3f", ROCAUC)), 
            color = "white", size = 7, fontface = "bold", nudge_y = 0.15) + 
  
  # 2. The 95% CI: Slightly smaller, regular weight, and shifted slightly DOWN
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

# Plot B: PR AUC
pr_plot <- ggplot(auc_data_pr, aes(x = Test, y = Model, fill = PRAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 
  # 1. The Mean AUC: Large, bold, and shifted slightly UP
  geom_text(aes(label = sprintf("%.3f", PRAUC)), 
            color = "white", size = 7, fontface = "bold", nudge_y = 0.15) + 
  
  # 2. The 95% CI: Slightly smaller, regular weight, and shifted slightly DOWN
  geom_text(aes(label = sprintf("[%.3f - %.3f]", Lower_CI, Upper_CI)), 
            color = "white", size = 4.5, nudge_y = -0.15) + 
  # Note: You may need to change limits = c(0.5, 1) to c(0, 1) if your PR AUC scores drop below 0.5
  scale_fill_gradient(low = "#0072B2", high = "#D55E00", name = "PR AUC", limits = c(0.1, 0.5)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "Test Dataset", y = "") + # Y-axis label removed to avoid clutter
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

# ==============================================================================
# 4. COMBINE AND EXPORT WITH PATCHWORK
# ==============================================================================
combined_plot <- roc_plot + pr_plot + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 24, face = "bold"))

# Width is set to 16 to accommodate both 8-width plots side-by-side
ggsave("combined_child_vs_adult_heatmaps.png", plot = combined_plot, width = 18, height = 8, dpi = 300)

##### FOREST PLOTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

# Define a clean, colorblind-friendly palette for the test datasets
test_colors <- c("Adults" = "#0072B2", "Children" = "#009E73", "Pooled" = "#7F7F7F")

# Plot A: ROC AUC
roc_plot <- ggplot(auc_data_roc, aes(x = ROCAUC, y = Model, color = Test)) + 
  # Add the chance baseline FIRST so it sits behind the data
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkred", linewidth = 1, alpha = 0.7) +
  
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), 
                width = 0.3, linewidth = 1.2, position = position_dodge(width = 0.6)) + 
  geom_point(size = 5, position = position_dodge(width = 0.6)) + 
  
  # Set the wider x-axis range
  scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, by = 0.1)) +
  
  scale_color_manual(values = test_colors, name = "Test Dataset") +
  labs(x = "ROC AUC Score", y = "Model Trained On") +
  theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 26),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 18),
    axis.text.y = element_text(colour = "black", size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    # Emphasize the horizontal grid lines to separate the y-axis groups
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5), 
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10) 
  )

# Plot B: PR AUC
pr_plot <- ggplot(auc_data_pr, aes(x = PRAUC, y = Model, color = Test)) + 
  # Add the chance baseline (~0.15 average based on prevalence)
  geom_vline(xintercept = 0.15, linetype = "dashed", color = "darkred", linewidth = 1, alpha = 0.7) +
  
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), 
                width = 0.3, linewidth = 1.2, position = position_dodge(width = 0.6)) + 
  geom_point(size = 5, position = position_dodge(width = 0.6)) + 
  
  # Set the wider x-axis range
  scale_x_continuous(limits = c(0.0, 0.6), breaks = seq(0.0, 0.6, by = 0.1)) +
  
  scale_color_manual(values = test_colors, name = "Test Dataset") +
  labs(x = "PR AUC Score", y = "") + # Y-axis label removed to avoid clutter
  theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 26),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 18),
    axis.text.y = element_blank(), # Hide y-axis text since it shares with Plot A
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10) 
  )

# ==============================================================================
# 4. COMBINE AND EXPORT WITH PATCHWORK
# ==============================================================================
# Combine plots and collect the legend at the bottom so it isn't duplicated
combined_plot <- roc_plot + pr_plot + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 24, face = "bold"),
    legend.position = "bottom"
  )

# Save the plot
ggsave("combined_child_vs_adult_forest_plots.png", plot = combined_plot, width = 16, height = 7, dpi = 300)

# =======================================================
# STATISTICAL VALIDATION (PAIRED T-TESTS)
# =======================================================

# # Helper function to print formatted p-values
# get_pval <- function(model_A_scores, model_B_scores) {
#   # We use paired = TRUE because the AUCs come from the exact same folds
#   res <- t.test(model_A_scores, model_B_scores, paired = TRUE)
#   return(res$p.value)
# }
# 
# # Extract the raw AUC vectors (each is length 10)
# # Structure: results_array[Model, Test, Fold]
# 
# dimnames(results_array)[[1]] <- c("All", "Children", "Adults") 
# dimnames(results_array)[[2]] <- c("All", "Children", "Adults")
# 
# # 1. TEST ON CHILDREN: Child Model vs. Adult Model
# # Question: Do we really need a child-specific model?
# auc_child_model_on_child <- results_array["Children", "Children", ]
# auc_adult_model_on_child <- results_array["Adults", "Children", ]
# 
# p_val_child_vs_adult <- get_pval(auc_child_model_on_child, auc_adult_model_on_child)
# 
# # 2. TEST ON CHILDREN: Child Model vs. All Model
# # Question: Is the specific model better than the general model?
# auc_all_model_on_child   <- results_array["All", "Children", ]
# 
# p_val_child_vs_all <- get_pval(auc_child_model_on_child, auc_all_model_on_child)
# 
# # 3. TEST ON ADULTS: Adult Model vs. Child Model
# # Question: Does the child model fail on adults?
# auc_adult_model_on_adult <- results_array["Adults", "Adults", ]
# auc_child_model_on_adult <- results_array["Children", "Adults", ]
# 
# p_val_adult_vs_child <- get_pval(auc_adult_model_on_adult, auc_child_model_on_adult)
# 
# # --- PRINT RESULTS ---
# cat("\n--- Statistical Validation (Paired T-test, k=10) ---\n")
# cat("Testing on Child Data:\n")
# cat(sprintf("Child Model vs Adult Model: p = %.5f %s\n", 
#             p_val_child_vs_adult, ifelse(p_val_child_vs_adult < 0.05, "*Sig*", "")))
# cat(sprintf("Child Model vs All Model:   p = %.5f %s\n", 
#             p_val_child_vs_all, ifelse(p_val_child_vs_all < 0.05, "*Sig*", "")))
# 
# cat("\nTesting on Adult Data:\n")
# cat(sprintf("Adult Model vs Child Model: p = %.5f %s\n", 
#             p_val_adult_vs_child, ifelse(p_val_adult_vs_child < 0.05, "*Sig*", "")))
