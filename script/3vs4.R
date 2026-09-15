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

set.seed(123)
APPLY_FILTER <- FALSE
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
data4elec <- read.csv("~/spindles/adults_data_2025.csv")
samp_freq <- ADULT_SAMP_FREQ

# 1. Setup: Sort and assign unique ID
data4elec <- data4elec %>% 
  arrange(patient_id, startSample) %>% 
  mutate(row_id = row_number())

# 2. Convert to data.table for interval matching
dt <- as.data.table(data4elec)
setkey(dt, patient_id, startSample, endSample)

# 3. Find Overlaps (Strict: Any overlap > 0 counts)
# 'type = "any"' ensures if events touch even slightly, they are matched
overlaps <- foverlaps(dt, dt, type = "any", which = FALSE)

# 4. Filter for valid pairs (remove self-matches)
valid_links <- overlaps[row_id < i.row_id, .(row_id, i.row_id)]

# 5. Graph Clustering
# Group all connected events (A touches B, B touches C -> all one cluster)
g <- graph_from_data_frame(valid_links, directed = FALSE, vertices = data4elec$row_id)
comps <- components(g)

if (APPLY_FILTER) {
# 6. Apply Filter
  df_filtered <- data4elec %>%
    mutate(cluster_id = comps$membership) %>%  # Assign Cluster IDs
    group_by(patient_id, cluster_id) %>%       # Group by Patient and Cluster
    slice_sample(n = 1) %>%                    # Randomly keep 1 per group
    ungroup() %>%
    select(-cluster_id, -row_id)               # Remove helper columns
} else {
  df_filtered <- data4elec
}

# View result
cat("Original rows:", nrow(data4elec), "\n")
cat("Filtered rows:", nrow(df_filtered), "\n")

data4elec <- df_filtered
data4elec$coupling_label <- as.factor(data4elec$coupling_label)
data4elec <- subset(data4elec, !grepl("O", channel))
data4elec$channel <- ifelse(grepl("3", data4elec$channel), "3", data4elec$channel)
data4elec$channel <- ifelse(grepl("4", data4elec$channel), "4", data4elec$channel)
data4elec <- tibble::rowid_to_column(data4elec, "ID")
data4elec$refrPeriod <- na.aggregate(data4elec$refrPeriod, by = data4elec$patient_id)
data4elec <- subset(data4elec, select= -numCycles)
data4elec$peakLoc <- data4elec$peakLoc / samp_freq

# Number of folds for cross-validation
k <- 10
# folds <- createFolds(data4elec$patient_id, k = k, list = TRUE, returnTrain = FALSE)
folds <- rsample::group_vfold_cv(data4elec, group = patient_id, v = k)

# Create an empty matrix for AUC values
final_auc_mat <- matrix(0, nrow = 3, ncol = 3)

tuning_results_adults <- read_csv("~/spindles/tuning_results_adults.csv")
tune_results <- tuning_results_adults

max_score_idx = which.max(tune_results$`as.numeric(score[1])`)
best_params <- tune_results[max_score_idx,-(ncol(tune_results)-1)]
best_params_list <- as.list(best_params[,-ncol(best_params)])

all_indices <- seq_len(nrow(data4elec))
# Create an empty dataframe to hold all predictions
all_oof_preds <- data.frame()

for (i in 1:k) {
  cat("Processing fold", i, "\n")
  
  # Get the specific split object for the i-th fold
  current_split <- folds$splits[[i]]
  
  # Get train and test indices from the split object
  train_indices <- current_split$in_id
  test_indices <- setdiff(all_indices, train_indices)
  
  # The rest of your code remains the same
  data4model <- data4elec[ , !names(data4elec) %in% c("patient_id", "spindle_idx", "ID", "detSample", "startSample", "endSample")]
  train4elec <- data4model[train_indices, ]
  test4elec <- data4model[test_indices, ]
  
  three_train <- train4elec[train4elec$channel == "3", ]
  four_train <-train4elec[train4elec$channel == "4", ]
  min_samples <- min(nrow(train4elec), nrow(three_train), nrow(four_train))
  all_inds <- createDataPartition(train4elec$coupling_label, p = min_samples /  nrow(train4elec), list = FALSE, times = 1)
  all_train <- train4elec[all_inds,]
  three_inds <- createDataPartition(three_train$coupling_label, p = min_samples /  nrow(three_train), list = FALSE, times = 1)
  three_train <- three_train[three_inds,]
  four_inds <- createDataPartition(four_train$coupling_label, p = min_samples /  nrow(four_train), list = FALSE, times = 1)
  four_train <- four_train[four_inds,]
  three_test <- test4elec[test4elec$channel == "3", ]
  four_test <- test4elec[test4elec$channel == "4", ]
  
  # Keep track of patient IDs for the test sets
  test_patient_ids_all <- data4elec$patient_id[test_indices]
  test_patient_ids_3 <- test_patient_ids_all[test4elec$channel == "3"]
  test_patient_ids_4 <- test_patient_ids_all[test4elec$channel == "4"]
  
  ## Classification model - XGBoost
  
  all_train_mat <- as.matrix(all_train[,sapply(all_train, is.numeric)])
  three_train_mat <- as.matrix(three_train[,sapply(three_train, is.numeric)])
  four_train_mat <- as.matrix(four_train[,sapply(four_train, is.numeric)])
  
  all_test_mat <- as.matrix(test4elec[,sapply(test4elec, is.numeric)])
  three_test_mat <- as.matrix(three_test[,sapply(three_test, is.numeric)])
  four_test_mat <- as.matrix(four_test[,sapply(four_test, is.numeric)])
  
  xgb_all_train <-  xgb.DMatrix(data = all_train_mat, label = as.numeric(all_train$coupling_label) -1)
  xgb_three_train <- xgb.DMatrix(data = three_train_mat, label = as.numeric(three_train$coupling_label) -1)
  xgb_four_train <- xgb.DMatrix(data = four_train_mat, label = as.numeric(four_train$coupling_label) -1)
  xgb_all_test <-  xgb.DMatrix(data = all_test_mat, label = as.numeric(test4elec$coupling_label) -1)
  xgb_three_test <- xgb.DMatrix(data = three_test_mat, label = as.numeric(three_test$coupling_label) -1)
  xgb_four_test <- xgb.DMatrix(data = four_test_mat, label = as.numeric(four_test$coupling_label) -1)
  
  imbalance_weight_all <-  sum(all_train$coupling_label == 0) / sum(all_train$coupling_label == 1)
  imbalance_weight_three <- sum(three_train$coupling_label == 0) / sum(three_train$coupling_label == 1)
  imbalance_weight_four <- sum(four_train$coupling_label == 0) / sum(four_train$coupling_label == 1)
  
    ## Train model 
  xgb_model_all <- xgb.train(
    params = best_params_list,
    data = xgb_all_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_all
  )
  
  xgb_model_three <- xgb.train(
    params = best_params_list,
    data = xgb_three_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_three
  )
  xgb_model_four <- xgb.train(
    params = best_params_list,
    data = xgb_four_train,
    nrounds = best_params$`as.numeric(score[2])`,
    verbose = 1,
    scale_pos_weight = imbalance_weight_four
  )
  
  
  # Define models and datasets in lists for looping
  model_names <- c("All", "3", "4")
  models <- list(xgb_model_all, xgb_model_three, xgb_model_four)
  
  datasets <- list(
    list(name = "All", data = xgb_all_test, labels = test4elec$coupling_label, pids = test_patient_ids_all),
    list(name = "3", data = xgb_three_test, labels = three_test$coupling_label, pids = test_patient_ids_3),
    list(name = "4", data = xgb_four_test, labels = four_test$coupling_label, pids = test_patient_ids_4)
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
    }
  }
}

# ==============================================================================
# CALCULATE PER-SUBJECT METRICS
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

# 2. Identify and print skipped patients
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

# Optional: Save these per-subject metrics for your records
saveRDS(subject_metrics, file = "adults_3vs4_subject_metrics.rds")

# ==============================================================================
# 1. HELPER FUNCTION: Clean data and rename to Left/Right
# ==============================================================================
process_auc_data_3vs4 <- function(metrics_df, metric_col) {
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
  
  # Rename to "Pooled", "Central", "Frontal" and enforce factor order
  desired_order <- c("Left", "Right", "Pooled")
  
  summary_df <- summary_df %>%
    mutate(
      Model = case_when(
        Model == "3"   ~ "Left",
        Model == "4"   ~ "Right",
        Model == "All" ~ "Pooled",
        TRUE ~ as.character(Model)
      ),
      Test = case_when(
        Test == "3"   ~ "Left",
        Test == "4"   ~ "Right",
        Test == "All" ~ "Pooled",
        TRUE ~ as.character(Test)
      )
    ) %>%
    mutate(
      # rev() puts Pooled at the top of the Y-axis. 
      Model = factor(Model, levels = desired_order), 
      Test  = factor(Test, levels = desired_order)
    )
  
  return(summary_df)
}

# ==============================================================================
# 2. GENERATE PLOT DATA
# ==============================================================================

subject_metrics <- readRDS("adults_3vs4_subject_metrics.rds")

auc_data_roc <- process_auc_data_3vs4(subject_metrics, "ROCAUC")
auc_data_pr  <- process_auc_data_3vs4(subject_metrics, "PRAUC")

# ==============================================================================
# 3. BUILD INDIVIDUAL PLOTS (With Dual-Layer Text)
# ==============================================================================

# --- Plot A: ROC AUC ---
roc_plot <- ggplot(auc_data_roc, aes(x = Test, y = Model, fill = ROCAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 
  
  # Dual-Layer Text trick to keep leading zeros while using large fonts
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

# --- Plot B: PR AUC ---
pr_plot <- ggplot(auc_data_pr, aes(x = Test, y = Model, fill = PRAUC)) + 
  geom_tile(color = "white", linewidth = 1) + 
  
  # Dual-Layer Text trick
  geom_text(aes(label = sprintf("%.3f", PRAUC)), 
            color = "white", size = 7, fontface = "bold", nudge_y = 0.15) + 
  geom_text(aes(label = sprintf("[%.3f - %.3f]", Lower_CI, Upper_CI)), 
            color = "white", size = 4.5, nudge_y = -0.15) + 
  
  # Adjust limits to c(0, 1) if baseline PR AUC drops below 0.5
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

# ==============================================================================
# 4. COMBINE AND EXPORT WITH PATCHWORK
# ==============================================================================
combined_plot <- roc_plot + pr_plot + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 28, face = "bold"))

# Save at manuscript dimensions (Wide enough for 2 plots + tall enough for stacked text)
ggsave("combined_adults_lateral_heatmaps.png", plot = combined_plot, width = 18, height = 8, dpi = 300)



