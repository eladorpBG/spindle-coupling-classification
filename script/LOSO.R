library(xgboost)
library(data.table)  
library(pROC)        
library(PRROC)
library(ggplot2)
library(SHAPforxgboost)
library(zoo)
library(tidyverse)
library(patchwork)

set.seed(123)
ADULT_SAMP_FREQ <- 200
CHILD_SAMP_FREQ <- 256

## 1. LOAD AND POOL DATA
children_data <- read.csv("~/spindles/children_data_2025.csv")
children_data$Cohort <- "Child"
children_data$peakLoc <- children_data$peakLoc / CHILD_SAMP_FREQ

adults_data <- read.csv("~/spindles/adults_data_2025.csv")
adults_data$Cohort <- "Adult"
adults_data$peakLoc <- adults_data$peakLoc / ADULT_SAMP_FREQ

# Combine into one pooled dataset
data <- rbind(children_data, adults_data)

# Preprocessing
data$coupling_label <- as.factor(data$coupling_label)
data <- subset(data, !grepl("O", channel))
data$channel <- sub("-M[12]$", "", data$channel)
data <- tibble::rowid_to_column(data, "ID")
data$refrPeriod <- na.aggregate(data$refrPeriod, by = data$patient_id)
data <- subset(data, select = -numCycles)

## 2. SET PARAMETERS
tuning_results <- read_csv("~/spindles/tuning_results_adults.csv") # UPDATE THIS IF NEEDED
max_score_idx <- which.max(tuning_results$`as.numeric(score[1])`)
best_params <- tuning_results[max_score_idx, -(ncol(tuning_results)-1)]
params <- as.list(best_params[, -ncol(best_params)])
n_rounds <- best_params$`as.numeric(score[2])`

## 3. INITIALIZE TRACKING DATAFRAME
# Using a data.frame is much cleaner for stratified analysis than separate lists
loso_results <- data.frame()
all_shap_values <- list()

unique_patients <- unique(data$patient_id)
i <- 0
for (patient in unique_patients) {
  i <- i + 1
  
  # Split data into training and validation sets
  train_data <- data[data$patient_id != patient,]
  val_data <- data[data$patient_id == patient,]
  
  # Identify if this patient is an Adult or Child
  patient_cohort <- unique(val_data$Cohort)
  
  imbalance_weight <- sum(train_data$coupling_label == 0) / sum(train_data$coupling_label == 1)
  
  # Setup Matrices
  cols_to_remove <- c("ID", "patient_id", "spindle_idx", "detSample", "startSample", "endSample", "Cohort")
  
  pat_train_mat <- as.matrix(train_data[, sapply(train_data, is.numeric)])
  pat_train_mat <- pat_train_mat[, !(colnames(pat_train_mat) %in% cols_to_remove)]
  
  val_mat <- as.matrix(val_data[, sapply(val_data, is.numeric)])
  val_mat <- val_mat[, !(colnames(val_mat) %in% cols_to_remove)]
  
  # Convert to xgb.DMatrix
  dtrain <- xgb.DMatrix(data = pat_train_mat, label = as.numeric(train_data$coupling_label) - 1)
  dval   <- xgb.DMatrix(data = val_mat, label = as.numeric(val_data$coupling_label) - 1)
  
  # Train model
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = n_rounds,
    verbose = 0,
    scale_pos_weight = imbalance_weight
  )
  
  # Predict and Calculate Metrics
  pred <- predict(model, newdata = dval)
  
  auc <- pROC::roc(val_data$coupling_label, pred, quiet = TRUE)$auc
  pr_result <- PRROC::pr.curve(scores.class0 = pred, weights.class0 = (as.numeric(val_data$coupling_label) - 1))
  pr_auc <- pr_result$auc.integral
  
  prev <- sum(val_data$coupling_label == 1) / nrow(val_data)
  pr_auc_prev_ratio <- if(prev > 0) pr_auc / prev else NA
  
  # Store results in the dataframe
  loso_results <- bind_rows(loso_results, data.frame(
    patient_id = patient,
    Cohort = patient_cohort,
    ROC_AUC = as.numeric(auc),
    PR_AUC = pr_auc,
    PRAUC_Prev_Ratio = pr_auc_prev_ratio
  ))
  
  # SHAP
  # shap_result_fold <- shap.prep(xgb_model = model, X_train = val_mat)
  # all_shap_values[[as.character(patient)]] <- shap_result_fold
  
  cat(sprintf("Processed patient %d/%d (%s): ROC AUC = %.3f\n", i, length(unique_patients), patient_cohort, auc))
}

# Save raw results
write.csv(loso_results, "pooled_loso_metrics.csv", row.names = FALSE)

# ==============================================================================
# AGGREGATE RESULTS (WITH 95% CONFIDENCE INTERVALS)
# ==============================================================================

# Helper functions to safely extract 95% CI bounds
calc_lower_ci <- function(x) {
  x <- na.omit(x)
  if(length(x) < 2) return(NA_real_)
  t.test(x)$conf.int[1]
}

calc_upper_ci <- function(x) {
  x <- na.omit(x)
  if(length(x) < 2) return(NA_real_)
  t.test(x)$conf.int[2]
}

# Calculate overall pooled metrics
overall_summary <- loso_results %>%
  summarize(
    Cohort = "Pooled (All)",
    N = n(),
    Avg_ROC = mean(ROC_AUC, na.rm = TRUE), 
    Lower_CI_ROC = calc_lower_ci(ROC_AUC),
    Upper_CI_ROC = calc_upper_ci(ROC_AUC),
    
    Avg_PR = mean(PR_AUC, na.rm = TRUE), 
    Lower_CI_PR = calc_lower_ci(PR_AUC),
    Upper_CI_PR = calc_upper_ci(PR_AUC)
  )

# Calculate cohort-specific metrics
cohort_summary <- loso_results %>%
  group_by(Cohort) %>%
  summarize(
    N = n(),
    Avg_ROC = mean(ROC_AUC, na.rm = TRUE), 
    Lower_CI_ROC = calc_lower_ci(ROC_AUC),
    Upper_CI_ROC = calc_upper_ci(ROC_AUC),
    
    Avg_PR = mean(PR_AUC, na.rm = TRUE), 
    Lower_CI_PR = calc_lower_ci(PR_AUC),
    Upper_CI_PR = calc_upper_ci(PR_AUC),
    .groups = "drop"
  )

# Combine and print the final table
final_summary <- bind_rows(cohort_summary, overall_summary)
print(final_summary)

write.csv(final_summary, "pooled_loso_summary_stats.csv", row.names = FALSE)


# ==============================================================================
# PLOTTING
# ==============================================================================

loso_results <- read.csv("pooled_loso_metrics.csv")

# 1. Calculate the overall prevalence per cohort for the PR AUC baseline
prevalence_baselines <- data %>%
  group_by(Cohort) %>%
  summarize(chance_level = sum(coupling_label == 1) / n(), .groups = 'drop')

# Colors: Adults = #e31a1c, Children = #1f78b4, Purple (Pooled) = #984ea3
cohort_colors <- c("Child" = "#1f78b4", "Adult" = "#e31a1c")

# Plot A: ROC AUC Histogram (Purple baseline at 0.5)
roc_plot <- ggplot(loso_results, aes(x = ROC_AUC)) +
  geom_histogram(aes(fill = Cohort), color = "black", bins = 15, position = "identity", alpha = 0.7) +
  geom_vline(aes(xintercept = 0.5, linetype = "Chance Level"), color = "#984ea3", linewidth = 1.2) +
  scale_fill_manual(values = cohort_colors) +
  scale_linetype_manual(name = NULL, values = c("Chance Level" = "dashed")) + 
  
  # FIX 1: Force the linetype legend key to be black, regardless of the plot's line color
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  
  labs(x = "ROC AUC", y = "Frequency", fill = "Cohort") +
  theme_minimal(base_size = 20) + 
  theme(            
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 16, color = "black")
    # legend.position removed from here so Patchwork handles it globally
  )

# Plot B: PR AUC Histogram (Cohort-specific baselines)
pr_plot <- ggplot(loso_results, aes(x = PR_AUC)) +
  geom_histogram(aes(fill = Cohort), color = "black", bins = 15, position = "identity", alpha = 0.7) +
  geom_vline(data = prevalence_baselines, 
             aes(xintercept = chance_level, color = Cohort), 
             linetype = "dashed", linewidth = 1.2, show.legend = FALSE) +
  scale_fill_manual(values = cohort_colors) +
  scale_color_manual(values = cohort_colors, guide = "none") + 
  labs(x = "PR AUC", y = "Frequency", fill = "Cohort") +
  theme_minimal(base_size = 20) + 
  theme(            
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 16, color = "black")
    # legend.position removed from here so Patchwork handles it globally
  )

# Combine using patchwork
combined_plot <- roc_plot + pr_plot + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 24, face = "bold"),
    
    # FIX 2: Move the global collected legend to the bottom so it doesn't clash with A/B tags
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.width = unit(2.8, "lines"),
    legend.key.height = unit(2.8, "lines")
  )

ggsave("combined_auc_pooled_loso.png", plot = combined_plot, width = 16, height = 7, dpi = 300)

