# Run pooled leave-one-subject-out XGBoost evaluation and plot cohort metrics.

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

## Load and preprocess data ----

children_data <- read.csv("~/spindles/children_data_2025.csv")
children_data$Cohort <- "Child"
children_data$peakLoc <- children_data$peakLoc / CHILD_SAMP_FREQ

adults_data <- read.csv("~/spindles/adults_data_2025.csv")
adults_data$Cohort <- "Adult"
adults_data$peakLoc <- adults_data$peakLoc / ADULT_SAMP_FREQ

data <- rbind(children_data, adults_data)

data$coupling_label <- as.factor(data$coupling_label)
data <- subset(data, !grepl("O", channel))
data$channel <- sub("-M[12]$", "", data$channel)
data <- tibble::rowid_to_column(data, "ID")
data$refrPeriod <- na.aggregate(data$refrPeriod, by = data$patient_id)
data <- subset(data, select = -numCycles)

## Model settings ----

# Pooled evaluation uses the adult tuning results.
tuning_results <- read_csv("~/spindles/tuning_results_adults.csv")
max_score_idx <- which.max(tuning_results$`as.numeric(score[1])`)
best_params <- tuning_results[max_score_idx, -(ncol(tuning_results)-1)]
params <- as.list(best_params[, -ncol(best_params)])
n_rounds <- best_params$`as.numeric(score[2])`

## Subject-level evaluation ----

loso_results <- data.frame()
all_shap_values <- list()

unique_patients <- unique(data$patient_id)
i <- 0
for (patient in unique_patients) {
  i <- i + 1

  # Hold out every row from the current patient.
  train_data <- data[data$patient_id != patient,]
  val_data <- data[data$patient_id == patient,]

  patient_cohort <- unique(val_data$Cohort)

  imbalance_weight <- sum(train_data$coupling_label == 0) / sum(train_data$coupling_label == 1)

  # Exclude identifiers and cohort metadata from the model matrix.
  cols_to_remove <- c("ID", "patient_id", "spindle_idx", "detSample", "startSample", "endSample", "Cohort")

  pat_train_mat <- as.matrix(train_data[, sapply(train_data, is.numeric)])
  pat_train_mat <- pat_train_mat[, !(colnames(pat_train_mat) %in% cols_to_remove)]

  val_mat <- as.matrix(val_data[, sapply(val_data, is.numeric)])
  val_mat <- val_mat[, !(colnames(val_mat) %in% cols_to_remove)]

  dtrain <- xgb.DMatrix(data = pat_train_mat, label = as.numeric(train_data$coupling_label) - 1)
  dval   <- xgb.DMatrix(data = val_mat, label = as.numeric(val_data$coupling_label) - 1)

  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = n_rounds,
    verbose = 0,
    scale_pos_weight = imbalance_weight
  )

  pred <- predict(model, newdata = dval)

  # ROC and PR AUC are calculated for each held-out patient.
  auc <- pROC::roc(val_data$coupling_label, pred, quiet = TRUE)$auc
  pr_result <- PRROC::pr.curve(scores.class0 = pred, weights.class0 = (as.numeric(val_data$coupling_label) - 1))
  pr_auc <- pr_result$auc.integral

  prev <- sum(val_data$coupling_label == 1) / nrow(val_data)
  pr_auc_prev_ratio <- if(prev > 0) pr_auc / prev else NA

  loso_results <- bind_rows(loso_results, data.frame(
    patient_id = patient,
    Cohort = patient_cohort,
    ROC_AUC = as.numeric(auc),
    PR_AUC = pr_auc,
    PRAUC_Prev_Ratio = pr_auc_prev_ratio
  ))

  cat(sprintf("Processed patient %d/%d (%s): ROC AUC = %.3f\n", i, length(unique_patients), patient_cohort, auc))
}

write.csv(loso_results, "pooled_loso_metrics.csv", row.names = FALSE)

## Summarize metrics ----

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

final_summary <- bind_rows(cohort_summary, overall_summary)
print(final_summary)

write.csv(final_summary, "pooled_loso_summary_stats.csv", row.names = FALSE)

## Plot metrics ----

loso_results <- read.csv("pooled_loso_metrics.csv")

# PR AUC chance baselines use prevalence within each cohort.
prevalence_baselines <- data %>%
  group_by(Cohort) %>%
  summarize(chance_level = sum(coupling_label == 1) / n(), .groups = 'drop')

cohort_colors <- c("Child" = "#1f78b4", "Adult" = "#e31a1c")

# The ROC reference line marks chance performance at 0.5.
roc_plot <- ggplot(loso_results, aes(x = ROC_AUC)) +
  geom_histogram(aes(fill = Cohort), color = "black", bins = 15, position = "identity", alpha = 0.7) +
  geom_vline(aes(xintercept = 0.5, linetype = "Chance Level"), color = "#984ea3", linewidth = 1.2) +
  scale_fill_manual(values = cohort_colors) +
  scale_linetype_manual(name = NULL, values = c("Chance Level" = "dashed")) + 

  guides(linetype = guide_legend(override.aes = list(color = "black"))) +

  labs(x = "ROC AUC", y = "Frequency", fill = "Cohort") +
  theme_minimal(base_size = 20) + 
  theme(            
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 16, color = "black")
  )

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
  )

combined_plot <- roc_plot + pr_plot + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 24, face = "bold"),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.width = unit(2.8, "lines"),
    legend.key.height = unit(2.8, "lines")
  )

ggsave("combined_auc_pooled_loso.png", plot = combined_plot, width = 16, height = 7, dpi = 300)
