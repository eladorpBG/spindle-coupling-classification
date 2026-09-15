# 1. Point directly to the python binary in your env
# (This path exists based on your previous logs)
python_path <- "/home/eladorp/.conda/envs/ebm_env/bin/python"

# 2. Tell R exactly where python is
Sys.setenv(RETICULATE_PYTHON = python_path)

# 3. THE TRICK: Tell R this is the environment root
# This prevents it from searching for 'conda' to activate it
Sys.setenv(RETICULATE_PYTHON_ENV = "/home/eladorp/.conda/envs/ebm_env")

# 4. Load reticulate
library(reticulate)

# 5. Verify it worked (it should show 'ebm_env' and the path above)
py_config()

# 6. NOW load your package
library(ebm)
library(xgboost)
library(data.table)  # For efficient data handling
library(pROC)        # For computing AUC
library(PRROC)
library(ggplot2)
library(SHAPforxgboost)
library(zoo)
library(tidyverse)

set.seed(123)

## Load and Clean Data
adults_data_2025 <- read.csv("~/spindles/adults_data_2025.csv")
data <- adults_data_2025
data$coupling_label <- as.factor(data$coupling_label) 
data <- subset(data, !grepl("O", channel))
data$channel <- sub("-M[12]$", "", data$channel)
data <- tibble::rowid_to_column(data, "ID")
data$refrPeriod <- na.aggregate(data$refrPeriod, by = data$coupling_label)
data <- subset(data, select= -numCycles)
data <- subset(data, select= -corrCoef)

# ---------------------------------------------------------
# PART 1: LOPO-CV for Performance Metrics (Unchanged)
# ---------------------------------------------------------
lopo_auc <- list()
lopo_prauc <- list()
unique_patients <- unique(data$patient_id)

print("Starting LOPO-CV...")

for (patient in unique_patients) {
  train_data <- data[data$patient_id != patient,]
  val_data <- data[data$patient_id == patient,]
  
  # Clean columns for formula
  cols_to_remove <- c("ID", "patient_id", "spindle_idx", "detSample", "startSample", "endSample")
  train_clean <- train_data[, !(colnames(train_data) %in% cols_to_remove)]
  val_clean   <- val_data[, !(colnames(val_data) %in% cols_to_remove)]
  
  # Train EBM
  model <- ebm(
    coupling_label ~ ., 
    data = train_clean,
    interactions = 10,
    random_state = 42
  )
  
  # Predict
  pred_probs <- predict(model, newdata = val_clean, type = "response")
  if(is.matrix(pred_probs)) pred <- pred_probs[,2] else pred <- pred_probs
  
  # Metrics
  val_y_num <- as.numeric(val_data$coupling_label) - 1
  auc <- roc(val_data$coupling_label, pred, quiet = TRUE)$auc
  pr_result <- pr.curve(scores.class0 = pred, weights.class0 = val_y_num, curve = TRUE)
  
  lopo_auc[[as.character(patient)]] <- auc
  lopo_prauc[[as.character(patient)]] <- pr_result$auc.integral
  
  print(paste("Patient:", patient, "- AUC:", round(auc, 3)))
}

print(paste("Average AUC:", mean(unlist(lopo_auc))))

# ---------------------------------------------------------
# PART 2: Global Interpretability (The "Glassbox" View)
# ---------------------------------------------------------
# We train ONE model on ALL data to inspect what the model learned globally.
# This replaces the SHAP plots.

print("Training Global Model for Interpretation...")

# Prepare full dataset
cols_to_remove <- c("ID", "patient_id", "spindle_idx", "detSample", "startSample", "endSample")
full_data_clean <- data[, !(colnames(data) %in% cols_to_remove)]

global_model <- ebm(
  coupling_label ~ ., 
  data = full_data_clean,
  interactions = 10,
  random_state = 42
)

# --- A. Get Feature Importance (Mean Absolute Score) ---
# We calculate importance by looking at the average magnitude of terms across the dataset
all_terms <- predict(global_model, newdata = full_data_clean, type = "terms")
term_importance <- colMeans(abs(all_terms))
term_importance <- sort(term_importance, decreasing = TRUE)
top_features <- names(term_importance)[1:4] # Pick top 4 features to plot

print("Top 4 Important Features:")
print(top_features)

# --- B. Plot the Shape Functions (The EBM Curves) ---
# Function to plot the exact EBM lookup table for a feature
plot_ebm_shape <- function(model, data, feature_name) {
  
  # 1. Create a synthetic grid for this feature
  if(is.numeric(data[[feature_name]])) {
    # If numeric, scan from min to max
    grid_vals <- seq(min(data[[feature_name]], na.rm=TRUE), 
                     max(data[[feature_name]], na.rm=TRUE), length.out = 200)
  } else {
    # If categorical, take unique levels
    grid_vals <- unique(data[[feature_name]])
  }
  
  # 2. Create a dummy dataframe (other columns don't matter for additive terms!)
  dummy_df <- data[1:length(grid_vals), ] 
  # We fill it with the grid values for our target feature
  dummy_df[[feature_name]] <- grid_vals
  
  # 3. Predict TERMS (not probability)
  # This gives us the additive score contribution f(x)
  terms <- predict(model, newdata = dummy_df, type = "terms")
  
  # 4. Extract the specific column for this feature
  # Note: EBM might rename interaction terms, but main effects usually match
  if(feature_name %in% colnames(terms)) {
    score_vals <- terms[, feature_name]
    
    plot_df <- data.frame(FeatureValue = grid_vals, Score = score_vals)
    
    # 5. Plot
    p <- ggplot(plot_df, aes(x = FeatureValue, y = Score)) +
      geom_line(color = "darkblue", size = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      labs(title = paste("EBM Shape:", feature_name),
           y = "Contribution to Log-Odds (Score)",
           x = feature_name) +
      theme_minimal()
    
    # Add density rug (histogram at bottom) to show where data actually exists
    # We sample the real data for the rug to keep it fast
    real_data_sample <- data[[feature_name]][sample(1:nrow(data), min(500, nrow(data)))]
    p <- p + geom_rug(data = data.frame(v = real_data_sample), aes(x = v), inherit.aes = F, alpha = 0.2)
    
    return(p)
  } else {
    warning(paste("Could not find term for", feature_name))
    return(NULL)
  }
}

# --- C. Save Plots for Top Features ---
for (feat in top_features) {
  # Skip interaction terms (which contain " & ") for simple plotting
  if(!grepl(" & ", feat)) {
    p <- plot_ebm_shape(global_model, full_data_clean, feat)
    if(!is.null(p)) {
      ggsave(paste0("ebm_shape_", feat, ".png"), p, width = 6, height = 4)
      print(paste("Saved plot for:", feat))
    }
  }
}