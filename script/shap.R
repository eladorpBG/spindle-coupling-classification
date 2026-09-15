library(xgboost)
library(data.table)
library(pROC)  
library(PRROC)
library(ggplot2)
library(SHAPforxgboost)
library(zoo)
library(tidyverse)
library(data.table)

aggregate_cohort_shap <- function(shap_list, cohort) {
  fold_tables <- lapply(names(shap_list), function(fold_name) {
    tbl <- shap_list[[fold_name]][[cohort]]
    if (is.null(tbl)) return(NULL)
    tbl <- as.data.table(tbl)
    tbl[, ID := paste(fold_name, ID, sep = "_")]
    tbl
  })
  combined <- rbindlist(Filter(Negate(is.null), fold_tables))
  
  # Global recompute so colors and ordering reflect the pooled data
  combined[, mean_value := mean(abs(value)), by = variable]
  
  # --- THE FIX: Use Quantile Clipping for the Color Scale ---
  combined[, stdfvalue := {
    # 1. Find the 1st and 99th percentiles to ignore extreme outliers
    q_vals <- quantile(rfvalue, probs = c(0.01, 0.99), na.rm = TRUE)
    min_val <- q_vals[1]
    max_val <- q_vals[2]
    
    # 2. Clip the raw values to these boundaries just for color scaling
    clipped_rf <- pmin(pmax(rfvalue, min_val), max_val)
    
    # 3. Calculate the standard min-max scaling using the clipped range
    rng_diff <- max_val - min_val
    if (rng_diff == 0) 0 else (clipped_rf - min_val) / rng_diff
  }, by = variable]
  
  combined
}

# ==============================================================================
# 1. HELPER FUNCTION: Process and Plot SHAP (WITH MASSIVE FONTS)
# ==============================================================================
create_shap_plot <- function(shap_data, plot_title) {
  shap_data <- as.data.frame(shap_data)
  
  # Rename Variables
  variable_renames <- c(
    "n"               = "Spindle Count",
    "Spindle Density" = "Spindle Density",
    "droppingSlope"   = "Dropping Slope",
    "raisingSlope"    = "Rising Slope",
    "peakAmp"         = "Peak Amplitude",
    "peakFreq"        = "Peak Frequency",
    "sigmaPower"      = "Sigma Power",
    "numBumps"        = "Number of Bumps",
    "refrPeriod"      = "Refractory Period",
    "corrCoef"        = "Correlation Coefficient",
    "freqGradient"    = "Frequency Gradient",
    "peakLoc"         = "Peak Time",
    "duration"        = "Duration",
    "energy"          = "Energy",
    "symmetry"        = "Symmetry",
    "fano"            = "Fano Factor"
  )
  
  shap_data$variable <- as.character(shap_data$variable) 
  for (old_name in names(variable_renames)) {
    shap_data$variable[shap_data$variable == old_name] <- variable_renames[old_name]
  }
  
  # Calculate Importance & Force Order
  importance_scores <- shap_data %>%
    group_by(variable) %>%
    summarise(mean_abs_shap = mean(abs(value))) %>%
    arrange(mean_abs_shap) 
  
  correct_order <- importance_scores$variable
  shap_data$variable <- factor(shap_data$variable, levels = rev(correct_order))
  
  # Generate Base Plot
  p <- shap.plot.summary(shap_data)
  
  # Dilute Points
  set.seed(123)
  diluted_points <- p$data %>%  
    group_by(variable) %>%
    sample_frac(0.1) %>%        
    ungroup()
  
  p$layers[[2]]$data <- diluted_points
  p$layers[[2]]$aes_params$alpha <- 0.2
  p$layers[[2]]$position$width <- 0.25
  
  # Increase Text Size
  # ggplot text size is in mm. 5 is roughly a 14pt font, 6 is roughly 17pt.
  p$layers[[3]]$aes_params$size <- 5.5 
  p <- p + 
    ggtitle(plot_title) +
    theme_bw(base_size = 20) + # Raises the baseline size for the whole plot
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold", color = "black"), # Larger variable names
      axis.text.x = element_text(size = 16, color = "black"),
      axis.title.x = element_text(size = 20, face = "bold", margin = ggplot2::margin(t = 10)),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      panel.grid.minor = element_blank() # Cleans up the background slightly
    )
  
  return(p)
}

# ==============================================================================
# 2. GENERATE AND COMBINE PLOTS
# ==============================================================================

all_shap_values <- readRDS("ChildvsAdults_shap_values.rds")
plot_children <- create_shap_plot(aggregate_cohort_shap(all_shap_values, "Children"), "Children")
plot_adults   <- create_shap_plot(aggregate_cohort_shap(all_shap_values, "Adults"),   "Adults")

combined_shap <- plot_adults /  plot_children + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') & 
  theme(
    plot.tag = element_text(size = 30, face = "bold"), 
    legend.position = "bottom"
  )

# Width 12 gives the x-axis room to stretch; Height 16 gives both plots room to stack
ggsave("combined_shap_summary.png", plot = combined_shap, width = 12, height = 16, dpi = 300)
