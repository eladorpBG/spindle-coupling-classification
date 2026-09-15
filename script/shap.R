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

# ### Dependence plots
# 
# # 1. Create a new directory to keep things organized
# output_dir <- "dependence_plots_children"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# 
# # 2. Get a list of all unique variables
# all_variables <- unique(children_shap$variable)
# 
# # 3. Loop through each variable and save the plot
# for (var_name in all_variables) {
#   
#   # A. Filter data for the current variable
#   plot_data <- children_shap[children_shap$variable == var_name, ]
#   
#   # B. Create the plot (Using GAM for speed)
#   p <- ggplot(plot_data, aes(x = rfvalue, y = value)) +
#     geom_point(alpha = 0.1, size = 0.8, color = "#2b8cbe") + # Low alpha for density
#     geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "red") + 
#     # FORCE Y-AXIS RANGE [-1.5, 1.5]
#     coord_cartesian(ylim = c(-3, 3)) +
#     labs(
#       title = paste("Dependence Plot:", var_name),
#       x = paste(var_name, "(Raw Value)"),
#       y = "SHAP Value (Impact on Model)"
#     ) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       axis.title = element_text(size = 12)
#     )
#   
#   # C. Construct a safe filename (remove spaces/special chars)
#   safe_name <- gsub(" ", "_", var_name) # Replace space with underscore
#   filename <- file.path(output_dir, paste0("shap_dep_", safe_name, ".png"))
#   
#   # D. Save
#   ggsave(filename, plot = p, width = 6, height = 4)
#   
#   message(paste("Saved:", filename))
# }
# 
# message("All plots generated successfully in folder: ", output_dir)
# 
# ##### Interactions plots
# 
# # 1. Setup Output Directory
# output_dir <- "interaction_plots_children"
# if (!dir.exists(output_dir)) dir.create(output_dir)
# 
# # 2. Define the Pairs you want to plot
# # Format: c("X_Variable", "Color_Variable")
# # Note: We use the RENAMED variable names (e.g., "Peak Frequency", not "peakFreq")
# pairs_to_plot <- list(
#   c("Fano Factor", "Correlation Coefficient"), # 1. Top 2 variables
#   c("Peak Frequency", "Peak Amplitude"),       # 2. peakFreq - peakAmp
#   c("Peak Frequency", "Duration"),             # 3. peakFreq - duration
#   c("Peak Amplitude", "Duration")              # 4. peakAmp - duration
# )
# 
# # 3. Safe Plotting Loop
# for (pair in pairs_to_plot) {
#   
#   x_name <- pair[1]
#   col_name <- pair[2]
#   
#   # A. Check if variables exist
#   if (!x_name %in% children_shap$variable || !col_name %in% children_shap$variable) {
#     message("Skipping: ", x_name, " - ", col_name)
#     next
#   }
#   
#   # B. PREPARE DATA (The Fix)
#   # We filter for the two variables and create a 'row_id' to guarantee alignment
#   pair_data <- children_shap %>%
#     filter(variable %in% c(x_name, col_name)) %>%
#     group_by(variable) %>%
#     mutate(row_id = row_number()) %>% # <--- CRITICAL STEP: Creates unique index 1..N
#     ungroup() %>%
#     select(row_id, variable, rfvalue, value)
#   
#   # C. PIVOT WIDER 
#   # This converts the data so we have columns for X and Color side-by-side
#   # Result: row_id | rfvalue_X | value_X | rfvalue_Color | value_Color
#   wide_data <- pair_data %>%
#     pivot_wider(
#       id_cols = row_id,
#       names_from = variable,
#       values_from = c(rfvalue, value),
#       names_sep = "__" # Separator to handle spaces in names
#     )
#   
#   # D. Extract vectors dynamically using the new column names
#   # Since names have spaces (e.g. "rfvalue__Fano Factor"), we use [[ ]] to select them
#   x_col_name   <- paste0("rfvalue__", x_name)
#   y_col_name   <- paste0("value__",   x_name)
#   col_col_name <- paste0("rfvalue__", col_name)
#   
#   # Remove rows where any needed value is NA (just in case of uneven lengths)
#   plot_ready <- wide_data %>%
#     filter(!is.na(!!sym(x_col_name)) & !is.na(!!sym(col_col_name)))
#   
#   # E. Plot
#   p <- ggplot(plot_ready, aes(x = .data[[x_col_name]], 
#                               y = .data[[y_col_name]], 
#                               color = .data[[col_col_name]])) +
#     
#     geom_point(alpha = 0.6, size = 1.2) + 
#     scale_color_viridis_c(option = "C", name = col_name) +
#     
#     # Fix Range
#     coord_cartesian(ylim = c(-1.5, 1.5)) +
#     
#     labs(
#       title = paste(x_name, "vs.", col_name),
#       x = paste(x_name, "(Raw Value)"),
#       y = paste("SHAP Value for", x_name)
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom")
#   
#   # F. Save
#   safe_filename <- paste0(gsub(" ", "", x_name), "_vs_", gsub(" ", "", col_name), ".png")
#   ggsave(file.path(output_dir, safe_filename), plot = p, width = 6, height = 5)
#   
#   message("Saved: ", safe_filename)
# }
# 
# message("All plots saved in: ", output_dir)
