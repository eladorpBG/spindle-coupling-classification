library(ggcorrplot)
library(patchwork)
library(dplyr)
library(zoo)

set.seed(123)


# --- A. VARIABLE RENAMING (Format: "Old Name" = "New Pretty Name") ---
variable_renames <- c(
  "n"               = "Spindle Count",
  "Spindle Density" = "Spindle Density", # Created in script
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

# Define labels as character strings using plotmath notation
unit_mapping <- c(
  "Spindle Density" = "Spindle~Density~(SP/min)",
  "n"               = "Count~(n)",
  "peakFreq"        = "Peak~Frequency~(Hz)",
  "duration"        = "Duration~(sec)",
  "energy"          = "Energy~(a.u.)",
  
  # Plotmath syntax for Greek letters and superscripts:
  "sigmaPower"      = "'Sigma'~Power~(mu*V^2/Hz)", 
  "peakAmp"         = "Peak~Amplitude~(mu*V)",
  "raisingSlope"    = "Rising~Slope~(mu*V/s)",
  "droppingSlope"   = "Dropping~Slope~(mu*V/s)",
  
  # Standard units
  "refrPeriod"      = "Refractory~Period~(sec)",
  "peakLoc"         = "Peak~Time~(sec)",
  "freqGradient"    = "Frequency~Gradient~(ms/cycle)",
  
  # Unitless variables (Labels only)
  "numBumps"        = "Number~of~Bumps",
  "corrCoef"        = "Correlation~Coefficient",
  "symmetry"        = "Symmetry~Index",
  "fano"            = "Fano~Factor"
)

# --- C. CUSTOM LABELLER ---
# This looks up the *Original Name* to find the *Pretty Name* AND *Unit*
label_pretty_with_units <- function(original_names) {
  # 1. Get the pretty name (default to original if not defined)
  pretty <- variable_renames[original_names]
  pretty[is.na(pretty)] <- original_names[is.na(pretty)]
  
  # 2. Get the unit
  units <- unit_mapping[original_names]
  units[is.na(units)] <- ""
  
  # 3. Combine: "Pretty Name\n(Unit)" or just "Pretty Name"
  ifelse(units == "", 
         pretty, 
         paste0(pretty, "\n(", units, ")"))
}

get_stars <- function(p) {
  stars <- rep("", length(p))
  stars[p <= 0.05] <- "*"
  stars[p <= 0.01] <- "**"
  stars[p <= 0.001] <- "***"
  stars[is.na(p)] <- "" # Handles the blank half of the plot safely
  return(stars)
}

process_data <- function(raw_data, fixed_variable_order) {
  data <- raw_data
  data$coupling_label <- as.factor(data$coupling_label)
  data <- subset(data, !grepl("O", channel))
  data$channel <- sub("-M[12]$", "", data$channel)
  data <- tibble::rowid_to_column(data, "ID")
  data$refrPeriod <- na.aggregate(data$refrPeriod, by = data$coupling_label)
  data <- subset(data, select = -numCycles)
  
  # 1. Filter out only the numeric variables and FORCE base R data.frame
  raw_numeric_data <- data %>%
    select(any_of(names(variable_renames))) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.data.frame() # <-- This prevents tibble/dplyr interference
  
  # 2. Rename columns
  current_names <- names(raw_numeric_data)
  pretty_names <- variable_renames[current_names]
  pretty_names[is.na(pretty_names)] <- current_names[is.na(pretty_names)]
  names(raw_numeric_data) <- pretty_names
  
  # 3. Calculate matrices (explicitly using ggcorrplot to avoid package conflicts)
  corr_matrix <- cor(raw_numeric_data, use = "everything")
  p_matrix <- ggcorrplot::cor_pmat(raw_numeric_data) # <-- Forced ggcorrplot function
  
  # 4. Force to basic R matrices
  corr_matrix <- as.matrix(corr_matrix)
  p_matrix <- as.matrix(p_matrix)
  
  # Ensure the row/col names match exactly just to be safe
  dimnames(p_matrix) <- dimnames(corr_matrix)
  
  # 5. Order both matrices using your fixed variable order
  corr_matrix <- corr_matrix[fixed_variable_order, fixed_variable_order]
  p_matrix <- p_matrix[fixed_variable_order, fixed_variable_order]
  
  # Return both matrices as a list
  return(list(cor = corr_matrix, p = p_matrix))
}


## Load data
children_data_2025 <- read.csv("~/spindles/children_data_2025.csv")
adults_data_2025 <- read.csv("~/spindles/adults_data_2025.csv")

# Process data
children_results <- process_data(children_data_2025, fixed_variable_order)
adults_results <- process_data(adults_data_2025, fixed_variable_order)

# Extract correlation matrices
children_corr_mat <- children_results$cor
adults_corr_mat <- adults_results$cor

# Extract p-value matrices
children_p_mat <- children_results$p
adults_p_mat <- adults_results$p

# 3. Generate the Cleaned Plot
# --- CHILDREN PLOT ---

png("children_correlation_plot.png", width = 3000, height = 3000, res = 300)

p_corr_children <- ggcorrplot(children_corr_mat, 
                     p.mat = children_p_mat,          
                     sig.level = 0.05,                
                     method = "square",       
                     type = "lower",          
                     hc.order = FALSE,         
                     insig = "blank",                 
                     lab = FALSE,              
                     colors = c("#1f78b4", "white", "#e31a1c"), 
                     ggtheme = theme_bw(base_size = 14)) + 
  
  # --- NEW: Adds asterisks to the significant tiles ---
  geom_text(aes(label = get_stars(pvalue)), size = 8, color = "black", vjust = 0.75) +
  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL) 

print(p_corr_children)
dev.off()


# --- ADULTS PLOT ---
png("adults_correlation_plot.png", width = 3000, height = 3000, res = 300)

p_corr_adults <- ggcorrplot(adults_corr_mat, 
                     p.mat = adults_p_mat,            
                     sig.level = 0.05,                
                     method = "square",       
                     type = "lower",          
                     hc.order = FALSE,         
                     insig = "blank",                 
                     lab = FALSE,              
                     colors = c("#1f78b4", "white", "#e31a1c"), 
                     ggtheme = theme_bw(base_size = 14)) + 
  
  # --- NEW: Adds asterisks to the significant tiles ---
  geom_text(aes(label = get_stars(pvalue)), size = 8, color = "black", vjust = 0.75) +
  
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL) 

print(p_corr_adults)
dev.off()

### Combine plots
message("Combining plots with patchwork...")

# 1. Combine the plots
combined_figure <- p_corr_adults + p_corr_children +
  
  # 2. Merge the redundant color bar legends into one
  plot_layout(guides = "collect") +
  
  # 3. Add the A/B tags
  plot_annotation(tag_levels = "A") &
  
  # 4. Master Theme Overrides (applies to the whole combined graphic)
  theme(
    plot.tag = element_text(size = 32, face = "bold"), # Made slightly bigger for the huge canvas
    
    # --- LEGEND ENLARGEMENT ---
    legend.title = element_text(size = 24, face = "bold"),  # Bigger "Corr" text
    legend.text = element_text(size = 20),                  # Bigger numbers (1.0, 0.5...)
    legend.key.height = unit(3, "cm"),                      # Makes the color bar much TALLER
    legend.key.width = unit(1.5, "cm")                      # Makes the color bar THICKER
  )
  

# 5. Export the final figure
# Notice we double the width (from 3000 to 6000) because they are side-by-side!
png("manuscript_combined_correlations.png", width = 6000, height = 2500, res = 300)
print(combined_figure)
dev.off()

message("Done! Saved combined figure to 'manuscript_combined_correlations.png'")


