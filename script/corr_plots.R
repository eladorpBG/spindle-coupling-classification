# Plot adult and child feature correlations and combine the manuscript panels.

library(ggcorrplot)
library(patchwork)
library(dplyr)
library(zoo)

set.seed(123)

## Feature labels ----

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

unit_mapping <- c(
  "Spindle Density" = "Spindle~Density~(SP/min)",
  "n"               = "Count~(n)",
  "peakFreq"        = "Peak~Frequency~(Hz)",
  "duration"        = "Duration~(sec)",
  "energy"          = "Energy~(a.u.)",

  "sigmaPower"      = "'Sigma'~Power~(mu*V^2/Hz)", 
  "peakAmp"         = "Peak~Amplitude~(mu*V)",
  "raisingSlope"    = "Rising~Slope~(mu*V/s)",
  "droppingSlope"   = "Dropping~Slope~(mu*V/s)",

  "refrPeriod"      = "Refractory~Period~(sec)",
  "peakLoc"         = "Peak~Time~(sec)",
  "freqGradient"    = "Frequency~Gradient~(ms/cycle)",

  "numBumps"        = "Number~of~Bumps",
  "corrCoef"        = "Correlation~Coefficient",
  "symmetry"        = "Symmetry~Index",
  "fano"            = "Fano~Factor"
)

## Plot helpers ----

label_pretty_with_units <- function(original_names) {
  pretty <- variable_renames[original_names]
  pretty[is.na(pretty)] <- original_names[is.na(pretty)]

  units <- unit_mapping[original_names]
  units[is.na(units)] <- ""

  ifelse(units == "", 
         pretty, 
         paste0(pretty, "\n(", units, ")"))
}

get_stars <- function(p) {
  stars <- rep("", length(p))
  stars[p <= 0.05] <- "*"
  stars[p <= 0.01] <- "**"
  stars[p <= 0.001] <- "***"
  stars[is.na(p)] <- ""
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

  # Use a base data frame for correlation functions.
  raw_numeric_data <- data %>%
    select(any_of(names(variable_renames))) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.data.frame()

  current_names <- names(raw_numeric_data)
  pretty_names <- variable_renames[current_names]
  pretty_names[is.na(pretty_names)] <- current_names[is.na(pretty_names)]
  names(raw_numeric_data) <- pretty_names

  corr_matrix <- cor(raw_numeric_data, use = "everything")
  # Call cor_pmat explicitly to avoid package conflicts.
  p_matrix <- ggcorrplot::cor_pmat(raw_numeric_data)

  corr_matrix <- as.matrix(corr_matrix)
  p_matrix <- as.matrix(p_matrix)

  # Align p-value labels with the correlation matrix.
  dimnames(p_matrix) <- dimnames(corr_matrix)

  corr_matrix <- corr_matrix[fixed_variable_order, fixed_variable_order]
  p_matrix <- p_matrix[fixed_variable_order, fixed_variable_order]

  return(list(cor = corr_matrix, p = p_matrix))
}

## Load and process data ----

children_data_2025 <- read.csv("~/spindles/children_data_2025.csv")
adults_data_2025 <- read.csv("~/spindles/adults_data_2025.csv")

children_results <- process_data(children_data_2025, fixed_variable_order)
adults_results <- process_data(adults_data_2025, fixed_variable_order)

children_corr_mat <- children_results$cor
adults_corr_mat <- adults_results$cor

children_p_mat <- children_results$p
adults_p_mat <- adults_results$p

## Child correlations ----

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

## Adult correlations ----

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

message("Combining plots with patchwork...")

## Combine and export plots ----

combined_figure <- p_corr_adults + p_corr_children +

  plot_layout(guides = "collect") +

  plot_annotation(tag_levels = "A") &

  theme(
    plot.tag = element_text(size = 32, face = "bold"),

    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(1.5, "cm")
  )

# Use a wide canvas for the side-by-side panels.
png("manuscript_combined_correlations.png", width = 6000, height = 2500, res = 300)
print(combined_figure)
dev.off()

message("Done! Saved combined figure to 'manuscript_combined_correlations.png'")
