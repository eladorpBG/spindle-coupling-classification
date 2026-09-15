library(tidyverse)
library(stringr)

# ==============================================================================
# 1. SETUP: RENAMING, UNITS & HELPERS
# ==============================================================================

# --- A. VARIABLE RENAMING (Format: "Old Name" = "New Pretty Name") ---
unit_mapping <- c(
  "Spindle Density" = "atop(Spindle~Density, '(SP/min)')",
  "n"               = "atop(Count, '(n)')",
  "peakFreq"        = "atop(Peak~Frequency, '(Hz)')",
  "duration"        = "atop(Duration, '(sec)')",
  "energy"          = "atop(Energy, '(a.u.)')",
  
  # Plotmath syntax for Greek letters and superscripts 
  # (No string quotes around the bottom row so the math symbols render correctly)
  "sigmaPower"      = "atop('Sigma'~Power, (mu*V^2/Hz))", 
  "peakAmp"         = "atop(Peak~Amplitude, (mu*V))",
  "raisingSlope"    = "atop(Rising~Slope, (mu*V/s))",
  "droppingSlope"   = "atop(Dropping~Slope, (mu*V/s))",
  
  # Standard units (Wrapped in single quotes for safe string parsing)
  "refrPeriod"      = "atop(Refractory~Period, '(sec)')",
  "freqGradient"    = "atop(Frequency~Gradient, '(ms/cycle)')",
  "peakLoc"         = "atop(Peak~Time, '(sec from NREM start)')",
  "symmetry"        = "atop(Symmetry~Index, '(0-1; 0.5 denotes middle)')",
  
  # Unitless variables (Using an empty string for the second row to keep box heights uniform)
  "numBumps"        = "atop(Number~of~Bumps, '')",
  "corrCoef"        = "atop(Correlation~Coefficient, '')",
  "fano"            = "atop(Fano~Factor, '')"
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

# --- D. HELPER FUNCTIONS ---
extract_mean <- function(x) {
  if (grepl("\\(", x)) {
    as.numeric(sub("^(.*) \\(.*", "\\1", x))
  } else {
    as.numeric(gsub("[^0-9\\.]", "", x))
  }
}

# Function specific for 'n' (Single integer -> "Mean (0.00)")
calc_density_from_count <- function(x, divisor=720) {
  val <- as.numeric(gsub("[^0-9\\.]", "", x))
  if (is.na(val)) return(NA)
  sprintf("%.2f (0.00)", val / divisor) # Fake SD of 0 for plotting
}

# Function for regular "Mean (SD)" strings
scale_mean_sd_string <- function(x, divisor) {
  val_clean <- gsub("[^0-9\\.]", " ", x)
  parts <- as.numeric(unlist(strsplit(trimws(val_clean), "\\s+")))
  if(length(parts) < 2) return(NA)
  sprintf("%.2f (%.2f)", parts[1]/divisor, parts[2]/divisor)
}

# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

message("Processing data...")

# 1. Clean variable names
final_table_with_base <- final_table %>%
  mutate(base_variable = trimws(sub(" \\(mean \\(SD\\)\\)$", "", variable)))

# 2. Find 'n' and Calculate Spindle Density
n_row_index <- which(final_table_with_base$base_variable == "n")

if (length(n_row_index) > 0) {
  message("Found 'n' variable. Calculating Spindle Density...")
  n_row <- final_table_with_base[n_row_index, ]
  
  spindle_density_row <- n_row %>%
    mutate(
      base_variable = "Spindle Density",
      variable = "Spindle Density (mean (SD))",
      adults_uncoupled_mean_std   = sapply(adults_uncoupled_mean_std, calc_density_from_count),
      adults_coupled_mean_std     = sapply(adults_coupled_mean_std, calc_density_from_count),
      children_uncoupled_mean_std = sapply(children_uncoupled_mean_std, calc_density_from_count),
      children_coupled_mean_std   = sapply(children_coupled_mean_std, calc_density_from_count)
    )
  final_table_with_base <- bind_rows(final_table_with_base, spindle_density_row)
}

# ==============================================================================
# 3. ORDERING
# ==============================================================================

# Calculate diffs for sorting
mean_diff_adults_df <- final_table_with_base %>%
  filter(base_variable %in% vars_continuous | base_variable == "Spindle Density") %>% 
  select(base_variable, adults_uncoupled_mean_std, adults_coupled_mean_std) %>% 
  mutate(
    mean_uncoupled = sapply(adults_uncoupled_mean_std, extract_mean),
    mean_coupled = sapply(adults_coupled_mean_std, extract_mean),
    Group = "Adults"
  ) %>%
  mutate(Mean_Difference = mean_coupled - mean_uncoupled) %>%
  select(base_variable, Mean_Difference, Group)

# Define sort order (Largest Abs Diff First)
order_vars <- mean_diff_adults_df %>%
  arrange(desc(abs(Mean_Difference))) %>% 
  pull(base_variable) %>% 
  unique()

# ==============================================================================
# 4. PREPARE PLOT DATA
# ==============================================================================

# Combine Adults and Children
plot_data_children <- final_table_with_base %>% 
  filter(base_variable %in% vars_continuous | base_variable == "Spindle Density") %>% 
  mutate(
    mean_uncoupled = sapply(children_uncoupled_mean_std, extract_mean),
    mean_coupled = sapply(children_coupled_mean_std, extract_mean),
    SMD = as.numeric(children_smd),
    Group = "Children"
  ) %>% select(base_variable, mean_uncoupled, mean_coupled, SMD, Group)

plot_data_adults <- final_table_with_base %>% 
  filter(base_variable %in% vars_continuous | base_variable == "Spindle Density") %>% 
  mutate(
    mean_uncoupled = sapply(adults_uncoupled_mean_std, extract_mean),
    mean_coupled = sapply(adults_coupled_mean_std, extract_mean),
    SMD = as.numeric(adults_smd),
    Group = "Adults"
  ) %>% select(base_variable, mean_uncoupled, mean_coupled, SMD, Group)

combined_data <- bind_rows(plot_data_children, plot_data_adults) %>%
  pivot_longer(cols = c("mean_uncoupled", "mean_coupled"), names_to = "Mean_Type", values_to = "Mean_Value") %>%
  filter(!is.na(Mean_Value)) %>%
  filter(base_variable != "Spindle Density" & base_variable != "n") %>%
  mutate(Highlight_SMD = ifelse(SMD > 0.1, "SMD > 0.1", "SMD <= 0.1"))



# Apply Factor Order (Must match the ORIGINAL names)
combined_data$base_variable <- factor(combined_data$base_variable, levels = order_vars)

# ==============================================================================
# 5. GENERATE PLOT
# ==============================================================================

# ADD THIS: Create invisible boundary points to force the 'symmetry' axis to 0-1
symmetry_boundaries <- data.frame(
  base_variable = factor(c("symmetry", "symmetry"), levels = order_vars),
  Mean_Value = c(0.45, 0.55),
  Group = "Adults",             # Dummy value to satisfy the plot's global aesthetics
  Mean_Type = "mean_coupled"    # Dummy value to satisfy the plot's global aesthetics
)

png("combined_mean_values_plot.png", width = 1500, height = 1000)

p <- ggplot(combined_data, aes(x = Mean_Value, y = Group, color = Group, shape = Mean_Type, group = interaction(base_variable, Group))) +
  geom_line(aes(alpha = Highlight_SMD), linewidth = 1.5) +
  geom_point(aes(alpha = Highlight_SMD), size = 5) +
  # geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey", linewidth = 1.2) +
  geom_blank(data = symmetry_boundaries, aes(x = Mean_Value, y = Group)) +
  # --- USE THE PRETTY LABELLER HERE ---
  facet_wrap(~ base_variable, scales = "free_x", ncol = 3, 
             labeller = as_labeller(unit_mapping, default = label_parsed)) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  
  scale_color_manual(values = c("Adults" ="#e31a1c", "Children" ="#1f78b4")) +
  scale_shape_manual(values = c("mean_uncoupled" = 1, "mean_coupled" = 19), 
                     labels = c("Coupled", "Uncoupled")) +
  scale_alpha_manual(values = c("SMD > 0.1" = 1.0, "SMD <= 0.1" = 0.25),
                     na.translate = FALSE) +
  labs(title = NULL, subtitle = NULL, 
       x = NULL, y = "Variable", 
       color = "Group", shape = "Mean Type", alpha = "SMD Significance") +
  theme_bw(base_size = 25) +
  theme(legend.position = "bottom",
        # legend.box = "vertical",
        strip.background = element_rect(fill="grey90"), 
        strip.text = element_text(face = "bold", size = 25, lineheight = 1.1))

print(p)
dev.off()
message("Done! Saved to 'combined_mean_values_plot.png'")

