library(ggpubr)
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
library(dplyr)
library(tidyr)

set.seed(0)
ADULT_SAMP_FREQ <- 200
CHILD_SAMP_FREQ <- 256

# ==========================================
# 1. LOAD & CLEAN DATA
# ==========================================

## Process Adults Data
adults_data <- adults_data_2025
adults_data$coupling_label <- as.factor(adults_data$coupling_label)
adults_data <- tibble::rowid_to_column(adults_data, "ID")
adults_data$refrPeriod <- na.aggregate(adults_data$refrPeriod, by = adults_data$coupling_label)
adults_data <- subset(adults_data, select= -numCycles)
adults_data <- subset(adults_data, !grepl("O", channel))
adults_data$channel <- sub("-M[12]$", "", adults_data$channel)
adults_data$peakLoc <- adults_data$peakLoc / ADULT_SAMP_FREQ

## Process Children Data
children_data <- children_data_2025 
children_data$coupling_label <- as.factor(children_data$coupling_label)
children_data <- tibble::rowid_to_column(children_data, "ID")
children_data$refrPeriod <- na.aggregate(children_data$refrPeriod, by = children_data$coupling_label)
children_data <- subset(children_data, select= -numCycles)
children_data <- subset(children_data, !grepl("O", channel))
children_data$channel <- sub("-M[12]$", "", children_data$channel)
children_data$peakLoc <- children_data$peakLoc / CHILD_SAMP_FREQ

# ==========================================
# 2. TABLE ONE GENERATION
# ==========================================

# Adults Table
adults_data_table <- adults_data[, 7:(ncol(adults_data))]
vars <- colnames(adults_data_table)
adults_tab_one <- CreateTableOne(vars = vars, data = adults_data_table, strata = "coupling_label")
adults_table_matrix <- print(adults_tab_one, smd = TRUE, quote = FALSE, noSpaces = TRUE)

# Children Table
children_data_table <- children_data[, 7:(ncol(children_data))]
children_tab_one <- CreateTableOne(vars = vars, data = children_data_table, strata = "coupling_label")
children_table_matrix <- print(children_tab_one, smd = TRUE, quote = FALSE, noSpaces = TRUE)

# ==========================================
# 3. IDENTIFY CONTINUOUS VARIABLES
# ==========================================
vars_continuous_adults <- names(adults_data_table)[sapply(adults_data_table, is.numeric)]
vars_continuous_children <- names(children_data_table)[sapply(children_data_table, is.numeric)]

# Use intersection of continuous variables present in both datasets
vars_continuous <- intersect(vars_continuous_adults, vars_continuous_children)
vars_continuous <- vars_continuous[vars_continuous != "coupling_label"] 

# ==========================================
# 4. SUMMARY DATA FRAMES & MERGING
# ==========================================
adults_df <- as.data.frame(adults_table_matrix)
children_df <- as.data.frame(children_table_matrix)

adults_summary <- data.frame(
  variable = rownames(adults_df),
  adults_uncoupled_mean_std = adults_df$`0`,
  adults_coupled_mean_std = adults_df$`1`,
  adults_smd = adults_df$SMD,
  stringsAsFactors = FALSE
)

children_summary <- data.frame(
  variable = rownames(children_df),
  children_uncoupled_mean_std = children_df$`0`,
  children_coupled_mean_std = children_df$`1`,
  children_smd = children_df$SMD,
  stringsAsFactors = FALSE
)

# Merge the two summary data frames
combined_summary_table <- merge(children_summary, adults_summary, by = "variable", all = TRUE)
combined_summary_table <- combined_summary_table[!grepl("coupling_label =", combined_summary_table$variable), ]

# Re-order columns
final_table <- combined_summary_table[, c(
  "variable",
  "children_uncoupled_mean_std",
  "children_coupled_mean_std",
  "children_smd",
  "adults_uncoupled_mean_std",
  "adults_coupled_mean_std",
  "adults_smd"
)]

message("Saving combined summary table to 'combined_summary_table.csv'")
write.csv(final_table, "combined_summary_table.csv", row.names = FALSE)
print(head(final_table, 10))

# ==========================================
# 5. COMBINED SMD PLOT (RESTORED)
# ==========================================
message("Generating combined SMD plot...")

smd_adults <- ExtractSmd(adults_tab_one)
smd_children <- ExtractSmd(children_tab_one)

dataPlot_adults <- data.frame(variable = rownames(smd_adults), SMD = as.numeric(smd_adults), Group = "Adults")
dataPlot_children <- data.frame(variable = rownames(smd_children), SMD = as.numeric(smd_children), Group = "Children")

dataPlot_combined <- rbind(dataPlot_children, dataPlot_adults)

var_to_remove <- rownames(ExtractSmd(adults_tab_one))[16]
if (!is.null(var_to_remove) && !is.na(var_to_remove)) {
  dataPlot_combined <- dataPlot_combined[dataPlot_combined$variable != var_to_remove, ]
}

adult_smds_for_ordering <- dataPlot_combined[dataPlot_combined$Group == "Adults", ]
varNames <- as.character(adult_smds_for_ordering$variable)[order(adult_smds_for_ordering$SMD)]
dataPlot_combined$variable <- factor(dataPlot_combined$variable, levels = varNames)

png("combined_smd_plot.png", width = 800, height = 1000)
combined_plot <- ggplot(data = dataPlot_combined,
                        mapping = aes(x = variable, y = SMD, group = Group, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) + 
  geom_hline(yintercept = 0.1, color = "black", linetype = "dashed", linewidth = 0.5) + 
  coord_flip() +
  theme_bw(base_size = 18) +
  labs(
    title = "Standardized Mean Differences (SMD)",
    subtitle = "Comparing Adults vs. Children Datasets",
    x = "Variable", y = "SMD"
  ) +
  theme(legend.position = "bottom", legend.title.align = 0.5,
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
print(combined_plot)
dev.off()

# ==========================================
# 6. MEAN VALUES (DUMBBELL) PLOT (RESTORED)
# ==========================================
extract_mean <- function(x) { as.numeric(sub("^(.*) \\(.*", "\\1", x)) }

mean_diff_adults_df <- final_table %>%
  filter(variable %in% vars_continuous) %>% 
  select(variable, adults_uncoupled_mean_std, adults_coupled_mean_std) %>%
  mutate(
    mean_uncoupled = as.numeric(sapply(adults_uncoupled_mean_std, extract_mean)) - 1,
    mean_coupled = as.numeric(sapply(adults_coupled_mean_std, extract_mean)) - 1,
    Group = "Adults"
  ) %>%
  mutate(Mean_Difference = mean_coupled - mean_uncoupled)

order_vars_by_adult_mean_diff <- mean_diff_adults_df %>%
  arrange(abs(Mean_Difference)) %>%
  pull(variable) %>%
  unique()

mean_vals_adults_df <- final_table %>%
  filter(variable %in% vars_continuous) %>%
  select(variable, adults_uncoupled_mean_std, adults_coupled_mean_std, adults_smd) %>%
  mutate(
    mean_uncoupled = sapply(adults_uncoupled_mean_std, extract_mean),
    mean_coupled = sapply(adults_coupled_mean_std, extract_mean),
    SMD = as.numeric(adults_smd),
    Group = "Adults"
  ) %>% select(variable, mean_uncoupled, mean_coupled, SMD, Group)

mean_vals_children_df <- final_table %>%
  filter(variable %in% vars_continuous) %>%
  select(variable, children_uncoupled_mean_std, children_coupled_mean_std, children_smd) %>%
  mutate(
    mean_uncoupled = as.numeric(sapply(children_uncoupled_mean_std, extract_mean)) - 1,
    mean_coupled = as.numeric(sapply(children_coupled_mean_std, extract_mean)) - 1,
    SMD = as.numeric(children_smd),
    Group = "Children"
  ) %>% select(variable, mean_uncoupled, mean_coupled, SMD, Group)

combined_mean_vals_plot_data <- rbind(mean_vals_children_df, mean_vals_adults_df) %>%
  pivot_longer(cols = c("mean_uncoupled", "mean_coupled"), names_to = "Mean_Type", values_to = "Mean_Value") %>%
  filter(!is.na(Mean_Value) & !is.na(SMD)) %>%
  mutate(Highlight_SMD = ifelse(SMD > 0.1, "SMD > 0.1", "SMD <= 0.1"))

combined_mean_vals_plot_data$variable <- factor(combined_mean_vals_plot_data$variable, levels = order_vars_by_adult_mean_diff)

png("combined_mean_values_plot.png", width = 1000, height = 1000)
mean_vals_plot <- ggplot(data = combined_mean_vals_plot_data,
                         aes(x = Mean_Value, y = variable, color = Group, shape = Mean_Type,
                             group = interaction(variable, Group))) +
  geom_line(aes(alpha = Highlight_SMD), position = position_dodge(width = 0.2), linewidth = 0.8) +
  geom_point(aes(alpha = Highlight_SMD), position = position_dodge(width = 0.2), size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey") +
  scale_color_manual(values = c("Adults" = "#1f78b4", "Children" = "#e31a1c")) +
  scale_shape_manual(values = c("mean_uncoupled" = 1, "mean_coupled" = 19), labels = c("Coupled", "Uncoupled")) +
  scale_alpha_manual(values = c("SMD > 0.1" = 1.0, "SMD <= 0.1" = 0.25)) +
  labs(title = "Mean Values for Coupled vs. Uncoupled Groups", x = "Mean Value", y = "Variable") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
print(mean_vals_plot)
dev.off()

# ==========================================
# 7. SPINDLE DENSITY PLOTS (RESTORED)
# ==========================================
adults_data_with_group <- adults_data %>% mutate(Group = "Adults")
children_data_with_group <- children_data %>% mutate(Group = "Children")
combined_spindle_data <- rbind(children_data_with_group, adults_data_with_group)

subject_density_data <- combined_spindle_data %>%
  group_by(Group, patient_id, channel) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(nesting(Group, patient_id), channel, fill = list(count = 0)) %>%
  mutate(Density = count / 180)

summary_density_data <- subject_density_data %>%
  group_by(Group, channel) %>%
  summarise(
    mean_density = mean(Density),
    sd_density = sd(Density),
    n_subjects = n(),
    sem = sd_density / sqrt(n_subjects),
    .groups = 'drop'
  ) %>%
  mutate(
    upper_limit = mean_density + sd_density,
    lower_limit = pmax(0, mean_density - sd_density) 
  )

png("spindle_density_per_channel_avg.png", width = 1000, height = 800)
spindle_density_plot <- ggplot(summary_density_data, aes(x = channel, y = mean_density)) +
  geom_col(aes(fill = Group)) +
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = 0.25, color = "black") +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  labs(title = "Mean Spindle Density per Channel (Averaged by Subject)", x = "Channel", y = "Mean Spindle Density (spindles/min)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
print(spindle_density_plot)
dev.off()