##############################################
#             OSTEOCYTE NB MATURATION       #
##############################################
# Date: Thu Jul 04 11:13:17 2024 -------------
# Updated by: Daniela M. Roth


library(readxl)
library(ggplot2)
library(dplyr)
library(here)
library(car)   # For Levene's Test
library(FSA)   # For Kruskal-Wallis post-hoc test
library(tidyr)
library(ggrepel)

# Set file path 
file_path <- here("postnatal_nb/raw-data/nb_osteocyte_maturation_measurements.xlsx")


# Dmp1:Sost ratio ---------------------------------------------------------

# Read the data
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype, Sost, Dmp1, or ID is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(Sost) & !is.na(Dmp1) & !is.na(ID))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")),
         ID = factor(ID))

# Calculate Dmp1 to Sost ratio
filtered_data <- filtered_data %>%
  mutate(Dmp1_to_Sost = Dmp1 / Sost)

# Check normality for Dmp1 to Sost ratio
shapiro_test_Dmp1_to_Sost <- shapiro.test(filtered_data$Dmp1_to_Sost)
cat("\nShapiro-Wilk Normality Test Results for Dmp1 to Sost ratio:\n")
print(shapiro_test_Dmp1_to_Sost)

# Check homogeneity of variances for Dmp1 to Sost ratio
levene_test_Dmp1_to_Sost <- leveneTest(Dmp1_to_Sost ~ genotype, data = filtered_data)
cat("\nLevene's Test for Homogeneity of Variance Results for Dmp1 to Sost ratio:\n")
print(levene_test_Dmp1_to_Sost)

# Perform ANOVA or Kruskal-Wallis test for Dmp1 to Sost ratio
if (!is.null(shapiro_test_Dmp1_to_Sost$p.value) && !is.null(levene_test_Dmp1_to_Sost$`Pr(>F)`[1])) {
  if (shapiro_test_Dmp1_to_Sost$p.value > 0.05 & levene_test_Dmp1_to_Sost$`Pr(>F)`[1] > 0.05) {
    # ANOVA if assumptions are met
    anova_result_Dmp1_to_Sost <- aov(Dmp1_to_Sost ~ genotype, data = filtered_data)
    cat("\nANOVA Results for Dmp1 to Sost ratio:\n")
    print(summary(anova_result_Dmp1_to_Sost))
    
    # Tukey's HSD test for pairwise comparisons
    tukey_result_Dmp1_to_Sost <- TukeyHSD(anova_result_Dmp1_to_Sost)
    cat("\nTukey HSD Post-Hoc Test Results for Dmp1 to Sost ratio:\n")
    print(tukey_result_Dmp1_to_Sost)
  } else {
    # Kruskal-Wallis test if assumptions are not met
    kruskal_test_Dmp1_to_Sost <- kruskal.test(Dmp1_to_Sost ~ genotype, data = filtered_data)
    cat("\nKruskal-Wallis Test Results for Dmp1 to Sost ratio:\n")
    print(kruskal_test_Dmp1_to_Sost)
    
    # Pairwise Wilcoxon tests with p-value adjustment if Kruskal-Wallis test is significant
    if (kruskal_test_Dmp1_to_Sost$p.value < 0.05) {
      pairwise_wilcox_Dmp1_to_Sost <- pairwise.wilcox.test(filtered_data$Dmp1_to_Sost, filtered_data$genotype, p.adjust.method = "BH")
      cat("\nPairwise Wilcoxon Test Results with BH Adjustment for Dmp1 to Sost ratio:\n")
      print(pairwise_wilcox_Dmp1_to_Sost)
    }
  }
} else {
  cat("\nError: Unable to perform tests for Dmp1 to Sost ratio due to invalid test results.\n")
}

# Example: Adjusting for zero division
filtered_data$Dmp1_to_Sost <- with(filtered_data, ifelse(Dmp1 == 0 & Sost == 0, 0, Dmp1 / Sost))

# Define custom colors for genotypes
genotype_colors <- c("ctrl" = "dodgerblue", "Bmp2 cko" = "red1", "Bmp7 cko" = "goldenrod1")

# Plotting Dmp1:Sost ratio
bar_plot <- ggplot(filtered_data, aes(x = genotype, y = Dmp1_to_Sost, fill = genotype, color = genotype)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(), width = 0.6, alpha = 0.4) +
  geom_jitter(aes(y = Dmp1_to_Sost), position = position_jitter(width = 0.2), size = 3, alpha = 1) +
  labs(
    title = "Dmp1:Sost Nasal Bone",
    x = "Genotype",
    y = "Dmp1:Sost",
    fill = "Genotype",
    color = "Genotype"
  ) +
  scale_fill_manual(values = genotype_colors) +
  scale_color_manual(values = genotype_colors) +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold")
  )

print(bar_plot)

# Save plot
ggsave(
  here("postnatal_nb/figs/Dmp1_to_Sost_ratio_barplot.pdf"),
  plot = bar_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("postnatal_nb/figs/Dmp1_to_Sost_ratio_barplot.png"),
  plot = bar_plot,
  width = 7.5,
  height = 6
)


# Osteocyte maturation ----------------------------------------------------

# Filter out rows where genotype, Sost, Dmp1, Dmp1_Sost, DAPI, or side is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(Sost) & !is.na(Dmp1) & !is.na(Dmp1_Sost) & !is.na(DAPI) & !is.na(ID) & !is.na(side))

# Convert genotype and ID to factors
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")),
         ID = factor(ID),
         side = factor(side))

# Prepare raw counts for each marker
filtered_data <- filtered_data %>%
  mutate(
    Dmp1_count = Dmp1,
    Dmp1_Sost_count = Dmp1_Sost,
    Sost_count = Sost,
    neg_count = DAPI - (Dmp1 + Dmp1_Sost + Sost)
  )

# Calculate proportions and percentages normalized to DAPI
filtered_data <- filtered_data %>%
  mutate(
    Dmp1_pos = (Dmp1 / DAPI) * 100,
    Dmp1_Sost_pos = (Dmp1_Sost / DAPI) * 100,
    Sost_pos = (Sost / DAPI) * 100,
    neg = 100 - (Dmp1_pos + Dmp1_Sost_pos + Sost_pos)
  )

# Prepare data for plotting (convert to long format)
filtered_data_long <- filtered_data %>%
  pivot_longer(cols = c(Dmp1_pos, Dmp1_Sost_pos, Sost_pos, neg), names_to = "Marker", values_to = "Percentage") %>%
  mutate(ID_side = interaction(ID, side))  # Create the interaction column

# Add the raw counts to the long format data
filtered_data_long <- filtered_data_long %>%
  mutate(Count = case_when(
    Marker == "Dmp1_pos" ~ Dmp1_count,
    Marker == "Dmp1_Sost_pos" ~ Dmp1_Sost_count,
    Marker == "Sost_pos" ~ Sost_count,
    Marker == "neg" ~ neg_count
  ))

# Define the order of markers (bottom to top)
marker_order <- c("neg", "Dmp1_pos", "Dmp1_Sost_pos", "Sost_pos")

# Define custom colors for markers
marker_colors <- c(
  "Dmp1_pos" = "magenta",
  "Dmp1_Sost_pos" = "lightpink",
  "Sost_pos" = "floralwhite",
  "neg" = "lightblue1"
)

# Create stacked bar plot by measurement (ID) and genotype
stacked_bar_plot <- ggplot(filtered_data_long, aes(x = ID_side, y = Percentage, fill = factor(Marker, levels = marker_order))) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_label(data = filtered_data_long %>% filter(Count > 0),
             aes(label = Count, fill = factor(Marker, levels = marker_order)),
             position = position_stack(vjust = 0.5),
             size = 2,
             fontface = "bold",
             inherit.aes = TRUE,
             color = "black") +  # Add labels conditionally
  facet_wrap(~ genotype, scales = "free_x", nrow = 1) +
  labs(
    title = "Osteocyte maturation",
    x = "ID.side",
    y = "Proportion of nuclei",
    fill = "Marker"
  ) +
  scale_fill_manual(values = marker_colors, labels = c(
    "Dmp1_pos" = "Dmp1_pos",
    "Dmp1_Sost_pos" = "Dmp1_Sost_pos",
    "Sost_pos" = "Sost_pos",
    "neg" = "DAPI_only"
  )) +
  theme_minimal() +
  theme(text = element_text(size = 12),  # Adjust text size
    panel.grid = element_blank(),  # Remove grid lines
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.title.x = element_blank(),  # Remove x-axis title
    plot.title = element_text(size = 12, face = "bold"),  # Adjust plot title size
    panel.background = element_blank(),
    strip.text.x = element_text(size = 12, margin = margin(b = -0.5)),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
  )

print(stacked_bar_plot)

# Save the plot with adjusted dimensions
ggsave(
  here("postnatal_nb/figs/Osteocyte_maturation_stacked_barplot.pdf"),
  plot = stacked_bar_plot,
  width = 8,  # Adjusted width
  height = 5.5  # Adjusted height
)
ggsave(
  here("postnatal_nb/figs/Osteocyte_maturation_stacked_barplot.png"),
  plot = stacked_bar_plot,
  width = 9,  # Adjusted width
  height = 5.5  # Adjusted height
)


# Total nuclei ------------------------------------------------------------

# Scatterplot with individual points and box plot for Total DAPI counts by genotype
scatter_box_plot <- ggplot(filtered_data, aes(x = genotype, y = DAPI, color = genotype)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6, fill = "white", color = "dodgerblue") +  # Box plot with fill and border colors
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 1) +  # Scatterplot with jittered points
  labs(
    title = "Total number of nuclei",
    x = "Genotype",
    y = "DAPI Counts",
    color = "Genotype"
  ) +
  scale_color_manual(values = genotype_colors) +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold")
  )

print(scatter_box_plot)

# Save the plot
ggsave(
  here("postnatal_nb/figs/Total_DAPI_counts_scatter_box_plot.pdf"),
  plot = scatter_box_plot,
  width = 7.5,
  height = 6
)
ggsave(
  here("postnatal_nb/figs/Total_DAPI_counts_scatter_box_plot.png"),
  plot = scatter_box_plot,
  width = 4.5,
  height = 4
)

# Calculate summary statistics for DAPI counts by genotype
dapi_summary_stats <- filtered_data %>%
  group_by(genotype) %>%
  summarise(
    mean_dapi = mean(DAPI, na.rm = TRUE),
    sd_dapi = sd(DAPI, na.rm = TRUE),
    median_dapi = median(DAPI, na.rm = TRUE),
    min_dapi = min(DAPI, na.rm = TRUE),
    max_dapi = max(DAPI, na.rm = TRUE),
    .groups = 'drop'
  )

# Print the summary statistics
print(dapi_summary_stats)

# Perform Kruskal-Wallis test to compare DAPI counts between genotypes
kruskal_test_dapi <- kruskal.test(DAPI ~ genotype, data = filtered_data)
cat("\nKruskal-Wallis Test Results for DAPI counts:\n")
print(kruskal_test_dapi)

# If Kruskal-Wallis test is significant, perform post-hoc pairwise comparisons
if (kruskal_test_dapi$p.value < 0.05) {
  pairwise_wilcox_dapi <- pairwise.wilcox.test(filtered_data$DAPI, filtered_data$genotype, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon Test Results with BH Adjustment for DAPI counts:\n")
  print(pairwise_wilcox_dapi)
}



# Proportion summary ------------------------------------------------------

# Calculate proportions of Dmp1 and Sost positive cells normalized to DAPI counts
filtered_data <- filtered_data %>%
  mutate(
    Dmp1_pos_proportion = (Dmp1 / DAPI) * 100,
    Sost_pos_proportion = (Sost / DAPI) * 100
  )

# Calculate summary statistics for Dmp1 and Sost proportions by genotype
proportion_summary_stats <- filtered_data %>%
  group_by(genotype) %>%
  summarise(
    mean_Dmp1_proportion = mean(Dmp1_pos_proportion, na.rm = TRUE),
    sd_Dmp1_proportion = sd(Dmp1_pos_proportion, na.rm = TRUE),
    mean_Sost_proportion = mean(Sost_pos_proportion, na.rm = TRUE),
    sd_Sost_proportion = sd(Sost_pos_proportion, na.rm = TRUE),
    .groups = 'drop'
  )

# Print the summary statistics
print(proportion_summary_stats)

# Perform Kruskal-Wallis test to compare Dmp1 proportions between genotypes
kruskal_test_Dmp1 <- kruskal.test(Dmp1_pos_proportion ~ genotype, data = filtered_data)
cat("\nKruskal-Wallis Test Results for Dmp1 proportions:\n")
print(kruskal_test_Dmp1)

# Perform Kruskal-Wallis test to compare Sost proportions between genotypes
kruskal_test_Sost <- kruskal.test(Sost_pos_proportion ~ genotype, data = filtered_data)
cat("\nKruskal-Wallis Test Results for Sost proportions:\n")
print(kruskal_test_Sost)

# If Kruskal-Wallis tests are significant, perform post-hoc pairwise comparisons
if (kruskal_test_Dmp1$p.value < 0.05) {
  pairwise_wilcox_Dmp1 <- pairwise.wilcox.test(filtered_data$Dmp1_pos_proportion, filtered_data$genotype, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon Test Results with BH Adjustment for Dmp1 proportions:\n")
  print(pairwise_wilcox_Dmp1)
}

if (kruskal_test_Sost$p.value < 0.05) {
  pairwise_wilcox_Sost <- pairwise.wilcox.test(filtered_data$Sost_pos_proportion, filtered_data$genotype, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon Test Results with BH Adjustment for Sost proportions:\n")
  print(pairwise_wilcox_Sost)
}
