# ## ######################################## ## #
#                 SCORING COLX POSITIVITY        #
# ## ######################################## ## #
# Date: Sat Jul 06 21:14:19 2024 ------------------
# Updated by: Daniela M. Roth

# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# Set file path 
file_path <- here("postnatal_cb/raw-data/colX_measurements_SOS.xlsx")

# Read the Excel file
data <- read_excel(file_path, sheet = "Sheet1")

# Filter out rows where genotype, DAPI_area, or ColX_area is NA
filtered_data <- data %>%
  filter(!is.na(genotype) & !is.na(DAPI_area) & !is.na(ColX_area))

# Convert genotype to factor
filtered_data <- filtered_data %>%
  mutate(genotype = factor(genotype, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")))

# Calculate proportion of ColX_area to DAPI_area
filtered_data <- filtered_data %>%
  mutate(proportion_ColX_to_DAPI = ColX_area / DAPI_area)

# Summary statistics by genotype
summary_stats <- filtered_data %>%
  group_by(genotype) %>%
  summarize(
    mean_proportion = mean(proportion_ColX_to_DAPI),
    sd_proportion = sd(proportion_ColX_to_DAPI)
  )

# Plotting the proportion of ColX_area to DAPI_area by genotype
plot <- ggplot(filtered_data, aes(x = genotype, y = proportion_ColX_to_DAPI, fill = genotype)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 1) +
  labs(
    title = "Chondrocyte hypertrophy",
    x = "Genotype",
    y = "ColX_area:DAPI_area",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c("ctrl" = "dodgerblue", "Bmp2 cko" = "red1", "Bmp7 cko" = "goldenrod1"),
                    name = "Genotype",
                    labels = c("ctrl" = "ctrl", "Bmp2 cko" = "Bmp2 cko", "Bmp7 cko" = "Bmp7 cko")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(color = "white"),
    panel.grid.major = element_line(color = "floralwhite"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "vertical"
  )

print(plot)

# Save the plot
ggsave(
  here("postnatal_cb/figs/Proportion_ColX_to_DAPI_plot.pdf"),
  plot = plot,
  width = 10,
  height = 7
)
ggsave(
  here("postnatal_cb/figs/Proportion_ColX_to_DAPI_plot.png"),
  plot = plot,
  width = 3,
  height = 7
)

# Open a sink connection to save the statistical results to a file
sink(here("postnatal_cb/figs/statistical_results_ColX.txt"))

# Statistical tests
shapiro_test <- filtered_data %>%
  group_by(genotype) %>%
  summarise(shapiro_p_value = shapiro.test(proportion_ColX_to_DAPI)$p.value)

cat("Shapiro-Wilk Normality Test Results:\n")
print(shapiro_test)

# Check normality assumption
if (all(shapiro_test$shapiro_p_value > 0.05)) {
  # If data is normal, perform ANOVA
  anova_result <- aov(proportion_ColX_to_DAPI ~ genotype, data = filtered_data)
  cat("\nANOVA Results for Proportion of ColX_area to DAPI_area:\n")
  print(summary(anova_result))
  
  # Post-hoc Tukey's HSD test
  tukey_result <- TukeyHSD(anova_result)
  cat("\nTukey HSD Post-Hoc Test Results:\n")
  print(tukey_result)
} else {
  # If data is not normal, perform Kruskal-Wallis test
  kruskal_test <- kruskal.test(proportion_ColX_to_DAPI ~ genotype, data = filtered_data)
  cat("\nKruskal-Wallis Test Results for Proportion of ColX_area to DAPI_area:\n")
  print(kruskal_test)
  
  # Pairwise Wilcoxon tests with p-value adjustment
  pairwise_wilcox <- pairwise.wilcox.test(filtered_data$proportion_ColX_to_DAPI, filtered_data$genotype, p.adjust.method = "BH")
  cat("\nPairwise Wilcoxon Test Results with BH Adjustment:\n")
  print(pairwise_wilcox)
}

# Calculate Cohen's d for Bmp2 cko vs ctrl
mean_diff_bc <- mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp2 cko"]) - 
  mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "ctrl"])
sd_pool_bc <- sqrt((sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp2 cko"])^2 + 
                      sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "ctrl"])^2) / 2)
cohens_d_bc <- mean_diff_bc / sd_pool_bc

# Calculate Cohen's d for Bmp7 cko vs ctrl
mean_diff_b7c <- mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp7 cko"]) - 
  mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "ctrl"])
sd_pool_b7c <- sqrt((sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp7 cko"])^2 + 
                       sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "ctrl"])^2) / 2)
cohens_d_b7c <- mean_diff_b7c / sd_pool_b7c

# Calculate Cohen's d for Bmp2 cko vs Bmp7 cko
mean_diff_b2b7 <- mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp2 cko"]) - 
  mean(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp7 cko"])
sd_pool_b2b7 <- sqrt((sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp2 cko"])^2 + 
                        sd(filtered_data$proportion_ColX_to_DAPI[filtered_data$genotype == "Bmp7 cko"])^2) / 2)
cohens_d_b2b7 <- mean_diff_b2b7 / sd_pool_b2b7

cat("Cohen's d for Bmp2 cko vs ctrl:", cohens_d_bc, "\n")
cat("Cohen's d for Bmp7 cko vs ctrl:", cohens_d_b7c, "\n")
cat("Cohen's d for Bmp2 cko vs Bmp7 cko:", cohens_d_b2b7, "\n")

# Close the sink connection
sink()
