# ## ######################################## ## #
#              SCORING NB TRABECULATION          #
# ## ######################################## ## #

# Date: Wed Jul 03 20:11:33 2024 ------------------
# Updated by: Daniela M. Roth

library(tidyverse)  # for data manipulation and ggplot
library(readxl)     # for reading Excel files
library(here)       # for setting the working directory
library(FSA)        # for statistical tests

# Set the working directory
nb.scoring <- here::here("postnatal_nb")

# Read the Excel file
data <- read_excel(here::here(nb.scoring, "raw-data/nb_trabeculation_measurements.xlsx"))

# Remove 'K' prefix from ID column and convert to numeric
data$ID <- as.numeric(sub("K", "", data$ID))


# Summarize data by ID and genotype
summarised_data <- data %>%
  group_by(ID, genotype) %>%
  summarise(
    total_bone = mean(total_bone, na.rm = TRUE),
    num_trab = mean(num_trab, na.rm = TRUE),
    avg_trab = mean(avg_trab, na.rm = TRUE),
    total_trab = mean(total_trab, na.rm = TRUE),
    num_measurements = n()  # Count the number of measurements per unique ID
  ) %>%
  ungroup()

# Calculate normalized average trabecular area
summarised_data <- summarised_data %>%
  mutate(normalized_avg_trab = avg_trab / total_bone)

# Set the order of the genotype factor levels
summarised_data$genotype <- factor(summarised_data$genotype, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko"))


# Get the number of unique IDs per genotype
unique_ids_per_genotype <- summarised_data %>%
  group_by(genotype) %>%
  summarise(unique_ids = n_distinct(ID), .groups = 'drop')

# Print the summary of unique IDs per genotype
print(unique_ids_per_genotype)

# Get descriptive statistics for each genotype
descriptive_stats <- summarised_data %>%
  group_by(genotype) %>%
  summarise(
    mean_avg_trab = mean(avg_trab, na.rm = TRUE),
    sd_avg_trab = sd(avg_trab, na.rm = TRUE),
    median_avg_trab = median(avg_trab, na.rm = TRUE),
    min_avg_trab = min(avg_trab, na.rm = TRUE),
    max_avg_trab = max(avg_trab, na.rm = TRUE),
    mean_total_trab = mean(total_trab, na.rm = TRUE),
    sd_total_trab = sd(total_trab, na.rm = TRUE),
    median_total_trab = median(total_trab, na.rm = TRUE),
    min_total_trab = min(total_trab, na.rm = TRUE),
    max_total_trab = max(total_trab, na.rm = TRUE),
    mean_num_measurements = mean(num_measurements, na.rm = TRUE),
    sd_num_measurements = sd(num_measurements, na.rm = TRUE),
    median_num_measurements = median(num_measurements, na.rm = TRUE),
    min_num_measurements = min(num_measurements, na.rm = TRUE),
    max_num_measurements = max(num_measurements, na.rm = TRUE),
    .groups = 'drop'
  )

# Print the descriptive statistics
print(descriptive_stats)

# Get the number of measurements per unique ID
measurements_per_id <- data %>%
  group_by(ID) %>%
  summarise(num_measurements = n(), .groups = 'drop')

# Print the summary of measurements per ID
print(measurements_per_id)

# Conduct Kruskal-Wallis test to compare genotypes for avg_trab
kruskal_result <- kruskal.test(avg_trab ~ genotype, data = data)
print(kruskal_result)

# If Kruskal-Wallis test is significant, perform post-hoc tests (Dunn's test)
if (kruskal_result$p.value < 0.05) {
  posthoc <- dunnTest(data$avg_trab, g = data$genotype, method = "bonferroni")
  print(posthoc)
}

# Create "figs" directory if it doesn't exist
figs_dir <- here::here(nb.scoring, "figs")
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir)
}

# Define a custom minimal theme for the plots
custom_theme <- theme_minimal() +
  theme(
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    panel.grid.major = element_blank(),
    plot.background = element_rect(color = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 9),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "vertical",
    axis.text.y = element_blank(),
    plot.title = element_text(size = 9, face = "bold")
  )

# Visualization: Violin plot with boxplot to visualize genotype differences in average trabeculae area
plot_avg_trab <- ggplot(summarised_data, aes(x = genotype, y = avg_trab, fill = genotype)) +
  geom_violin(trim = FALSE, alpha = 0.2) +  # Keep orientation default (vertical)
  geom_boxplot(width = 0.1, outlier.size = 1, alpha = 0.5) +
  labs(x = "genotype", y = "avg trabecular area (px)",
       title = "Average trabecular area (AT)") +
  custom_theme +
  scale_fill_manual(values = c("ctrl" = "dodgerblue", "Bmp2 cko" = "red1", "Bmp7 cko" = "goldenrod1")) +  # Set violin plot colors
  coord_flip()  # Flip coordinates to horizontal orientation

# Save the plot as a PNG file
ggsave(here::here(figs_dir, "01_avg_trab_boxplot.png"), plot = plot_avg_trab, width = 2.5, height = 5, units = "in")

# Display the plot
print(plot_avg_trab)

# Visualization: Violin plot with boxplot to visualize genotype differences in total trabeculae area
plot_total_trab <- ggplot(summarised_data, aes(x = genotype, y = total_trab, fill = genotype)) +
  geom_violin(trim = FALSE, alpha = 0.2) +  # Keep orientation default (vertical)
  geom_boxplot(width = 0.1, outlier.size = 1, alpha = 0.5) +
  labs(x = "genotype", y = "total trabecular area",
       title = "Total trabecular area") +
  custom_theme +
  scale_fill_manual(values = c("ctrl" = "dodgerblue", "Bmp2 cko" = "red1", "Bmp7 cko" = "goldenrod1")) +  # Set violin plot colors
  coord_flip()  # Flip coordinates to horizontal orientation

# Save the plot as a PNG file
ggsave(here::here(figs_dir, "02_total_trab_boxplot.png"), plot = plot_total_trab, width = 8, height = 6, units = "in")

# Display the plot
print(plot_total_trab)

# Visualization: Violin plot with boxplot to visualize genotype differences in normalized average trabecular area
plot_normalized_avg_trab <- ggplot(summarised_data, aes(x = genotype, y = normalized_avg_trab, fill = genotype)) +
  geom_violin(trim = FALSE, alpha = 0.2) +  # Keep orientation default (vertical)
  geom_boxplot(width = 0.1, outlier.size = 1, alpha = 0.5) +
  labs(x = "genotype", y = "avg trabecular area/bone area",
       title = "Normalized average trabecular area") +
  custom_theme +
  scale_fill_manual(values = c("ctrl" = "dodgerblue", "Bmp2 cko" = "red1", "Bmp7 cko" = "goldenrod1")) +  # Set violin plot colors
  coord_flip()  # Flip coordinates to horizontal orientation

# Save the plot as a PNG file
ggsave(here::here(figs_dir, "03_normalized_avg_trab_boxplot.png"), plot = plot_normalized_avg_trab, width = 8, height = 6, units = "in")

# Display the plot
print(plot_normalized_avg_trab)

