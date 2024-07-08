# ## ######################################## ## #
#         SCORING SPHENOID MALFORMATION          #
# ## ######################################## ## #

# Date: Sat Jul 06 21:46:04 2024 ------------------
# Updated by: Daniela M. Roth

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(here)
library(umap)

# Load the data from Excel
data <- read_excel(here::here("postnatal_cb/raw-data/embryonic_sphenoid_scoring.xlsx"))

# Convert relevant columns to factors if they are not already factors
data <- mutate(data,
               genotype = factor(genotype),
               missing_psph = factor(missing_psph),
               cleft = factor(cleft),
               frag_basisph = factor(frag_basisph),
               basisph_notch = factor(basisph_notch))

# Extract binary phenotype columns using dplyr::select
binary_data <- data %>%
  dplyr::select(missing_psph, cleft, frag_basisph, basisph_notch)

# Convert binary variables to numeric and convert list back to data frame
binary_data <- as.data.frame(lapply(binary_data, as.numeric))

# Function to calculate unique combination of binary phenotypes
get_combinations <- function(row) {
  paste(row, collapse = "")
}

# Apply function to each row to get phenotype combinations
phenotype_combinations <- apply(binary_data, 1, get_combinations)

# Adjust UMAP parameters for better separation
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.15

# Apply UMAP for dimensionality reduction with custom parameters
umap_result <- umap(binary_data, config = umap_config)

# Create dataframe for plotting
plot_data <- data.frame(
  genotype = data$genotype,
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  combination = phenotype_combinations
)

# Define custom colors
custom_colors <- c(
  "Bmp7 cko" = "black",
  "ctrl" = "dodgerblue",
  "Bmp2 ncko" = "goldenrod1",
  "Bmp7 ko" = "mediumseagreen",
  "Bmp2hetBmp7ncko" = "pink",
  "Bmp2hetBmp7hetncko" = "powderblue",
  "Bmp2Bmp7ncko" = "red1",
  "Bmp2flBmp7hetncko" = "olivedrab2",
  "Bmp7 ncko" = "violet"
)

# Plotting UMAP-like plot with ellipses and custom colors
plot <- ggplot(plot_data, aes(UMAP1, UMAP2, color = genotype)) +
  geom_point(size = 3) +  # Adjust point size as needed
  stat_ellipse(aes(fill = genotype), geom = "polygon", alpha = 0.15, color = NA) +  # Ellipses for clusters
  scale_color_manual(name = "Genotype", values = custom_colors) +
  scale_fill_manual(name = "Genotype", values = custom_colors) +
  labs(title = "Sphenoid observations") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = "white")
  )

ggsave(
  here("postnatal_cb/figs/UMAP_sphenoid_observations.png"),
  plot = plot,
  width = 9,
  height = 4.5
)

# Save the statistical results to a file
sink(here("postnatal_cb/figs/sphenoid_stats.txt"))

# Perform Kruskal-Wallis test for UMAP1 by genotype
kruskal_result_UMAP1 <- kruskal.test(UMAP1 ~ genotype, data = plot_data)
cat("Kruskal-Wallis Test Results for UMAP1:\n")
print(kruskal_result_UMAP1)

# Perform Kruskal-Wallis test for UMAP2 by genotype
kruskal_result_UMAP2 <- kruskal.test(UMAP2 ~ genotype, data = plot_data)
cat("Kruskal-Wallis Test Results for UMAP2:\n")
print(kruskal_result_UMAP2)

# Pairwise Wilcoxon rank-sum tests for UMAP1
pairwise_results_UMAP1 <- pairwise.wilcox.test(plot_data$UMAP1, plot_data$genotype, p.adjust.method = "bonferroni")
cat("Pairwise Wilcoxon Test Results for UMAP1:\n")
print(pairwise_results_UMAP1)

# Pairwise Wilcoxon rank-sum tests for UMAP2
pairwise_results_UMAP2 <- pairwise.wilcox.test(plot_data$UMAP2, plot_data$genotype, p.adjust.method = "bonferroni")
cat("Pairwise Wilcoxon Test Results for UMAP2:\n")
print(pairwise_results_UMAP2)

sink()



# Count the number of unique biological replicates per genotype
biological_replicates <- data %>%
  group_by(genotype) %>%
  summarise(count = n_distinct(ID))

# Print the number of biological replicates per genotype
print(biological_replicates)

library(dplyr)
library(readr)
library(here)

# Load the CSV file
file_path <- here("postnatal_cb/figs/phenotype_stats.csv")
phenotype_stats <- read_csv(file_path)

# Calculate the percentages
phenotype_stats <- phenotype_stats %>%
  mutate(
    missing_psph_percent = (missing_psph_count / count) * 100,
    cleft_percent = (cleft_count / count) * 100,
    frag_basisph_percent = (frag_basisph_count / count) * 100,
    basisph_notch_percent = (basisph_notch_count / count) * 100
  )

# Save the results
write_csv(phenotype_stats, here("postnatal_cb/figs/phenotype_stats_with_percentages.csv"))

# Print the results to check
print(phenotype_stats)
