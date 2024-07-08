# ## ######################################## ## #
#                 DATA VISUALIZATION             #
# ## ######################################## ## #
# Date: Mon Jun 03 16:05:52 2024 ------------------
# Updated by: Daniela M. Roth



# Synchondrosis fusion scoring --------------------------------------------

# Measurements are collected in excel spreadsheet titled "cb_scoring"
library(here)
dir.create(here(ct.scoring <- here("postnatal_cb","figs")), recursive = TRUE) 

source(here("postnatal_cb/docs/packages.R"))

# Read the data
df <- read.csv(here("postnatal_cb/raw-data/cb_scoring.csv"))

# Remove rows with any missing values
df_clean <- na.omit(df)

# Ensure relevant columns are converted to factors
df_clean$genotype <- factor(df_clean$genotype, levels = c("Bmp2 ctrl", "Bmp2 cko", "Bmp7 ctrl", "Bmp7 cko"))
df_clean$ISS <- factor(df_clean$ISS, levels = rev(c(0, 1, 2, 3)))
df_clean$SOS <- factor(df_clean$SOS, levels = rev(c(0, 1, 2, 3)))
df_clean$sex <- factor(df_clean$sex, levels = c("F", "M"))

# Combine "Bmp2 ctrl" and "Bmp7 ctrl" into a single category "ctrl"
df_clean <- df_clean %>%
  mutate(genotype_grouped = case_when(
    genotype == "Bmp2 ctrl" | genotype == "Bmp7 ctrl" ~ "ctrl",
    TRUE ~ as.character(genotype)
  )) %>%
  mutate(genotype_grouped = factor(genotype_grouped, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")))

# Define custom colors for cb scores
score_colors <- c("firebrick3", "lightsalmon", "floralwhite")

# Ensure there are rows left in the dataframe
if (nrow(df_clean) > 0) {
  # Filter data to exclude PT10 age category
  df_clean_filtered <- df_clean %>%
    filter(age != "PT10")
  
  # Stacked bar plot for genotype with SOS categories
  plot_sos <- ggplot(df_clean_filtered, aes(x = genotype_grouped, fill = SOS)) +
    geom_bar(position = position_fill(), stat = "count") +
    scale_fill_manual(values = score_colors) +
    facet_grid(age ~ ifelse(sex == "F", "female", "male"), scales = "free_y") +  # labeling sexes as female and male
    labs(x = "Genotype", y = "Frequency") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), labels = scales::percent(seq(0, 1, 0.5))) +  # Custom breaks and labels for y-axis
    geom_hline(yintercept = 0, color = "navy") +  # navy line at 0%
    theme(panel.background = element_rect(fill = "white"),  # white background
          panel.grid.major = element_line(color = "#F0F8FF"),  # light blue grid lines
          strip.background = element_rect(fill = "navy"),  # navy background for facet wrap labels
          strip.text = element_text(color = "floralwhite", face = "bold")) +  # floralwhite bold text for facet wrap labels
    geom_text(stat='count', aes(label=after_stat(count)), position=position_fill(vjust=0.5), color="black", size = 3)
  
  # Stacked bar plot for genotype with ISS categories
  plot_iss <- ggplot(df_clean_filtered, aes(x = genotype_grouped, fill = ISS)) +
    geom_bar(position = position_fill(), stat = "count") +
    scale_fill_manual(values = score_colors) +
    facet_grid(age ~ ifelse(sex == "F", "female", "male"), scales = "free_y") +  # labeling sexes as female and male
    labs(title = "Cranial base fusion scoring", x = "", y = "Frequency") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), labels = scales::percent(seq(0, 1, 0.5))) +  # Custom breaks and labels for y-axis
    geom_hline(yintercept = 0, color = "navy") +  # navy line at 0%
    theme(panel.background = element_rect(fill = "white"),  # white background
          panel.grid.major = element_line(color = "#F0F8FF"),  # light blue grid lines
          strip.background = element_rect(fill = "navy"),  # navy background for facet wrap labels
          strip.text = element_text(color = "floralwhite", face = "bold")) +  # floralwhite bold text for facet wrap labels
    geom_text(stat='count', aes(label=after_stat(count)), position=position_fill(vjust=0.5), color="black", size = 3)
  
  # Arrange plots vertically
  plot_grid(plot_iss, plot_sos, ncol = 1)
  
  # Save the plot
  ggsave(here(ct.scoring, "01a_frequency-of-ISS-SOS_scoring_by-sex.png"), width = 5.5, height = 4.5)
} else {
  print("No data available after removing rows with missing values.")
}

# Filter data to exclude PT10 age category for plots without splitting by sex
df_clean_filtered_no_sex <- df_clean %>%
  filter(age != "PT10")
# Stacked bar plot for genotype with SOS categories (without splitting by sex)
plot_sos_no_sex <- ggplot(df_clean_filtered_no_sex, aes(x = genotype_grouped, fill = SOS)) +
  geom_bar(position = position_fill(), stat = "count") +
  scale_fill_manual(values = score_colors) +
  facet_grid(age ~ ., scales = "free_y") +  # no splitting by sex
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = scales::percent(seq(0, 1, 0.5))) +  # Custom breaks and labels for y-axis
  geom_hline(yintercept = 0, color = "navy") +  # navy line at 0%
  theme(panel.background = element_rect(fill = "white"),  # white background
        panel.grid.minor = element_line(color = "#F0F8FF"),  # light blue grid lines
        panel.grid.major = element_line(color = "#F0F8FF"),  # light blue grid lines
        strip.background = element_rect(fill = "navy"),  # navy background for facet wrap labels
        axis.text = element_text(size = 7),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(color = "floralwhite", face = "bold"),
        text = element_text(size = 7)) +  # floralwhite bold text for facet wrap labels
  geom_text(stat='count', aes(label=after_stat(count)), position=position_fill(vjust=0.5), color="navy", size = 2)

# Stacked bar plot for genotype with ISS categories (without splitting by sex)
plot_iss_no_sex <- ggplot(df_clean_filtered_no_sex, aes(x = genotype_grouped, fill = ISS)) +
  geom_bar(position = position_fill(), stat = "count") +
  scale_fill_manual(values = score_colors) +
  facet_grid(age ~ ., scales = "free_y") +  # no splitting by sex
  labs(title = "Cranial base fusion scoring",x = "Genotype", y = "Frequency") +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = scales::percent(seq(0, 1, 0.5))) +  # Custom breaks and labels for y-axis
  geom_hline(yintercept = 0, color = "navy") +  # navy line at 0%
  theme(panel.background = element_rect(fill = "white"),  # white background
        panel.grid.minor = element_line(color = "#F0F8FF"),  # light blue grid lines
        panel.grid.major = element_line(color = "#F0F8FF"),  # light blue grid lines
        strip.background = element_rect(fill = "navy"),  # navy background for facet wrap labels
        axis.text = element_text(size = 7),
        strip.text = element_text(color = "floralwhite", face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        text = element_text(size = 7)) +  # floralwhite bold text for facet wrap labels
  geom_text(stat='count', aes(label=after_stat(count)), position=position_fill(vjust=0.5), color="navy", size = 2)


# Arrange plots horizontally
plot_grid(plot_iss_no_sex, plot_sos_no_sex, ncol = 2)

# Save the plot
ggsave(here(ct.scoring, "01b_frequency-of-ISS-SOS_scoring_no_sex_split.pdf"), width = 5, height = 4)

ggsave(here(ct.scoring, "01b_frequency-of-ISS-SOS_scoring_no_sex_split.png"), width = 4.5, height = 7)


# Statistical analysis ----------------------------------------------------

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(here)
library(ggpubr)

# Read the data
df <- read.csv(here("postnatal_cb/raw-data/cb_scoring.csv"))

# Ensure the age column is correctly interpreted and group PT5 and PT10 as PT
df$age <- factor(df$age, levels = c("PT5", "PT10", "AT"))
df <- df %>%
  mutate(age_grouped = ifelse(age %in% c("PT5", "PT10"), "PT", as.character(age)))

# Remove rows with any missing values in ISS or SOS
df_clean <- df %>%
  filter(!is.na(ISS) & ISS != "" & !is.na(SOS) & SOS != "")

# Combine "Bmp2 ctrl" and "Bmp7 ctrl" into a single category "ctrl"
df_clean <- df_clean %>%
  mutate(genotype_grouped = case_when(
    genotype == "Bmp2 ctrl" | genotype == "Bmp7 ctrl" ~ "ctrl",
    TRUE ~ as.character(genotype)
  )) %>%
  mutate(genotype_grouped = factor(genotype_grouped, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")))

# Check the unique values in the age column
cat("Unique values in the age_grouped column after cleaning:\n")
print(unique(df_clean$age_grouped))

# Filter data to exclude rows with any missing values for age_grouped
df_clean_filtered_no_sex <- df_clean %>%
  filter(!is.na(age_grouped) & age_grouped != "")

# Check the number of rows with age "PT" and "AT"
cat("Number of rows with age 'PT':", nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "PT")), "\n")
cat("Number of rows with age 'AT':", nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "AT")), "\n")

# Perform Fisher's exact test for ISS stratified by age
if (nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "PT")) > 0) {
  fisher_test_iss_pt <- fisher.test(table(df_clean_filtered_no_sex %>% filter(age_grouped == "PT") %>% select(genotype_grouped, ISS)))
  cat("Fisher's Exact Test Results for ISS (PT):\n")
  print(fisher_test_iss_pt)
} else {
  cat("No data available for ISS (PT)\n")
}

# Perform Fisher's exact test for SOS stratified by age
if (nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "PT")) > 0) {
  fisher_test_sos_pt <- fisher.test(table(df_clean_filtered_no_sex %>% filter(age_grouped == "PT") %>% select(genotype_grouped, SOS)))
  cat("Fisher's Exact Test Results for SOS (PT):\n")
  print(fisher_test_sos_pt)
} else {
  cat("No data available for SOS (PT)\n")
}

# Perform Fisher's exact test for ISS stratified by age for AT
if (nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "AT")) > 0) {
  fisher_test_iss_at <- fisher.test(table(df_clean_filtered_no_sex %>% filter(age_grouped == "AT") %>% select(genotype_grouped, ISS)))
  cat("Fisher's Exact Test Results for ISS (AT):\n")
  print(fisher_test_iss_at)
} else {
  cat("No data available for ISS (AT)\n")
}

# Perform Fisher's exact test for SOS stratified by age for AT
if (nrow(df_clean_filtered_no_sex %>% filter(age_grouped == "AT")) > 0) {
  fisher_test_sos_at <- fisher.test(table(df_clean_filtered_no_sex %>% filter(age_grouped == "AT") %>% select(genotype_grouped, SOS)))
  cat("Fisher's Exact Test Results for SOS (AT):\n")
  print(fisher_test_sos_at)
} else {
  cat("No data available for SOS (AT)\n")
}

# Count the number of samples per genotype and age group
sample_sizes <- df_clean_filtered_no_sex %>%
  group_by(genotype_grouped, age_grouped) %>%
  summarise(sample_size = n(), .groups = 'drop')

# Print the sample sizes
print(sample_sizes)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(here)
library(ggpubr)

# Set up the sink to capture output
sink(here("postnatal_cb", "figs", "freq_ISS-SOS_sex_statistical_analysis_results.txt"))

# Read the data
df <- read.csv(here("postnatal_cb", "raw-data", "cb_scoring.csv"))

# Clean and prepare data
df_clean <- df %>%
  filter(!is.na(ISS), ISS != "", !is.na(SOS), SOS != "") %>%
  mutate(
    age = factor(age, levels = c("PT5", "PT10", "AT")),
    age_grouped = ifelse(age %in% c("PT5", "PT10"), "PT", as.character(age)),
    genotype_grouped = case_when(
      genotype == "Bmp2 ctrl" | genotype == "Bmp7 ctrl" ~ "ctrl",
      TRUE ~ as.character(genotype)
    ),
    genotype_grouped = factor(genotype_grouped, levels = c("ctrl", "Bmp2 cko", "Bmp7 cko")),
    sex = factor(sex, levels = c("F", "M"))
  ) %>%
  filter(!is.na(age_grouped) & age_grouped != "")

# Display sample sizes per group to verify sufficient data
cat("Data count by genotype, age, and sex:\n")
df_clean %>%
  group_by(genotype_grouped, age_grouped, sex) %>%
  summarise(count = n(), .groups = 'drop') %>%
  print()

# Function to perform the tests
perform_test <- function(data, score) {
  table_score <- table(data$sex, data[[score]])
  if (any(table_score < 5)) {
    test_result <- fisher.test(table_score)
  } else {
    test_result <- chisq.test(table_score)
  }
  return(test_result)
}

# Perform statistical analysis for each genotype and age group
results_list <- list()
for (genotype in unique(df_clean$genotype_grouped)) {
  for (age in unique(df_clean$age_grouped)) {
    data_subset <- df_clean %>%
      filter(genotype_grouped == genotype, age_grouped == age)
    if (nrow(data_subset) > 1 && length(unique(data_subset$sex)) > 1) {
      results_list[[paste(genotype, age, "ISS", sep = "_")]] <- perform_test(data_subset, "ISS")
      results_list[[paste(genotype, age, "SOS", sep = "_")]] <- perform_test(data_subset, "SOS")
    } else {
      results_list[[paste(genotype, age, "ISS", sep = "_")]] <- NA
      results_list[[paste(genotype, age, "SOS", sep = "_")]] <- NA
    }
  }
}

# Print results
print("Statistical test results by genotype and age group:")
print(results_list)

# Close the sink to stop diverting output to the file
sink()
