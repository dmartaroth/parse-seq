# ## ######################################## ## #
#                 DATA VISUALIZATION             #
# ## ######################################## ## #
# Date: Mon Jun 03 16:05:52 2024 ------------------
# Updated by: Daniela M. Roth

# Data visualization for semi-quantitative measures (not scRNA-seq)

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
  ggsave(here(ct.scoring, "01a_frequency-of-ISS-SOS_scoring_by-sex.pdf"), width = 7.2, height = 8.1)
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
        axis.text = element_text(size = 5),
        strip.text = element_text(color = "floralwhite", face = "bold"),
        text = element_text(size = 6)) +  # floralwhite bold text for facet wrap labels
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
        axis.text = element_text(size = 5),
        strip.text = element_text(color = "floralwhite", face = "bold"),
        text = element_text(size = 6)) +  # floralwhite bold text for facet wrap labels
  geom_text(stat='count', aes(label=after_stat(count)), position=position_fill(vjust=0.5), color="navy", size = 2)


# Arrange plots horizontally
plot_grid(plot_iss_no_sex, plot_sos_no_sex, ncol = 2)

# Save the plot
ggsave(here(ct.scoring, "01b_frequency-of-ISS-SOS_scoring_no_sex_split.pdf"), width = 5, height = 4)
