# Load the necessary packages
library(ggplot2)
library(gridExtra)
library(dplyr)

# Custom colors
my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A")

# Custom theme for barplot
custom_theme_bar <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(color = "white", fill = "white"),
      panel.grid.major = element_blank(),           # No major gridlines
      panel.grid.minor = element_blank(),           # No minor gridlines
      legend.position = "none",                     # No legend
      axis.ticks = element_line(color = "black"),
      axis.ticks.x = element_line(color = NA),
      axis.title = element_text(color = "black",size = 8),
      axis.text = element_text(size = 8),           # Size of axis text
      axis.line.y = element_line(colour = NA), # Set y-axis line color
      axis.line.x = element_line(colour = NA),  
      axis.text.y = element_text(colour = "black"), # Set y-axis text color
      axis.ticks.length = unit(0.2, "cm"),          # Shorten tick length
      panel.border = element_blank(),   
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),         # No vertical gridlines
      text = element_text(colour = "black", size = 8),                        # Regular text
      plot.title = element_text(face = "bold", size = 8,hjust = 0),     # Bold plot title
      # plot.margin = margin(20, 20, 20, 20),         # Adjust plot margins
      plot.caption = element_text(color = "black")  # Caption color and position
    )
}

# Custom theme for density plot
custom_theme_density <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),                        # No vertical gridlines
      
      # Axis
      axis.ticks = element_line(color = "black"),                  # Color of axis ticks
      axis.line.y = element_line(colour = "black"),                # Color of y-axis line
      axis.line.x = element_line(colour = NA),                     # Remove x-axis line
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(size = 10),                         # Size of axis text
      axis.text.y = element_text(colour = "black"),                # Color of y-axis text
      axis.ticks.length = unit(0.2, "cm"),                         # Shorten tick length
      
      # Panel
      panel.border = element_rect(color = "black", fill = NA),     # Panel border
      
      # Text
      plot.title = element_text(face = "bold", vjust = 0.5),       # Bold plot title, centered vertically
      text = element_text(),                                       # Regular text
      
      # Legend and Caption
      legend.position = "none",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0)      # Caption color and position
    )
}

# Functions
#-------------------------------
prepro.plots <- function(data_1, data_2, output_dir) {
  plot_number <- 0  # Starting plot number
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create subdirectory "01_Preprocessing"
  preprocessing_dir <- file.path(output_dir, "01_Preprocessing")
  if (!dir.exists(preprocessing_dir)) {
    dir.create(preprocessing_dir)
  }
  
  # Plot Cell Counts
  plot_cell_counts <- function(data_1, data_2, plot_number) {
    data_1 <- as.data.frame(data_1@meta.data)
    data_2 <- as.data.frame(data_2@meta.data)
    plot <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>%
      ggplot(aes(x = sample, fill = sample)) +
      geom_bar(colour = "black", position = "dodge") +
      labs(title = "Cell Counts")+
      custom_theme_bar() +
      facet_wrap(~ Genotype, scales = "free_x") +
      scale_fill_manual(values = my_colors)
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Cell_Counts.png", plot_number)), width = 3.5, height = 4.5, plot)
  }
  
  # Plot UMI Counts
  plot_umi_counts <- function(data_1, data_2, plot_number) {
    data_1 <- as.data.frame(data_1@meta.data)
    data_2 <- as.data.frame(data_2@meta.data)
    plot <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>%
      ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      labs(title = "UMI Counts")
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_UMI_Counts.png", plot_number)), plot)
  }
  
  # Plot Genes Detected per Cell
  plot_genes_detected_per_cell <- function(data_1, data_2, plot_number) {
    data_1 <- as.data.frame(data_1@meta.data)
    data_2 <- as.data.frame(data_2@meta.data)
    p1 <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>%
      ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      geom_vline(xintercept = 300) +
      custom_theme_density() +
      labs(title = "Density Plot")
    
    p2 <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>% 
      ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
      geom_boxplot() + 
      custom_theme_density() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ggtitle("NCells vs NGenes")
    
    combined_plot <- plot_grid(p1, p2, rel_widths = c(2, 1.5), scale = c(1, 0.8), labels = "AUTO")
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Genes_Detected_Per_Cell.png", plot_number)), combined_plot)
  }
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_cell_counts(data_1, data_2, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_umi_counts(data_1, data_2, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_genes_detected_per_cell(data_1, data_2, plot_number)
  
  # Print message
  cat("Plots saved in", preprocessing_dir, "\n")
}

# Example usage:
# prepro.plots(object_1, object_2, "output_directory")