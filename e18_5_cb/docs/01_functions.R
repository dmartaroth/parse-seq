# Load the necessary packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(crayon)

# Custom colors
my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A")
error <- red $ bold
note1 <- red $ bold
note2 <- black $ bold
note3 <- blue $ bold

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
      panel.border = element_rect(color = NA, fill = NA),     # Panel border
      
      # Text
      plot.title = element_text(face = "bold", vjust = 0.5),       # Bold plot title, centered vertically
      text = element_text(),                                       # Regular text
      
      # Legend and Caption
      legend.position = "right",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0)      # Caption color and position
    )
}

# Custom theme for scatter plots
custom_theme_scatter <- function() {
  theme_minimal() +
    theme(
      # Background and Gridlines
      plot.background = element_rect(fill = "white", color = NA),  # White background
      panel.background = element_rect(fill = "white"),            # Panel background
      panel.grid.major = element_blank(),                          # No major gridlines
      panel.grid.minor = element_blank(),                          # No minor gridlines
      
      # Axis
      axis.line = element_line(color = "black"),                   # Color of axis lines
      axis.title = element_text(color = "black"),                  # Color of axis titles
      axis.text = element_text(color = "black", size = 8),         # Size and color of axis text
      
      # Legend and Caption
      legend.position = "none",                                    # No legend
      plot.caption = element_text(color = "black", hjust = 0.5),    # Caption color and position
      
      # Remove axis ticks
      axis.ticks = element_blank()                                 # No axis ticks
    )
}

# Custom violin plot theme
custom_theme_violin <- function() {
  theme_minimal() +
    theme(
      legend.position = "none",  # Remove legend
      panel.background = element_rect(fill = "white"),  # White background for panels
      panel.grid.major.y = element_line(colour = "grey98"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),  # Color of axis lines
      axis.title = element_text(color = "black"),  # Color of axis titles
      axis.text = element_text(color = "black"),  # Color of axis text
      axis.ticks = element_line(color = "black"),  # Color of axis ticks
      text = element_text(color = "black"),  # Color of text
      plot.title = element_text(face = "bold", hjust = 0.5),  # Bold plot title, centered
      plot.caption = element_text(hjust = 0)  # Caption position
    )
}
# Functions
#-------------------------------

# Custom function to format labels in scientific notation if digits are more than 4
format_labels <- function(x) {
  ifelse(abs(x) >= 10000, format(x, scientific = TRUE), format(x))
}

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
    
    cat(note1("Cell counts represents the number of unique cellular barcodes detected.\n"))
    
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
    
    cat(note2("UMI counts per cell should be above 500. If between 500 and 1000 counts, the cells may have benefited from deeper sequencing. The x intercept of this plot is set at 500.\n"))
    
    plot <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>%
      ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "UMI Counts", fill = "Genotype", color = "Genotype")
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_UMI_Counts.png", plot_number)), width = 5, height = 4, plot)
  }
  
  
  # Plot Genes Detected per Cell
  plot_genes_detected_per_cell <- function(data_1, data_2, plot_number) {
    data_1 <- as.data.frame(data_1@meta.data)
    data_2 <- as.data.frame(data_2@meta.data)
    
    p1 <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>%
      ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident)) + 
      geom_density(alpha = 0.7) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500) +
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(title = "UMI Counts", fill = "Genotype", color = "Genotype")
    
    p2 <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>% 
      ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
      geom_boxplot() + 
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ggtitle("NCells vs NGenes")
    
    combined_plot <- plot_grid(p1, p2, rel_widths = c(2, 1.5), scale = c(1, 0.8), labels = "AUTO")
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Genes_Detected_Per_Cell.png", plot_number)),width = 7, height = 3, combined_plot)
  }
  
  # Plot Mitochondrial Content
  plot_mitochondrial_content <- function(data_1, data_2, plot_number) {
    
    # Calculate mitochondrial content for data_1
    mt_content_1 <- PercentageFeatureSet(data_1, pattern = "^mt-")
    mt_ratio_1 <- mt_content_1 / 100
    
    # Calculate mitochondrial content for data_2
    mt_content_2 <- PercentageFeatureSet(data_2, pattern = "^mt-")
    mt_ratio_2 <- mt_content_2 / 100
    
    # Add mitochondrial content to data frames
    data_1$percent.mt <- mt_ratio_1
    data_2$percent.mt <- mt_ratio_2
    
    # Add mitoRatio to meta.data
    data_1$mitoRatio <- mt_ratio_1
    data_2$mitoRatio <- mt_ratio_2
    
    # Plot mitochondrial content
    mitoratio_plot <- bind_rows(data_1 = data_1@meta.data, data_2 = data_2@meta.data, .id = "Genotype") %>%
      ggplot(aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point(aes(color=mitoRatio)) +
      scale_colour_gradient(low = "gray90", high = "black") +
      geom_smooth(method=lm) +   
      scale_x_log10() +
      scale_y_log10() +
      theme_classic() +
      geom_vline(xintercept = 250) +
      geom_hline(yintercept = 250) +
      facet_wrap(~sample)
    
    # Save plot
    ggsave(
      filename = file.path(preprocessing_dir, sprintf("%02d_Mitochondrial_Content.png", plot_number)),
      mitoratio_plot,
      width = 6,
      height = 3
    )
  }
  
  # Function to create QC violin plots 
  generate_QC_violin_plots <- function(output_dir, data_1, data_2, plot_number) {
    cat((green("\nPlotting QC Violin plots \n")))
    
    # Calculate mitochondrial content for data_1
    mt_content_1 <- PercentageFeatureSet(data_1, pattern = "^mt-")
    mt_ratio_1 <- mt_content_1 / 100
    
    # Calculate mitochondrial content for data_2
    mt_content_2 <- PercentageFeatureSet(data_2, pattern = "^mt-")
    mt_ratio_2 <- mt_content_2 / 100
    
    # Add mitochondrial content to data frames
    data_1$percent.mt <- mt_ratio_1
    data_2$percent.mt <- mt_ratio_2
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # Combine data frames for plotting
    plot_data <- bind_rows(
      data_1 = data_1@meta.data,
      data_2 = data_2@meta.data,
      .id = "sample"
    )
    
    # Function to generate violin plots
    # Function to generate violin plots with dots and without legend
    library(gridExtra)
    
    # Function to generate violin plots with dots and without legend
    generate_violin_plot <- function(data, feature, title) {
      ggplot(data, aes(x = sample, y = !!sym(feature), fill = sample)) +
        geom_violin(alpha = 0.8) +
        geom_point(position = position_jitter(width = 0.2), size = 0.1, alpha = 0.4, aes(color = sample)) +  # Add dots with jitter, adjust transparency, and set color aesthetic
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = my_colors) +  
        scale_y_continuous(labels = format_labels) +  # Apply custom label formatting for y-axis
        custom_theme_violin() +
        ggtitle(title)
    }
    
    # Features to plot
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    
    # Create and save violin plots for each feature arranged in a grid
    pdf(file.path(output_dir, sprintf("%02d_b_QC_violin_plots.pdf", plot_number)), width=10, height=4, useDingbats=FALSE)  # Adjust width and height as needed
    grid.arrange(
      generate_violin_plot(plot_data, features[1], paste(features[1])),
      generate_violin_plot(plot_data, features[2], paste(features[2])),
      generate_violin_plot(plot_data, features[3], paste(features[3])),
      nrow = 1, ncol = 3
    )
    invisible(dev.off())  # Close the PDF device
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
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  plot_mitochondrial_content(data_1, data_2, plot_number)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  suppressWarnings(generate_QC_violin_plots(preprocessing_dir, data_1, data_2, plot_number))
  
  
  # Print message
  cat("Plots saved in", preprocessing_dir, "\n")
}
