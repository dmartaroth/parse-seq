# ## ######################################## ## #
#                    FUNCTIONS                   #
# ## ######################################## ## #

# Custom function to format labels in scientific notation if digits are more than 4
format_labels <- function(x) {
  ifelse(abs(x) >= 10000, format(x, scientific = TRUE), format(x))
}



# prepro.plots ------------------------------------------------------------


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
  
  # Function to generate QC violin plots 
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
    
    # Function to generate violin plots with dynamically focused zoom
    generate_violin_plot <- function(data, feature, title) {
      # Main violin plot
      violin_plot <- ggplot(data, aes(x = sample, y = !!sym(feature), fill = sample)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +  # Adjust width of violins and position them closer
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.4, aes(color = sample)) +  # Add dots with jitter, adjust transparency, and set color aesthetic
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = my_colors) +  
        scale_y_continuous(labels = format_labels) +  # Apply custom label formatting for y-axis
        custom_theme_violin() +
        ggtitle(title)
      
      # Dynamically determine zoom range based on data distribution
      density <- density(data[[feature]])
      max_density <- max(density$y)
      density_threshold <- 0.3 * max_density  # 50% of maximum density
      density_peaks <- which(density$y >= density_threshold)  # Find regions with density above threshold
      if (length(density_peaks) > 1) {
        peak_values <- density$x[density_peaks]
        upper_border <- max(peak_values)  # Focus on the highest peak
      } else {
        upper_border <- density$x[density_peaks]
      }
      buffer <- 0.15 * (max(data[[feature]]) - min(data[[feature]]))  # 15% buffer
      zoomed_in_range <- c(0, upper_border + buffer)  # Focus on upper border with a buffer
      
      # Darker colors for points
      darker_colors <- alpha(my_colors, 0.8)  # Adjust alpha for darker points
      
      # Zoomed-in violin plot
      zoomed_in_plot <- ggplot(data, aes(x = sample, y = !!sym(feature), fill = sample)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.6, aes(color = sample)) +  # Adjust point size and alpha for darker points
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = darker_colors) +  # Use darker colors for points
        scale_y_continuous(labels = format_labels, limits = zoomed_in_range) +  # Adjust y-axis limits
        custom_theme_violin() +
        ggtitle(paste("zoom", title))
      
      # Red rounded rectangular outline to indicate zoom region on the original plot
      rect_plot <- geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = min(zoomed_in_range), ymax = max(zoomed_in_range)), 
                             color = "black", fill = NA, linetype = "dashed", size = 0.1, 
                             inherit.aes = FALSE, 
                             lineend = "round", 
                             lwd = 0.4, lty = 2, 
                             dash = c(3, 5))
      
      # Add the rectangle to the main violin plot
      violin_plot <- violin_plot + rect_plot
      
      return(list(violin_plot, zoomed_in_plot))
    }
    
    # Features to plot
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    
    # Create and save violin plots for each feature arranged in a grid
    pdf(file.path(output_dir, sprintf("%02d_b_QC_violin_plots.pdf", plot_number)), width=6, height=3, useDingbats=FALSE)  # Adjust width and height as needed
    for (i in 1:length(features)) {
      plots <- generate_violin_plot(plot_data, features[i], paste(features[i]))
      grid.arrange(plots[[1]], plots[[2]], nrow = 1, ncol = 2)
    }
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

# Add Percentage of Mitochondrial Gene Expression to Seurat Object Metadata
# 
# This function calculates the percentage of mitochondrial gene expression for each cell in a Seurat object and adds it to the metadata. It provides a convenient way to include mitochondrial content information in the Seurat object, which is useful for quality control and downstream analysis. The function takes a Seurat object as input and optionally allows specifying a pattern for mitochondrial genes. After execution, the `percent.mito` column is added to the metadata of the Seurat object, facilitating further analysis and visualization of mitochondrial content.
add_percent_mito <- function(seurat_object, pattern = "^mt-") {
  # Calculate mitochondrial content for the Seurat object
  mt_content <- PercentageFeatureSet(seurat_object, pattern = pattern)
  percent_mito <- mt_content / 100
  
  # Add percent.mito to the metadata of the Seurat object
  seurat_object$percent.mito <- percent_mito
  
  return(seurat_object)
}


# QC plots after filtering ------------------------------------------------------------


filt.plots <- function(data_1, data_2, output_dir) {
  plot_number <- 0  # Starting plot number
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Create subdirectory "01_Preprocessing"
  preprocessing_dir <- file.path(output_dir, "02_QC-after-filtering")
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
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Cell_Counts_filtered.png", plot_number)), width = 3.5, height = 4.5, plot)
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
    
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_UMI_Counts_filtered.png", plot_number)), width = 5, height = 4, plot)
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
      labs(title = "Filtered UMI Counts", fill = "Genotype", color = "Genotype")
    
    p2 <- bind_rows(data_1 = data_1, data_2 = data_2, .id = "Genotype") %>% 
      ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
      geom_boxplot() + 
      custom_theme_density() +
      scale_fill_manual(values = my_colors) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      ggtitle("Filtered NCells vs NGenes")
    
    combined_plot <- plot_grid(p1, p2, rel_widths = c(2, 1.5), scale = c(1, 0.8), labels = "AUTO")
    ggsave(filename = file.path(preprocessing_dir, sprintf("%02d_Genes_Detected_Per_Cell_Filtered.png", plot_number)),width = 7, height = 3, combined_plot)
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
      filename = file.path(preprocessing_dir, sprintf("%02d_Mitochondrial_Content_Filtered.png", plot_number)),
      mitoratio_plot,
      width = 6,
      height = 3
    )
  }
  
  # Function to generate QC violin plots 
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
    
    # Function to generate violin plots with dynamically focused zoom
    generate_violin_plot <- function(data, feature, title) {
      # Main violin plot
      violin_plot <- ggplot(data, aes(x = sample, y = !!sym(feature), fill = sample)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +  # Adjust width of violins and position them closer
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.4, aes(color = sample)) +  # Add dots with jitter, adjust transparency, and set color aesthetic
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = my_colors) +  
        scale_y_continuous(labels = format_labels) +  # Apply custom label formatting for y-axis
        custom_theme_violin() +
        ggtitle(title)
      
      # Dynamically determine zoom range based on data distribution
      density <- density(data[[feature]])
      max_density <- max(density$y)
      density_threshold <- 0.3 * max_density  # 50% of maximum density
      density_peaks <- which(density$y >= density_threshold)  # Find regions with density above threshold
      if (length(density_peaks) > 1) {
        peak_values <- density$x[density_peaks]
        upper_border <- max(peak_values)  # Focus on the highest peak
      } else {
        upper_border <- density$x[density_peaks]
      }
      buffer <- 0.15 * (max(data[[feature]]) - min(data[[feature]]))  # 15% buffer
      zoomed_in_range <- c(0, upper_border + buffer)  # Focus on upper border with a buffer
      
      # Darker colors for points
      darker_colors <- alpha(my_colors, 0.8)  # Adjust alpha for darker points
      
      # Zoomed-in violin plot
      zoomed_in_plot <- ggplot(data, aes(x = sample, y = !!sym(feature), fill = sample)) +
        geom_violin(alpha = 0.8, position = position_dodge(width = 0.6)) +
        geom_point(position = position_jitter(width = 0.3), size = 0.3, alpha = 0.6, aes(color = sample)) +  # Adjust point size and alpha for darker points
        scale_fill_manual(values = my_colors) +
        scale_color_manual(values = darker_colors) +  # Use darker colors for points
        scale_y_continuous(labels = format_labels, limits = zoomed_in_range) +  # Adjust y-axis limits
        custom_theme_violin() +
        ggtitle(paste("zoom", title))
      
      # Red rounded rectangular outline to indicate zoom region on the original plot
      rect_plot <- geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = min(zoomed_in_range), ymax = max(zoomed_in_range)), 
                             color = "black", fill = NA, linetype = "dashed", size = 0.1, 
                             inherit.aes = FALSE, 
                             lineend = "round", 
                             lwd = 0.4, lty = 2, 
                             dash = c(3, 5))
      
      # Add the rectangle to the main violin plot
      violin_plot <- violin_plot + rect_plot
      
      return(list(violin_plot, zoomed_in_plot))
    }
    
    # Features to plot
    features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    
    # Create and save violin plots for each feature arranged in a grid
    pdf(file.path(output_dir, sprintf("%02d_b_QC_violin_plots_filtered.pdf", plot_number)), width=6, height=3, useDingbats=FALSE)  # Adjust width and height as needed
    for (i in 1:length(features)) {
      plots <- generate_violin_plot(plot_data, features[i], paste(features[i]))
      grid.arrange(plots[[1]], plots[[2]], nrow = 1, ncol = 2)
    }
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



# Gene-level filtering ----------------------------------------------------

# The presence of genes with zero counts can dramatically reduce avg expression for a cell; may remove from data with these transformations:
#   
# 1.  Remove genes with zero expression in **all** cells
# 2.  Perform filtering by prevalence, in which genes expressed in only a handful of cells are not particularly meaningful and still bring down averages for all cells not expressed in
# 
# In other words, keep only genes which are expressed in 10 or more cells. This should be adjusted depending on data.

genelvlfilt <- function(object) {
  
  # Output a logical vector for every gene on whether the more than zero counts per cell
  counts <- GetAssayData(object = object, assay = "RNA", layer = "counts") # extract counts
  nonzero <- counts > 0
  
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  
  # Only keeping those genes expressed in more than 10 cells
  filtered_counts <- counts[keep_genes, ]
  
  # Create a new Seurat object with filtered counts
  filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = object@meta.data)
  
  return(filtered_seurat)
}


# Define a function to find top markers for each cluster
findTopMarkers <- function(obj, cluster_names) {
  top_markers <- lapply(cluster_names, function(cluster) {
    response <- FindMarkers(obj, 
                            ident.1 = paste(cluster, "Bmp2_ctrl", sep = "_"), 
                            ident.2 = paste(cluster, "Bmp2_ncko", sep = "_"), 
                            verbose = FALSE)
    top3 <- response %>%
      slice_max(n = 3, order_by = avg_log2FC)
    return(top3)
  })
  names(top_markers) <- cluster_names
  return(top_markers)
}


# Define a function to extract and order top markers
extractTopMarkers <- function(all_top_markers) {
  top_markers <- lapply(all_top_markers, function(markers) {
    markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
    return(markers)
  })
  return(top_markers)
}


# Define a function to collect all top markers
collectAllTopMarkers <- function(all_top_markers) {
  # Extract and order top markers
  ordered_top_markers <- extractTopMarkers(all_top_markers)
  
  # Initialize an empty character vector to store top markers
  all_top_markers_combined <- character()
  
  # Collect top markers for each cluster
  for (cluster_name in names(ordered_top_markers)) {
    cluster_top_markers <- rownames(ordered_top_markers[[cluster_name]])
    all_top_markers_combined <- c(all_top_markers_combined, cluster_top_markers)
  }
  
  # Remove duplicates
  all_top_markers_combined <- unique(all_top_markers_combined)
  
  return(all_top_markers_combined)
}

# Define a function to collect top positive markers
collectTopPosMarkers <- function(all_top_pos_markers) {
  # Extract and order top markers
  ordered_top_pos_markers <- extractTopPosMarkers(all_top_pos_markers)
  
  # Initialize an empty character vector to store top markers
  all_top_pos_markers_combined <- character()
  
  # Collect top markers for each cluster
  for (cluster_name in names(ordered_top_pos_markers)) {
    cluster_top_pos_markers <- rownames(ordered_top_pos_markers[[cluster_name]])
    all_top_pos_markers_combined <- c(all_top_pos_markers_combined, cluster_top_pos_markers)
  }
  
  # Remove duplicates
  all_top_pos_markers_combined <- unique(all_top_pos_markers_combined)
  
  return(all_top_pos_markers_combined)
}  

# Define a function to find top positive markers for each cluster 
findTopPosMarkers <- function(obj, cluster_names) {
  top_pos_markers <- lapply(cluster_names, function(cluster) {
    response <- FindMarkers(obj, 
                            ident.1 = paste(cluster, "Bmp2_ctrl", sep = "_"), 
                            ident.2 = paste(cluster, "Bmp2_ncko", sep = "_"),
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
    top3 <- response %>%
      slice_max(n = 3, order_by = avg_log2FC)
    return(top3)
  })
  names(top_pos_markers) <- cluster_names
  return(top_pos_markers)
}


# Define a function to extract and order top positive markers
extractTopPosMarkers <- function(all_top_pos_markers) {
  top_pos_markers <- lapply(all_top_pos_markers, function(markers) {
    markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
    return(markers)
  })
  return(top_pos_markers)
}


# Multi feature plot function --------------------------------------------

multi_feature_plot <- function(seurat_object, features, reduction = "umap", na_cutoff = 0, colors_use) {
  FeaturePlot_scCustom(
    seurat_object = seurat_object,
    reduction = reduction,
    split.by = "genotype",
    na_cutoff = na_cutoff,
    features = features,
    colors_use = colors_use
  )
}

# Convenience function to create and save multi_feature_plot
convenient_multi_feature_plot <- function(seurat_object = obj, features, colors_use, name, number = plot_number, width = 6, height = 9, dir = comparison_dir) {
  
  # Create the multi-feature plot
  plot <- multi_feature_plot(seurat_object = seurat_object, features = features, colors_use = colors_use)
  
  # Increment plot number
  number <- number + 1
  
  # Save the plot
  ggsave(filename = file.path(dir, sprintf("%02d_%s_by_genotype_featureplots.png", number, name)), 
         width = width, height = height, plot)
  
  # Update plot_number in the global environment
  assign("plot_number", number, envir = .GlobalEnv)
  
  return(plot)
}


