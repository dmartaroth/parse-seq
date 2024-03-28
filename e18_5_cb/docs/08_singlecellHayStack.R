# ## ######################################## ## #
#                     schaystack                 #
# ## ######################################## ## #

# Date: Tue Mar 26 15:45:27 2024 ------------------

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

install.packages("singleCellHaystack")
library(singleCellHaystack)


# Create subdirectory "09_singleCellHaystack"
schaystack.dir <- file.path(figs, "09_singleCellHaystack")
if (!dir.exists(schaystack.dir)) {
  dir.create(schaystack.dir)
}

plot_number <- 0

set.seed(123)


obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))
obj$cell_types <- Idents(obj)

res.pc20 <- haystack(x = obj, coord = "integrated.cca",assay = "RNA", slot = "data")

show_result_haystack(res.haystack = res.pc20, n = 10)

p <- SCpubr::do_FeaturePlot(sample = obj,
                       features = "Ccdc187", split.by = "genotype")


convenient_save_plot(p, name = "Ccdc187_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = obj,
                            features = "Csf2rb", split.by = "genotype")


convenient_save_plot(p, name = "Csf2rb_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = obj,
                            features = "Dcx", split.by = "genotype")


convenient_save_plot(p, name = "Dcx_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = obj,
                            features = "Hecw1", split.by = "genotype")


convenient_save_plot(p, name = "Hecw1_featureplot", dir = schaystack.dir,height = 6, width = 8)



# Apply singlecellHayStack to osteochondro cluster only -------------------

# Subset chondro clusters
osteochondro.clusters <- c("chondro.1","chondro.2","chondro.3","chondro.4","mes.1","mes.2","mes.3")
osteochondro_subset <- subset(obj, idents = osteochondro.clusters)

(plot <- DimPlot(osteochondro_subset, reduction = "umap", label = FALSE,split.by = "genotype",cols = osteochondro.colors))


res.pc20 <- haystack(x = osteochondro_subset, coord = "integrated.cca",assay = "RNA", slot = "data")

results <- show_result_haystack(res.haystack = res.pc20, n = 50)
write.csv(results, file = here(schaystack.dir,"haystack_result.csv"))

p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Col1a1", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)

convenient_save_plot(p, name = "Col1a1_featureplot", dir = schaystack.dir,height = 6, width = 8)

p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Ptn", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Ptn_featureplot", dir = schaystack.dir,height = 6, width = 8)

p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Cdh11", split.by = "genotype",use_viridis = TRUE
                       ,viridis.palette = "cividis", viridis.direction = 1,
                       na.value = "floralwhite", plot_cell_borders = FALSE, 
                       order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Cdh11_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Col3a1", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Col3a1_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Col5a2", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Col5a2_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Ptprd", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Ptprd_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Cped1", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Cped1_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Postn", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Postn_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Cacna1g", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Cacna1g_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Fbn2", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Fbn2_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Aspn", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Aspn_featureplot", dir = schaystack.dir,height = 6, width = 8)


p <- SCpubr::do_FeaturePlot(sample = osteochondro_subset,
                            features = "Thbs2", split.by = "genotype",use_viridis = TRUE
                            ,viridis.palette = "cividis", viridis.direction = 1,
                            na.value = "floralwhite", plot_cell_borders = FALSE, 
                            order = TRUE,pt.size = 1.5)


convenient_save_plot(p, name = "Thbs2_featureplot", dir = schaystack.dir,height = 6, width = 8)


# Violin plot of top 10 hits showing differences between chondro clusters as a whole and osteo as a whole

top10 <- c("Col1a1","Ptn","Cdh11","Col3a1","Col5a2","Ptprd","Cped1","Postn","Cacna1g","Fbn2")


selected_cells_chondro <- names(obj$cell_types[obj$cell_types %in% c("chondro.1", "chondro.2", "chondro.3", "chondro.4")])
selected_cells_osteo <- names(obj$cell_types[obj$cell_types %in% c("mes.1", "mes.2", "mes.3")])

# Fetch data for both chondro and osteo groups
data_chondro <- FetchData(obj,
                          vars = c(top10, "genotype"),
                          cells = selected_cells_chondro,
                          layer = "data")

data_osteo <- FetchData(obj,
                        vars = c(top10, "genotype"),
                        cells = selected_cells_osteo,
                        layer = "data")

# Filter data to include only the top 10 genes
data_chondro <- data_chondro[, c(top10, "genotype")]
data_osteo <- data_osteo[, c(top10, "genotype")]

# Reshape data into long format
long_data_chondro <- reshape2::melt(data_chondro)
long_data_osteo <- reshape2::melt(data_osteo)

# Combine long-format data for each group
combined_data <- rbind(transform(long_data_chondro, cluster = "Chondro"),
                       transform(long_data_osteo, cluster = "Osteo"))


# Plot combined data using a violin plot
ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
  geom_violin(position = position_dodge(width = 0), width = 10) +  # Adjust width as needed
  geom_jitter(size = 0.5, position = position_dodge2(width = 5), alpha = 0.9) +
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~ variable + cluster, scales = "fixed", nrow = 1, strip.position = "bottom") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
  theme(panel.spacing = unit(0, "lines"))  # Adjust spacing between panels




# Automated with a function
library(ggplot2)
library(reshape2)



# Define a function to generate a violin plot for grouped chondro clusters and osteo clusters
generate_violin_plot <- function(obj, top_genes, chondro_cells, osteo_cells, my_colors) {
  # Fetch data for chondro and osteo groups
  data_chondro <- FetchData(obj, vars = c(top_genes, "genotype"), cells = chondro_cells, layer = "data")
  data_osteo <- FetchData(obj, vars = c(top_genes, "genotype"), cells = osteo_cells, layer = "data")
  
  # Filter data to include only the top genes
  data_chondro <- data_chondro[, c(top_genes, "genotype")]
  data_osteo <- data_osteo[, c(top_genes, "genotype")]
  
  # Reshape data into long format
  long_data_chondro <- reshape2::melt(data_chondro)
  long_data_osteo <- reshape2::melt(data_osteo)
  
  # Combine long-format data for each group
  combined_data <- rbind(transform(long_data_chondro, cluster = "Chondro"),
                         transform(long_data_osteo, cluster = "Osteo"))
  
  # Generate the violin plot
  plot <- ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
    geom_violin(position = position_dodge(width = 0), width = 10) +
    geom_jitter(size = 0.5, position = position_dodge2(width = 5), alpha = 0.8) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white")) +
    facet_wrap(~ variable + cluster, scales = "fixed", nrow = 1, strip.position = "bottom") +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
    theme(panel.spacing = unit(0, "lines"))
  
  return(plot)
}

# Example usage of the function
# Replace 'obj', 'top_genes', 'chondro_cells', 'osteo_cells', and 'my_colors' with actual data and variables
# plot <- generate_violin_plot(obj, top10, selected_cells_chondro, selected_cells_osteo, my_colors)

# get top genes from haystack result file
top_genes_number <- 40 # change as needed
increment <- 4
# Extract top gene names from first column of haystack result, called results
top_genes <- rownames(results)[1:top_genes_number]
print(top_genes)

# Loop through gene sets and generate plots
for (i in seq(1, top_genes_number, by = increment)) {
  start_index <- i
  end_index <- min(i + increment - 1, top_genes_number)
  genes <- rownames(results)[start_index:end_index]
  
  # Generate the violin plot for the current gene set
  vlnplot <- generate_violin_plot(obj, genes, selected_cells_chondro, selected_cells_osteo, my_colors = my_colors)
  
  # Display or save the plot as needed
  print(vlnplot)  # Display the plot
  
  # Save the plot using convenient_save_plot function
  convenient_save_plot(vlnplot, name = paste0("tophaystack_violin_plot_chondro-vs-osteo_", start_index, "_to_", end_index), dir = schaystack.dir, height = 7, width = 12, file_format = "pdf")
}





# Define a function to generate the violin plot for multiple cell types with iterative gene selection
generate_violin_plot_multiple <- function(obj, top_genes_number, increment, cell_types_list, my_colors, file_format = "pdf") {
  # Extract top gene names from first column of haystack result, called results
  top_genes <- rownames(results)[1:top_genes_number]
  
  # Initialize a list to store individual plots
  plots <- list()
  
  # Loop through gene sets and generate plots for each cell type
  for (i in seq(1, top_genes_number, by = increment)) {
    start_index <- i
    end_index <- min(i + increment - 1, top_genes_number)
    genes <- rownames(results)[start_index:end_index]
    
    # Initialize data frame to store combined data for all cell types
    combined_data <- data.frame()
    
    # Iterate over each cell type in the list
    for (cell_type in cell_types_list) {
      # Fetch data for the current cell type
      cell_cells <- names(obj$cell_types[obj$cell_types %in% cell_type])
      data_cell <- FetchData(obj, vars = c(genes, "genotype"), cells = cell_cells, layer = "data")
      
      # Filter data to include only the top genes
      data_cell <- data_cell[, c(genes, "genotype")]
      
      # Reshape data into long format
      long_data_cell <- reshape2::melt(data_cell)
      
      # Add cell type information to the data
      long_data_cell$cluster <- cell_type
      
      # Combine data for the current cell type with the combined data frame
      combined_data <- rbind(combined_data, long_data_cell)
    }
    
    # Generate the violin plot
    plot <- ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
      geom_violin(position = position_dodge(width = 0), width = 10) +
      geom_jitter(size = 0.5, position = position_dodge2(width = 5), alpha = 0.9) +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            panel.grid = element_blank(),
            plot.background = element_rect(fill = "white")) +
      facet_wrap(~ variable + cluster, scales = "fixed", nrow = 1, strip.position = "bottom") +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
      theme(panel.spacing = unit(0, "lines"))
    
    # Save the plot using convenient_save_plot function
    convenient_save_plot(plot, name = paste0("tophaystack_violin_plot_chondro1-4_vs_mes1-3_", start_index, "_to_", end_index), dir = schaystack.dir, height = 6, width = 16, file_format = file_format)
    
    # Add the plot to the list of plots
    plots[[paste0(start_index, "_to_", end_index)]] <- plot
  }
  
  return(plots)  # Return the list of plots
}

# Example usage of the function with multiple cell types and iterative gene selection
top_genes_number <- 30  # change as needed
increment <- 3
cell_types_list <- c("chondro.1", "chondro.2", "chondro.3", "chondro.4", "mes.1", "mes.2", "mes.3")
plots <- generate_violin_plot_multiple(obj, top_genes_number, increment, cell_types_list, my_colors)



# Define a function to generate FeaturePlots with iterative gene selection
generate_feature_plots <- function(sample, top_genes_number, increment, gene_names, genotype_col, dir) {
  # Extract top gene names based on top_genes_number and increment
  gene_sets <- split(gene_names, ceiling(seq_along(gene_names) / increment))
  
  # Loop through gene sets and generate FeaturePlots
  for (i in seq_along(gene_sets)) {
    gene_set <- gene_sets[[i]]
    
    # Generate FeaturePlot for the current gene set
    p <- SCpubr::do_FeaturePlot(
      sample = sample,
      features = gene_set,
      split.by = genotype_col,
      font.size = 7,
      use_viridis = TRUE,
      viridis.palette = "cividis",
      viridis.direction = 1,
      na.value = "floralwhite",
      plot_cell_borders = FALSE,
      order = TRUE,
      pt.size = 1,
      legend.position = "right"
    )
    
    # Save the plot using convenient_save_plot function with incremental naming
    convenient_save_plot(
      p,
      name = paste0("tophaystack_FeaturePlots_", i * increment - increment + 1, "_to_", i * increment),
      dir = dir,
      height = 4,
      width = 18,
      file_format = "pdf"
    )
  }
}

# Example usage of the function
top_genes_number <- 24  # change as needed
increment <- 4
gene_names <- rownames(results)[1:top_genes_number]
genotype_col <- "genotype"

# Replace 'sample' with your actual sample data and 'schaystack_dir' with your directory path
generate_feature_plots(sample = osteochondro_subset, top_genes_number, increment, gene_names, genotype_col, dir = schaystack.dir)
