# Annotate cells using monocle
# devtools::install_github("satijalab/seurat-wrappers", force = TRUE)
library(monocle3)

cds <- SeuratWrappers::as.cell_data_set(readRDS(file = paste0(data.output,"/integrated_filtered_",control,"_",mutant,"_",sample,".Rds")))

cds <- preprocess_cds(cds, num_dim = 100)
plot <- plot_pc_variance_explained(cds) +custom_theme_scatter()

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_PC_variance_explained.png", plot_number)), width = 6, height = 5, plot)


cds <- reduce_dimension(cds)

plot_cells(cds)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
(plot <- plot_cells(cds, genes = c("Gli1","Sp7","Col2a1","Col10a1"))+
    theme_feature_plot() + scale_color_gradient(low = "darkseagreen1", high = "darkgreen", na.value = "cornsilk"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_Feature_plots_Gli1_Sp7_Col2a1_Col10a1.png", plot_number)), width = 6, height = 5, plot)


# Group cells into clusters -----------------------------------------------

cds <- cluster_cells(cds, resolution=4e-4)
plot <- plot_cells(cds)+ umap_theme()+ scale_color_manual(values = pastel_palette)  

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_res_4e-4.png", plot_number)), width = 5, height = 4, plot)


plot <- plot_cells(cds, color_cells_by = "partition", group_cells_by = "partition") + umap_theme() + scale_color_manual(values = pastel_palette)

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_partitions.png", plot_number)), width = 6, height = 5, plot)

# Cells are now grouped by "cluster" or "partition". Coloring will be by cluster
# by default unless partition is specified
head(partitions(cds))
head(clusters(cds))
# Find markers -----------------------------------------------------
marker_test_res <- top_markers(cds,
                               group_cells_by="cluster",
                               reduction_method = "UMAP",
                               marker_sig_test = TRUE,
                               reference_cells=1000,
                               cores=8)

# Explore markers
head(marker_test_res)
dim(marker_test_res)
length(which(duplicated(marker_test_res$gene_id)))
duplicate_markers <- names(which(table(marker_test_res$gene_id) > 1))
head(marker_test_res[marker_test_res$gene_id %in% duplicate_markers,])
unique_markers <- marker_test_res[!(marker_test_res$gene_id %in% duplicate_markers),]
head(unique_markers)
rm(marker_test_res, duplicate_markers)

top_specific_markers <- unique_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  arrange(desc(specificity), .by_group = TRUE) %>%
  dplyr::slice(1:3) %>%
  pull(gene_id)

(plot <- plot_genes_by_group(cds,
                             top_specific_markers,
                             group_cells_by="cluster",
                             ordering_type="maximal_on_diag",
                             max.size=3) + scale_color_gradient(low = "cornsilk", high = "red"))


# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_dotplot_clusters_unique_markers.png", plot_number)), width = 5, height = 5, plot)


# Defined marker list visualization ---------------------------------------

my_genes <-c("Prrx1",  "Gli1",  "Axin2",   "Sox9",   "Six2",  "Col1a1",   "Alpl",  "Bglap",
             "Runx2",  "Crabp1",   "Sp7",  "Phex",   "Dmp1",   "Sost",   "Col2a1",
             "Col10a1",  "Acan",   "Chad",   "Ctsk",   "Acp5")


(  plot <- plot_genes_by_group(
  cds,
  my_genes,
  group_cells_by = "cluster",
  norm_method = "size_only",
  ordering_type = "cluster_row_col",
  max.size = 5, 
) + scale_color_gradient(low = "goldenrod1", high = "red"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_dotplot_clusters_osteochondrodifferentiation.png", plot_number)), width =5, height = 3.5, plot)

colData(cds)$monocle_clusters <- as.character(clusters(cds))
colData(cds)$monocle_partitions <- as.character(partitions(cds))

# Visualize marker genes
(plot <- plot_cells(cds, genes = c("Gli1","Sp7","Col2a1","Col10a1"))+
    theme_feature_plot() + scale_color_gradient(low = "darkseagreen1", high = "darkgreen", na.value = "cornsilk"))


# Save Monocle data to return to Seurat -----------------------------------


# Extract expression matrix and cell metadata from Monocle object
expression_matrix <- counts(cds)
cell_metadata <- as.data.frame(colData(cds))
umap_coords <- reducedDims(cds)$UMAP

# Check column names in cell_metadata
print(colnames(cell_metadata))

# Create Seurat object
clus.data <- CreateSeuratObject(counts = expression_matrix, project = "monocleclusters")

# Assign metadata to Seurat object
clus.data@meta.data <- cell_metadata
# Create new columns in meta.data of Seurat object for UMAP coordinates
clus.data@meta.data$UMAP_1 <- umap_coords[, 1]
clus.data@meta.data$UMAP_2 <- umap_coords[, 2]
clus.data
colnames(clus.data@meta.data)

DimPlot(clus.data, reduction = "UMAP_1:UMAP_2", split.by = "genotype")+
  umap_theme() + scale_color_manual(values = pastel_palette) 
