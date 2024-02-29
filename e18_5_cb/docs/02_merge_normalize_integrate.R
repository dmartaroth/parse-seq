# ## ######################################## ## #
#           NORMALIZE FILTERED OBJECTS           #
# ## ######################################## ## #

# This script normalizes filtered objects prior to integration
# Date: Thu Feb 29 12:28:45 2024 ------------------
# Daniela M. Roth

# Packages, directories, and functions ------------------------------------

library(here)
source(here("e18_5_cb","docs","packages.R")) # load packages
source(here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here("e18_5_cb","docs","functions.R")) # load functions
source(here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "03_Integration"
integration_dir <- file.path(figs, "03_Integration")
if (!dir.exists(integration_dir)) {
  dir.create(integration_dir)
}
plot_number <- 0  # Starting plot number

# Load objects ------------------------------------------------------------

ctrl <- readRDS(file = paste0(data.output,"/filtered_",control,"_",sample,".Rds"))
ncko <- readRDS(file = paste0(data.output,"/filtered_",mutant,"_",sample,".Rds"))


# Merge objects -----------------------------------------------------------

ctrl$genotype <- control # assign genotype labels
ncko$genotype <- mutant

# Merge seurat objects
merged_sample <- merge(ctrl, y = ncko, add.cell.ids = c(control, mutant), project = sample)
merged_sample <- JoinLayers(merged_sample)

# Split RNA measurements into 2 layers: one for control cells, one for mutant cells
merged_sample[["RNA"]] <- split(merged_sample[["RNA"]], f = merged_sample$genotype)
merged_sample


# Perform analysis without integration ------------------------------------

# Standard analysis workflow
merged_sample <- NormalizeData(merged_sample)
merged_sample <- FindVariableFeatures(merged_sample)
merged_sample <- ScaleData(merged_sample)
merged_sample <- RunPCA(merged_sample)

merged_sample <- FindNeighbors(merged_sample, dims = 1:30, reduction = "pca")
merged_sample <- FindClusters(merged_sample, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_sample <- RunUMAP(merged_sample, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
plot <- DimPlot(merged_sample, reduction = "umap.unintegrated", group.by = c("genotype", "seurat_clusters")) +
  umap_theme() + scale_color_manual(values = pastel_palette)   

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_unintegrated_umap.png", plot_number)), width = 8, height = 3, plot)



# Perform integration -----------------------------------------------------


merged_sample <-
  IntegrateLayers(
    object = merged_sample,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_sample[["RNA"]] <- JoinLayers(merged_sample[["RNA"]]) # re-join layers after integration

merged_sample <- FindNeighbors(merged_sample, reduction = "integrated.cca", dims = 1:30)
merged_sample <- FindClusters(merged_sample, resolution = 0.6)

merged_sample <- RunUMAP(merged_sample, dims = 1:30, reduction = "integrated.cca")

plot <- DimPlot(merged_sample, reduction = "umap", group.by = c("genotype", "seurat_clusters")) +
  umap_theme() + scale_color_manual(values = pastel_palette) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated-cca_umap.png", plot_number)), width = 8, height = 3, plot)

# Visualize conditions side by side
plot <- DimPlot(merged_sample, reduction = "umap", split.by = "genotype")+
  umap_theme() + scale_color_manual(values = pastel_palette) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated-cca_umap_split.png", plot_number)), width = 7, height = 3, plot)


# Save integrated filtered object -----------------------------------------

# Save filtered cells -----------------------------------------------------
saveRDS(merged_sample, file = paste0(data.output,"/integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))



