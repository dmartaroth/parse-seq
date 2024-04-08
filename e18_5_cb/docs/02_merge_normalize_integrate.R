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


# Merge control with controltopup -----------------------------------------

# Load objects ------------------------------------------------------------

ctrl <- readRDS(file = paste0(data.output,"/filtered_",control,"_",sample,".Rds"))
ctrltopup <- readRDS(file =  paste0(data.output,"/filtered_",control,"_",sample,"_topup",".Rds"))


# Merge control objects ---------------------------------------------------

ctrl$replicate <- "1" # assign replicate labels
ctrltopup$replicate <- "2"                  
                     
# Merge seurat objects
merged_control <- merge(ctrl, y = ctrltopup, add.cell.ids = c("1", "2"), project = control)
merged_control <- JoinLayers(merged_control)                   
                     
# Split RNA measurements into 2 layers: one for each replicate
merged_control[["RNA"]] <- split(merged_control[["RNA"]], f = merged_control$replicate)
merged_control


# Perform analysis without integration ------------------------------------

# Standard analysis workflow
merged_control <- NormalizeData(merged_control)
merged_control <- FindVariableFeatures(merged_control)
merged_control <- ScaleData(merged_control)
merged_control <- RunPCA(merged_control)

merged_control <- FindNeighbors(merged_control, dims = 1:30, reduction = "pca")
merged_control <- FindClusters(merged_control, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_control <- RunUMAP(merged_control, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
plot <- DimPlot(merged_control, reduction = "umap.unintegrated", group.by = c("replicate", "seurat_clusters")) +
  umap_theme() + scale_color_manual(values = pastel_palette)   

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_unintegrated_umap.png", plot_number)), width = 8, height = 3, plot)



# Perform integration -----------------------------------------------------

merged_control <-
  IntegrateLayers(
    object = merged_control,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_control[["RNA"]] <- JoinLayers(merged_control[["RNA"]]) # re-join layers after integration

merged_control <- FindNeighbors(merged_control, reduction = "integrated.cca", dims = 1:30)

merged_control <- FindClusters(merged_control, verbose=FALSE,
                              resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)



(plot <- clustree(merged_control, prefix = "RNA_snn_res.") )

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated_UMAP_clustree.png", plot_number)), width = 8, height = 9, plot)

merged_control <- RunUMAP(merged_control, dims = 1:30, reduction = "integrated.cca")

# Low resolution clustering best at 0.1 for this control
(plot <- DimPlot(merged_control, reduction = "umap", group.by = c("replicate", "RNA_snn_res.0.1"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_lowres_integrated-cca_umap_res.0.1.png", plot_number)), width = 8, height = 3, plot)


# Higher res clustering good at 0.4
(plot <- DimPlot(merged_control, reduction = "umap", group.by = c("RNA_snn_res.0.4"))+
    umap_theme() + scale_color_manual(values = pastel_palette))
(plot <- plot + theme(legend.position = "bottom"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_highres_integrated-cca_umap_res.0.4.png", plot_number)), width = 4, height =5, plot)


# Visualize conditions side by side
plot <- DimPlot(merged_control, reduction = "umap", split.by = "replicate", group.by = "RNA_snn_res.0.4")+
  umap_theme() + scale_color_manual(values = pastel_palette) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated-cca_umap_split.png", plot_number)), width = 7, height = 3, plot)


# Save integrated filtered object -----------------------------------------

# Save filtered cells -----------------------------------------------------
saveRDS(merged_control, file = paste0(data.output,"/integrated_filtered_",control,"_",sample,"_mergedreplicates",".Rds"))



# Merge ncko with nckotopup -----------------------------------------

# Load objects ------------------------------------------------------------

ncko <- readRDS(file = paste0(data.output,"/filtered_",mutant,"_",sample,".Rds"))
nckotopup <- readRDS(file =  paste0(data.output,"/filtered_",mutant,"_",sample,"_topup",".Rds"))


# Merge ncko objects ---------------------------------------------------

ncko$replicate <- "1" # assign replicate labels
nckotopup$replicate <- "2"                  

# Merge seurat objects
merged_ncko <- merge(ncko, y = nckotopup, add.cell.ids = c("1", "2"), project = mutant)
merged_ncko <- JoinLayers(merged_ncko)                   

# Split RNA measurements into 2 layers: one for each replicate
merged_ncko[["RNA"]] <- split(merged_ncko[["RNA"]], f = merged_ncko$replicate)
merged_ncko


# Perform analysis without integration ------------------------------------

# Standard analysis workflow
merged_ncko <- NormalizeData(merged_ncko)
merged_ncko <- FindVariableFeatures(merged_ncko)
merged_ncko <- ScaleData(merged_ncko)
merged_ncko <- RunPCA(merged_ncko)

merged_ncko <- FindNeighbors(merged_ncko, dims = 1:30, reduction = "pca")
merged_ncko <- FindClusters(merged_ncko, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_ncko <- RunUMAP(merged_ncko, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
plot <- DimPlot(merged_ncko, reduction = "umap.unintegrated", group.by = c("replicate", "seurat_clusters")) +
  umap_theme() + scale_color_manual(values = pastel_palette)   

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_unintegrated_umap.png", plot_number)), width = 8, height = 3, plot)



# Perform integration -----------------------------------------------------

merged_ncko <-
  IntegrateLayers(
    object = merged_ncko,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_ncko[["RNA"]] <- JoinLayers(merged_ncko[["RNA"]]) # re-join layers after integration

merged_ncko <- FindNeighbors(merged_ncko, reduction = "integrated.cca", dims = 1:30)

merged_ncko <- FindClusters(merged_ncko, verbose=FALSE,
                               resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)



(plot <- clustree(merged_ncko, prefix = "RNA_snn_res.") )

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated_UMAP_clustree.png", plot_number)), width = 8, height = 9, plot)

merged_ncko <- RunUMAP(merged_ncko, dims = 1:30, reduction = "integrated.cca")

# Low resolution clustering best at 0.1 for this ncko
(plot <- DimPlot(merged_ncko, reduction = "umap", group.by = c("replicate", "RNA_snn_res.0.1"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_lowres_integrated-cca_umap_res.0.1.png", plot_number)), width = 8, height = 3, plot)


# Higher res clustering good at 0.4
(plot <- DimPlot(merged_ncko, reduction = "umap", group.by = c("RNA_snn_res.0.4"))+
    umap_theme() + scale_color_manual(values = pastel_palette))
(plot <- plot + theme(legend.position = "bottom"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_highres_integrated-cca_umap_res.0.4.png", plot_number)), width = 4, height =5, plot)


# Visualize conditions side by side
plot <- DimPlot(merged_ncko, reduction = "umap", split.by = "replicate", group.by = "RNA_snn_res.0.4")+
  umap_theme() + scale_color_manual(values = pastel_palette) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated-cca_umap_split.png", plot_number)), width = 7, height = 3, plot)


# Save integrated filtered object -----------------------------------------

# Save filtered cells -----------------------------------------------------
saveRDS(merged_ncko, file = paste0(data.output,"/integrated_filtered_",mutant,"_",sample,"_mergedreplicates",".Rds"))





# Merge integrated control with integrated mutant -------------------------

ctrl <- readRDS(file = paste0(data.output,"/integrated_filtered_",control,"_",sample,"_mergedreplicates",".Rds"))
ncko <- readRDS(file = paste0(data.output,"/integrated_filtered_",mutant,"_",sample,"_mergedreplicates",".Rds"))
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

merged_sample <- FindClusters(merged_sample, verbose=FALSE,
  resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)



(plot <- clustree(merged_sample, prefix = "RNA_snn_res.") )

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated_UMAP_clustree.png", plot_number)), width = 8, height = 9, plot)

merged_sample <- RunUMAP(merged_sample, dims = 1:30, reduction = "integrated.cca")

# Low resolution clustering best at 0.1 for this sample
(plot <- DimPlot(merged_sample, reduction = "umap", group.by = c("genotype", "RNA_snn_res.0.1"))+
  umap_theme() + scale_color_manual(values = pastel_palette)) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_lowres_integrated-cca_umap_res.0.1.png", plot_number)), width = 8, height = 3, plot)


# Higher res clustering good at 0.3
(plot <- DimPlot(merged_sample, reduction = "umap", group.by = c("RNA_snn_res.0.3"))+
  umap_theme() + scale_color_manual(values = pastel_palette))
(plot <- plot + theme(legend.position = "bottom"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_highres_integrated-cca_umap_res.0.3.png", plot_number)), width = 4, height =5, plot)


# Visualize conditions side by side
plot <- DimPlot(merged_sample, reduction = "umap", split.by = "genotype", group.by = "RNA_snn_res.0.3")+
  umap_theme() + scale_color_manual(values = pastel_palette) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(integration_dir, sprintf("%02d_integrated-cca_umap_split.png", plot_number)), width = 7, height = 3, plot)


# Save integrated filtered object -----------------------------------------

# Save filtered cells -----------------------------------------------------
saveRDS(merged_sample, file = paste0(data.output,"/integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))



