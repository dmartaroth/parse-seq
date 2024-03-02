# ## ######################################## ## #
#               COMPARATIVE ANALYSIS             #
# ## ######################################## ## #

library(here)
source(here("e18_5_cb","docs","packages.R")) # load packages
source(here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here("e18_5_cb","docs","functions.R")) # load functions
source(here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "05_Analysis"
analysis_dir <- file.path(figs, "05_Analysis")
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir)
}

plot_number <- 0  # Starting plot number


# Load annotated integrated object ----------------------------------------

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))

my_colors =  c("#E6B0C2","#FADBD8","#FFB5B5","pink3", "thistle1",
               "#ABEBC6",  "powderblue","red2", 
               "#B7A4DB",  "#76448A", "#F1948A", "thistle3","#2E86C1", 
               "#424949","#9A7D0A","#1C7F82", "steelblue","#7EBDC2",
               "#F4D03F","#C7CC8F", "#1B4F72","#CB4335",
               "darkgreen", "#873600", "#4A235A", "#F1C41F",
               "red2")
(plot <- DimPlot(obj, reduction = "umap", label = FALSE)+
    umap_theme() + scale_color_manual(values = my_colors))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_UMAP_annotated.png", plot_number)), width = 6, height = 4, plot)


# Identify conserved cell type markers ------------------------------------
# Create an empty dataframe to store all markers
all_markers <- data.frame()

# Loop through each cluster
for (cluster_id in levels(obj)) {
  # Find conserved markers for the current cluster
  cluster_markers <- FindConservedMarkers(obj, ident.1 = cluster_id, grouping.var = "genotype", verbose = FALSE)
  
  # Add cluster ID column
  cluster_markers$Cluster <- cluster_id
  
  # Append cluster markers to all_markers dataframe
  all_markers <- rbind(all_markers, cluster_markers)
}

# Export the combined markers dataframe to a CSV file
write.csv(all_markers, file = here(data.output, "annotated_cluster_markers_all.csv"), row.names=TRUE)



# Identify DEG across conditions ------------------------------------------
obj$cell_types <- Idents(obj)

# Which genes change in different conditions for cells of the same type?
obj$celltype_genotype <- paste(obj$cell_types, obj$genotype, sep = "_")
Idents(obj) <- "celltype_genotype"

Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))


splitcolors =  c("#E6B0C2","#E6B0C9","#FADBD8","#FADBD1","#FFB5B5","#FFB5B9","pink3","pink4", "thistle1",
                 "thistle2",
               "#ABEBC6", "#ABEBC9", "powderblue","powderblue","red2","red2", 
               "#B7A4DB",   "#B7A4DB",  "#76448A","#76448A", "#F1948A", "#F1948A", "thistle3", "thistle3",
               "#2E86C1","#2E86C1", 
               "#424949","#424949","#9A7D0A","#9A7D0A","#1C7F82","#1C7F82", "steelblue", "steelblue",
               "#7EBDC2","#7EBDC2", "#F4D03F","#F4D03F","#C7CC8F","#C7CC8F",
               "#1B4F72", "#1B4F72","#CB4335","#CB4335", "darkgreen",
               "darkgreen", "#873600", "#873600", "#4A235A","#4A235A", "#F1C41F",
               "#F1C41F","red2","red2")

(plot <- DimPlot(obj, reduction = "umap", label = FALSE,split.by = "genotype")+
    umap_theme() + scale_color_manual(values = splitcolors)) + theme(legend.position = "bottom")

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_UMAP_celltype_genotype_split.png", plot_number)), width = 12, height = 7, plot)


# Define your cluster names
cluster_names <- c("chondro.1", "chondro.2", "chondro.3", "chondro.4",
                   "endocr", "epith", "mes", "musc", "myeloid",
                   "neu.1", "neu.2", "neu.3", "neu.4", "neu.5",
                   "olf", "osteo.1", "osteo.2", "osteocl",
                   "unknown.1", "unknown.2", "unknown.3", "vasc")

# Call the function to find top markers for all clusters
all_top_pos_markers <- findTopPosMarkers(obj, cluster_names)

# Access top markers for a specific cluster
head(all_top_pos_markers$chondro.3)

# Call the function to extract and order top markers
ordered_top_pos_markers <- extractTopPosMarkers(all_top_pos_markers)

sink(here(data.output,"top_pos_markers_output.txt"))
# Print top markers for each cluster
for (cluster_name in cluster_names) {
  cat("Top markers for", cluster_name, "ordered by avg_log2FC:\n")
  print(ordered_top_pos_markers[[cluster_name]])
  cat("\n")
}

sink()

topposmarkers <- collectTopPosMarkers(all_top_pos_markers)



plots <- VlnPlot(obj, features = topposmarkers, split.by = "genotype", group.by = "cell_types",
                 pt.size = 0.1, combine = FALSE, split.plot = TRUE)
p <- wrap_plots(plots = plots, ncol = 3)

p <- p+ theme(plot.background = element_rect(fill = "white"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_top-pos-response-markers_vln_split.png", plot_number)), width = 49, height = 49, p)



# chondro.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Xist","Cntn3","Peak1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.1_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# chondro.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Xist","Brip1","Smpd3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.2_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# chondro.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("4932411K12Rik","Xist","Cemip"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.3_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# chondro.4
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Xist","Vit","Gm26992"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.4_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# mes
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Mfap5","Xist","Cfh"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)


# osteo.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Xist","Ibsp","Hist1h4a"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_osteo.1_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# osteo.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Xist","Kcnd2","Cfh"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_osteo.2_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)


# osteoclasts
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Tmem267","Sdk1","Col12a1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_osteocl_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)

# unknown.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Peak1","Cdk8","Col1a2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_unknown.1_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)


# unknown.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Cdk8","mt-Rnr2","Ahnak"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_unknown.2_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)


# unknown.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Nkain3","Zbtb26","Eif4a3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_unknown.3_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)


# vasc
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Stab1","Peak1","Mrps6"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_vasc_DEG_featureplots.png", plot_number)), width = 4, height = 7, plot)



