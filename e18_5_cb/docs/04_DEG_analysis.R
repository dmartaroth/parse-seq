# ## ######################################## ## #
#          COMPARATIVE ANALYSIS (DEG)            #
# ## ######################################## ## #

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "05_DEG_Analysis"
analysis_dir <- file.path(figs, "05_DEG_Analysis")
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir)
}

plot_number <- 0  # Starting plot number


# Load annotated integrated object ----------------------------------------

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))

my_colors =  c("#E6B0C2","#FADBD8","#FFB5B5","thistle1",
               "#424949",
               "#ABEBC6", "#1C7F82",
               "#7A8D0A","#C7CC8F", "darkolivegreen3",
               "powderblue", "#7EBDC2", "#2E86C1", 
               "pink3", "#F1C41F", "#B7A4DB",  "#76448A","darkseagreen",
               "#F1948A", "#CB4335")
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


splitcolors =  c("#E6B0C2","#E6B0C2","#FADBD8","#FADBD8","#FFB5B5","#FFB5B5",
                 "thistle1","thistle1","#424949",
                 "#424949","#ABEBC6",
                 "#ABEBC6", "#1C7F82", "#1C7F82",
                 "#7A8D0A", "#7A8D0A","#C7CC8F", "#C7CC8F", "darkolivegreen3",
                 "darkolivegreen3","powderblue", "powderblue", "#7EBDC2","#7EBDC2", 
                 "#2E86C1", "#2E86C1","pink3", "pink3", "#F1C41F","#F1C41F",
                 "#B7A4DB","#B7A4DB",  "#76448A","#76448A","darkseagreen",
                 "darkseagreen","#F1948A","#F1948A", "#CB4335", "#CB4335")

(plot <- DimPlot(obj, reduction = "umap", label = FALSE,split.by = "genotype")+
    umap_theme() + scale_color_manual(values = splitcolors)) + theme(legend.position = "bottom")

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_UMAP_celltype_genotype_split.png", plot_number)), width = 12, height = 5, plot)


# Define your cluster names
cluster_names <- c("chondro.1", "chondro.2", "chondro.3", "chondro.4",
                   "div", "epith.1","epith.2","imm.1","imm.2","imm.3",
                   "mes.1","mes.2","mes.3", "neu.1", "neu.2", "neu.3", "neu.4", 
                   "neu.5","vasc")

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
                 pt.size = 0.1, combine = FALSE, split.plot = TRUE,cols = my_colors)
p <- wrap_plots(plots = plots, ncol = 3)

p <- p+ theme(plot.background = element_rect(fill = "white"))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_top-pos-response-markers_vln_split.png", plot_number)), width = 49, height = 49, p)

# Xist is top for many clusters, but this is a sex-linked difference. Will plot next top 3 markers where present.

# chondro.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Cntn3","Peak1","Rtl3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# chondro.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Gm26992","Brip1", "Spp1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# chondro.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("4932411K12Rik","Cemip", "Steap4"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# chondro.4
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Vit","Gm26992", "Cntn3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.4_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# div
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Thbs2","Pdzrn4","Fat3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_div_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



# epith.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Gh","Chgb","Nrg1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_epith.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# epith.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Lmntd1","Wdr49",'Kndc1'),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_epith.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# imm.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Ccl7","Cmpk2", "Ifit2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_imm.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# imm.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Ibsp","Hist1h4a", "Eya2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_imm.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# imm.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Pi15","Pard3b","Tns1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_imm.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



# mes.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Tec","Tollip","Rpl10"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# mes.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Mfap5","Cfh","Htra3"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# mes.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Cfh","Ston2","Prss23"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# neu.1
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Pomc","Kcnj6","Gm45846"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



# neu.2
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Ptpro","Rfx4","Pitx2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# neu.3
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Gh","Nr3c2","Pou1f1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# neu.4
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Scn7a","Tmem108","Pi15"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.4_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# neu.5
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Sostdc1","Scn7a","Ugt8a"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.5_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# vasc
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Myct1","Mmrn2","Clec14a"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_vasc_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



