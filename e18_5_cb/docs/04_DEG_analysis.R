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
               "#ABEBC6","darkolivegreen3", "#1C7F82",
               "powderblue", "#7EBDC2", "#2E86C1", "#1B4F72",
               "#7A8D0A","#C7CC8F", 
               "pink3", "#F1C41F", "#B7A4DB",  "#76448A","darkseagreen",
               
               "coral3", "#CB4335", "thistle3",   "#F1948A",
               
               "#4A235A", "steelblue","red2", 
               "#F4D03F", "red2")
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
                 "thistle1","thistle1",
                 "#424949","#424949",
                 "#ABEBC6","#ABEBC6","darkolivegreen3","darkolivegreen3", "#1C7F82",
                 "#1C7F82",
                 "powderblue",                 "powderblue", "#7EBDC2",  "#7EBDC2", 
                 "#2E86C1","#2E86C1", "#1B4F72", "#1B4F72",
                 "#7A8D0A","#7A8D0A","#C7CC8F","#C7CC8F", 
                 "pink3","pink3", "#F1C41F", "#F1C41F", "#B7A4DB", "#B7A4DB",
                 "#76448A","#76448A","darkseagreen","darkseagreen",
                 "coral3","coral3", "#CB4335", "#CB4335", "thistle3", "thistle3",
                 "#F1948A","#F1948A",
                 "#4A235A","#4A235A", "steelblue", "steelblue")

(plot <- DimPlot(obj, reduction = "umap", label = FALSE,split.by = "genotype")+
    umap_theme() + scale_color_manual(values = splitcolors)) + theme(legend.position = "bottom")

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_UMAP_celltype_genotype_split.png", plot_number)), width = 12, height = 5, plot)


# Define your cluster names
cluster_names <- c("chondro.1", "chondro.2", "chondro.3",
                   "ciliated", "epith","imm.1","imm.2",
                   "mes.1","mes.2","mes.3", "mes.4","mito.1","mito.2",
                   "neu.1", "neu.2", "neu.3", "neu.4", 
                   "neu.5","neu.6","vasc.1","vasc.2")

# Call the function to find top markers for all clusters
all_top_pos_markers <- findTopPosMarkers(obj, cluster_names)

# Access top markers for a specific cluster
head(all_top_pos_markers$chondro.1)

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
      features = c("Gm26992","Cntn3","Elmo1"),
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
head(all_top_pos_markers$chondro.2)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Ndnf","Brinp1", "Cd55"),
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
head(all_top_pos_markers$chondro.3)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("4932411K12Rik","Cemip", "Brip1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_chondro.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# ciliated
head(all_top_pos_markers$ciliated)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Lrrc36","Slc9a3r1", "Dnah7a"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_ciliated_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# epith
head(all_top_pos_markers$epith)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Gh","Chgb","Nnat"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_epith_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# imm.1
head(all_top_pos_markers$imm.1)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Col1a2","Ror1", "Cdk8"),
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
head(all_top_pos_markers$imm.2)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Eif2b5","Car13", "Akap6"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_imm.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



# mes.1
head(all_top_pos_markers$mes.1)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Tec","1700029H14Rik","Gli1"),
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
head(all_top_pos_markers$mes.2)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Mfap5","Cfh","Il1rapl1"),
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
head(all_top_pos_markers$mes.3)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Gdf10","Cyp1b1","Kcnd2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes.3_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# mes.4
head(all_top_pos_markers$mes.4)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Cfh","Chrdl1","Kcnd2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mes.4_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)



# mito.1
head(all_top_pos_markers$mito.1)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Vit","Cntn3","Gm26992"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mito.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# mito.2
head(all_top_pos_markers$mito.2)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Ibsp","Mmp13","Spp1"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_mito.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# neu.1
head(all_top_pos_markers$neu.1)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Scn7a","Ak5","Epha6"),
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
head(all_top_pos_markers$neu.2)

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Rab34","Mgat4c","Phykpl"),
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
head(all_top_pos_markers$neu.3)

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Pomc","Ptpro","Pitx2"),
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
head(all_top_pos_markers$neu.4)

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Pomc","Gpc5","Gh"),
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
head(all_top_pos_markers$neu.5)

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Pomc","Gm45455","Rit2"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.5_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# neu.6
head(all_top_pos_markers$neu.6)

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Scn7a","Entpd2","Gm26992"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_neu.6_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

# vasc.1
head(all_top_pos_markers$vasc.1)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Peak1","Syne1","Filip1l"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_vasc.1_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)


# vasc.2
head(all_top_pos_markers$vasc.2)
(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      split.by = "genotype",
      na_cutoff = 0,
      features = c("Dync1i1","Cyp26b1","Fam78b"),
      colors_use = c(
        "grey95",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(analysis_dir, sprintf("%02d_vasc.2_DEG_featureplots.png", plot_number)), width = 6, height = 7, plot)

