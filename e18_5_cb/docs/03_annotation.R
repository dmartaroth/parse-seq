# ## ######################################## ## #
#                 CLUSTER ANNOTATION             #
# ## ######################################## ## #

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "04_Clustering"
clustering_dir <- file.path(figs, "04_Clustering")
if (!dir.exists(clustering_dir)) {
  dir.create(clustering_dir)
}

plot_number <- 0  # Starting plot number


# Cluster marker identification -------------------------------------------

clustered_renamed <- readRDS(file = paste0(data.output,"/integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))

## Low resolution cluster annotation ---------------------------------------
lowres.clusters <- "RNA_snn_res.0.1"
(plot <- DimPlot(clustered_renamed, reduction = "umap", group.by = lowres.clusters)+
  umap_theme() + scale_color_manual(values = pastel_palette)) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_%s.png", plot_number, lowres.clusters)), width = 3, height = 3, plot)


Idents(clustered_renamed) <- clustered_renamed$RNA_snn_res.0.1

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(clustered_renamed,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "0.1res_all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "0.1res_top50_markers_pct.csv"))


(plot <- DotPlot(
  object = clustered_renamed,
  features =   top_5,
  scale.by = "size",
  dot.scale = 10,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_DotPlot_top5_markers_%s.png", plot_number, lowres.clusters)), width = 16, height = 4, plot)


# List of genes for feature plots
progenitors <- c("Axin2","Gli1","Prrx1","Six2")
osteogenic <- c("Crabp1","Runx2","Sp7","Dmp1")
chondrogenic <- c("Col2a1","Acan","Mgp","Sox9")
hypertrophy <- c("Runx2","Col10a1","Mmp13","Fgfr1")
osteoclasts <- c("Ctsk","Mmp9","Pheta1","Cd44")
vascular <- c("Mcam","Vwf","Pecam1","Pdgfrb")
myeloid_lymphocyte <- c("Pou2f2","Il1rl1","Gata2")
neurons_gli1 <- c("Neurod1","Cplx3","Otx2","Gfra3","Sox10","Foxd3")
erythrocytes <- c("Hba-a1","Hba-a2","Hbb-bs","Gypa","Gybp","Alas2","Klf1","Slc25a37","Slc2a1")
smoothmuscle <- c("Acta2","Tagln","Myh11","Des")

# Progenitors

(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_progenitor_featureplots.png", plot_number)), width = 6, height = 5, plot)


(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_osteogenic_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Chondrogenic
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_chondrogenic_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Hypertrophy
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = hypertrophy,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_hypertrophy_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Osteoclasts
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_osteoclasts_featureplots.png", plot_number)), width = 6, height = 5, plot)


# Vascular
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_vascular_featureplots.png", plot_number)), width = 6, height = 5, plot)


# Myeloid and lymphocyte
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_myeloid_lymphocyte_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Neurons and glial cells
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_neurons_glia_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Erythrocytes
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_erythrocytes_featureplots.png", plot_number)), width = 6, height = 5, plot)

# Smooth muscle
(plot <-
    FeaturePlot_scCustom(
      seurat_object = clustered_renamed,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_smoothmuscle_featureplots.png", plot_number)), width = 6, height = 5, plot)



# Based on these plots, identify potential clusters at low resolution

# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "0.1res_cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"0.1res_cluster_annotation.csv"))

# Rename clusters
clustered_renamed_updated <- Rename_Clusters(seurat_object = clustered_renamed, new_idents = annotation_info$new_cluster_idents, meta_col_name = "lowres_annotation")

(plot <- DimPlot(clustered_renamed_updated, reduction = "umap",label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = my_colors)+
    umap_theme()) 



# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_lowres_annotated_%s.png", plot_number, lowres.clusters)), width = 3, height = 3, plot)



## High res cluster annotation ---------------------------------------------
highres.clusters <- "RNA_snn_res.0.4"

(plot <- DimPlot(clustered_renamed, reduction = "umap", label = TRUE,group.by = highres.clusters)+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_%s.png", plot_number, highres.clusters)), width = 5, height = 3, plot)


Idents(clustered_renamed) <- clustered_renamed$RNA_snn_res.0.4

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(clustered_renamed,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "0.4res_all_markers_pct.csv"))

top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "0.4res_top50_markers_pct.csv"))

# Extract the top N marker genes per cluster
top_3 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 3, named_vector = FALSE,make_unique = TRUE,rank_by = "avg_log2FC")
write.csv(top_3, file = here(data.output, "0.4res_top3_markers_pct.csv"))



(plot <- DotPlot( object = clustered_renamed,
  features =   top_3,
  scale.by = "size",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = TRUE
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_DotPlot_top3_markers_%s.png", plot_number, highres.clusters)), width = 17, height = 7, plot)


# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "0.4res_cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"0.4res_cluster_annotation.csv"))

Idents(clustered_renamed_updated) <- clustered_renamed_updated$RNA_snn_res.0.4
# Rename clusters
annotated_integrated <- Rename_Clusters(seurat_object = clustered_renamed_updated, new_idents = annotation_info$new_cluster_idents)

Idents(annotated_integrated) <- factor(x = Idents(annotated_integrated), levels = sort(levels(annotated_integrated)))

my_colors =  c("#E6B0C2","#FADBD8","#FFB5B5","thistle1",
               "#424949",
               "#ABEBC6",  "powderblue",
                "#76448A",
                "#9A7D0A","#C7CC8F", 
               "pink3", "#F1C41F", "#B7A4DB", 
               "#2E86C1", "#1C7F82",
                "#F1948A", "thistle3",  "darkgreen", "#873600", "red2", 
               
               "#4A235A", 
               
               "steelblue","#7EBDC2",
               "#F4D03F","#1B4F72","#CB4335",
             
             "red2")
(plot <- DimPlot(annotated_integrated, reduction = "umap", label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = my_colors)+
    umap_theme() )

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(clustering_dir, sprintf("%02d_UMAP_highres_annotated_%s.png", plot_number, highres.clusters)), width = 6, height = 4, plot)


# Save annotated clustered object -----------------------------------------
saveRDS(annotated_integrated, file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))



