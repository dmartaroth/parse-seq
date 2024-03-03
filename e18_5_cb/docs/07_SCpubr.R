# ## ######################################## ## #
#                      SCpubr                   #
# ## ######################################## ## #

# Vignette for SCpubr with application to this project


library(here)
source(here("e18_5_cb","docs","packages.R")) # load packages
source(here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here("e18_5_cb","docs","functions.R")) # load functions
source(here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "08_SCpubr"
SCpubr.dir <- file.path(figs, "08_SCpubr")
if (!dir.exists(SCpubr.dir)) {
  dir.create(SCpubr.dir)
}

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))
obj$cell_types <- Idents(obj)

plot_number <- 0
# Installation ------------------------------------------------------------
# 
# install.packages("SCpubr")
# 
# # Install CRAN packages.
# cran_packages <- c("assertthat",
#                    "circlize",
#                    "colorspace",
#                    "dplyr",
#                    "ggbeeswarm",
#                    "ggdist",
#                    "ggExtra",
#                    "ggnewscale",
#                    "ggplot2",
#                    "ggplotify",
#                    "ggrastr",
#                    "ggrepel",
#                    "ggridges",
#                    "ggsignif",
#                    "graphics",
#                    "magrittr",
#                    "patchwork",
#                    "pheatmap",
#                    "plyr",
#                    "rlang",
#                    "scales",
#                    "scattermore",
#                    "Seurat",
#                    "tibble",
#                    "tidyr",
#                    "forcats",
#                    "Matrix",
#                    "purrr",
#                    "stringr",
#                    "svglite",
#                    "viridis")
# 
# install.packages(cran_packages)
# 
# # Install bioconductor packages.
# bioconductor_packages <- c("AUCell",
#                            "ComplexHeatmap",
#                            "clusterProfiler",
#                            "enrichplot",
#                            "infercnv",
#                            "Nebulosa",
#                            "UCell")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(bioconductor_packages)
# 
# # Install github packages.
# github_packages <- c("ggsankey",
#                      "liana",
#                      "monocle3")
# 
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# 
# remotes::install_github(github_packages)
# 


# Dim plots ---------------------------------------------------------------

# Seurat's DimPlot.
p1 <- Seurat::DimPlot(obj)

# SCpubr's DimPlot.
p2 <- SCpubr::do_DimPlot(sample = obj)

p <- p1 | p2
p

convenient_save_plot(p, name = "seurat_vs_scpubr_dimplot", dir = SCpubr.dir,height = 8, width = 11)


## Modifying axes behavior -------------------------------------------------

# Example using PCA reduction.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         reduction = "pca")

# Example using a non-canonical set of dimensions.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         dims = c(2, 1))

p <- p1 | p2
p

convenient_save_plot(p, name = "DimPlot_SCpubr_axis mod", dir = SCpubr.dir,height = 8, width = 12)

# You can cancel this behavior by setting plot.axes = TRUE

# Bring back the Axes.
p <- SCpubr::do_DimPlot(sample = obj,
                        plot.axes = TRUE)
p


## Change colors -----------------------------------------------------------
colors <-  c("chondro.1" ="#E6B0C2","chondro.2" = "#FADBD8",
             "chondro.3" ="#FFB5B5","chondro.4" = "pink3","endocr"=  "thistle1",
                          "epith"= "#ABEBC6", "mes"= "powderblue",
                          "musc"= "red2", "myeloid"="#B7A4DB","neu.1"=  "#76448A",
                          "neu.2"=  "#F1948A","neu.3"= "thistle3",
                          "neu.4"= "#2E86C1", "neu.5"="#424949","olf"= "#9A7D0A",
                          "osteo.1"= "#1C7F82","osteo.2"= "steelblue", 
                          "osteocl"=  "#7EBDC2", "unknown.1"= "#F4D03F",
                          "unknown.2"=  "#C7CC8F","unknown.3"= "#1B4F72",
                          "vasc"= "#CB4335")


## Modify legend appearance ------------------------------------------------
# Modify the number of columns in the legend.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.ncol = 2)

# Modify the number of rows in the legend.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.nrow = 3)

p <- p1 | p2 
p

# Fill the legend by column.
p1 <- SCpubr::do_DimPlot(sample = obj,
                         legend.byrow = FALSE)

# Fill the legend by rows.
p2 <- SCpubr::do_DimPlot(sample = obj,
                         legend.byrow = TRUE)

p <- p1 | p2 
p

# Put labels on top of the clusters.
p <- SCpubr::do_DimPlot(sample = obj, 
                        label = TRUE)
p


# Labels as text
p <- SCpubr::do_DimPlot(sample = obj, 
                        label = TRUE,
                        label.box = FALSE)
p



# Change the color of the label text.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black")

# Change the color of the text.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black",
                         label.box = FALSE)
p <- p1 | p2
p


# Change the size of the label text.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black",
                         label.size = 4)

# Change the size of the text.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black",
                         label.box = FALSE,
                         label.size = 4)
p <- p1 | p2
p



# Repel the labels.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black",
                         repel = TRUE)

# Repel the text.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         label.color = "black",
                         label.box = FALSE,
                         repel = TRUE)
p <- p1 | p2
p





## Legend behavior ---------------------------------------------------------
# Remove the legend from the plot.
p <- SCpubr::do_DimPlot(sample = obj, 
                        label = TRUE, 
                        legend.position = "none")
p

# Top
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.position = "top")

p2 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.position = "bottom")

p3 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.position = "left")

p4 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.position = "right")

p <- (p1 | p2) / (p3 | p4)
p



# Add a legend title.
p <- SCpubr::do_DimPlot(sample = obj, 
                        legend.title = "Annotated clusters")
p



# Top
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.title = "My clusters",
                         legend.title.position = "top")

# Bottom
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.title = "My clusters",
                         legend.title.position = "bottom")

# Left
p3 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.title = "My clusters",
                         legend.title.position = "left")

# Right
p4 <- SCpubr::do_DimPlot(sample = obj, 
                         legend.title = "My clusters",
                         legend.title.position = "right")

p <- (p1 | p2) / (p3 | p4)
p

# Incorporating these:
p <- SCpubr::do_DimPlot(sample = obj, 
                        label = TRUE,
                        colors.use = colors,
                        label.color = "white",
                        label.fill = NULL,
                        repel = TRUE,
                        legend.title = "Annotated clusters")
p


convenient_save_plot(p, name = "DimPlot_SCpubr", dir = SCpubr.dir,height = 12, width = 10)



## Change the order of plotting --------------------------------------------
# Regular SCpubr DimPlot.
p1 <- SCpubr::do_DimPlot(sample = obj,
                         reduction = "pca",
                         plot.title = "Normal DimPlot",
                         legend.position = "none")

# Using order with one value and shuffle = TRUE.
p2 <- SCpubr::do_DimPlot(sample = obj,
                         shuffle = TRUE,
                         order = "5",
                         reduction = "pca",
                         plot.title = "shuffle = TRUE",
                         legend.position = "none")

# Using order with one value and shuffle = FALSE.
p3 <- SCpubr::do_DimPlot(sample = obj,
                         shuffle = FALSE,
                         order = "5",
                         reduction = "pca",
                         plot.title = "shuffle = FALSE",
                         legend.position = "none")

# Using order with all values.
p4 <- SCpubr::do_DimPlot(sample = obj,
                         shuffle = FALSE,
                         order = c("5", "8", "4",
                                   "9", "3", "1",
                                   "6", "0", "7", "2"),
                         reduction = "pca",
                         plot.title = "shuffle = FALSE all identities",
                         legend.position = "none")

p <- (p1 | p2) / (p3 | p4)
p


## Highlighting cells ------------------------------------------------------
cells.use <- WhichCells(obj, idents = "chondro.3")

# Compare Seurat and SCpubr way of highlighting cells.
p1 <- Seurat::DimPlot(obj, 
                      cells.highlight = cells.use)

p2 <- SCpubr::do_DimPlot(sample = obj,
                         cells.highlight = cells.use)

p <- p1 | p2
p



# Change color of highlighted and non-highlighted cells.
p <- SCpubr::do_DimPlot(sample = obj, 
                        cells.highlight = cells.use,
                        colors.use = "black",
                        na.value = "grey90")
p



# Increase the size of the highlighted cells.
p <- SCpubr::do_DimPlot(sample = obj, 
                        cells.highlight = cells.use, 
                        sizes.highlight = 1)
p



# Using cells.highlight.
p1 <- SCpubr::do_DimPlot(sample = obj, 
                         cells.highlight = cells.use)

# Using idents.highlight.
p2 <- SCpubr::do_DimPlot(sample = obj, 
                         idents.highlight = c("chondro.1"))

# Using both.
p3 <- SCpubr::do_DimPlot(sample = obj, 
                         cells.highlight = cells.use, 
                         idents.highlight = c("chondro.1"))

p <- p1 | p2 | p3
p



## Restrict the identities displayed ---------------------------------------
# Subset desired identities in a DimPlot.
p <- SCpubr::do_DimPlot(sample = obj[, obj$cell_types %in% c("chondro.1","chondro.2", "chondro.3")],colors.use = colors)

p


convenient_save_plot(p, name = "Selective_DimPlot_chondro.1_2_3_SCpubr", dir = SCpubr.dir,height = 10, width = 10)

# The approach above loses UMAP silhouette
# Select identities with idents.keep.
p1 <- SCpubr::do_DimPlot(sample = obj,
                         idents.keep = c("chondro.1","chondro.2", "chondro.3"))

# Also, non-selected cells's color can be modified.
p2 <- SCpubr::do_DimPlot(sample = obj,
                         idents.keep = c("chondro.1","chondro.2", "chondro.3"),
                         na.value = "grey50")
p <- p1 | p2
p


## Group by another metadata variable --------------------------------------

# Group by another metadata variable.
p1 <- SCpubr::do_DimPlot(obj, 
                         group.by = "lowres_Idents",
                         legend.title = "low res clusters")

p2 <- SCpubr::do_DimPlot(obj, 
                         group.by = "cell_types",
                         legend.position = "right",colors.use = colors)

p <- p1 | p2
p


convenient_save_plot(p, name = "lowres_highres_clustering_SCpubr", dir = SCpubr.dir,height = 9, width = 12)


## Splitting by a category -------------------------------------------------
# Seurat's DimPlot using split.by
p <- Seurat::DimPlot(obj, 
                     split.by = "cell_types", 
                     ncol = 6,cols = colors)

p

# SCpubr's DimPlot using split.by
p <- SCpubr::do_DimPlot(obj, 
                        split.by = "cell_types", 
                        ncol = 5, 
                        legend.position = "none",
                        font.size = 12, colors.use = colors)

p


convenient_save_plot(p, name = "split_clusters_SCpubr", dir = SCpubr.dir,height = 20, width = 20)


# Using split.by and restricting the number of output plots with idents.keep.
p <- SCpubr::do_DimPlot(obj, 
                        split.by = "cell_types", 
                        ncol = 3, 
                        idents.keep = c("chondro.1","chondro.2","chondro.3"),
                        legend.position = "none",
                        font.size = 12)

p



## Group by a variable but split by another --------------------------------

# Using split.by and group.by in combination.


p <- SCpubr::do_DimPlot(obj, 
                        group.by = "cell_types",
                        split.by = "lowres_Idents", 
                        font.size = 12)

p






# Feature Plots -----------------------------------------------------------

# Seurat's Feature Plot.
p1 <- Seurat::FeaturePlot(obj, 
                          features = "Acan")

# SCpubr's Feature Plot.
p2 <- SCpubr::do_FeaturePlot(sample = obj,
                             features = "Acan")

p <- p1 | p2
p

# Use case with PCA embedding.
p1 <- SCpubr::do_FeaturePlot(sample = obj, 
                             features = "Acan",
                             plot.title = "Plotting PCA coordinates",
                             reduction = "pca")

# Use case with non-canonical dimensions.                             
p2 <- SCpubr::do_FeaturePlot(sample = obj, 
                             features = "Acan",
                             plot.title = "Plotting UMAP coordinates in a different order",
                             dims = c(2, 1))

p <- p1 | p2
p


## Multiple features -------------------------------------------------------

p <- SCpubr::do_FeaturePlot(obj, features = c("nCount_RNA", 
                                                 "nFeature_RNA", 
                                                 "percent.mito", 
                                                 "Runx2"), 
                            plot.title = "A collection of features", 
                            ncol = 2)

p


## Working with subsets of cells -------------------------------------------

# Useful if a certain subset of cells drives the ends of the color scales
# Let's say the unknown clusters are responsible
# Remove them with the following command completely

cells.plot <- colnames(obj[, !(obj$cell_types %in% c("unknown.1","unknown.2","unknown.3"))])

p <- SCpubr::do_FeaturePlot(obj[, cells.plot], 
                            features = c("Ankrd11"))

p

# But lose silhouette of the plot

# Using cells.highlight parameter to select the cells we want to include in the plot.
p <- SCpubr::do_FeaturePlot(sample = obj, 
                            cells.highlight = cells.plot, 
                            features = c("Ankrd11"))

p



# Selecting given identitites to include in the plot.
p <- SCpubr::do_FeaturePlot(sample = obj, 
                            idents.highlight = levels(obj)[!(levels(obj) %in% c("unknown.1","unknown.2","unknown.3","neu.1","neu.2","neu.3","neu.4","neu.5","endocrine"))], 
                            features = c("Col10a1"),
                            split.by = "genotype")

p



## Splitting the FeaturePlot by a variable ---------------------------------

# Group clusters into three values for visualization purposes.
obj$grouped_clusters <- as.character(obj$cell_types)
obj$grouped_clusters[obj$grouped_clusters %in% c("chondro.1","chondro.2","chondro.3","chondro.4")] <- "chondro"
obj$grouped_clusters[obj$grouped_clusters %in% c("mes","osteo.1","osteo.2","musc")] <- "osteo"
obj$grouped_clusters[obj$grouped_clusters %in% c("neu.1","neu.2","neu.3","neu.4","neu.5","endocr","olf")] <- "neu"
obj$grouped_clusters[obj$grouped_clusters %in% c("myeloid","osteocl")] <- "imm"
obj$grouped_clusters[obj$grouped_clusters %in% c("unknown.1","unknown.2","unknown.3")] <- "unknown"

# Seurat Feature Plot using split.by.
p <- Seurat::FeaturePlot(obj, 
                         features = "Acan", 
                         split.by = "grouped_clusters",
                         ncol = 2)
p


# SCpubr Feature Plot using split.by
p <- SCpubr::do_FeaturePlot(sample = obj, 
                            features = "Acan", 
                            split.by = "grouped_clusters")

p



## Subset the color scale to a min and max ---------------------------------

# Use min.cutoff and max.cutoff.
p1 <- SCpubr::do_FeaturePlot(obj, 
                             features = c("Col10a1"))

p2 <- SCpubr::do_FeaturePlot(obj, 
                             features = c("Col10a1"),
                             min.cutoff = 1,
                             max.cutoff = 2)
p <- p1 | p2
p




## Apply symmetrical color scales ------------------------------------------

# Enforce two-end symmetrical color scale.
p1 <- SCpubr::do_FeaturePlot(obj, 
                             features = "PC_1",
                             enforce_symmetry = FALSE)

p2 <- SCpubr::do_FeaturePlot(obj, 
                             features = "PC_1",
                             enforce_symmetry = TRUE)

p <- p1 | p2
p



# Nebulosa plots ----------------------------------------------------------

p <- SCpubr::do_NebulosaPlot(sample = obj, 
                             features = "Col10a1")
p


convenient_save_plot(p, name = "Col10a1_density_SCpubr", dir = SCpubr.dir,height = 7, width = 7)


# Use alongside featureplot to see expression and density of surrounding cells
p1 <- SCpubr::do_FeaturePlot(sample = obj, 
                             features = "Col10a1")  + 
  facet_grid(.~obj$genotype)

p2 <- SCpubr::do_NebulosaPlot(sample = obj, 
                              features = "Col10a1") + 
  facet_grid(.~obj$genotype)
p <- p1 | p2
p

convenient_save_plot(p, name = "Col2a1_featureplot_density_bygenotype_SCpubr", dir = SCpubr.dir,height = 6, width = 15)


# Multiple features at the same time
p <- SCpubr::do_NebulosaPlot(obj, 
                             features = c("Col2a1", "Col10a1"))
p




## Compute joint densities -------------------------------------------------

p <- SCpubr::do_NebulosaPlot(sample = obj, 
                             features = c("Col10a1","Sp7"), 
                             joint = TRUE)+ 
  facet_grid(.~obj$genotype)
p 

convenient_save_plot(p, name = "Col10a1_Sp7_joint_density_bygenotype_SCpubr", dir = SCpubr.dir,height = 6, width = 16)


features.use <- c("Col10a1", "Sp7")

p <- SCpubr::do_NebulosaPlot(sample = obj, 
                             features = features.use, 
                             joint = TRUE, 
                             return_only_joint = TRUE,
                             plot.title = "Joint density Col10a1+ Sp7+")+ 
  facet_grid(.~obj$genotype)

p




convenient_save_plot(p, name = "Col10a1_Sp7_joint_density_bygenotype_SCpubr", dir = SCpubr.dir,height = 6, width = 10)


# Bee Swarm plots ---------------------------------------------------------

# Rank cells in a given variable, which must be continuous
# Cells are grouped into another variable of interest and displayed in a scatter plot

# Compute enrichment.
gene_list <- c("Col10a1", "Col2a1")
obj <- Seurat::AddModuleScore(obj, features = list(gene_list), name = "testing_list")

# Rank the enrichment scores.
obj$rank <- rank(obj$testing_list1)

# Visualize the two distribution.
p1 <- SCpubr::do_ViolinPlot(sample = obj,
                            feature = "testing_list1",
                            group.by = "genotype")

p2 <- obj@meta.data %>% # Extract metadata
  dplyr::mutate(cell_name = rownames(obj@meta.data)) %>% # Get the cell names.
  dplyr::select(cell_name, rank) %>% # Select the columns to plot.
  dplyr::arrange(rank) %>% # Reorder the rows.
  dplyr::mutate(cell_name = factor(cell_name, levels = cell_name)) %>% # Convert to factor for plotting.
  ggplot2::ggplot(mapping = ggplot2::aes(x = cell_name, y = rank)) +
  ggplot2::geom_point() +
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank()) + 
  ggplot2::xlab("Cell name") + 
  ggplot2::ylab("Rank")

p <- p1 | p2
p

# further subset distribution by another variable
p <- obj@meta.data %>% # Extract metadata
  dplyr::mutate(cell_name = rownames(obj@meta.data)) %>% # Get the cell names.
  dplyr::select(cell_name, rank, grouped_clusters) %>% # Select the columns to plot.
  dplyr::arrange(rank) %>% # Reorder the rows.
  dplyr::mutate(cell_name = factor(cell_name, levels = cell_name)) %>% # Convert to factor for plotting.
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$cell_name, y = .data$rank)) +
  ggplot2::geom_point() +
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank()) + 
  ggplot2::xlab("Cell name") + 
  ggplot2::ylab("Rank") + 
  ggplot2::facet_wrap("grouped_clusters", ncol = 4)
p

# Now we see not all clusters have the same distribution or ranks
# Bee Swarm Plot aims to plot these distributions in a nice way


## Using categorical variables ---------------------------------------------

p1 <- SCpubr::do_DimPlot(obj, 
                         reduction = "pca", group.by = "grouped_clusters",
                         label = TRUE, 
                         legend.position = "none", 
                         dims = c(1, 2)) 
p2 <- SCpubr::do_DimPlot(obj, 
                         reduction = "pca",group.by = "grouped_clusters", 
                         label = TRUE, 
                         legend.position = "none",
                         dims = c(3, 4)) 

p <- p1 | p2
p

# Already see that some clusters separate on PC_1 from the rest

p1 <- SCpubr::do_DimPlot(sample = obj, 
                         reduction = "pca", 
                         group.by = "cell_types",
                         label = TRUE, 
                         repel = TRUE,
                         label.fill = NULL,
                         legend.position = "none",
                         dims = c(1, 2), colors.use = colors)

p2 <- SCpubr::do_DimPlot(sample = obj, 
                         reduction = "pca", 
                         group.by = "cell_types",
                         label = TRUE, 
                         repel = TRUE,
                         label.fill = NULL,
                         legend.position = "none",
                         dims = c(3, 4), colors.use = colors) 

p3 <- SCpubr::do_BeeSwarmPlot(sample = obj, 
                              feature_to_rank = "PC_1", 
                              group.by = "cell_types", 
                              continuous_feature = FALSE,
                              legend.position = "none", colors.use = colors)

p4 <- SCpubr::do_BeeSwarmPlot(sample = obj, 
                              feature_to_rank = "PC_2",
                              group.by = "cell_types", 
                              continuous_feature = FALSE,
                              legend.position = "none",colors.use = colors)

p <- (p1 | p3) / (p2 | p4)
p




## Using continuous variables ----------------------------------------------


