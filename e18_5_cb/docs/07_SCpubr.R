# ## ######################################## ## #
#                      SCpubr                   #
# ## ######################################## ## #

# Vignette for SCpubr with application to this project


library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

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
cran_packages <- c("assertthat",
                   "circlize",
                   "colorspace",
                   "dplyr",
                   "ggbeeswarm",
                   "ggdist",
                   "ggExtra",
                   "ggnewscale",
                   "ggplot2",
                   "ggplotify",
                   "ggrastr",
                   "ggrepel",
                   "ggridges",
                   "ggsignif",
                   "graphics",
                   "magrittr",
                   "patchwork",
                   "pheatmap",
                   "plyr",
                   "rlang",
                   "scales",
                   "scattermore",
                   "Seurat",
                   "tibble",
                   "tidyr",
                   "forcats",
                   "Matrix",
                   "purrr",
                   "stringr",
                   "svglite",
                   "viridis")

lapply(cran_packages, library, character.only = TRUE)
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

# Set up list of a genes to compute enrichment. 
genes.use <- c("Mmp13","Col10a1")

# Compute enrichment and rename the output.
obj <- Seurat::AddModuleScore(obj, 
                                 features = genes.use, 
                                 name = "Hypertrophic_chondrocytes")
obj$Hypertrophic_chondrocytes <- obj$Hypertrophic_chondrocytes1
obj$Hypertrophic_chondrocytes1 <- NULL

p1 <- SCpubr::do_DimPlot(sample = obj, 
                         label = TRUE, 
                         legend.position = "none")

p2 <- SCpubr::do_FeaturePlot(sample = obj, 
                             features = "Hypertrophic_chondrocytes",
                             legend.title = "Hypertrophic chondrocytes") 

p3 <- SCpubr::do_BeeSwarmPlot(sample = obj, 
                              feature_to_rank = "Hypertrophic_chondrocytes", 
                              group.by = "cell_types", 
                              continuous_feature = TRUE,
                              legend.title = "Hypertrophic chondrocytes")
p <- p1 | p2 | p3
p


SCpubr::do_BeeSwarmPlot(sample = obj, 
                        feature_to_rank = "Hypertrophic_chondrocytes", 
                        group.by = "cell_types",
                        continuous_feature = TRUE,
                        legend.title = "Hypertrophic chondrocytes")+ 
  facet_grid(.~obj$genotype)




convenient_save_plot(p, name = "Col10a1_Mmp13_Ranking_bygenotype_SCpubr", dir = SCpubr.dir,height = 6, width = 10)


# Violin plots ------------------------------------------------------------
interestinggenes <- c("Col10a1","Mmp13","Alpl","Tnfsf11","Igf1","Fgfr3")

p <- Seurat::VlnPlot(obj, 
                     features = "Mmp13", cols = iss.colors)
p



# "Surgically" add the alpha parameter in the ggplot2 object.
p$layers[[2]]$aes_params$alpha <- 0.05
p


# Basic violin plot.
p <- SCpubr::do_ViolinPlot(sample = obj, 
                           features = "Mmp13",colors.use = colors)
p


# Remove the box plots.
p <- SCpubr::do_ViolinPlot(sample = obj, 
                           features = "Mmp13",
                           plot_boxplot = FALSE,colors.use = colors)
p


# Rotate x axis labels.
p <- SCpubr::do_ViolinPlot(sample = obj, 
                            features = "Mmp13",
                           colors.use = colors,
                           rotate_x_axis_labels = 45)
p



## Violin plots as a means of QC -------------------------------------------
# Add horizontal lines.
p <- SCpubr::do_ViolinPlot(sample = obj, 
                           features = "percent.mito", 
                           colors.use = colors,
                           y_cut = 0.06)
p



## Modifying aesthetics ----------------------------------------------------
# Increase line width.
p1 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "percent.mito",
                            legend.position = "none")

p2 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "percent.mito",
                            line_width = 2)

p <- p1 / p2
p


# Decrease boxplot width.
p1 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "percent.mito",
                            colors.use = colors,
                            legend.position = "none")

p2 <- SCpubr::do_ViolinPlot(sample = obj,
                            features = "percent.mito",
                            colors.use = colors,
                            boxplot_width =  0.1)

p <- p1 / p2
p



## Force the same limits on different violin plots -------------------------
# Share the same Y axis.
p <- SCpubr::do_ViolinPlot(sample = obj,
                           features = c("nCount_RNA", "nFeature_RNA"),
                           ncol = 1,
                           colors.use = colors,
                           share.y.lims = TRUE,
                           legend.position = "none")
p



## Split by another variable -----------------------------------------------

# Split violin plots.


genotypecolors <- c("Bmp2_ctrl" = "#FFB6C1","Bmp2_ncko"= "#ADD8E6")

p<- SCpubr::do_ViolinPlot(sample = obj,
                          features = "Mmp13",
                          split.by = "genotype",
                          plot_boxplot = FALSE,
                          legend.position = "bottom")

p


# Ridge plots -------------------------------------------------------------


# Compute the most basic ridge plot.
p <- SCpubr::do_RidgePlot(sample = obj,
                          feature = "percent.mito")
p


# Use continuous color scale.
p1 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "percent.mito",
                           continuous_scale = TRUE,
                           viridis.direction = 1,
                           legend.position = "none")

p2 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "percent.mito",
                           continuous_scale = TRUE,
                           viridis.direction = -1)

p <- p1 / p2
p



# Plot quantiles of the distribution --------------------------------------

## Draw quantiles of the distribution.
p1 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "nFeature_RNA",
                           group.by = "genotype",
                           continuous_scale = TRUE,
                           compute_quantiles = TRUE,
                           compute_custom_quantiles = TRUE,
                           legend.position = "right")

p2 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "nFeature_RNA",
                           group.by = "genotype",
                           continuous_scale = TRUE,
                           compute_quantiles = TRUE,
                           compute_custom_quantiles = TRUE,
                           quantiles = c(0.1, 0.5, 0.75),
                           legend.position = "right")

p <- p1 / p2
p



## Compute probability tails -----------------------------------------------

# Draw probability tails.
p1 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "percent.mito",
                           group.by = "genotype",
                           legend.position = "right",
                           continuous_scale = TRUE,
                           compute_quantiles = TRUE,
                           compute_distribution_tails = TRUE)

p2 <- SCpubr::do_RidgePlot(sample = obj,
                           feature = "percent.mito",
                           group.by = "genotype",
                           legend.position = "right",
                           continuous_scale = TRUE,
                           compute_quantiles = TRUE,
                           compute_distribution_tails = TRUE,
                           prob_tails = 0.3)

p <- p1 / p2
p



## Compute probability densities -------------------------------------------

# Draw probability tails.
p <- SCpubr::do_RidgePlot(sample = obj,
                          feature = "percent.mito",
                          group.by = "genotype",
                          legend.position = "right",
                          continuous_scale = TRUE,
                          compute_quantiles = TRUE,
                          color_by_probabilities = TRUE)
p


# Dot plots ---------------------------------------------------------------

# Seurat's dot plot.
p <- Seurat::DotPlot(obj, 
                     features = interestinggenes)
p


# SCpubr's dot plot.
p <- SCpubr::do_DotPlot(sample = obj, 
                        use_viridis = TRUE,
                        viridis.palette = "cividis",
                        features = interestinggenes)
p


# Query multiple genes as a named list
genes <- list("Proliferation" = c("Sox9","Ptch1","Fgfr3","Igf1"),
              "Differentiation" = c("Col2a1","Acan","Sox5","Runx2"),
              "Maturation" = c("Col10a1","Vegfa","Pthlh","Mmp13"))

p <- SCpubr::do_DotPlot(sample = obj,  
                        use_viridis = TRUE,
                        viridis.palette = "E",
                        viridis.direction = 1,
                        features = genes)
p



## Clustering the identities -----------------------------------------------

p1 <- SCpubr::do_DotPlot(sample = obj,  
                         use_viridis = TRUE,
                         viridis.palette = "E",
                         viridis.direction = 1,
                         features = genes,
                         plot.title = "Not clustered",
                         legend.position = "right")

p2 <- SCpubr::do_DotPlot(sample = obj,  
                         use_viridis = TRUE,
                         viridis.palette = "E",
                         viridis.direction = 1, 
                         features = genes, 
                         cluster =  TRUE, 
                         plot.title = "Clustered",
                         legend.position = "none")

p <- p1 / p2
p


convenient_save_plot(p, name = "dotplot_chondrocyte_state_SCpubr", dir = SCpubr.dir,height = 12, width = 7)



## Inverting the axes ------------------------------------------------------

p1 <- SCpubr::do_DotPlot(sample = obj,  
                         use_viridis = TRUE,
                         viridis.palette = "E",
                         viridis.direction = 1,
                         features = genes,
                         plot.title = "Not clustered",
                         legend.position = "right",
                         flip = TRUE)

p2 <- SCpubr::do_DotPlot(sample = obj,  
                         use_viridis = TRUE,
                         viridis.palette = "E",
                         viridis.direction = 1,
                         features = genes,
                         plot.title = "Clustered",
                         legend.position = "none",
                         flip = TRUE,
                         rotate_x_axis_labels = 45)
p <- p1 | p2
p


convenient_save_plot(p, name = "dotplot_chondrocyte_state_flipped_SCpubr", dir = SCpubr.dir,height = 7, width = 14)

# Use to characterize cell states between control and mutant
SCpubr::do_DotPlot(sample = obj,  
                   use_viridis = TRUE,
                   group.by =c("genotype"),
                   viridis.palette = "E",
                   viridis.direction = 1, 
                   scale = FALSE,
                   dot.scale = 10,
                   features = genes, 
                   cluster =  TRUE, 
                   plot.title = "Chondrocyte differentiation",
                   legend.position = "right")

# Bar plots ---------------------------------------------------------------

# Basic bar plot, horizontal.
p1 <- SCpubr::do_BarPlot(sample = obj, 
                         group.by = "cell_types", 
                         colors.use = colors,
                         legend.position = "none", 
                         plot.title = "Number of cells per cluster")

# Basic bar plot, vertical.
p2 <- SCpubr::do_BarPlot(sample = obj, 
                         group.by = "cell_types", 
                         colors.use = colors,
                         legend.position = "none", 
                         plot.title = "Number of cells per cluster", 
                         flip = TRUE)
p <- p1 | p2
p



## Grouping by a second variable -------------------------------------------
# Split by a second variable.
p1 <- SCpubr::do_BarPlot(sample = obj, 
                         group.by = "cell_types", 
                         split.by = "genotype",
                         colors.use = colors,
                         plot.title = "Number of cells per cluster in each sample",
                         position = "stack",
                         legend.position = "left",
                         legend.ncol = 1) + 
  theme(plot.title = element_text(size = 10), legend.title = element_text(size = 10))

p2 <- SCpubr::do_BarPlot(sample = obj, 
                         group.by = "genotype",  
                         split.by = "cell_types",
                         colors.use = genotypecolors,
                         plot.title = "Number of cells per sample in each cluster",
                         position = "stack",
                         legend.position = "bottom")+ 
  theme(plot.title = element_text(size = 10), legend.title = element_text(size = 10))
p <- p1 | p2
p


# Position stack and fill with and without split.by.
p1 <- SCpubr::do_BarPlot(sample = obj, 
                         colors.use = colors,
                         group.by = "cell_types", 
                         plot.title = "Without split.by - position = stack",
                         position = "stack",
                         flip = FALSE,
                         legend.position = "none")

p2 <- SCpubr::do_BarPlot(sample = obj, 
                         colors.use = colors,
                         group.by = "cell_types", 
                         plot.title = "Without split.by - position = fill",
                         flip = FALSE,
                         legend.position = "none")

p3 <- SCpubr::do_BarPlot(sample = obj, 
                         colors.use = colors,
                         group.by = "cell_types", 
                         split.by = "grouped_clusters",
                         plot.title = "With split.by - position = stack",
                         position = "stack",
                         flip = FALSE,legend.position = "none")

p4 <- SCpubr::do_BarPlot(sample = obj, 
                         colors.use = colors,
                         group.by = "cell_types", 
                         split.by = "grouped_clusters",
                         plot.title = "With split.by - position = fill",
                         position = "fill",
                         flip = FALSE,
                         legend.position = "right")
p <- (p1 | p2) / (p3 | p4)
p 


# Box plots ---------------------------------------------------------------

# Basic box plot.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA")
p


# Use custom grouping.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = genotypecolors,
                        feature = "nCount_RNA",
                        group.by = "genotype")
p

# Flip the box plot.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA",
                        flip = TRUE)
p


## Modify aesthetic style --------------------------------------------------
# Use silhouette style.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA",
                        use_silhouette = TRUE)
p


## Reorder by mean values --------------------------------------------------
# Order by mean values.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA",
                        order = TRUE)
p

p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "percent.mito",
                        order = TRUE)
p

# Apply second grouping.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = genotypecolors,
                       feature = "nCount_RNA",
                        split.by = "genotype")
p



## Apply statistical tests to compare groups -------------------------------
# Apply statistical tests.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA",
                        use_test = TRUE,
                        comparisons = list(c("chondro.2","chondro.3"),
                                           c("osteo.1","osteo.2"),
                                           c("chondro.1","osteo.1")))
p


# Apply statistical tests and show the p-value.
p <- SCpubr::do_BoxPlot(sample = obj,
                        colors.use = colors,
                        feature = "nCount_RNA",
                        use_test = TRUE,
                        comparisons = list(c("chondro.2","chondro.3"),
                                           c("osteo.1","osteo.2"),
                                           c("chondro.1","osteo.1")),
                        map_signif_level = FALSE)
p


# Geyser plots ------------------------------------------------------------

# Geyser plot with categorical color scale.
p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = colors,
                            feature = "nCount_RNA",
                            scale_type = "categorical")

# Geyser plot with continuous color scale.
p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = colors,
                            feature = "nCount_RNA",
                            scale_type = "continuous")


p <- p1 / p2
p


# Geyser plot with categorical color scale without ordering by mean.
p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = colors,
                            feature = "nCount_RNA",
                            scale_type = "categorical",
                            order = FALSE)

# Geyser plot with continuous color scale without ordering by mean.
p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = colors,
                            feature = "nCount_RNA",
                            scale_type = "continuous",
                            order = FALSE)


p <- p1 / p2
p



## Plotting symmetrical scales ---------------------------------------------
# For plotting a continuous variable that spans both positive and negative values

# Geyser plot with continuous color scale.
p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "PC_1",
                            scale_type = "continuous",
                            enforce_symmetry = FALSE)

# Geyser plot with continuous and symmetrical color scale.
p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "PC_1",
                            scale_type = "continuous",
                            enforce_symmetry = TRUE)


p <- p1 / p2
p


# Select the groups displayed on the X axis -------------------------------

# Geyser plot with categorical color scale default X axis grouping.
p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = colors,
                            features = "Alpl",
                            scale_type = "categorical",
                            group.by = NULL,
                            xlab = "cell_types")

# Geyser plot with categorical color scale and custom grouping.
p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Alpl",
                            scale_type = "categorical",
                            group.by = "genotype",
                            colors.use = genotypecolors,
                            xlab = "Control and mutant individual")


p <- p1 / p2
p


## Split the plots by groups -----------------------------------------------
# Geyser plot with categorical color scale split by seurat clusters.
p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = genotypecolors,
                            features = "Alpl",
                            scale_type = "categorical",
                            group.by = "genotype",
                            split.by = "cell_types")

# Geyser plot with continuous color scale split by seurat clusters.
p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            colors.use = genotypecolors,
                            features = "Alpl",
                            scale_type = "continuous",
                            group.by = "genotype",
                            split.by = "cell_types")


p <- p1 / p2
p



# Geyser plot with different jitter.
p0 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.01,
                            legend.position = "none")

p1 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.1,
                            legend.position = "none")

p2 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.2,
                            legend.position = "none")

p3 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.3,
                            legend.position = "none")

p4 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.4,
                            legend.position = "none")

p5 <- SCpubr::do_GeyserPlot(sample = obj,
                            features = "Xist",
                            colors.use = colors,
                            scale_type = "continuous",
                            enforce_symmetry = TRUE,
                            jitter = 0.49)


p <- p0 / p1 / p2 / p3 / p4 / p5
p


# Alluvial plots ----------------------------------------------------------

# Used to visualize how cells "flow" from a given group to another
# For ex how cells from each cluster distribute across the individ datasets

# Compute basic sankey plot.
p <- SCpubr::do_AlluvialPlot(sample = obj,
                             first_group = "genotype",
                             last_group = "cell_types",
                             colors.use = colors,
                             legend.position = "bottom")

p+theme(text = element_text(size = 5))



convenient_save_plot(p, name = "alluvial_genotype_to_celltype_SCpubr", dir = SCpubr.dir,height = 16, width = 7)


# Compute basic sankey plot.
p <- SCpubr::do_AlluvialPlot(sample = obj,
                             first_group = "genotype",
                             last_group = "cell_types",
                             colors.use = colors,
                             flip = TRUE,
                             legend.position = "top")

p


convenient_save_plot(p, name = "alluvial_genotype_to_celltype_flipped_SCpubr", dir = SCpubr.dir,height = 8, width = 17)



## Add more groups ---------------------------------------------------------

# Compute basic sankey plot.
p <- SCpubr::do_AlluvialPlot(sample = obj,
                             first_group = "genotype",
                             middle_groups = "grouped_clusters",
                             last_group = "cell_types",
                             colors.use = colors,legend.position = "bottom")

p


convenient_save_plot(p, name = "alluvial_genotype_to_groupedcluster_to_celltype_SCpubr", dir = SCpubr.dir,height = 16, width = 7)


# Control overplotting.
p1 <- SCpubr::do_AlluvialPlot(sample = obj,
                              first_group = "genotype",
                              middle_groups = "grouped_clusters",
                              last_group = "cell_types",
                              colors.use = colors,legend.position = "none",
                              use_labels = FALSE)

p2 <- SCpubr::do_AlluvialPlot(sample = obj,
                              first_group = "genotype",
                              middle_groups = "grouped_clusters",
                              last_group = "cell_types",
                              colors.use = colors,legend.position = "none",
                              use_labels = TRUE)

p3 <- SCpubr::do_AlluvialPlot(sample = obj,
                              first_group = "genotype",
                              middle_groups = "grouped_clusters",
                              last_group = "cell_types",
                              colors.use = colors,legend.position = "none",
                              use_labels = FALSE,
                              repel = TRUE)

p4 <- SCpubr::do_AlluvialPlot(sample = obj,
                              first_group = "genotype",
                              middle_groups = "grouped_clusters",
                              last_group = "cell_types",
                              colors.use = colors,legend.position = "none",
                              use_labels = TRUE,
                              repel = TRUE)

p <- (p1 | p2) / (p3 | p4)
p


convenient_save_plot(p, name = "alluvial_genotype_to_celltype_relabeled_SCpubr", dir = SCpubr.dir,height = 18, width = 18)



# Color by another column.
p <- SCpubr::do_AlluvialPlot(sample = obj,
                             first_group = "genotype",
                             middle_groups = "grouped_clusters",
                             last_group = "cell_types",
                             fill.by = "genotype",
                             use_labels = TRUE,
                             repel = TRUE)

p


# Use custom colors for borders and fill
p <- SCpubr::do_AlluvialPlot(sample = obj,
                             first_group = "genotype",
                             middle_groups = "grouped_clusters",
                             last_group = "cell_types",
                             fill.by = "genotype",
                             use_labels = TRUE,
                             repel = TRUE,
                             stratum.color = "goldenrod2",
                             alluvium.color = "black",
                             stratum.fill = "floralwhite")

p


# Chord diagram plots -----------------------------------------------------


p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",
                                 big.gap = 40)

p

p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",
                                 small.gap = 5)

p


## Control alignment of the diagram ----------------------------------------
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",
                                 big.gap = 40,
                                 alignment = "horizontal")

p


## Control the directions of the links -------------------------------------

# We need to set direction.type to diffHeight only as arrows are, by nature, directional.
# Without any direction
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",colors.from = colors,
                                 colors.to = genotypecolors,
                                 directional = 0,
                                 direction.type = "diffHeight")

p


# With links going "from" to "to"
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",colors.from = colors,
                                 colors.to = genotypecolors,
                                 directional = 1)

p


# With links going "to" to "from"



p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",colors.from = colors,
                                 colors.to = genotypecolors,
                                 directional = -1)

p

# With links in both directions

p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",colors.from = colors,
                                 colors.to = genotypecolors,
                                 directional = 2,
                                 direction.type = "diffHeight")

p



## Add padding to the labels -----------------------------------------------

p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "genotype",
                                 colors.from = colors,
                                 colors.to = genotypecolors,
                                 padding_labels = 8)

p


## Scale the nodes ---------------------------------------------------------

p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 to = "cell_types",
                                 from = "genotype",
                                 colors.to = colors,
                                 colors.from = genotypecolors,
                                 scale = TRUE,
                                 padding_labels = 8)

p


## Self linking -----------------------------------------------------------
# Prevent self linking.
obj$cell_types2 <- obj$cell_types
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "cell_types2",
                                 self.link = 1,
                                 scale = TRUE)

p

# With self linking.
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 from = "cell_types",
                                 to = "cell_types2",
                                 self.link = 2,
                                 scale = TRUE)

p



## Control appearance of arrows --------------------------------------------

# set triangle arrows
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 to = "cell_types",
                                 from = "genotype",
                                 colors.to = colors,
                                 colors.from = genotypecolors,
                                 padding_labels = 8,
                                 link.arr.type = "triangle")

p

# Set big arrows
p <- SCpubr::do_ChordDiagramPlot(sample = obj,
                                 to = "cell_types",
                                 from = "genotype",
                                 colors.to = colors,
                                 colors.from = genotypecolors,
                                 padding_labels = 8,
                                 link.arr.type = "big.arrow")

p





# Volcano plots -----------------------------------------------------------
Idents(obj) <- obj$genotype
de_genes <- Seurat::FindMarkers(obj,ident.1 = "Bmp2_ctrl",ident.2 = "Bmp2_ncko",
                                logfc.threshold = 0.25,only.pos = FALSE)


# Generate a volcano plot.
p <- SCpubr::do_VolcanoPlot(sample = obj, de_genes = de_genes)
p


# Modify cutoffs
p <- SCpubr::do_VolcanoPlot(sample = obj, de_genes = de_genes,
                            pval_cutoff = 1e-100,
                            FC_cutoff = 5)
p


# Modify number of gene tags

p <- SCpubr::do_VolcanoPlot(sample = obj, de_genes = de_genes,
                            n_genes = 15)
p



# Group-wise DE analysis plots --------------------------------------------

obj <- obj

# Set identities
Seurat::Idents(obj) <- obj$genotype

# Compute DE genes and transform to a tibble
de_genes <- tibble::tibble(Seurat::FindAllMarkers(object = obj,min.pct = 0.25, min.cells.feature = 5))


# Default output.
p <- SCpubr::do_GroupwiseDEPlot(sample = obj,
                                de_genes = de_genes,
                                group.by = "grouped_clusters",
                                min.cutoff = 1,
                                top_genes = 50,
                                legend.position = "bottom",
                                font.size = 7,
                                legend.length = 4,
                                border.color = "white"
                                )

p

# Add other rows of gene expression

p <- SCpubr::do_GroupwiseDEPlot(sample = obj,
                                de_genes = de_genes,
                                group.by = c("grouped_clusters","genotype"),
                                min.cutoff = 1,
                                top_genes = 50,
                                legend.position = "bottom",
                                font.size = 7,
                                legend.length = 4,
                                border.color = "white"
)

p


# Modify color scales

p <- SCpubr::do_GroupwiseDEPlot(sample = obj,
                                de_genes = de_genes,
                                group.by = c("grouped_clusters","genotype"),
                                min.cutoff = 1,
                                top_genes = 50,
                                use_viridis = TRUE,
                                viridis.palette.expression = "F",
                                viridis.palette.logfc = "mako",
                                viridis.palette.pvalue = "D",
                                legend.position = "bottom",
                                font.size = 7,
                                legend.length = 4,
                                border.color = "white"
)

p




# Term enrichment plots ---------------------------------------------------

# Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to Human genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})

# Choose the dataset to query against.
dbs_use <- c("GO_Biological_Process_2021", 
             "GO_Cellular_Component_2021")

# List of genes to use as input.
genes <- c(genes.use)

# Retrieve the enriched terms.
enriched_terms <- enrichR::enrichr(genes, dbs_use)

# Default plot.
p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
p


# Increased number of terms.
p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                   nterms = 15)
p


# Control the length of the terms.
p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                    nterms = 15)
p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                    nterms = 15,
                                    nchar_wrap = 30)
p <- p1 / p2
p


# Modify font size of the terms.
p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                    text_labels_size = 6)

p <- p1 / p2
p


# Expression heatmaps -----------------------------------------------------


## Single grouping variable ------------------------------------------------
genes <- list("Proliferation" = c("Sox9","Ptch1","Fgfr3","Igf1"),
              "Differentiation" = c("Col2a1","Acan","Sox5","Runx2"),
              "Maturation" = c("Col10a1","Vegfa","Pthlh","Mmp13"))

# Default parameters.
p <- SCpubr::do_ExpressionHeatmap(sample = obj,
                                  features = genes,
                                  group.by = c("cell_types","genotype"),
                                  use_viridis = TRUE,
                                  viridis.direction = -1,
                                  flip = TRUE)
p



# Enrichment score heatmaps -----------------------------------------------
# Define list of genes.
genes

# Default parameters.
p <- SCpubr::do_EnrichmentHeatmap(sample = obj,
                                  group.by = c("cell_types","genotype","grouped_clusters"),
                                  input_gene_list = genes,
                                  viridis.direction = -1,
                                  legend.length = 5,
                                  legend.position = "right",
                                  flip = FALSE)
p

# Other scorings
p <- SCpubr::do_EnrichmentHeatmap(sample = obj,
                                  group.by = c("cell_types","genotype","grouped_clusters"),
                                  input_gene_list = genes,
                                  viridis.direction = -1,
                                  legend.length = 5,
                                  legend.position = "right",
                                  flip = FALSE,
                                  flavor = "UCell")
p


p <- SCpubr::do_EnrichmentHeatmap(sample = obj,
                                  group.by = c("cell_types","genotype","grouped_clusters"),
                                  input_gene_list = genes,
                                  viridis.direction = -1,
                                  legend.length = 5,
                                  legend.position = "right",
                                  flip = FALSE,
                                  flavor = "AUCell")
p



# Correlation matrix heatmaps ---------------------------------------------


## Using Highly Variable Genes ---------------------------------------------

# Default values.
p <- SCpubr::do_CorrelationPlot(sample = obj,
                                group.by = "cell_types")
p
# Would be useful for clustering analysis


# Cellular state plots ----------------------------------------------------


gene_set <- list("A" = Seurat::VariableFeatures(obj)[1:10],
                 "B" = Seurat::VariableFeatures(obj)[11:20],
                 "C" = Seurat::VariableFeatures(obj)[21:30],
                 "D" = Seurat::VariableFeatures(obj)[31:40])

# 2 Variables
p <- SCpubr::do_CellularStatesPlot(sample = obj,
                                   input_gene_list = gene_set,
                                   x1 = "A",
                                   y1 = "B",
                                   enforce_symmetry = TRUE,
                                   colors.use = genotypecolors)
p

# 3 variables
p <- SCpubr::do_CellularStatesPlot(sample = obj,
                                   input_gene_list = gene_set,
                                   x1 = "A",
                                   y1 = "B",
                                   x2 = "C",
                                   enforce_symmetry = TRUE, colors.use = genotypecolors)
p

# 4 variable plots
p <- SCpubr::do_CellularStatesPlot(sample = obj,
                                   input_gene_list = gene_set,
                                   x1 = "A",
                                   y1 = "C",
                                   x2 = "B",
                                   y2 = "D",
                                   enforce_symmetry = TRUE,
                                   colors.use = genotypecolors)
p


## Continuous features -----------------------------------------------------

# Plot continuous features.
out <- SCpubr::do_CellularStatesPlot(sample = obj,
                                     input_gene_list = gene_set,
                                     colors.use = genotypecolors,
                                     x1 = "A",
                                     y1 = "C",
                                     x2 = "B",
                                     y2 = "D",
                                     plot_cell_borders = TRUE,
                                     use_viridis = TRUE,
                                     enforce_symmetry = TRUE,
                                     plot_features = TRUE,
                                     features = c("PC_1", "nFeature_RNA"))
p <- out$main | out$PC_1 | out$nFeature_RNA
p



# Plot enrichment scores for the input gene lists.
out <- SCpubr::do_CellularStatesPlot(sample = obj,
                                     input_gene_list = gene_set,
                                     colors.use = genotypecolors,
                                     x1 = "A",
                                     y1 = "C",
                                     x2 = "B",
                                     y2 = "D",
                                     plot_cell_borders = TRUE,
                                     use_viridis = TRUE,
                                     enforce_symmetry = TRUE,
                                     plot_enrichment_scores = TRUE)
layout <- "AABC
           AADE"
p <- patchwork::wrap_plots(A = out$main,
                           B = out$A,
                           C = out$B,
                           D = out$C,
                           E = out$D,
                           design = layout)
p


# Pathway Activity Plots --------------------------------------------------

# Define your sample and assay.
assay <- "RNA"

# Retrieve prior knowledge network.
network <- decoupleR::get_progeny(organism = "mouse")

# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(obj@assays$RNA$data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)


# General heatmap.
out <- SCpubr::do_PathwayActivityPlot(sample = obj,
                                      activities = activities)
p <- out$heatmaps$average_scores
p

