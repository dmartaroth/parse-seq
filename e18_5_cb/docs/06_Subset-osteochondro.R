# ## ######################################## ## #
#                 CHONDROCYTE SUBSET             #
# ## ######################################## ## #


library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "07_Subset_chondro"
subset_osteochondro_dir <- file.path(figs, "07_Subset_osteochondro")
if (!dir.exists(subset_osteochondro_dir)) {
  dir.create(subset_osteochondro_dir)
}

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))

# For visualization of chondrogenic genes, would be useful to look in scope of
# chondrocyte clusters only

# Subset chondro clusters
osteochondro.clusters <- c("chondro.1","chondro.2","chondro.3","chondro.4","mes.1","mes.2","mes.3")
osteochondro_subset <- subset(obj, idents = osteochondro.clusters)

(plot <- DimPlot(osteochondro_subset, reduction = "umap", label = FALSE,split.by = "genotype",cols = osteochondro.colors))


plot_number <- 0  # Starting plot number

# Cell states -------------------------------------------------------------

# Proliferation stage



prolif.chondro <- c("Sox9","Ptch1","Fgfr3","Igf1")
convenient_multi_feature_plot(seurat_object = osteochondro_subset,features = prolif.chondro, colors_use = violet.gradient, name = "proliferation_stage",
                              dir = subset_osteochondro_dir)


# Differentiation stage
differ.chondro <- c("Col2a1","Acan","Sox5","Runx2")
convenient_multi_feature_plot(seurat_object = osteochondro_subset,features = differ.chondro, 
                              colors_use = violet.gradient, name = "differentiation_stage",
                              dir = subset_osteochondro_dir)

# Maturation stage

mat.stage <- c("Col10a1","Vegf","Pthrp","Mmp13")
convenient_multi_feature_plot(seurat_object = osteochondro_subset,features = mat.stage, colors_use = violet.gradient,
                              name = "maturation_stage", height = 7, dir = subset_osteochondro_dir)

# Hypertrophy stage

hypertr.stage <- c("Runx2","Col10a1","Mmp13","Vegf")
convenient_multi_feature_plot(seurat_object = osteochondro_subset, features = hypertr.stage, 
                              colors_use = violet.gradient, name = "hypertrophy_stage",
                              dir = subset_osteochondro_dir, height = 10)






# Your defined list of genes
chondro.diff <- c("Sox9", "Ptch1", "Fgfr3", "Igf1", "Col2a1", "Acan", "Sox5", 
                  "Runx2", "Col10a1", "Vegfa", "Pthlh", "Mmp13")





vln_plot <- VlnPlot(
  object = osteochondro_subset,
  cols = my_colors,
  features = chondro.diff,
  split.by = "genotype",      # Split plots by genotype
  pt.size = 0.1,              # Point size
  combine = TRUE,
  split.plot = TRUE
)
vln_plot

convenient_save_plot(vln_plot, "chondro.diff_osteochondro-subset_vlnplot_by-genotype",width = 15, height = 12, dir = subset_osteochondro_dir)




genes <- c("Cfh","Tns1","Cntn3","Eya2","Fat3","Vit","Peak1","Brip1","Cemip","Syne1","Ibsp","Col1a2","Htra3")


vln_plot <- VlnPlot(
  object = osteochondro_subset,
  cols = my_colors,
  features = genes,
  split.by = "genotype",      # Split plots by genotype
  pt.size = 0.1,              # Point size
  combine = TRUE,
  split.plot = TRUE
)
vln_plot

convenient_save_plot(vln_plot, "selected_degs_osteochondro-subset_vlnplot_by-genotype",width = 13, height = 10, dir = subset_osteochondro_dir)


smads <- c("Smad1","Smad5","Smad9")


vln_plot <- VlnPlot(
  object = osteochondro_subset,
  cols = my_colors,
  features = smads,
  split.by = "genotype",      # Split plots by genotype
  pt.size = 0.1,              # Point size
  combine = TRUE,
  split.plot = TRUE
)
vln_plot

convenient_save_plot(vln_plot, "smads_osteochondro-subset_vlnplot_by-genotype",width = 8, height = 4, dir = subset_osteochondro_dir)


