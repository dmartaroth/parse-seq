# ## ######################################## ## #
#                 CHONDROCYTE SUBSET             #
# ## ######################################## ## #


library(here)
source(here("e18_5_cb","docs","packages.R")) # load packages
source(here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here("e18_5_cb","docs","functions.R")) # load functions
source(here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "07_Subset_chondro"
subset_chondro_dir <- file.path(figs, "07_Subset_chondro")
if (!dir.exists(subset_chondro_dir)) {
  dir.create(subset_chondro_dir)
}

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))

# For visualization of chondrogenic genes, would be useful to look in scope of
# chondrocyte clusters only

# Subset chondro clusters
chondrocyte.clusters <- c("chondro.1","chondro.2","chondro.3","chondro.4")
chondro_subset <- subset(obj, idents = chondrocyte.clusters)

(plot <- DimPlot(chondro_subset, reduction = "umap", label = FALSE,split.by = "genotype",cols = iss.colors))


plot_number <- 0  # Starting plot number

# Cell states -------------------------------------------------------------

# Proliferation stage



prolif.chondro <- c("Sox9","Ptch1","Fgfr3","Igf1")
convenient_multi_feature_plot(seurat_object = chondro_subset,features = prolif.chondro, colors_use = violet.gradient, name = "proliferation_stage",
                              dir = subset_chondro_dir)


# Differentiation stage
differ.chondro <- c("Col2a1","Acan","Sox5","Runx2")
convenient_multi_feature_plot(seurat_object = chondro_subset,features = differ.chondro, 
                              colors_use = violet.gradient, name = "differentiation_stage",
                              dir = subset_chondro_dir)

# Maturation stage

mat.stage <- c("Col10a1","Vegf","Pthrp","Mmp13")
convenient_multi_feature_plot(seurat_object = chondro_subset,features = mat.stage, colors_use = violet.gradient,
                              name = "maturation_stage", height = 7, dir = subset_chondro_dir)

# Hypertrophy stage

hypertr.stage <- c("Runx2","Col10a1","Mmp13","Vegf")
convenient_multi_feature_plot(seurat_object = chondro_subset, features = hypertr.stage, 
                              colors_use = violet.gradient, name = "hypertrophy_stage",
                              dir = subset_chondro_dir, height = 10)

