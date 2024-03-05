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

