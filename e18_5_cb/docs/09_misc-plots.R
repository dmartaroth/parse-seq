# ## ######################################## ## #
#                     MISC PLOTS                 #
# ## ######################################## ## #

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "10_misc-plots"
miscplot.dir <- file.path(figs, "10_misc-plots")
if (!dir.exists(miscplot.dir)) {
  dir.create(miscplot.dir)
}

plot_number <- 0

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))
obj$cell_types <- Idents(obj)
obj$cluster <- obj$cell_types


# BMPs --------------------------------------------------------------------


genes <- c("Bmp3","Bmp4","Bmp6","Bmp7")

selected_cells <- names(obj$cell_types[obj$cell_types %in% c("chondro.1", "chondro.2", "chondro.3")])


# Fetch data for both chondro and osteo groups
data_cells <- FetchData(obj,
                          vars = c(genes, "genotype","cell_types"),
                          cells = selected_cells,
                          layer = "data")


# Filter data to include only the top 10 genes
data_cells <- data_cells[, c(genes, "genotype","cell_types")]


# Reshape data into long format
long_data_cells <- reshape2::melt(data_cells)
combined_data <- long_data_cells



# Plot combined data using a violin plot
p <- ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
  geom_violin(position = position_dodge(width = 0), width = 10) +  # Adjust width as needed
  geom_jitter(size = 0, position = position_dodge2(width = 5), alpha = 0.9) +
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~ variable + cell_types, scales = "fixed", nrow = 1, strip.position = "bottom") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
  theme(panel.spacing = unit(0, "lines"))  # Adjust spacing between panels

convenient_save_plot(p, name = "Bmp_violinplot_chondro1-3", dir = miscplot.dir,height = 3.5, width = 8)



# TGFs --------------------------------------------------------------------


genes <- c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1")

selected_cells <- names(obj$cell_types[obj$cell_types %in% c("chondro.1", "chondro.2", "chondro.3")])


# Fetch data for both chondro and osteo groups
data_cells <- FetchData(obj,
                        vars = c(genes, "genotype","cell_types"),
                        cells = selected_cells,
                        layer = "data")


# Filter data to include only the top 10 genes
data_cells <- data_cells[, c(genes, "genotype","cell_types")]


# Reshape data into long format
long_data_cells <- reshape2::melt(data_cells)
combined_data <- long_data_cells



# Plot combined data using a violin plot
p <- ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
  geom_violin(position = position_dodge(width = 0), width = 10) +  # Adjust width as needed
  geom_jitter(size = 0, position = position_dodge2(width = 5), alpha = 0.9) +
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~ variable + cell_types, scales = "fixed", nrow = 1, strip.position = "bottom") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
  theme(panel.spacing = unit(0, "lines"))  # Adjust spacing between panels

convenient_save_plot(p, name = "Tgfb_violinplot_chondro1-3", dir = miscplot.dir,height = 3.5, width = 8)


#
chondro.clusters <- c("chondro.1","chondro.2","chondro.3")
chondro_subset <- subset(obj, idents = chondro.clusters)

Bmps <- c("Bmp3","Bmp4","Bmp6","Bmp7")
convenient_multi_feature_plot(seurat_object = chondro_subset, features = Bmps, 
                              colors_use = violet.gradient, name = "Bone Morphogenetic Proteins",
                              dir = miscplot.dir, height = 10)


TGFBR <- c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1")
convenient_multi_feature_plot(seurat_object = chondro_subset, features = TGFBR, 
                              colors_use = green.gradient, name = "TGFB",
                              dir = miscplot.dir, height = 10)

attachment <- c("Ptn","Postn","Cacna1g","Svil")

convenient_multi_feature_plot(seurat_object = chondro_subset, features = attachment, 
                              colors_use = muted.blue.gradient, name = "Cell attachment & spreading",
                              dir = miscplot.dir, height = 10)


(p <- SCpubr::do_FeaturePlot(
  sample = chondro_subset,
  features = attachment,
  split.by = "genotype",
  font.size = 7,
  use_viridis = TRUE,
  viridis.palette = "magma",
  viridis.direction = 1,
  na.value = "floralwhite",
  plot_cell_borders = FALSE,
  order = TRUE,
  pt.size = 1,
  legend.position = "right")
)
convenient_save_plot(p, name = "featureplots_attachment_osteochondrosubset", dir = miscplot.dir,height = 4, width = 15)



# matrix --------------------------------------------------------------------


genes <- c("Col1a1","Col1a2","Col2a1","Col3a1","Col4a1","Col5a2","Col9a1","Col9a2","Col11a2")
# 
# (plot <- DotPlot( object = chondro_subset,
#                   features =   genes,
#                   cols = c("hotpink1", "cadetblue2"),
#                   scale.by = "radius",
#                   dot.scale = 8,
#                   scale = TRUE,
#                   split.by = "genotype",
#                   cluster.idents = FALSE
# ) +   custom_dotplot_theme() +RotatedAxis())

p <- VlnPlot(chondro_subset, features = genes, split.by = "genotype", 
        group.by = "cell_types",pt.size = 0, combine = TRUE, cols = c("#FFB6C1", "#ADD8E6"),
        ncol = 3)

convenient_save_plot(p, name = "collagens_osteochondrosubset", dir = miscplot.dir,height = 4, width = 8)


# others --------------------------------------------------------------------


genes <- c("Cdh11","Ptn","Cped1","Postn","Cacna1g","Zeb2","Thbs2","Igfbp5","Svil")
# 
# (plot <- DotPlot( object = chondro_subset,
#                   features =   genes,
#                   cols = c("hotpink1", "cadetblue2"),
#                   scale.by = "radius",
#                   dot.scale = 8,
#                   scale = TRUE,
#                   split.by = "genotype",
#                   cluster.idents = FALSE
# ) +   custom_dotplot_theme() +RotatedAxis())

p <- VlnPlot(chondro_subset, features = genes, split.by = "genotype", 
             group.by = "cell_types",pt.size = 0, combine = TRUE, cols = c("#FFB6C1", "#ADD8E6"),
             ncol = 9)
p
convenient_save_plot(p, name = "cell-attachment_chondrosubset", dir = miscplot.dir,height = 2.8, width = 14)

