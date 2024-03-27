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







# Subset chondro.3 and chondro.2 clusters
wnts <- c("Lrp5", "Axin2", "Nkd1", "Pax9", "Sfrp2")

obj$cell_types <- Idents(obj)

selected_cells_chondro3 <- names(obj$cell_types[obj$cell_types == "chondro.3"])
selected_cells_chondro2 <- names(obj$cell_types[obj$cell_types == "chondro.2"])

data_chondro3 <- FetchData(obj,
                           vars = c(wnts, "genotype"),
                           cells = selected_cells_chondro3,
                           layer = "data")

data_chondro2 <- FetchData(obj,
                           vars = c(wnts, "genotype"),
                           cells = selected_cells_chondro2,
                           layer = "data")

# Filter data to include only WNTs genes
data_chondro3 <- data_chondro3[, c(wnts, "genotype")]
data_chondro2 <- data_chondro2[, c(wnts, "genotype")]

long_data_chondro3 <- reshape2::melt(data_chondro3)
long_data_chondro2 <- reshape2::melt(data_chondro2)

ggplot() +
  geom_violin(data = long_data_chondro3, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.4),
              position = position_dodge(width = 0), width = 2) +  # Adjust width as needed
  geom_violin(data = long_data_chondro2, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.4),
              position = position_dodge(width = 0), width = 2) +  # Adjust width as needed
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "italic"),
        panel.grid = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  theme(panel.spacing = unit(0.5, "lines"))  # Adjust spacing between panels




chondro3 <- ggplot() +
  geom_violin(data = long_data_chondro3, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.4),
              position = position_dodge(width = 0), width = 2) +  # Adjust width as needed
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "italic"),
        panel.grid = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  theme(panel.spacing = unit(0.5, "lines"))  # Adjust spacing between panels



chondro2 <- ggplot() +
  geom_violin(data = long_data_chondro2, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.4),
              position = position_dodge(width = 0), width = 2) +  # Adjust width as needed
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "italic"),
        panel.grid = element_blank()) +
  coord_flip() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  theme(panel.spacing = unit(0.5, "lines"))  # Adjust spacing between panels





combined_data <- rbind(transform(long_data_chondro3, cluster = "chondro.3"),
                       transform(long_data_chondro2, cluster = "chondro.2"))


dotcols <- c("hotpink3","dodgerblue")
# Plot
combined_plot <- ggplot(combined_data, aes(x = variable, y = value, fill = genotype, color = genotype, alpha = 0.2)) +
  geom_violin(position = position_dodge(width = 0), width = 10) +  # Adjust width as needed
  geom_jitter(size = 0.5, position = position_dodge2(width = 5), alpha = 0.9)+
  scale_fill_manual(values = my_colors) +  # Adjust colors as needed
  scale_color_manual(values = my_colors) +  # Adjust outline colors as needed
  labs(x = "Genes", y = "Expression", fill = "Genotype", color = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~ variable + cluster, scales = "fixed", nrow = 1, strip.position = "bottom") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "italic", hjust = 0)) +
  theme(panel.spacing = unit(0, "lines"))  # Adjust spacing between panels

# Print the combined plot
print(combined_plot)

