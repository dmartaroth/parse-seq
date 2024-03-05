# ## ######################################## ## #
#               CLUSTER COMPARISON               #
# ## ######################################## ## #

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Create subdirectory "06_Cluster_Comparison"
comparison_dir <- file.path(figs, "06_Cluster_Comparison")
if (!dir.exists(comparison_dir)) {
  dir.create(comparison_dir)
}

plot_number <- 0  # Starting plot number

# This script is for characterization of known cell types between genotypes


# Load annotated integrated object ----------------------------------------

obj <- readRDS(file = paste0(data.output,"/annotated_integrated_filtered_",control,"_",mutant,"_",sample,".Rds"))
obj$cell_types <- Idents(obj)
Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))

(plot <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 3,label.box = TRUE,cols = iss.colors)+
    umap_theme() )

# Genes to use ------------------------------------------------------------

# List of genes for cell state feature plots
progenitors <- c("Axin2","Gli1","Prrx1","Six2")
osteogenic <- c("Crabp1","Runx2","Sp7","Dmp1")
chondrogenic <- c("Col2a1","Acan","Mgp","Sox9")
hypertrophy <- c("Runx2","Col10a1","Mmp13","Fgfr1")
osteoclasts <- c("Ctsk","Mmp9","Pheta1","Cd44")
vascular <- c("Mcam","Vwf","Pecam1","Pdgfrb")
myeloid_lymphocyte <- c("Pou2f2","Il1rl1","Gata2")
neurons_glia <- c("Neurod1","Cplx3","Otx2","Gfra3","Sox10","Foxd3")
erythrocytes <- c("Hba-a1","Hba-a2","Hbb-bs","Gypa","Gybp","Alas2","Klf1","Slc25a37","Slc2a1")
smoothmuscle <- c("Acta2","Tagln","Myh11","Des")


# Cell states --------------------------------------------------------------


## Progenitor cells --------------------------------------------------------

(plot <- multi_feature_plot(obj, features = progenitors, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_progenitors_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## Osteogenic cells --------------------------------------------------------

(plot <- multi_feature_plot(obj, features = osteogenic, colors_use = muted.blue.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_osteogenic_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## Chondrogenic cells ------------------------------------------------------

(plot <- multi_feature_plot(obj, features = chondrogenic, colors_use = violet.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_chondrogenic_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)

prolif.chondro <- c("Sox9","Ptch1","Fgfr3","Igf1")
convenient_multi_feature_plot(features = prolif.chondro, colors_use = violet.gradient, name = "proliferation_stage")

differ.chondro <- c("Col2a1","Acan","Sox5","Runx2")
convenient_multi_feature_plot(features = differ.chondro, 
                              colors_use = violet.gradient, name = "differentiation_stage")


## Hypertrophic cells ------------------------------------------------------

(plot <- multi_feature_plot(obj, features = hypertrophy, colors_use = purple.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_hypertrophic_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)

mat.stage <- c("Col10a1","Vegf","Pthrp","Mmp13")
convenient_multi_feature_plot(features = mat.stage, colors_use = violet.gradient,
                              name = "maturation_stage", height = 5.5)



## Osteoclasts -------------------------------------------------------------

(plot <- multi_feature_plot(obj, features = osteoclasts, colors_use = turquoise.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_osteoclasts_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)



## Vascular cells ----------------------------------------------------------

(plot <- multi_feature_plot(obj, features = vascular, colors_use = red.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_vascular_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## immune ------------------------------------------------------------------

(plot <- multi_feature_plot(obj, features = myeloid_lymphocyte, colors_use = turquoise.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_immune_by_genotype_featureplots.png", plot_number)), width = 6, height = 7, plot)


## neural ------------------------------------------------------------------


(plot <- multi_feature_plot(obj, features = neurons_glia, colors_use = turquoise.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_neural_by_genotype_featureplots.png", plot_number)), width = 6, height = 7, plot)



## smooth muscle ------------------------------------------------------------------


(plot <- multi_feature_plot(obj, features = smoothmuscle, colors_use = turquoise.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_smooth_muscle_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


# Signaling networks ------------------------------------------------------


## Bone Morphogenetic Proteins ---------------------------------------------

BMP.a <- c("Bmp2", "Bmp3","Bmp4","Bmp5")
BMP.b <- c("Bmp6","Bmp7","Bmp8a","Bmp1")

(plot <- multi_feature_plot(obj, features = BMP.a, colors_use = red.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_Bmp2_Bmp3_Bmp4_Bmp5_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)

(plot <- multi_feature_plot(obj, features = BMP.b, colors_use = red.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_Bmp6_Bmp7_Bmp8a_Bmp1_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## Wnt signaling -----------------------------------------------------------

Wnt.a <- c("Wnt4","Wnt5a", "Wnt11","Axin2")
Wnt.b <- c("Dkk2","Dkk3","Notum","Sostdc1")
Wnt.c <- c("Sfrp1","Sfrp2","Sfrp4", "Sfrp5")
Wnt.d <- c("Lrp5","Lrp6","Tnks2","Pax9")
Wnt.e <- c("Ctnnb1","Apcdd1","Nkd1","Nkd2")

(plot <- multi_feature_plot(obj, features = Wnt.a, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_wnt_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


(plot <- multi_feature_plot(obj, features = Wnt.b, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_wnt_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


(plot <- multi_feature_plot(obj, features = Wnt.c, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_wnt_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)



(plot <- multi_feature_plot(obj, features = Wnt.d, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_wnt_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)



(plot <- multi_feature_plot(obj, features = Wnt.e, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_wnt_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## Fgf signaling -----------------------------------------------------------

Fgf.a <- c("Fgf2","Fgf18","Fgf9")
Fgfrs <- c("Fgfr1","Fgfr2","Fgfr3","Fgfr4")
Fgf.b <- c("Hras","Kras","Nras","Araf")
Fgf.c <- c("Braf","Raf1","Map2k1","Map2k2")
Fgf.d <- c("Mapk3","Mapk1","Pik3ca","Pik3cb")

(plot <- multi_feature_plot(obj, features = Fgf.a, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_fgf_by_genotype_featureplots.png", plot_number)), width = 6, height = 7, plot)


(plot <- multi_feature_plot(obj, features = Fgfrs, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_fgfrs_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)




(plot <- multi_feature_plot(obj, features = Fgf.b, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_fgf_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


(plot <- multi_feature_plot(obj, features = Fgf.c, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_fgf_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


(plot <- multi_feature_plot(obj, features = Fgf.d, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_fgf_by_genotype_featureplots.png", plot_number)), width = 6, height = 9, plot)


## Akt ---------------------------------------------------------------------

Akt.a <- c("Akt1","Akt2","Akt3")


(plot <- multi_feature_plot(obj, features = Akt.a, colors_use = green.gradient))

# Save plot
plot_number <- plot_number + 1
ggsave(filename = file.path(comparison_dir, sprintf("%02d_Akt_by_genotype_featureplots.png", plot_number)), width = 6, height = 7, plot)


## Delta/Notch signaling ---------------------------------------------------

notch <- c("Dll4","Jag1","Notch1","Hes1")
convenient_multi_feature_plot(features = notch, colors_use = muted.blue.gradient, name = "Delta_Notch")


## Platelet-Derived Growth Factor (PDGF) -----------------------------------

pdgf <- c("Pdgfb","Pdgfa","Pdgfc","Pdgfd")
convenient_multi_feature_plot(features = pdgf, colors_use = muted.blue.gradient, name = "Pdgf")


# TGFb --------------------------------------------------------------------

tgfb.1 <- c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1")
tgfb.2 <- c("Tgfbr2","Acvr1","Bmpr1a","Bmpr2")


convenient_multi_feature_plot(features = tgfb.1, colors_use = violet.gradient, name = "tgfb_1")
convenient_multi_feature_plot(features = tgfb.2, colors_use = violet.gradient, name = "tgfb_2")


# Biological processes ----------------------------------------------------

## Angiogenesis ------------------------------------------------------------

angio.1 <- c("Vegfa","Vegfb","Vegfc","Vegfd")
angio.2 <- c("Fgf2","Fgf1","Fgf18","Fgf9")
angio.3 <- c("Angpt1","Angpt2","Angpt3","Angpt4")

convenient_multi_feature_plot(features = angio.1, colors_use = red.gradient, name ="angio.1")
convenient_multi_feature_plot(features = angio.2, colors_use = red.gradient, name = "angio.fgfs")
convenient_multi_feature_plot(features = angio.3,colors_use = red.gradient,name = "angiopoietins", height = 7)


## Chondrocyte differentiation -------------------------------------------------------------

pretohypertr <- c("Sox9","Col10a1","Runx2","Mmp13")

convenient_multi_feature_plot(features = pretohypertr, colors_use = violet.gradient, name = "prehyp_to_hypertrophic")


## Bone remodeling ---------------------------------------------------------

remodeling <- c("Tnfsf11","Tnfrsf11b","Alpl","Csf1")
convenient_multi_feature_plot(features = remodeling, colors_use = turquoise.gradient, name = "remodeling")

