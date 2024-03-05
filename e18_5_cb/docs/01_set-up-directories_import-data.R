# ## ################################### ## #
#        LOAD DATA AND PREPROCESS           #
# ## ################################### ## #

# Downloaded Parse Biosciences pipeline output from UBC bioinformatics team for
# Bmp2 ctrl and ncko cranial base samples.
# Date: Tue Mar 05 09:53:11 2024 ------------------


# Packages, directories, and functions ------------------------------------

library(here)
source(here::here("e18_5_cb","docs","packages.R")) # load packages
source(here::here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here::here("e18_5_cb","docs","functions.R")) # load functions
source(here::here("e18_5_cb","docs","themes.R")) # load themes

# Import data -------------------------------------------------------------

## Control sample ----------------------------------------------------------
mat_path <- data.path.ctrl
mat <- ReadParseBio(mat_path)

table(rownames(mat) == "") # check to see if empty gene names are present, add name if so.
# returns FALSE

rownames(mat)[rownames(mat) == ""] <- "unknown"

cell_meta <-
  read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1) # read in cell meta data

(ctrl <-
  CreateSeuratObject(
    mat,
    min.features = 100,
    min.cells = 100,
    names.field = 0,
    meta.data = cell_meta
  ))

ctrl$original <- ctrl$orig.ident
ctrl$orig.ident <- NULL
head(x=ctrl[[]])
Idents(ctrl) <- "ctrl"
ctrl$orig.ident <- "ctrl"
head(x=ctrl[[]])

## Mutant sample ----------------------------------------------------------
mat_path <- data.path.ncko
mat <- ReadParseBio(mat_path)

table(rownames(mat) == "") # check to see if empty gene names are present, add name if so.
# returns FALSE

rownames(mat)[rownames(mat) == ""] <- "unknown"

cell_meta <-
  read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1) # read in cell meta data

(ncko <-
    CreateSeuratObject(
      mat,
      min.features = 100,
      min.cells = 100,
      names.field = 0,
      meta.data = cell_meta
    ))

ncko$original <- ncko$orig.ident
ncko$orig.ident <- NULL
head(x=ncko[[]])
Idents(ncko) <- "ncko"
ncko$orig.ident <- "ncko"
head(x=ncko[[]])

# Quality control ---------------------------------------------------------
prepro.plots(ctrl,ncko,figs)

ctrl <- add_percent_mito(ctrl)
ncko <- add_percent_mito(ncko)

# Filtering low quality cells ---------------------------------------------
# Select subset values from prepro.plots output, especially violin plots
# Input needed here

filtered_ctrl <- subset(x = ctrl,
                        subset = (nFeature_RNA > 200) & # value from seurat vignette
                          (nFeature_RNA < 5500) &
                          (nCount_RNA < 50000) &
                          (percent.mito < 0.15))

filtered_ncko <- subset(x = ncko,
                        subset = (nFeature_RNA > 200) &
                          (nFeature_RNA < 5500) &
                          (nCount_RNA < 70000) &
                          (percent.mito < 0.15))

# Gene-level filtering
filtered_ctrl <- genelvlfilt(filtered_ctrl)
filtered_ncko <- genelvlfilt(filtered_ncko)

filt.plots(filtered_ctrl,filtered_ncko,figs)


# Save filtered cells -----------------------------------------------------
saveRDS(filtered_ctrl, file = paste0(data.output,"/filtered_",control,"_",sample,".Rds"))
saveRDS(filtered_ncko, file = paste0(data.output,"/filtered_",mutant,"_",sample,".Rds"))
