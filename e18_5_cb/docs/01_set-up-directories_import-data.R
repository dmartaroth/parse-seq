# ## ################################### ## #
#        LOAD DATA AND PREPROCESS           #
# ## ################################### ## #

# Downloaded Parse Biosciences pipeline output from UBC bioinformatics team for
# Bmp2 ctrl and ncko cranial base samples.
# Date: Wed Feb 28 10:20:55 2024 ------------------


# Packages, directories, and functions ------------------------------------

library(here)
source(here("e18_5_cb","docs","packages.R")) # load packages
source(here("e18_5_cb","docs","directories.R")) # load file paths/directories
source(here("e18_5_cb","docs","functions.R")) # load functions
source(here("e18_5_cb","docs","themes.R")) # load themes

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


