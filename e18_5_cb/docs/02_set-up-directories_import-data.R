# ## ################################### ## #
#        LOAD DATA AND PREPROCESS           #
# ## ################################### ## #

# Downloaded Parse Biosciences pipeline output from UBC bioinformatics team for
# Bmp2 ctrl and ncko cranial base samples.
# Date: Sat Feb 24 21:10:29 2024 ------------------

# Libraries ---------------------------------------------------------------

library(here)
library(Seurat)
library(tidyverse)  # Includes ggplot2, dplyr, and others
library(patchwork)
library(gprofiler2)
library(clustree)
library(SingleR)
library(Matrix)
library(cowplot)
library(DropletUtils)
library(crayon)

# Define paths and directories --------------------------------------------
# Set variables for sample, control, and mutant
sample <- "e18_5_cb" # replace quoted text with type of sample
control <- "Bmp2_ctrl" # replace quoted text with genotype of sample
mutant <- "Bmp2_ncko"

# Create home.path directory
home.path <- here(sample)
dir.create(home.path, recursive = TRUE)

# Initialize an empty list to store directory paths
dir_paths <- list()

# Create directories for data-output, docs, figures, src, and raw-data
dirs <- c("docs", "src", "raw-data")
for (dir in dirs) {
  dir_path <- here::here(home.path, dir)
  dir.create(dir_path, recursive = TRUE)
  dir_paths[[dir]] <- dir_path
}
dir.create(data.output <- here(home.path, "data-output"),recursive = TRUE)
dir.create(figs <- here(home.path, "figures"),recursive = TRUE)

# Create subdirectories for control and mutant data
data.path <- here(home.path, "raw-data")
dir.create(data.path.ctrl <- here(data.path, paste(control, sep = "_", sample)))
dir.create(data.path.ncko <- here(data.path, paste(mutant, sep = "_", sample)))

# Once directories are created, move DGE-unfiltered contents to appropriate
# raw-data directories


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
