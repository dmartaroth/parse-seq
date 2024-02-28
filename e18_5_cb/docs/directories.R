# ## #################################### ## #
#                 DIRECTORIES                #
# ## #################################### ## #

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