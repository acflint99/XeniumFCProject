# 1. Initialize the project-specific library
# renv::load() automatically finds the .Rprofile and sets the libPaths
if (file.exists("renv/activate.R")) source("renv/activate.R")

library(here)
library(Seurat)
library(dplyr)

# Capture Slurm Array ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No Task ID provided. This script should be run via Slurm array.")
}
task_id <- as.numeric(args[1])

# Define your 15 samples explicitly to match the array index
samples <- c(
  #"GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G", "FB330_1_X_G",
  #"FB78_X_G", "GZFB5_X_G", "GZFB_12_X_G_1", "GZFB_12_X_G_2", "GZFB_12_X_G_3",
  "GZFB_12_X_G_4", "GZFB_12_X_G_5", "GZFB_1_X_G"#, "GZFB_9_X_G_1", "GZFB_9_X_G_2"
)

current_sample <- samples[task_id]

# Use here() to build the input path
# Structure: project_root/outputs/XeniumAldingerABT_RDS/...
input_path <- here("outputs", "XeniumAldingerABT_RDS", paste0(current_sample, "_Aldinger_annotated.rds"))

message("Loading: ", input_path)
if (!file.exists(input_path)) stop("File not found: ", input_path)

obj <- readRDS(input_path)

# Identify the clusters you want to pull out for re-analysis
# (Update these strings to match your exact 'annotated' cluster names)
target_clusters <- c("Glia", "GABA", "Purkinje", "OPC") 

# Subset the object
obj_subset <- subset(obj, idents = target_clusters)

# Define and create output directory using here()
output_dir <- here("outputs", "Xenium_VZ_Subsets")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_path <- file.path(output_dir, paste0(current_sample, "_VZsubset.rds"))
saveRDS(obj_subset, file = output_path)

message("Successfully saved subset for ", current_sample, " to ", output_path)