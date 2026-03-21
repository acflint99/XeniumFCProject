# ==============================================================================
# Script: Merge Xenium VZ Subsets for BANKSY
# Purpose: Merge subsets while preserving spatial centroids and stripping images
# Date: 2026-03-17
# ==============================================================================

# 1. Environment Setup
rm(list = ls())

library(here)
library(Seurat)
library(dplyr)
library(future)

# HPC optimization: Disable parallel for reading/merging to save overhead
plan("sequential") 
options(future.globals.maxSize = 400 * 1024^3)

# 2. Path Setup
subset_path <- here("outputs", "XenAld_VZ_Subsets_RDS")
subset_files <- list.files(subset_path, pattern = "\\.rds$", full.names = TRUE)

# Generate clean sample names from filenames
sample_names <- gsub("_VZsubset.rds", "", basename(subset_files))

# 3. Process Subsets (Extract Centroids & Strip Images)
# We use mapply to pass the sample names directly into the objects
subsets <- mapply(function(f, s_name) {
  message("Processing: ", s_name)
  obj <- readRDS(f)
  
  # A. Extract Centroids
  # Using a flexible grep to find 'x' and 'y' columns in the coordinates
  img_name <- names(obj@images)[1]
  coords <- GetTissueCoordinates(obj, image = img_name)
  
  # Handle both Seurat v4 (x, y) and v5 (x_centroid, y_centroid) naming
  x_col <- grep("x", colnames(coords), ignore.case = TRUE, value = TRUE)[1]
  y_col <- grep("y", colnames(coords), ignore.case = TRUE, value = TRUE)[1]
  
  obj$x_centroid <- coords[[x_col]]
  obj$y_centroid <- coords[[y_col]]
  
  # B. Set explicit sample ID for BANKSY grouping later
  obj$sample_id <- s_name
  
  # C. Strip heavy spatial overhead (boundaries and molecules)
  # This makes the merged object manageable on the HPC
  obj@images <- list() 
  
  return(obj)
}, subset_files, sample_names, SIMPLIFY = FALSE)

# 4. Merge into one object
# add.cell.ids ensures barcodes are unique across samples
merged_obj <- merge(
  x = subsets[[1]],
  y = subsets[-1],
  add.cell.ids = sample_names,
  project = "Xenium_VZ_Refinement"
)

# Standardize orig.ident for easier downstream grouping
merged_obj$orig.ident <- merged_obj$sample_id

# 5. Memory Cleanup
rm(subsets)
gc()

# 6. Save the Merged Object
output_dir <- here("outputs", "XenAld_VZ_Subsets_RDS", "Merged")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_path <- file.path(output_dir, "Xenium_Merged_VZSubsets_spatial_31726.rds")

message("Saving merged object to: ", save_path)

# compress = FALSE is crucial for fast I/O on HPC network storage
saveRDS(merged_obj, file = save_path, compress = FALSE)

message("Workflow Complete.")