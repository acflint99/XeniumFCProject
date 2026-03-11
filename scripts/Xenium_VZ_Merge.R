# Clear the environment
rm(list = ls())

library(here)
library(Seurat)
library(dplyr)
library(future)
plan("sequential") # This disables parallel workers for the current session

# Set limit to 100GB (100 * 1024^3 bytes)
options(future.globals.maxSize = 250 * 1024^3)

# 1. List all the subset files you just created
subset_path <- here("outputs", "Xenium_VZ_Subsets")
subset_files <- list.files(subset_path, pattern = "\\.rds$", full.names = TRUE)

# 2. Read them into a list
# This uses the filename to name the list elements for easier tracking
subsets <- lapply(subset_files, readRDS)
names(subsets) <- gsub("_VZsubset.rds", "", basename(subset_files))

# 3. Merge into one object
# add.cell.ids prepends the sample name to the barcodes to prevent duplicates
merged_obj <- merge(
  x = subsets[[1]],
  y = subsets[-1],
  add.cell.ids = names(subsets),
  project = "Xenium_VZ_Refinement"
)

# Use the 'Remove everything after the LAST underscore' logic
merged_obj$orig.ident <- gsub("_[^_]+$", "", colnames(merged_obj))



# 4. Clean up the list to free memory
rm(subsets)
gc()

# 5. Save the merged object
# Adding compress = FALSE makes this take 5 mins instead of 50 mins
saveRDS(merged_obj, 
        here("outputs", "Xenium_VZ_Subsets", "Xenium_Merged_VZSubsets.rds"),
        compress = FALSE)