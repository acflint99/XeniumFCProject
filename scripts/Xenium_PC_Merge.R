# Clear the environment
rm(list = ls())

library(here)
library(Seurat)
library(dplyr)
library(future)

plan("sequential") # This disables parallel workers for the current session

# Set limit to 400GB
options(future.globals.maxSize = 400 * 1024^3)

# 1. List all the subset files you just created
subset_path <- here("outputs", "XenAld_PC_Subsets_RDS")
subset_files <- list.files(subset_path, pattern = "\\.rds$", full.names = TRUE)

# 2. Read them into a list and strip spatial overhead
subsets <- lapply(subset_files, function(f) {
  obj <- readRDS(f)
  
  # Remove all image data (boundaries and molecule coordinates)
  obj@images <- list() 
  
  # Optional: Clear out unnecessary assays if any exist besides "RNA" or "Xenium"
  # obj[["ControlCodewords"]] <- NULL 
  
  return(obj)
})
names(subsets) <- gsub("_PCsubset.rds", "", basename(subset_files))

# 3. Merge into one object
# add.cell.ids prepends the sample name to the barcodes to prevent duplicates
  merged_obj <- merge(
    x = subsets[[1]],
    y = subsets[-1],
    add.cell.ids = names(subsets),
    project = "Xenium_PC_Refinement"
  )

# Use the 'Remove everything after the LAST underscore' logic
merged_obj$orig.ident <- gsub("_[^_]+$", "", colnames(merged_obj))

# 4. Clean up the list to free memory
rm(subsets)
gc()

# 5. Save the merged object
# Adding compress = FALSE makes this take 5 mins instead of 50 mins
saveRDS(merged_obj, 
        here("outputs", "XenAld_PC_Subsets_RDS", "Merged", "Xenium_Merged_PCSubsets_31626.rds"), #edit date if redo
        compress = FALSE)