library(Seurat)
library(here)
library(dplyr)
library(future)

# 1. Memory and Parallelization Setup
options(future.globals.maxSize = Inf)
plan("sequential")

# 2. Path Setup
out_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

input_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# --- 3. Create and Fill obj_list with PRE-MERGE RENAMING ---
obj_list <- list()

for (i in 1:length(files)) {
  f_name <- basename(files[i])
  message("Reading (", i, "/", length(files), "): ", f_name)
  
  tmp <- readRDS(files[i])
  
  # Slim down the object using Seurat v5 layers syntax
  tmp <- DietSeurat(tmp, layers = c("counts", "data"), dimreducs = NULL, graphs = NULL)
  
  # Set Sample ID
  sample_id <- gsub("_Ald_VZ_RL_QC_Subclusters\\.rds$", "", f_name)
  tmp$sample_id <- sample_id
  
  # --- THE FIX: Pre-Merge Renaming ---
  # Force the cell names to be globally unique right now. 
  # This perfectly syncs the RNA metadata AND the spatial FOV image.
  clean_names <- paste0(sample_id, "_", colnames(tmp))
  tmp <- RenameCells(tmp, new.names = clean_names)
  
  obj_list[[sample_id]] <- tmp
  rm(tmp)
  gc()
}
# ------------------------------------------

# 4. Consolidated destructive merge
message("Starting destructive chain merge...")

# Initialize with the first object and REMOVE it from the list immediately
merged_obj <- obj_list[[1]]
obj_list[[1]] <- NULL 
gc()

# Loop through the rest of the list
for (j in 1:length(obj_list)) {
  sample_name <- names(obj_list)[1]
  message("Merging: ", sample_name, " (Remaining in list: ", length(obj_list), ")")
  
  # --- THE FIX: Remove add.cell.ids ---
  # Because we already renamed them, we let Seurat merge them natively
  merged_obj <- merge(merged_obj, y = obj_list[[1]]) 
  
  # Delete the item we just merged and free RAM
  obj_list[[1]] <- NULL
  gc()
}

message("Merge complete. Final cell count: ", ncol(merged_obj))

# Optional Verification: Check that Seurat didn't add any extra underscores
message("Verification Check (First 5 cells):")
print(head(colnames(merged_obj), 5))

# 5. Save
save_path <- file.path(out_dir, "XenAld_VZ_RL_QC_Subclusters_Spatial_Merged_4-22-26.rds")
message("Saving to: ", save_path)

saveRDS(merged_obj, save_path, compress = FALSE)

# Final Cleanup
rm(obj_list)
gc()

message("Done! Process complete.")