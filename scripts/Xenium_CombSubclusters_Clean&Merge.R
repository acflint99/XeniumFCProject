# Clear the environment
rm(list = ls())

library(Seurat)
library(here)

options(future.globals.maxSize = +Inf)

# 1. Setup Paths
# Project root: /home/acflint/R/Projects/XeniumFCProject/
i_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
o_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS")

if (!dir.exists(o_path)) dir.create(o_path, recursive = TRUE)

files <- list.files(i_path, pattern = "\\.rds$", full.names = TRUE)

# 2. Individual Cleaning Loop
for (f in files) {
  # Extract Sample ID (e.g., "GZFB5_X_G" from "GZFB5_X_G_Ald_VZ_RL_QC_Subclusters.rds")
  s_id <- gsub("_Ald_VZ_RL_QC_Subclusters\\.rds", "", basename(f))
  message("Processing Sample: ", s_id)
  
  # Load
  tmp <- readRDS(f)
  
  # Strip the 3GB+ dense layer
  tmp[["Xenium"]]@layers$scale.data <- NULL
  
  # Remove QC assays
  tmp[["BlankCodeword"]] <- NULL
  tmp[["ControlCodeword"]] <- NULL
  tmp[["ControlProbe"]] <- NULL
  tmp[["GenomicControl"]] <- NULL
  
  # REMOVE ALL SPATIAL DATA
  tmp@images <- list()
  
  # Assign the clean Sample ID to the project name
  tmp@project.name <- s_id
  
  # Save intermediate clean file
  saveRDS(tmp, file.path(o_path, paste0("clean_", s_id, ".rds")), compress = FALSE)
  
  rm(tmp)
  gc()
}

# 3. Merging
clean_files <- list.files(o_path, pattern = "^clean_.*\\.rds$", full.names = TRUE)

# Read all cleaned objects into a list
obj_list <- lapply(clean_files, readRDS)

# Extract IDs again for the list names to ensure merge uses them correctly
names(obj_list) <- gsub("clean_|.rds", "", basename(clean_files))

message("Merging ", length(obj_list), " samples...")

# Merge into one master object
merged_obj <- merge(
  x = obj_list[[1]], 
  y = obj_list[2:length(obj_list)], 
  add.cell.ids = names(obj_list)
)

# 4. Final Save
saveRDS(merged_obj, here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS", "XenAld_VZ&RL_clean_merge_4-6-26.rds"), compress = FALSE)

message("Done! Final object contains IDs like: ", paste(head(names(obj_list), 2), collapse = ", "))