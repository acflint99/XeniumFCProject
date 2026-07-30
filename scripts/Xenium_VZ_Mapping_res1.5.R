# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(Seurat)
library(future)
library(future.apply)
library(here)

# 1. Setup Parallel Backend
# Increased globals limit because we're passing a lookup table
options(future.globals.maxSize = 200 * 1024^3) 
plan(multisession, workers = 8) 

# 2. Define Paths and "Master Key"
input_dir  <- here("outputs", "Xenium_AldingerABT_VZ_QC_Res1.5_RDS")
output_dir <- here("outputs", "Xenium_AldingerABT_VZsubclusters_Res1.5_RDS")
if(!dir.exists(output_dir)) dir.create(output_dir)

# FIX: Load the object from the path provided
message("Loading master subclustered object...")
master_obj <- readRDS(here("outputs", "XenAld_VZ_Subclusters_Res1.5_RDS", "Xenium_VZ_Res1.5_newSubclusters_4-3-26.rds"))

# PRE-EXTRACT LABELS
all_new_labels <- as.character(Idents(master_obj))
names(all_new_labels) <- colnames(master_obj)

# Update pattern to match your original whole objects
sample_files <- list.files(input_dir, pattern = "_Ald_VZ_QC\\.rds$", full.names = TRUE)

# 3. Run Parallel Processing
message("Starting mapping to whole objects...")

updated_status <- future_lapply(sample_files, function(f) {
  
  # UPDATED: Strip the specific suffix to get the clean sample name
  # This turns "SampleName_Aldinger_annotated.rds" -> "SampleName"
  s_name <- gsub("_Ald_VZ_QC\\.rds", "", basename(f))
  
  # Load the WHOLE object
  temp_obj <- readRDS(f)
  
  # BARCODE MATCHING
  original_barcodes <- colnames(temp_obj)
  integrated_style_barcodes <- paste0(s_name, "_", original_barcodes)
  
  # Extract labels 
  matched_labels <- all_new_labels[integrated_style_barcodes]
  
  # 2. Extract and Label
  # We use match() to ensure the order perfectly aligns with temp_obj barcodes
  # This avoids the "length mismatch" error
  match_idx <- match(integrated_style_barcodes, names(all_new_labels))
  relevant_labels <- all_new_labels[match_idx]
  names(relevant_labels) <- colnames(temp_obj) # Rename back to original style
  
  # 3. Safe Injection
  temp_obj <- AddMetaData(temp_obj, metadata = relevant_labels, col.name = "VZ_subcluster")
  
  # 4. Success Check
  match_count <- sum(!is.na(temp_obj$VZ_subcluster))
  
  # Save the updated whole object
  out_path <- file.path(output_dir, paste0(s_name, "_Ald_VZ_QC_Subclusters.rds"))
  saveRDS(temp_obj, file = out_path, compress = FALSE)
  
  # Clean up memory for the next loop
  rm(temp_obj, relevant_labels, match_idx)
  gc()
  
  return(paste0(s_name, ": Mapped ", match_count, " VZ cells."))
}, future.seed = TRUE)

plan(sequential)
message("Done! Refined VZ labels mapped back to whole tissue objects.")