# Clear the environment
rm(list = ls())

# 1. INITIALIZATION
source("renv/activate.R")
library(Seurat)
library(here)

# 2. Define Paths - FIXED: Removed the trailing space
input_dir  <- here("outputs", "XeniumBanksyRDS")
output_dir <- here("outputs", "XeniumAldingerABT_Banksy_VZsubcluster_RDS")
if(!dir.exists(output_dir)) dir.create(output_dir)

# Load Master Key
message("Loading master subclustered object...")
master_obj <- readRDS(here("outputs", "XenAld_VZ_Banksy_RDS", "Xenium_VZ_Banksy_Annotated.rds"))

# PRE-EXTRACT LABELS
all_new_labels <- as.character(Idents(master_obj))
names(all_new_labels) <- colnames(master_obj)

# Get the 3 files
sample_files <- list.files(input_dir, pattern = "_CB_QC_Bcluster\\.rds", full.names = TRUE)
message("Found ", length(sample_files), " samples.")

# 3. Process Sequentially (Safer for small sample sizes)
results <- lapply(sample_files, function(f) {
  
  s_name <- gsub("_CB_QC_Bcluster\\.rds", "", basename(f))
  message("Processing: ", s_name)
  
  # Load
  temp_obj <- readRDS(f)
  original_barcodes <- colnames(temp_obj)
  
  # BARCODE MATCHING (Handles the FB124_X_G_ style)
  # Try both common prefix styles to be safe
  test_bc_1 <- paste0(s_name, "_", original_barcodes) # Style: Sample__Barcode
  test_bc_2 <- paste0(s_name, original_barcodes)      # Style: Sample_Barcode
  
  if (sum(test_bc_1 %in% names(all_new_labels)) > sum(test_bc_2 %in% names(all_new_labels))) {
    integrated_style_barcodes <- test_bc_1
  } else {
    integrated_style_barcodes <- test_bc_2
  }
  
  # Extract and Label
  match_idx <- match(integrated_style_barcodes, names(all_new_labels))
  relevant_labels <- all_new_labels[match_idx]
  
  # RENAME back to original barcodes so Seurat accepts the metadata
  names(relevant_labels) <- original_barcodes
  
  # Safe Injection
  temp_obj <- AddMetaData(temp_obj, metadata = relevant_labels, col.name = "Banksy_VZ_Refined")
  
  # Success Check
  match_count <- sum(!is.na(temp_obj$Banksy_VZ_Refined))
  
  # Save
  out_path <- file.path(output_dir, paste0(s_name, "_Banksy_VZ.rds"))
  saveRDS(temp_obj, file = out_path, compress = FALSE)
  
  return(paste0(s_name, ": Mapped ", match_count, " refined labels."))
})

print(unlist(results))