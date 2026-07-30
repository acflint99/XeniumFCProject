# Clear the environment
rm(list = ls())

library(Seurat)
library(here)
library(dplyr)
library(ggplot2)
library(reticulate)

# --- 1. Python/Conda Setup ---
env_python <- "/home/acflint/.conda/envs/spatrack_env/bin/python"

# Scrub cluster pollution
Sys.unsetenv("PYTHONPATH")
Sys.unsetenv("PYTHONHOME")

# Set the environment variable BEFORE initializing Python
Sys.setenv(RETICULATE_PYTHON = env_python)

# Load isolated spatrack_env
ad <- import("anndata")
spt <- import("spaTrack")

message("Python environment successfully loaded from: ", env_python)

# --- 2. Load Data and Metadata ---
source(here("scripts", "color_palette.R"))

# Pointing to the clean, integrated object you just built!
merged_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS", "XenAld_VZ_RL_QC_Subclusters_Spatial_Merged_Integrated2_4-23-26.rds")
obj <- readRDS(merged_path)

target_clusters <- c(
  "VZPs", "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs",
  "GABA Progenitors", "Golgi Cells", "MLIs", "iCN", 
  "RG Progenitors", "BG", "Astrocytes/Ependyma", "Ependymal Cells", 
  "OPCs", "Cycling Cells"
)

# --- 3. Subset and Re-process ---
message("Subsetting VZ lineages...")
vz_obj <- subset(obj, comb_subcluster %in% target_clusters)

# Robust factor conversion
vz_obj$comb_subcluster <- factor(vz_obj$comb_subcluster)
vz_obj$comb_subcluster <- droplevels(vz_obj$comb_subcluster)

DefaultAssay(vz_obj) <- "Xenium"
message("Re-processing integrated subset...")
vz_obj <- NormalizeData(vz_obj) %>% FindVariableFeatures() %>% ScaleData()
vz_obj <- RunPCA(vz_obj, npcs = 30, verbose = FALSE)

# --- 4. Handle Sample ID from Slurm & Setup Dirs ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No sample ID provided by Slurm.")
current_sample <- args[1]

h5ad_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_SpaTrack_h5ad")
if (!dir.exists(h5ad_dir)) dir.create(h5ad_dir, recursive = TRUE)

# --- 5. Process the Specific Slice ---
message("Starting SpaTrack analysis for: ", current_sample)
vz_slice <- subset(vz_obj, sample_id == current_sample)

if (ncol(vz_slice) > 50) {
  counts_mat <- t(as.matrix(LayerData(vz_slice, assay = "Xenium", layer = "data")))
  meta_data <- vz_slice@meta.data
  pca_coords <- Embeddings(vz_slice, reduction = "pca")[, 1:30]
  
  # --- THE ULTIMATE BYPASS: Steal Coordinates from the Unmerged File ---
  message("Loading pristine unmerged object to extract coordinates...")
  raw_sample_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS", paste0(current_sample, "_Ald_VZ_RL_QC_Subclusters.rds"))
  
  if (!file.exists(raw_sample_path)) stop("Could not find the raw unmerged RDS for this sample!")
  
  raw_obj <- readRDS(raw_sample_path)
  
  # Extract coordinates directly from the raw object
  coords_df <- GetTissueCoordinates(raw_obj, which = "centroids")
  
  # Clean up memory immediately
  rm(raw_obj)
  gc()
  
  # Force both sets of names to their raw base barcodes for perfect matching
  slice_cells <- rownames(meta_data)
  slice_base_cells <- gsub(".*_", "", slice_cells)
  rownames(coords_df) <- gsub(".*_", "", coords_df$cell)
  
  # Safely map coordinates to our integrated slice
  spatial_coords <- as.matrix(coords_df[slice_base_cells, c("x", "y")])
  rownames(spatial_coords) <- slice_cells # Restore the full SampleID_Barcode names for Python
  
  # Filter any genuine NAs
  valid_cells <- slice_cells[complete.cases(spatial_coords) & complete.cases(pca_coords)]
  
  if (length(valid_cells) > 50) {
    message("Filtered out ", nrow(meta_data) - length(valid_cells), " genuine NA cells.")
    
    # Subset all matrices strictly to the valid cells
    counts_mat <- counts_mat[valid_cells, , drop = FALSE]
    meta_data <- meta_data[valid_cells, , drop = FALSE]
    spatial_coords <- spatial_coords[valid_cells, , drop = FALSE]
    pca_coords <- pca_coords[valid_cells, , drop = FALSE]
    
    # Bridge to Python
    adata <- ad$AnnData(X = counts_mat, obs = meta_data, 
                        obsm = dict(spatial = spatial_coords, X_pca = pca_coords))
    
    py$adata <- adata
    py$save_path <- file.path(h5ad_dir, paste0("SpaTrack_", current_sample, ".h5ad"))
    
    tryCatch({
      py_run_string("
import spaTrack as spt
import numpy as np

print('Calculating Optimal Transport Matrix...')
adata.obsm['X_spatial'] = adata.obsm['spatial']
adata.obsp['trans'] = spt.get_ot_matrix(adata, data_type='spatial')

start_cells = np.where(adata.obs['comb_subcluster'] == 'VZPs')[0].tolist()

if len(start_cells) > 0:
    print('Calculating Pseudotime and Velocities...')
    adata.obs['ptime'] = spt.get_ptime(adata, start_cells)
    adata.uns['E_grid'], adata.uns['V_grid'] = spt.get_velocity(adata, basis='spatial', n_neigh_pos=100)
else:
    print('Warning: No VZP cells found in this specific slice.')

adata.write(save_path)
print(f'Successfully saved {save_path}')
      ")
    }, error = function(e) { message("Python execution failed for ", current_sample, ": ", e) })
  } else {
    message("Skipping ", current_sample, ": Not enough valid cells after NA filtering.")
  }
} else {
  message("Skipping ", current_sample, ": Not enough starting cells in VZ lineage.")
}

message("Sample ", current_sample, " complete.")