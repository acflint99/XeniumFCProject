# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

library(Seurat)
library(here)
library(ggplot2)
library(reticulate)

options(future.globals.maxSize = 100 * 1024^3)

# --- 1. Python/Conda Setup ---
env_python <- "/home/acflint/.conda/envs/spatrack_env/bin/python"
Sys.unsetenv("PYTHONPATH")
Sys.unsetenv("PYTHONHOME")
Sys.setenv(RETICULATE_PYTHON = env_python)

ad <- import("anndata")
message("Python environment loaded.")

# --- 2. Load and Subset the Seurat Object ---
merged_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS", "XenAld_VZ_RL_QC_Subclusters_Spatial_Merged_Integrated2_4-23-26.rds")
message("Loading full integrated object...")
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

vz_obj$comb_subcluster <- factor(vz_obj$comb_subcluster)
vz_obj$comb_subcluster <- droplevels(vz_obj$comb_subcluster)

# Clean up memory
rm(obj)
gc()

message("Calculating zoomed-in UMAP for VZ lineages...")
vz_obj <- RunUMAP(vz_obj, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.vz")

vz_obj$VZ_ptime <- NA

# --- 4. Loop Through the h5ad Files and Import ptime (MEMORY SAFE) ---
h5ad_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_SpaTrack_h5ad")
h5ad_files <- list.files(h5ad_dir, pattern = "\\.h5ad$", full.names = TRUE)

message("Importing pseudotime metadata using backed mode...")

for (file_path in h5ad_files) {
  message("Processing: ", basename(file_path))
  
  # Read the metadata WITHOUT loading the massive matrices into RAM
  adata <- ad$read_h5ad(file_path, backed = "r")
  obs_df <- py_to_r(adata$obs)
  
  if ("ptime" %in% colnames(obs_df)) {
    valid_cells <- rownames(obs_df)
    vz_obj@meta.data[valid_cells, "VZ_ptime"] <- obs_df$ptime
    message(" -> Successfully mapped ptime for ", length(valid_cells), " cells.")
  } else {
    message(" -> No ptime found - skipping.")
  }
  
  rm(adata, obs_df)
  gc()
  py_run_string("import gc; gc.collect()")
}

# --- 5. Native Seurat Visualization ---
plot_out_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_VZ_SpaTrack_Plots")
if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)

# 5A. Plot Pseudotime on the ISOLATED VZ UMAP
message("Generating UMAP...")
p1 <- FeaturePlot(vz_obj, features = "VZ_ptime", reduction = "umap.vz", raster = TRUE) +
  scale_color_viridis_c(option = "magma", na.value = "lightgrey") +
  ggtitle("Ventricular Zone Lineage: Pseudotime Trajectory")

ggsave(file.path(plot_out_dir, "VZ_Lineage_UMAP_Pseudotime.png"), p1, width = 10, height = 10)

# --- 5B. Spatial Feature Plot (The Safe Way) ---
sample_to_plot <- "FB124_X_G"
message("Loading pristine raw object for spatial mapping: ", sample_to_plot)

# Load the pristine, unmerged sample object
raw_sample_path <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS", paste0(sample_to_plot, "_Ald_VZ_RL_QC_Subclusters.rds"))
raw_obj <- readRDS(raw_sample_path)

# Extract the pseudotime data just for this sample from our VZ object
cells_for_this_sample <- Cells(vz_obj)[vz_obj$sample_id == sample_to_plot]
ptime_data <- vz_obj@meta.data[cells_for_this_sample, "VZ_ptime", drop = FALSE]

# Strip the prefix so the rownames match the untouched raw_obj
rownames(ptime_data) <- gsub(paste0("^", sample_to_plot, "_"), "", rownames(ptime_data))

# Inject the matched pseudotime data into the pristine raw object
raw_obj <- AddMetaData(raw_obj, metadata = ptime_data)

message("Generating Spatial Plot...")
p2 <- SpatialFeaturePlot(raw_obj, features = "VZ_ptime", pt.size.factor = 1.5) +
  scale_fill_viridis_c(option = "magma", na.value = "black") +
  ggtitle(paste("Spatial Pseudotime:", sample_to_plot))

ggsave(file.path(plot_out_dir, paste0(sample_to_plot, "_Spatial_Pseudotime.png")), p2, width = 10, height = 10)

# Clean up
rm(raw_obj)
gc()

message("Spatial plotting complete!")