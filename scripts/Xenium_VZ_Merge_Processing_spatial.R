# ==============================================================================
# Script: BANKSY Spatial Clustering for Xenium VZ Refinement
# Purpose: Perform spatially-aware clustering on merged VZ subsets
# Date: 2026-03-17
# ==============================================================================

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(Banksy)
library(SingleCellExperiment) # Added for spatialCoords()
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# --- HELPER: MEMORY MONITOR ---
check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 2. DELAYED PARALLELIZATION 
plan("sequential") 
options(future.globals.maxSize = 200 * 1024^3)

check_mem("PIPELINE START")

# Ensure directories exist
plot_dir <- here("outputs", "XenAld_VZ_Banksy_Plots")
rds_dir  <- here("outputs", "XenAld_VZ_Banksy_RDS")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(rds_dir))  dir.create(rds_dir, recursive = TRUE)

# 3. LOAD DATA
merged_path <- here("outputs", "XenAld_VZ_Subsets_RDS", "Merged", "Xenium_Merged_VZSubsets_spatial_31726.rds")
obj <- readRDS(merged_path)
check_mem("DATA LOADED")

# 4. REVISED STABLE PRE-PROCESSING
message("Step 4a: Normalizing...")
obj <- NormalizeData(obj, verbose = TRUE)
gc()

message("Step 4b: Finding Variable Features...")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
gc()

message("Step 4c: Scaling Data...")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
gc()

message("Step 4d: Running PCA...")
obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
gc()

message("Initialization of Parallel Workers for BANKSY...")
plan("multisession", workers = 20) 

# 5. CONVERT & COMPUTE BANKSY
message("Joining Seurat v5 layers before SCE conversion...")
obj[["Xenium"]] <- JoinLayers(obj[["Xenium"]])

message("Converting to SingleCellExperiment for Banksy...")
sce <- as.SingleCellExperiment(obj, assay = "Xenium")

message("Splitting dataset by sample to prevent cross-sample spatial edges...")
# Extract sample names (converting to character prevents factor-indexing bugs)
sample_names <- as.character(unique(sce$orig.ident))

# Split into a list so Banksy only finds neighbors within the same physical tissue
sce_list <- lapply(sample_names, function(x) sce[, sce$orig.ident == x])

message("Computing Banksy Neighborhood Matrices per sample...")
sce_list <- lapply(sce_list, Banksy::computeBanksy, 
                   assay_name = "logcounts", 
                   coord_names = c("x_centroid", "y_centroid"), 
                   compute_agf = FALSE, 
                   k_geom = 15)

message("Rebinding samples for joint analysis...")
sce <- do.call(cbind, sce_list)

# Clear out the list immediately to free up RAM
rm(sce_list)
gc()

# 6. BANKSY PCA
message("Running Banksy PCA with multi-sample scaling...")
# FIX: Use 'group' to scale samples separately, mitigating batch effects before joint clustering
sce <- Banksy::runBanksyPCA(sce, 
                            use_agf = FALSE, 
                            lambda = 0.1,
                            group = "orig.ident")
gc()

# 7. BANKSY CLUSTERING
message("Running Banksy Clustering...")
sce <- Banksy::clusterBanksy(sce, 
                             use_agf = FALSE, 
                             lambda = 0.1, 
                             resolution = c(0.5, 0.8, 1))

# --- IMPORTANT: MOVE RESULTS BACK TO SEURAT ---
message("Syncing Banksy results back to Seurat object...")

pca_name <- grep("PCA_M0", reducedDimNames(sce), value = TRUE)[1]
banksy_pca_matrix <- reducedDim(sce, pca_name)

# CRITICAL FIX: `cbind` grouped our cells by sample, changing their order. 
# We MUST reorder the PCA matrix to match the original Seurat object 
# exactly before creating the DimReduc.
banksy_pca_matrix <- banksy_pca_matrix[colnames(obj), ]

obj[["banksy_pca"]] <- CreateDimReducObject(embeddings = banksy_pca_matrix, 
                                            key = "banksypca_", 
                                            assay = "Xenium")

banksy_metadata <- as.data.frame(colData(sce))

# Reorder metadata to match Seurat as well
banksy_metadata <- banksy_metadata[colnames(obj), , drop = FALSE]
clust_cols <- grep("clust", colnames(banksy_metadata), value = TRUE)

obj <- AddMetaData(obj, metadata = banksy_metadata[, clust_cols, drop = FALSE])

# Clean up to dump the remaining 15GB+ of SCE data
rm(sce)
gc()

# 7.5 RUN UMAP
message("Running UMAP on Banksy PCA...")
obj <- RunUMAP(obj, reduction = "banksy_pca", dims = 1:20, reduction.name = "umap_banksy", min.dist = 0.5, spread = 1)

# 8. PLOTTING & VISUALIZATION
res_list <- c(0.5, 0.8, 1)
lambda_val <- 0.1

for(res in res_list) {
  # FIX: Dynamically find the column to account for Banksy injecting the 'k' parameter
  search_pattern <- paste0("clust_M0_lam", lambda_val, ".*_res", res)
  res_col <- grep(search_pattern, colnames(obj@meta.data), value = TRUE)[1]
  
  if(is.na(res_col) || !res_col %in% colnames(obj@meta.data)) {
    message("Column matching ", search_pattern, " not found. Skipping plot for res ", res)
    next
  }
  
  message(Sys.time(), ": Generating UMAP for Banksy Res: ", res, " (Column: ", res_col, ")")
  
  p <- DimPlot(obj, 
               reduction = "umap_banksy", 
               group.by = res_col, 
               label = TRUE, 
               label.size = 5,
               label.box = TRUE,
               raster = TRUE, 
               pt.size = 0.3,
               alpha = 0.5) + 
    ggtitle(paste0("VZ Banksy UMAP (Lambda ", lambda_val, ") - Res ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename = file.path(plot_dir, paste0("XenAld_VZ_Banksy_UMAP_Res_", res, ".png")), 
         p, width = 10, height = 8, dpi = 300)
}

# 9. Orig Cluster CHECK PLOT
p_batch <- DimPlot(obj, reduction = "umap_banksy", group.by = "cluster_weighted", raster = TRUE) + 
  ggtitle("Banksy UMAP: Colored by Orig Cluster (weighted)")

ggsave(filename = file.path(plot_dir, "XenAld_VZ_Banksy_OrigClusterWeighted_UMAP.png"), 
       p_batch, width = 10, height = 8, dpi = 300)

# 9. BATCH CHECK PLOT
p_batch <- DimPlot(obj, reduction = "umap_banksy", group.by = "orig.ident", raster = TRUE) + 
  ggtitle("Banksy UMAP: Colored by Sample")

ggsave(filename = file.path(plot_dir, "XenAld_VZ_Banksy_Batch_UMAP.png"), 
       p_batch, width = 10, height = 8, dpi = 300)


# 10. SAVE FINAL RESULT
output_path <- file.path(rds_dir, "Xenium_VZ_Banksy_Integrated_lambda0.1_31826.rds")
message("Saving final object to: ", output_path)
saveRDS(obj, output_path, compress = FALSE)

plan("sequential")
check_mem("PIPELINE COMPLETE")


# 1. Define the specific Banksy resolution column you want to visualize
# (Update this to whichever resolution gave you those clean groups in your UMAP)
cluster_col <- "clust_M0_lam0.1_k50_res1" 

meta <- obj@meta.data

# 3. Get the names of the original samples from orig.ident
sample_names <- c("FB78_X_G", "FB330_1_X_G", "GZFB_12_X_G_1")
message("Found ", length(sample_names), " tissue sections in metadata.")

# 4. Loop through each sample and plot using pure ggplot2
for(sample in sample_names) {
  message("Generating Spatial Plot for: ", sample)
  
  # Subset the cells that belong to this specific tissue section
  plot_data <- meta[meta$orig.ident == sample, ]
  
  # Build the spatial plot
  p <- ggplot(plot_data, aes(x = x_centroid, y = y_centroid, color = .data[[cluster_col]])) +
    geom_point(size = 0.5, stroke = 0) + # stroke = 0 prevents borders from blurring the colors
    coord_fixed() +                      # CRITICAL: Keeps aspect ratio 1:1 so the tissue isn't stretched
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_blank(),       # Hide coordinate numbers (not visually useful)
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    ) +
    # Make the legend dots large enough to actually see the colors
    guides(color = guide_legend(override.aes = list(size = 4), title = "Cluster")) +
    ggtitle(paste("Banksy Spatial Mapping -", sample))
  
  # Save the plot
  save_path <- file.path(plot_dir, paste0("XenAld_Spatial_Manual_", sample, "_Res1.png"))
  ggsave(filename = save_path, plot = p, width = 10, height = 8, dpi = 300)
}