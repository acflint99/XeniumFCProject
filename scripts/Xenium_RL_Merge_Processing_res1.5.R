# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(harmony)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# --- HELPER: MEMORY MONITOR ---
# This prints the current RAM usage to your console at every step
check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 2. PARALLELIZATION SETUP
#plan(multisession, workers = 4) 
plan("sequential")
options(future.globals.maxSize = 200 * 1024^3)

check_mem("PIPELINE START")

# Ensure both directories exist
plot_dir <- here("outputs", "XenAld_RL_Res1.5_Plots")
rds_dir  <- here("outputs", "XenAld_RL_Res1.5_RDS")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(rds_dir))  dir.create(rds_dir, recursive = TRUE)

# 3. LOAD DATA
merged_path <- here("outputs", "XenAld_RL_Subsets_Res1.5_RDS", "Merged", "Xenium_Merged_RLSubsets.rds")
obj <- readRDS(merged_path)
check_mem("DATA LOADED")

# 4. STANDARD WORKFLOW
message("Normalizing and Finding Variable Features...")
obj <- NormalizeData(obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message("Scaling data...")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
check_mem("POST-SCALE")

message("Running PCA...")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# 5. BATCH CHECK (UNCORRECTED UMAP)
message("Generating uncorrected UMAP...")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, 
               reduction.name = "umap_uncorrected", verbose = TRUE)
check_mem("POST-UNCORRECTED UMAP")

# 6. RUN HARMONY
message("Running Harmony integration...")
obj <- RunHarmony(obj, group.by.vars = "orig.ident", 
                  dims.use = 1:30, reduction.save = "harmony", verbose = TRUE)
check_mem("POST-HARMONY")

# 7. POST-HARMONY REDUCTION & CLUSTERING
message("Running Integrated UMAP and Neighbors...")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, 
               reduction.name = "umap_harmony", verbose = TRUE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = TRUE)

# Define resolutions for the loop
res_list <- c(0.5, 0.8, 1)

# 8. GENERATE CLUSTERS (RUN IN BULK FOR SPEED)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters in parallel...")

# Running all resolutions at once is faster with future/parallelization
obj <- FindClusters(obj, resolution = res_list, verbose = FALSE)
check_mem("POST-BATCH-CLUSTERING")

# Now loop only for sorting and plotting
for(res in res_list) {
  # Dynamically build the column name so it ALWAYS matches
  res_col <- paste0(assay_prefix, "_snn_res.", res)
  
  # 1 & 2. Sort numerically
  clusters <- unique(na.omit(obj@meta.data[[res_col]]))
  numeric_order <- sort(as.numeric(as.character(clusters)))
  
  # 3. Re-assign (using res_col, not the hardcoded Xenium string)
  obj@meta.data[[res_col]] <- factor(as.character(obj@meta.data[[res_col]]), 
                                     levels = as.character(numeric_order))
  
  message(Sys.time(), ": Generating & Saving UMAP for Resolution: ", res)
  
  p <- DimPlot(obj, 
               reduction = "umap_harmony", 
               group.by = res_col, 
               label = TRUE, 
               label.size = 5,
               label.box = TRUE,
               raster = TRUE, 
               pt.size = 0.6,
               alpha = 0.8) + 
    ggtitle(paste0("Refined RL UMAP - Resolution ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  Cairo::CairoTIFF(
    filename = file.path(plot_dir, paste0("XenAld_RL_UMAP_Res_", res, ".tif")),
    width = 10,
    height = 8,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
  
  rm(p)
  gc()
}

# 8. VISUALIZATION (Using Raster to prevent RStudio lag)
message("Saving comparison plots...")
# raster = TRUE converts points to pixels; essential for >100k cells
p1 <- DimPlot(obj, reduction = "umap_uncorrected", group.by = "orig.ident", raster = TRUE) + 
  NoLegend() + ggtitle("Pre-Harmony")

p2 <- DimPlot(obj, reduction = "umap_harmony", group.by = "orig.ident", raster = TRUE) + 
  ggtitle("Post-Harmony")

Cairo::CairoTIFF(
  filename = file.path(plot_dir, "XenAld_RL_Batch_Comp_UMAP.tif"),
  width = 16,
  height = 7,
  units = "in",
  res = 600
)
print(p1 + p2)
grDevices::dev.off()

# 2. Generate the Plot
p_orig <- DimPlot(obj, 
                  reduction = "umap_harmony", 
                  group.by = "consensus_label",
                  label = TRUE, 
                  label.size = 4,
                  label.box = TRUE,       # Makes original labels easier to see
                  raster = TRUE, 
                  pt.size = 0.5, 
                  alpha = 0.8) +
  ggtitle("Xenium RL Integrated UMAP: Consensus Labels") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# 3. Save the Plot
Cairo::CairoTIFF(
  filename = file.path(plot_dir, "Xenium_RL_ConsensusLabel_UMAP.tif"),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p_orig)
grDevices::dev.off()

# 9. SAVE FINAL RESULT
# Final RDS Save
output_path <- file.path(rds_dir, "Xenium_RL_Res1.5_33026.rds")
saveRDS(obj, output_path, compress = FALSE)

check_mem("PIPELINE COMPLETE - FILE SAVED")
