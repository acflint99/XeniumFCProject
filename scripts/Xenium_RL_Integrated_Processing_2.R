# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# Load your new palette and order
source(here("scripts", "color_palette.R"))


# --- HELPER: MEMORY MONITOR ---
check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 2. PARALLELIZATION SETUP
plan("multicore", workers = 40) 
options(future.globals.maxSize = 150 * 1024^3) 

check_mem("PIPELINE START")

plot_dir <- here("outputs", "XenAld_RL_Integrated_rmOther_Plots")
rds_dir  <- here("outputs", "XenAld_RL_Integrated_rmOther_RDS")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(rds_dir))  dir.create(rds_dir, recursive = TRUE)

# 3. LOAD DATA
obj_path <- here("outputs", "XenAld_RL_Integrated_RDS", "Xenium_RL_Integrated_newSubclusters.rds")
obj <- readRDS(obj_path)
check_mem("DATA LOADED")

obj <- subset(obj, subset = RL_subcluster %in% c("Purkinje Cells", "RG Progenitors", "Astrocytes"), invert = TRUE)

# 4. STANDARD WORKFLOW
message("Normalizing and Finding Variable Features...")
obj <- NormalizeData(obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message("Scaling data...")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = TRUE)
check_mem("POST-SCALE")

message("Running PCA...")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# 5. UMAP & NEIGHBORS (Based on PCA, no Harmony)
message("Running UMAP and Neighbors on PCA...")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, 
               reduction.name = "umap", verbose = TRUE)

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, verbose = TRUE)
check_mem("POST-REDUCTION")

# 6. GENERATE CLUSTERS
res_list <- c(0.3, 0.5, 0.8)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters in parallel...")

obj <- FindClusters(obj, resolution = res_list, verbose = FALSE)
check_mem("POST-BATCH-CLUSTERING")

# 7. LOOP FOR SORTING AND PLOTTING
for(res in res_list) {
  res_col <- paste0(assay_prefix, "_snn_res.", res)
  
  # Sort numerically
  clusters <- unique(na.omit(obj@meta.data[[res_col]]))
  numeric_order <- sort(as.numeric(as.character(clusters)))
  
  obj@meta.data[[res_col]] <- factor(as.character(obj@meta.data[[res_col]]), 
                                     levels = as.character(numeric_order))
  
  message(Sys.time(), ": Generating & Saving UMAP for Resolution: ", res)
  
  p <- DimPlot(obj, 
               reduction = "umap", 
               group.by = res_col, 
               label = TRUE, 
               label.size = 5,
               label.box = TRUE,
               raster = TRUE, 
               pt.size = 0.6,
               alpha = 0.8) + 
    ggtitle(paste0("RL UMAP (Remove Other Subclusters) - Resolution ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename = file.path(plot_dir, paste0("XenAld_RL_rmOther_recluster_Res_", res, "UMAP.png")), 
         p, width = 10, height = 8, dpi = 300)
  
  rm(p)
  gc()
}

# 8. VISUALIZATION
message("Saving Subcluster UMAP Plot...")
p_batch <- DimPlot(obj, reduction = "umap", group.by = "RL_subcluster", raster = TRUE, label = TRUE, cols = rl_palette) + 
  ggtitle("RL UMAP (Remove Other Subclusters) - RL Subclusters") +
  theme_classic()

ggsave(filename = file.path(plot_dir, "XenAld_RL_rmOther_RLsubclusters_UMAP.png"), 
       p_batch, width = 10, height = 8, dpi = 300)

# Original Labels Plot
p_orig <- DimPlot(obj, 
                  reduction = "umap", 
                  group.by = "cluster_weighted", 
                  label = TRUE, 
                  label.size = 4,
                  label.box = TRUE,
                  raster = TRUE, 
                  pt.size = 0.5, 
                  alpha = 0.8) +
  ggtitle("RL Integrated UMAP: Original Labels (Remove Other Subclusters)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

ggsave(filename = file.path(plot_dir, "XenAld_RL_OrigLabels_rmOther_UMAP.png"), 
       p_orig, width = 12, height = 9, dpi = 300)

# 9. SAVE FINAL RESULT
output_path <- file.path(rds_dir, "Xenium_RL_Integrated_rmOther_31426.rds")
# Keeping compress = FALSE for faster I/O since you have 440GB RAM
saveRDS(obj, output_path, compress = FALSE)

plan("sequential")
message("Workflow Complete.")