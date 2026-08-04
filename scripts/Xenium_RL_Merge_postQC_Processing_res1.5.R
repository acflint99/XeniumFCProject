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

# Load your new palette and order
source(here("scripts", "color_palette.R"))

check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 1. PARALLELIZATION 
plan(multisession, workers = 4) 
options(future.globals.maxSize = 200 * 1024^3)

plot_path <- here("outputs", "XenAld_RL_postQC_Res1.5_Plots")
if(!dir.exists(plot_path)) dir.create(plot_path)

table_path <- here("outputs", "XenAld_RL_postQC_Res1.5_Tables")
if(!dir.exists(table_path)) dir.create(table_path, recursive = TRUE)

merged_path <- here("outputs", "XenAld_RL_Res1.5_RDS", "Xenium_RL_Res1.5_33026.rds") #check/change date!
obj <- readRDS(merged_path)

# This merges the 15 separate sample layers into one unified matrix
obj <- JoinLayers(obj)

# 2. SET THE IDENTITY
Idents(obj) <- "Xenium_snn_res.0.8" #check change resolution

# This removes low QC cells from the 'obj' variable entirely for the rest of the script
obj <- subset(obj, idents = "7", invert = TRUE) #check/change identity to remove

# 1. Re-run PCA on the subset
obj <- RunPCA(obj, verbose = FALSE, reduction.name = "pca_clean", npcs = 50)

# 2. Re-run Harmony (If you used it originally)
# You must re-integrate because the batch-effect vectors change when 31k cells leave
obj <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca_clean", reduction.save = "harmony_clean")

# 3. Re-run UMAP
obj <- RunUMAP(obj, reduction = "harmony_clean", dims = 1:50, n.neighbors = 50, reduction.name = "umap_clean")

obj <- FindNeighbors(obj, reduction = "harmony_clean", dims = 1:50, k.param = 30, verbose = TRUE)

# Define resolutions for the loop
res_list <- c(0.3, 0.5, 0.8)

# 8. GENERATE CLUSTERS (RUN IN BULK FOR SPEED)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters in parallel...")

# Running all resolutions at once is faster with future/parallelization
obj <- FindClusters(obj, resolution = res_list, verbose = TRUE)
check_mem("POST-BATCH-CLUSTERING")

p1 <- DimPlot(obj, 
              reduction = "umap_clean", 
              group.by = "cluster_weighted", 
              label = TRUE, 
              label.size = 5,
              label.box = TRUE,
              raster = TRUE, 
              pt.size = 0.6,
              alpha = 0.8) + 
  ggtitle(paste0("Refined RL UMAP - Original Clusters")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = file.path(plot_path, paste0("XenAld_RL_PostQC_OrigCluster_UMAP.tif")), 
       plot = p1, device = "tiff", width = 10, height = 8, dpi = 600, compression = "lzw")

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
               reduction = "umap_clean", 
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
  
  ggsave(filename = file.path(plot_path, paste0("XenAld_RL_PostQC_UMAP_Res_", res, ".tif")), 
         plot = p, device = "tiff", width = 10, height = 8, dpi = 600, compression = "lzw")
  
  rm(p)
  gc()
}

output_path <- here("outputs", "XenAld_RL_postQC_Res1.5_RDS")
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

saveRDS(obj, file.path(output_path, "Xenium_RL_postQC_Res1.5_4-2-26.rds"), compress = FALSE)