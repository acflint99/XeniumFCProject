# ==============================================================================
# Script: Combine_VZ_RL_Subclusters.R
# Purpose: Merge specific subcluster metadata for Xenium cerebellum samples
# ==============================================================================

rm(list = ls())

library(Seurat)
library(here)
library(dplyr)
library(ggplot2)

# Load your new palette and order
source(here("scripts", "color_palette.R"))

# 1. SETUP DIRECTORIES ---------------------------------------------------------
out_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

table_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Tables")
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)


# 2. DEFINE THE FUNCTION -------------------------------------------------------
process_xenium_sample <- function(sample_id) {
  
  filename <- paste0(sample_id, "_Ald_VZ_RL_QC_Subclusters.rds")
  path <- here("outputs", "Xenium_AldingerABT_VZ&RLsubclusters_QC_Res1.5_RDS", filename)
  
  if (!file.exists(path)) {
    warning("File missing for: ", sample_id)
    return(NULL)
  }
  
  message("--- Processing: ", sample_id, " ---")
  obj <- readRDS(path)
  
  # Combine columns: Priority VZ -> RL -> Weighted
  obj@meta.data <- obj@meta.data %>%
    mutate(
      comb_subcluster = coalesce(VZ_subcluster, RL_subcluster, cluster_weighted),
      # Set the factor levels to your 32-label list
      # intersect() ensures it only sets levels for labels actually present in this sample
      comb_subcluster = factor(comb_subcluster, levels = intersect(master_subcluster_order, unique(comb_subcluster)))
    )
  
  # Set Identity for downstream plotting
  Idents(obj) <- "comb_subcluster"
  
  # --- NEW: COUNTING AND TABLES -----------------------------------------------
  
  # 1. Print total cell count to console
  total_cells <- ncol(obj)
  message("Total cells in ", sample_id, ": ", total_cells)
  
  # 2. Calculate counts for cluster_weighted
  weighted_counts <- as.data.frame(table(obj$cluster_weighted))
  colnames(weighted_counts) <- c("Cluster_Weighted", "Cell_Count")
  
  # 3. Calculate counts for comb_subcluster
  comb_counts <- as.data.frame(table(obj$comb_subcluster))
  colnames(comb_counts) <- c("Combined_Subcluster", "Cell_Count")
  
  # 4. Save tables to CSV
  write.csv(weighted_counts, 
            file = file.path(table_dir, paste0(sample_id, "_counts_clusterweighted.csv")), 
            row.names = FALSE)
  
  write.csv(comb_counts, 
            file = file.path(table_dir, paste0(sample_id, "_counts_VZ&RLsubclusters.csv")), 
            row.names = FALSE)
  
  # ----------------------------------------------------------------------------
  
  # Save Plot
  p <- ImageDimPlot(obj, group.by = "comb_subcluster", cols = subcluster_palette, size = 0.6) + 
    ggtitle(paste(sample_id, "Combined Subclusters"))
  
  ggsave(file.path(plot_dir, paste0(sample_id, "_combSubcluster_GlobalSpatial0.6.tif")), 
         plot = p, device = "tiff", width = 14, height = 10, dpi = 600, bg = "black", compression = "lzw")
  
  p2 <- ImageDimPlot(obj, group.by = "comb_subcluster", cols = subcluster_palette, size = 0.85) + 
    ggtitle(paste(sample_id, "Combined Subclusters"))
  
  ggsave(file.path(plot_dir, paste0(sample_id, "_combSubcluster_GlobalSpatial0.85.tif")), 
         plot = p2, device = "tiff", width = 14, height = 10, dpi = 600, bg = "black", compression = "lzw")
  
  # Generate the plot
  # Note: Ensure 'cluster_colors' is defined in your environment beforehand
  p_wei <- ImageDimPlot(
    obj, 
    group.by = "cluster_weighted", 
    size = 0.85, 
    cols = cluster_colors
  ) + 
    ggtitle(paste(sample_id, "-Post QC (Weighted)"))
  
  # Save with a specific background color to ensure ggsave doesn't add white
  ggsave(
    filename = file.path(plot_dir, paste0(sample_id, "_postQC_GlobalSpatial_plot0.85.tif")), 
    plot = p_wei,
    device = "tiff",
    width = 12,
    height = 10,
    dpi = 600,
    bg = "black",  # Forces the saved file background to be black
    compression = "lzw"
  )
  
  # Generate the plot
  # Note: Ensure 'cluster_colors' is defined in your environment beforehand
  p_wei2 <- ImageDimPlot(
    obj, 
    group.by = "cluster_weighted", 
    size = 0.6, 
    cols = cluster_colors
  ) + 
    ggtitle(paste(sample_id, "-Post QC (Weighted)"))
  
  # Save with a specific background color to ensure ggsave doesn't add white
  ggsave(
    filename = file.path(plot_dir, paste0(sample_id, "_postQC_GlobalSpatial_plot0.6.tif")), 
    plot = p_wei2,
    device = "tiff",
    width = 12,
    height = 10,
    dpi = 600,
    bg = "black",  # Forces the saved file background to be black
    compression = "lzw"
  )
  
  # Save RDS (Uncompressed for speed)
  saveRDS(obj, file = file.path(out_dir, filename), compress = FALSE)
  
  return(TRUE) # Return true to indicate success
}


# 3. EXECUTION LOOP ------------------------------------------------------------
samples <- c(
  "GZFB4_X_G", "FB124_X_G", "FB198_X_G", "FB328_1_X_G",
  "FB330_1_X_G", "FB78_X_G", "GZFB5_X_G", "GZFB_12_X_G_1"
  #"GZFB_12_X_G_2", "GZFB_12_X_G_3", "GZFB_12_X_G_4", "GZFB_12_X_G_5",
  #"GZFB_1_X_G", "GZFB_9_X_G_1", "GZFB_9_X_G_2", "GZFB_9_X_G_3"
)

# Iterate through all samples
for (s in samples) {
  process_xenium_sample(s)
  gc() # Garbage collection to free up RAM after each sample
}

message("DONE: All samples processed.")