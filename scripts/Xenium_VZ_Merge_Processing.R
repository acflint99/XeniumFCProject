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

# 2. PARALLELIZATION SETUP (Optimized for 44 Cores / 440GB RAM)
# Using 'multicore' is essential for 50GB+ objects on Linux to share memory
plan("multicore", workers = 40) 
options(future.globals.maxSize = 150 * 1024^3) # 150GB limit for globals

check_mem("PIPELINE START")

# 3. LOAD DATA
merged_path <- here("outputs", "Xenium_VZ_Subsets", "Xenium_Merged_VZSubsets.rds")
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
res_list <- c(0.3, 0.5, 0.8)

# Fix the order for all your resolution columns
for(res in res_list) {
  col_name <- paste0("Xenium_snn_res.", res)
  
  # 1. Get the unique cluster names
  clusters <- unique(obj@meta.data[[col_name]])
  
  # 2. Sort them numerically (converts "10" to 10 so it knows where it goes)
  numeric_order <- clusters[order(as.numeric(as.character(clusters)))]
  
  # 3. Re-assign the column as a factor with this specific order
  obj@meta.data[[col_name]] <- factor(obj@meta.data[[col_name]], 
                                      levels = numeric_order)
}

# 8. GENERATE CLUSTERS AND EXPORT PLOTS
for(res in res_list) {
  assay_prefix <- DefaultAssay(obj)
  res_col <- paste0(assay_prefix, "_snn_res.", res)
  
  # message(Sys.time(), ": Calculating clusters for Resolution: ", res)
  # obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  # check_mem(paste0("POST-RES-", res))
  
  # Generate Plot with Rastering for performance
  message(Sys.time(), ": Generating UMAP plot for Resolution: ", res)
  p <- DimPlot(obj, 
               reduction = "umap_harmony", 
               group.by = res_col, 
               label = TRUE, 
               label.size = 5,
               label.box = TRUE,   # Adds a small box behind labels for readability
               raster = TRUE, 
               pt.size = 0.6,      # Increased for visibility
               alpha = 0.8) +      # Makes dots more solid
    ggtitle(paste0("Refined VZ UMAP - Resolution ", res)) +
    theme_classic() +            # Standard axes, no grid lines
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 3. Save the plot
  plot_name <- paste0("XenAld_VZ_UMAP_Res_", res, "_Refined.png")
  plot_path <- here("outputs", "XenAld_VZ_Plots", plot_name)
  
  message(Sys.time(), ": Saving plot to ", plot_path)
  ggsave(plot_path, p, width = 10, height = 8, dpi = 300)
  
  # Remove plot object and trigger garbage collection to keep RStudio snappy
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

ggsave(here("outputs", "XenAld_VZ_Plots", "XenAld_VZ_Batch_Comp_UMAP.png"), 
       p1 + p2, width = 16, height = 7, dpi = 300)

# 2. Generate the Plot
p_orig <- DimPlot(obj, 
                  reduction = "umap_harmony", 
                  group.by = "cluster_weighted", # Use your original label column here
                  label = TRUE, 
                  label.size = 4,
                  label.box = TRUE,       # Makes original labels easier to see
                  raster = TRUE, 
                  pt.size = 0.5, 
                  alpha = 0.8) +
  ggtitle("Xenium Aldinger VZ Integrated UMAP: Original Cluster Weighted Labels") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# 3. Save the Plot
ggsave(here("outputs", "XenAld_VZ_Plots", "XenAld_VZ_OrigClusterWeighted_UMAP.png"), 
       p_orig, width = 12, height = 9, dpi = 300)

# 9. SAVE FINAL RESULT
output_path <- here("outputs", "XenAld_VZ_RDS", "Xenium_VZ_Integrated.rds")
saveRDS(obj, output_path, compress = FALSE)
check_mem("PIPELINE COMPLETE - FILE SAVED")