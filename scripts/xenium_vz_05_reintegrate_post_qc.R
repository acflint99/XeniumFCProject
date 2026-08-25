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
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id

check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 1. PARALLELIZATION 
plan(multisession, workers = 7) 
options(future.globals.maxSize = 200 * 1024^3)

plot_path <- here("outputs", "xenium", "vz", "05_post_qc", "plots")
if(!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

merged_path <- here("outputs", "xenium", "vz", "03_integrated", "rds", "Xenium_VZ_Res1.5.rds")
removal_path <- here("outputs", "xenium", "vz", "04_qc", "tables", "XenAld_VZ_QC_removal_manifest.csv")
review_path <- here("outputs", "xenium", "vz", "04_qc", "tables", "XenAld_VZ_QC_review_manifest.csv")
if (!file.exists(removal_path)) stop("Reviewed VZ removal manifest not found: ", removal_path)
if (!file.exists(review_path)) stop("VZ QC review manifest not found: ", review_path)
removal_manifest <- read.csv(removal_path, stringsAsFactors = FALSE)
if (nrow(removal_manifest) != length(sample_ids) ||
    !setequal(removal_manifest$sample_id, sample_ids)) {
  stop("VZ removal manifest must contain exactly the configured 34 samples.")
}
removal_decisions <- unique(as.character(removal_manifest$clusters_removed))
if (length(removal_decisions) != 1L) stop("VZ removal manifest contains inconsistent decisions.")
clusters_to_remove <- if (identical(removal_decisions, "none")) {
  character()
} else strsplit(removal_decisions, "|", fixed = TRUE)[[1]]
obj <- readRDS(merged_path)

review <- read.csv(review_path, stringsAsFactors = FALSE)
merged_info <- file.info(merged_path)
review_is_current <- nrow(review) == 1L &&
  isTRUE(all.equal(as.numeric(review$merged_size), as.numeric(merged_info$size))) &&
  isTRUE(all.equal(as.numeric(review$merged_mtime), as.numeric(merged_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj)))
if (!review_is_current) stop("VZ merged object changed after QC review. Repeat the QC decision stage.")

# Join the manifest-defined sample layers into one unified matrix.
obj <- JoinLayers(obj)

# 2. SET THE IDENTITY
Idents(obj) <- "Xenium_snn_res.0.8"

# Apply the explicit decision recorded after QC review.
cells_before_qc <- ncol(obj)
if (length(clusters_to_remove)) {
  unknown_clusters <- setdiff(clusters_to_remove, levels(Idents(obj)))
  if (length(unknown_clusters)) stop("Reviewed VZ cluster(s) are absent: ", paste(unknown_clusters, collapse = ", "))
  obj <- subset(obj, idents = clusters_to_remove, invert = TRUE)
}
removed_cells <- cells_before_qc - ncol(obj)
if (removed_cells != sum(removal_manifest$removed_cells)) {
  stop("Merged VZ removal count does not match the reviewed per-sample removal manifest.")
}

# 1. Re-run PCA on the subset
obj <- RunPCA(obj, verbose = FALSE, reduction.name = "pca_clean", npcs = 50)

# 2. Re-run Harmony (If you used it originally)
# You must re-integrate because the batch-effect vectors change when 31k cells leave
obj <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca_clean", reduction.save = "harmony_clean")

# 3. Re-run UMAP
obj <- RunUMAP(obj, reduction = "harmony_clean", dims = 1:50, n.neighbors = 100, reduction.name = "umap_clean")

obj <- FindNeighbors(obj, reduction = "harmony_clean", dims = 1:50, k.param = 30, verbose = TRUE)

# Define resolutions for the loop
res_list <- c(0.3, 0.5, 0.8, 1)

# 8. GENERATE CLUSTERS (RUN IN BULK FOR SPEED)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters in parallel...")

# Running all resolutions at once is faster with future/parallelization
obj <- FindClusters(obj, resolution = res_list, verbose = TRUE)
check_mem("POST-BATCH-CLUSTERING")

p1 <- DimPlot(obj, 
              reduction = "umap_clean", 
              group.by = "consensus_label", 
              label = TRUE, 
              label.size = 5,
              label.box = TRUE,
              raster = TRUE, 
              pt.size = 0.6,
              alpha = 0.8) + 
  ggtitle(paste0("Refined VZ UMAP - Original Clusters")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_VZ_PostQC_OrigCluster_UMAP.tif")),
  width = 10,
  height = 8,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()

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
    ggtitle(paste0("Refined VZ UMAP - Resolution ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  Cairo::CairoTIFF(
    filename = file.path(plot_path, paste0("XenAld_VZ_PostQC_UMAP_Res_", res, ".tif")),
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

output_path <- here("outputs", "xenium", "vz", "05_post_qc", "rds")
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

saveRDS(obj, file.path(output_path, "Xenium_VZ_postQC_Res1.5_4-2-26.rds"), compress = FALSE)
