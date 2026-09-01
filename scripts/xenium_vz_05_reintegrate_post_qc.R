#!/usr/bin/env Rscript

# Reintegrate the reviewed VZ object after the explicit stage-4 QC decision.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_05_reintegrate_post_qc.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "vz", "05_post_qc", "plots")
rds_dir <- file.path(output_root, "xenium", "vz", "05_post_qc", "rds")
table_dir <- file.path(output_root, "xenium", "vz", "05_post_qc", "tables")
merged_path <- file.path(output_root, "xenium", "vz", "03_integrated", "rds", "Xenium_VZ_Res1.5.rds")
removal_path <- file.path(output_root, "xenium", "vz", "04_qc", "tables", "XenAld_VZ_QC_removal_manifest.csv")
review_path <- file.path(output_root, "xenium", "vz", "04_qc", "tables", "XenAld_VZ_QC_review_manifest.csv")
output_path <- file.path(rds_dir, "Xenium_VZ_postQC_Res1.5.rds")
provenance_path <- file.path(table_dir, "Xenium_VZ_postQC_manifest.csv")
res_list <- c(0.3, 0.5, 0.8, 1)
plot_paths <- c(
  file.path(plot_path, "XenAld_VZ_PostQC_OrigCluster_UMAP.tif"),
  file.path(plot_path, paste0("XenAld_VZ_PostQC_UMAP_Res_", res_list, ".tif"))
)
expected_outputs <- c(output_path, provenance_path, plot_paths)

required_removal_columns <- c(
  "sample_id", "input_cells", "removed_cells", "remaining_cells", "clusters_removed"
)
required_review_columns <- c(
  "merged_path", "merged_size", "merged_mtime", "cells", "cluster_column", "clusters", "random_seed"
)
removal_manifest <- if (file.exists(removal_path)) {
  read.csv(removal_path, stringsAsFactors = FALSE)
} else NULL
review <- if (file.exists(review_path)) {
  read.csv(review_path, stringsAsFactors = FALSE)
} else NULL
removal_ready <- !is.null(removal_manifest) &&
  all(required_removal_columns %in% names(removal_manifest)) &&
  nrow(removal_manifest) == length(sample_ids) &&
  identical(as.character(removal_manifest$sample_id), sample_ids) &&
  all(is.finite(suppressWarnings(as.numeric(removal_manifest$removed_cells))))
review_ready <- !is.null(review) &&
  all(required_review_columns %in% names(review)) &&
  nrow(review) == 1L &&
  identical(as.character(review$cluster_column), "Xenium_snn_res.0.8")

if (dry_run) {
  compact_dry_run(
    "VZ post-QC reintegration",
    inputs = c(merged_path, review_path, removal_path),
    outputs = expected_outputs,
    checks = c(review_manifest = review_ready, removal_manifest = removal_ready)
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(merged_path)) stop("Merged VZ input not found: ", merged_path)
if (!review_ready) stop("VZ QC review manifest is missing or invalid: ", review_path)
if (!removal_ready) stop("VZ QC removal manifest is missing or invalid: ", removal_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ post-QC outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(patchwork)
})
source(here("scripts", "color_palette.R"))

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(config$runtime$random_seed)

with_future_workers <- function(expr) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  if (workers > 1L) plan(multisession, workers = workers)
  force(expr)
}

check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(merged_path)
required_metadata <- c("orig.ident", "consensus_label", "PCW", "Xenium_snn_res.0.8")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("Merged VZ input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!setequal(unique(as.character(obj$orig.ident)), sample_ids)) {
  stop("Merged VZ input does not contain exactly the configured sample IDs.")
}

merged_info <- file.info(merged_path)
review_is_current <-
  isTRUE(all.equal(as.numeric(review$merged_size), as.numeric(merged_info$size))) &&
  isTRUE(all.equal(as.numeric(review$merged_mtime), as.numeric(merged_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj)))
if (!review_is_current) stop("VZ merged object changed after QC review. Repeat the QC decision stage.")

removal_decisions <- unique(as.character(removal_manifest$clusters_removed))
if (length(removal_decisions) != 1L) stop("VZ removal manifest contains inconsistent decisions.")
clusters_to_remove <- if (identical(removal_decisions, "none")) {
  character()
} else strsplit(removal_decisions, "|", fixed = TRUE)[[1]]

obj <- JoinLayers(obj)
Idents(obj) <- "Xenium_snn_res.0.8"
cells_before_qc <- ncol(obj)
if (length(clusters_to_remove)) {
  unknown_clusters <- setdiff(clusters_to_remove, levels(Idents(obj)))
  if (length(unknown_clusters)) stop("Reviewed VZ cluster(s) are absent: ", paste(unknown_clusters, collapse = ", "))
  obj <- subset(obj, idents = clusters_to_remove, invert = TRUE)
}
removed_cells <- cells_before_qc - ncol(obj)
if (removed_cells != sum(as.numeric(removal_manifest$removed_cells))) {
  stop("Merged VZ removal count does not match the reviewed per-sample removal manifest.")
}

# 1. Re-run PCA on the subset
obj <- RunPCA(
  obj, verbose = FALSE, reduction.name = "pca_clean", npcs = 50,
  seed.use = config$runtime$random_seed
)

# 2. Re-run Harmony (If you used it originally)
# You must re-integrate because the batch-effect vectors change when 31k cells leave
obj <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca_clean", reduction.save = "harmony_clean")

# 3. Re-run UMAP
obj <- RunUMAP(
  obj, reduction = "harmony_clean", dims = 1:50, n.neighbors = 100,
  reduction.name = "umap_clean", seed.use = config$runtime$random_seed
)

obj <- FindNeighbors(obj, reduction = "harmony_clean", dims = 1:50, k.param = 30, verbose = TRUE)

# Define resolutions for the loop
# 8. GENERATE CLUSTERS (RUN IN BULK FOR SPEED)
assay_prefix <- DefaultAssay(obj)
message(Sys.time(), ": Calculating all clusters with ", workers, " configured worker(s)...")

# Scope multisession export of the large object to clustering.
obj <- with_future_workers(
  FindClusters(obj, resolution = res_list, verbose = TRUE)
)
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

saveRDS(obj, output_path, compress = FALSE)
output_info <- file.info(output_path)
write.csv(
  data.frame(
    branch = "VZ",
    input_rds = merged_path,
    review_manifest = review_path,
    removal_manifest = removal_path,
    cells_before_qc = cells_before_qc,
    removed_cells = removed_cells,
    cells_after_qc = ncol(obj),
    cluster_column = "Xenium_snn_res.0.8",
    pca_dimensions = 50L,
    umap_neighbors = 100L,
    neighbor_k = 30L,
    resolutions = paste(res_list, collapse = "|"),
    random_seed = config$runtime$random_seed,
    workers = workers,
    output_rds = output_path,
    output_size = as.numeric(output_info$size),
    output_mtime = as.numeric(output_info$mtime),
    stringsAsFactors = FALSE
  ),
  provenance_path,
  row.names = FALSE
)
plan(sequential)
check_mem("PIPELINE COMPLETE - FILE SAVED")
