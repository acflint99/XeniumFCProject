#!/usr/bin/env Rscript

# Reuse or reintegrate the reviewed RL object after the explicit stage-4 QC decision.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_vz_rl_qc_routing.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_rl_05_reintegrate_post_qc.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "rl", "05_post_qc", "plots")
rds_dir <- file.path(output_root, "xenium", "rl", "05_post_qc", "rds")
table_dir <- file.path(output_root, "xenium", "rl", "05_post_qc", "tables")
merged_path <- file.path(output_root, "xenium", "rl", "03_integrated", "rds", "Xenium_RL_integrated.rds")
removal_path <- file.path(output_root, "xenium", "rl", "04_qc", "tables", "Xenium_RL_QC_removal_manifest.csv")
review_path <- file.path(output_root, "xenium", "rl", "04_qc", "tables", "Xenium_RL_QC_review_manifest.csv")
output_path <- file.path(rds_dir, "Xenium_RL_postQC.rds")
provenance_path <- file.path(table_dir, "Xenium_RL_postQC_manifest.csv")
full_res_list <- c(0.3, 0.5, 0.8)
reuse_res_list <- c(0.5, 0.8)
post_qc_dimensions <- suppressWarnings(as.integer(
  config$regional_subsets$post_qc_dimensions_observed
))
post_qc_umap_neighbors <- suppressWarnings(as.integer(
  config$regional_subsets$umap_neighbors_post_qc
))
post_qc_neighbor_k <- suppressWarnings(as.integer(
  config$regional_subsets$neighbor_k_post_qc
))
configured_post_qc_parameters <- c(
  dimensions = post_qc_dimensions,
  umap_neighbors = post_qc_umap_neighbors,
  neighbor_k = post_qc_neighbor_k
)
if (length(configured_post_qc_parameters) != 3L ||
    anyNA(configured_post_qc_parameters) ||
    any(!is.finite(configured_post_qc_parameters)) ||
    any(configured_post_qc_parameters < 1L)) {
  stop(
    "RL post-QC dimensions, UMAP neighbors, and graph k must be positive integers."
  )
}

required_review_columns <- c(
  "merged_path", "merged_size", "merged_mtime", "cells", "cluster_column", "clusters", "random_seed"
)
removal_manifest <- if (file.exists(removal_path)) {
  read.csv(removal_path, stringsAsFactors = FALSE)
} else NULL
review <- if (file.exists(review_path)) {
  read.csv(review_path, stringsAsFactors = FALSE)
} else NULL
removal_inspection <- inspect_regional_qc_manifest(removal_manifest, sample_ids)
removal_ready <- isTRUE(removal_inspection$valid)
review_ready <- !is.null(review) &&
  all(required_review_columns %in% names(review)) &&
  nrow(review) == 1L &&
  identical(as.character(review$cluster_column), "Xenium_snn_res.0.8")
zero_removal_fast_path <- removal_ready && isTRUE(removal_inspection$zero_removal)
res_list <- if (zero_removal_fast_path) reuse_res_list else full_res_list
plot_paths <- c(
  file.path(plot_path, "Xenium_RL_PostQC_OrigCluster_UMAP.tif"),
  file.path(plot_path, paste0("Xenium_RL_PostQC_UMAP_Res_", res_list, ".tif"))
)
expected_outputs <- c(output_path, provenance_path, plot_paths)
excluded_fast_path_outputs <- if (zero_removal_fast_path) {
  file.path(
    plot_path,
    paste0("Xenium_RL_PostQC_UMAP_Res_", setdiff(full_res_list, reuse_res_list), ".tif")
  )
} else {
  character()
}

if (dry_run) {
  compact_dry_run(
    if (zero_removal_fast_path) {
      "RL post-QC zero-removal reuse"
    } else {
      "RL post-QC reintegration"
    },
    inputs = c(merged_path, review_path, removal_path),
    outputs = expected_outputs,
    checks = c(
      review_manifest = review_ready,
      removal_manifest = removal_ready,
      no_excluded_full_outputs = !length(
        excluded_fast_path_outputs[file.exists(excluded_fast_path_outputs)]
      )
    )
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(merged_path)) stop("Merged RL input not found: ", merged_path)
if (!review_ready) stop("RL QC review manifest is missing or invalid: ", review_path)
if (!removal_ready) {
  stop(
    "RL QC removal manifest is missing or invalid: ", removal_path, "\n- ",
    paste(removal_inspection$errors, collapse = "\n- ")
  )
}
stale_excluded_outputs <- excluded_fast_path_outputs[file.exists(excluded_fast_path_outputs)]
if (length(stale_excluded_outputs)) {
  stop(
    "Refusing zero-removal reuse while excluded full-reintegration output(s) exist:\n- ",
    paste(stale_excluded_outputs, collapse = "\n- "),
    "\nArchive or remove them explicitly before running the fast path."
  )
}
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL post-QC outputs:\n- ",
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
required_metadata <- c(
  "orig.ident", "consensus_label", "regional_consensus_label",
  "regional_consensus_label_pooled", "regional_unknown_label",
  "PCW", "Xenium_snn_res.0.8"
)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("Merged RL input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!setequal(unique(as.character(obj$orig.ident)), sample_ids)) {
  stop("Merged RL input does not contain exactly the configured sample IDs.")
}
if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
  stop("Merged RL input has invalid sample-qualified Unknown labels.")
}

merged_info <- file.info(merged_path)
review_is_current <-
  isTRUE(all.equal(as.numeric(review$merged_size), as.numeric(merged_info$size))) &&
  isTRUE(all.equal(as.numeric(review$merged_mtime), as.numeric(merged_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj)))
if (!review_is_current) stop("RL merged object changed after QC review. Repeat the QC decision stage.")

removal_decisions <- unique(as.character(removal_manifest$clusters_removed))
if (length(removal_decisions) != 1L) stop("RL removal manifest contains inconsistent decisions.")
clusters_to_remove <- if (identical(removal_decisions, "none")) {
  character()
} else strsplit(removal_decisions, "|", fixed = TRUE)[[1]]

Idents(obj) <- "Xenium_snn_res.0.8"
cells_before_qc <- ncol(obj)
if (length(clusters_to_remove)) {
  unknown_clusters <- setdiff(clusters_to_remove, levels(Idents(obj)))
  if (length(unknown_clusters)) stop("Reviewed RL cluster(s) are absent: ", paste(unknown_clusters, collapse = ", "))
  obj <- subset(obj, idents = clusters_to_remove, invert = TRUE)
}
removed_cells <- cells_before_qc - ncol(obj)
if (removed_cells != sum(as.numeric(removal_manifest$removed_cells))) {
  stop("Merged RL removal count does not match the reviewed per-sample removal manifest.")
}

assay_prefix <- DefaultAssay(obj)
processing_mode <- if (zero_removal_fast_path) {
  "zero_removal_reuse_stage3"
} else {
  "full_post_removal_reintegration"
}
reused_reductions <- character()
reused_resolutions <- numeric()
computed_resolutions <- numeric()
pca_dimensions <- NA_integer_
umap_neighbors <- NA_integer_
neighbor_k <- NA_integer_

if (zero_removal_fast_path) {
  message(
    "Stage-4 manifest records zero removals across all ", length(sample_ids),
    " samples; reusing the reviewed stage-3 reductions, graph, and clusters."
  )
  required_stage3_reductions <- c("pca", "harmony", "umap_harmony")
  missing_stage3_reductions <- setdiff(required_stage3_reductions, Reductions(obj))
  if (length(missing_stage3_reductions)) {
    stop(
      "Zero-removal reuse requires stage-3 reduction(s): ",
      paste(missing_stage3_reductions, collapse = ", ")
    )
  }
  preserved_resolutions <- reuse_res_list
  preserved_columns <- paste0(assay_prefix, "_snn_res.", preserved_resolutions)
  missing_preserved_columns <- setdiff(preserved_columns, colnames(obj[[]]))
  if (length(missing_preserved_columns)) {
    stop(
      "Zero-removal reuse requires stage-3 cluster column(s): ",
      paste(missing_preserved_columns, collapse = ", ")
    )
  }
  graph_name <- paste0(assay_prefix, "_snn")
  if (!graph_name %in% Graphs(obj)) {
    stop("Zero-removal reuse requires the stage-3 SNN graph: ", graph_name)
  }
  preserved_cell_ids <- colnames(obj)
  preserved_cluster_values <- lapply(
    preserved_columns,
    function(column) as.character(obj[[column, drop = TRUE]])
  )
  names(preserved_cluster_values) <- preserved_columns

  stage3_umap <- Embeddings(obj, reduction = "umap_harmony")
  if (!identical(rownames(stage3_umap), colnames(obj))) {
    stop("Stage-3 UMAP cell order does not exactly match the RL object.")
  }
  obj[["umap_clean"]] <- CreateDimReducObject(
    embeddings = stage3_umap,
    key = "UMAPCLEAN_",
    assay = assay_prefix
  )
  reused_reductions <- required_stage3_reductions
  reused_resolutions <- preserved_resolutions
  pca_dimensions <- ncol(Embeddings(obj, reduction = "pca"))

  # Marker testing in stage 6 requires joined assay layers, but joining does
  # not change cell IDs, expression values, reductions, graphs, or clusters.
  obj <- JoinLayers(obj)
  if (!identical(colnames(obj), preserved_cell_ids)) {
    stop("Joining layers changed RL cell IDs or their order during zero-removal reuse.")
  }
  changed_cluster_columns <- preserved_columns[
    !vapply(
      preserved_columns,
      function(column) {
        identical(
          as.character(obj[[column, drop = TRUE]]),
          preserved_cluster_values[[column]]
        )
      },
      logical(1)
    )
  ]
  if (length(changed_cluster_columns)) {
    stop(
      "Joining layers changed preserved stage-3 cluster assignments: ",
      paste(changed_cluster_columns, collapse = ", ")
    )
  }
} else {
  obj <- JoinLayers(obj)

  # Recompute the latent space and graph because removing cells changes the
  # batch correction and neighbor relationships.
  obj <- RunPCA(
    obj,
    verbose = FALSE,
    reduction.name = "pca_clean",
    npcs = post_qc_dimensions,
    seed.use = config$runtime$random_seed
  )
  obj <- RunHarmony(
    obj,
    group.by.vars = "orig.ident",
    reduction = "pca_clean",
    dims.use = seq_len(post_qc_dimensions),
    reduction.save = "harmony_clean"
  )
  obj <- RunUMAP(
    obj,
    reduction = "harmony_clean",
    dims = seq_len(post_qc_dimensions),
    n.neighbors = post_qc_umap_neighbors,
    reduction.name = "umap_clean",
    seed.use = config$runtime$random_seed
  )
  obj <- FindNeighbors(
    obj,
    reduction = "harmony_clean",
    dims = seq_len(post_qc_dimensions),
    k.param = post_qc_neighbor_k,
    verbose = TRUE
  )

  message(Sys.time(), ": Calculating all clusters with ", workers, " configured worker(s)...")
  obj <- with_future_workers(
    FindClusters(obj, resolution = res_list, verbose = TRUE)
  )
  computed_resolutions <- full_res_list
  pca_dimensions <- post_qc_dimensions
  umap_neighbors <- post_qc_umap_neighbors
  neighbor_k <- post_qc_neighbor_k
  check_mem("POST-BATCH-CLUSTERING")
}

required_post_qc_columns <- paste0(assay_prefix, "_snn_res.", res_list)
missing_post_qc_columns <- setdiff(required_post_qc_columns, colnames(obj[[]]))
if (length(missing_post_qc_columns)) {
  stop(
    "RL post-QC object lacks required cluster column(s): ",
    paste(missing_post_qc_columns, collapse = ", ")
  )
}
if (!"umap_clean" %in% Reductions(obj)) {
  stop("RL post-QC object lacks the required umap_clean reduction.")
}

pooled_regional_levels <- order_regional_consensus_levels(
  obj$regional_consensus_label_pooled, celltype_order
)
obj$regional_consensus_label_pooled <- factor(
  as.character(obj$regional_consensus_label_pooled),
  levels = pooled_regional_levels
)
pooled_regional_colors <- resolve_regional_consensus_colors(pooled_regional_levels)

p1 <- DimPlot(obj,
              reduction = "umap_clean",
              group.by = "regional_consensus_label_pooled",
              cols = pooled_regional_colors,
              label = TRUE,
              label.size = 5,
              label.box = TRUE,
              raster = TRUE,
              pt.size = 0.6,
              alpha = 0.8) +
  ggtitle("Refined RL UMAP - Pooled input labels") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_PostQC_OrigCluster_UMAP.tif")),
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
    ggtitle(paste0("Refined RL UMAP - Resolution ", res)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  Cairo::CairoTIFF(
    filename = file.path(plot_path, paste0("Xenium_RL_PostQC_UMAP_Res_", res, ".tif")),
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
    branch = "RL",
    input_rds = merged_path,
    review_manifest = review_path,
    removal_manifest = removal_path,
    cells_before_qc = cells_before_qc,
    removed_cells = removed_cells,
    cells_after_qc = ncol(obj),
    cluster_column = "Xenium_snn_res.0.8",
    processing_mode = processing_mode,
    pca_dimensions = pca_dimensions,
    umap_neighbors = umap_neighbors,
    neighbor_k = neighbor_k,
    resolutions = paste(res_list, collapse = "|"),
    preferred_stage6_resolution = 0.5,
    reused_reductions = if (length(reused_reductions)) {
      paste(reused_reductions, collapse = "|")
    } else {
      "none"
    },
    reused_resolutions = if (length(reused_resolutions)) {
      paste(reused_resolutions, collapse = "|")
    } else {
      "none"
    },
    computed_resolutions = if (length(computed_resolutions)) {
      paste(computed_resolutions, collapse = "|")
    } else {
      "none"
    },
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
