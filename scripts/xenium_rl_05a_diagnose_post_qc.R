#!/usr/bin/env Rscript

# Non-destructive diagnostics for understanding why RL stage-3 and stage-5
# UMAPs and partitions differ after the reviewed stage-4 removal decision.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_vz_rl_qc_routing.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_rl_05a_diagnose_post_qc.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
stage3_path <- file.path(
  output_root, "xenium", "rl", "03_integrated", "rds", "Xenium_RL_integrated.rds"
)
removal_path <- file.path(
  output_root, "xenium", "rl", "04_qc", "tables", "Xenium_RL_QC_removal_manifest.csv"
)
stage4_marker_path <- file.path(
  output_root, "xenium", "rl", "04_qc", "tables",
  "Xenium_RL_RawSubcluster_top10_Markers_Res0.8.csv"
)
stage5_path <- file.path(
  output_root, "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC.rds"
)
stage5_manifest_path <- file.path(
  output_root, "xenium", "rl", "05_post_qc", "tables", "Xenium_RL_postQC_manifest.csv"
)
diagnostic_root <- file.path(output_root, "xenium", "rl", "05_post_qc", "diagnostics")
plot_dir <- file.path(diagnostic_root, "plots")
table_dir <- file.path(diagnostic_root, "tables")

focus_clusters <- c("6", "15")
dimensions <- 30L
umap_neighbors <- 30L
neighbor_k_values <- c(30L, 50L)
resolutions <- c(0.3, 0.5, 0.8)
cluster_column <- "Xenium_snn_res.0.8"
random_seed <- as.integer(config$runtime$random_seed)

retained_plot_name <- "RL_Diagnostic_01_Stage3_RetainedCells_OriginalUMAP_Res0.8.tif"
highlight_plot_names <- paste0(
  "RL_Diagnostic_02_Stage3_RemovedCluster_", focus_clusters, "_OriginalUMAP.tif"
)
graph_plot_names <- unlist(lapply(neighbor_k_values, function(k_value) {
  paste0(
    "RL_Diagnostic_04_Dims30_k", k_value, "_Res", resolutions, ".tif"
  )
}), use.names = FALSE)
table_names <- c(
  "RL_Diagnostic_03_RemovedCluster_Overview.csv",
  "RL_Diagnostic_03_PooledLabel_Composition.csv",
  "RL_Diagnostic_03_Sample_Composition.csv",
  "RL_Diagnostic_03_QC_Summary.csv",
  "RL_Diagnostic_03_Stage4_Top10_Markers.csv",
  "RL_Diagnostic_04_Cluster_Sizes_And_Centroids.csv",
  "RL_Diagnostic_04_k30_vs_k50_Contingency.csv",
  "RL_Diagnostic_Manifest.csv"
)
expected_outputs <- c(
  file.path(plot_dir, c(retained_plot_name, highlight_plot_names, graph_plot_names)),
  file.path(table_dir, table_names)
)
required_inputs <- c(
  stage3_path, removal_path, stage4_marker_path, stage5_path, stage5_manifest_path
)

removal_manifest <- if (file.exists(removal_path)) {
  read.csv(removal_path, stringsAsFactors = FALSE)
} else NULL
removal_inspection <- inspect_regional_qc_manifest(removal_manifest, sample_ids)

if (dry_run) {
  compact_dry_run(
    "RL stage-3 versus stage-5 diagnostic",
    inputs = required_inputs,
    outputs = expected_outputs,
    checks = c(
      sample_manifest_has_34_samples = length(sample_ids) == 34L,
      removal_manifest = isTRUE(removal_inspection$valid),
      nonzero_removal = isTRUE(removal_inspection$valid) &&
        !isTRUE(removal_inspection$zero_removal)
    )
  )
  quit(save = "no", status = 0L)
}

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing RL diagnostic input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
if (length(sample_ids) != 34L || anyDuplicated(sample_ids) || any(!nzchar(sample_ids))) {
  stop("config/samples.csv must contain exactly 34 unique, nonblank sample IDs.")
}
if (!isTRUE(removal_inspection$valid) || isTRUE(removal_inspection$zero_removal)) {
  stop(
    "RL removal manifest must be a valid, nonzero 34-sample decision:\n- ",
    paste(removal_inspection$errors, collapse = "\n- ")
  )
}
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL diagnostic outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(future)
})

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(random_seed)

check_mem <- function(step_label) {
  memory <- gc(full = TRUE)
  message("\n[", Sys.time(), "] --- ", step_label, " ---")
  message("Memory in use: ", round(sum(memory[, 2]), 1), " MB\n")
}

numeric_cluster_factor <- function(values) {
  values <- as.character(values)
  observed <- unique(values[!is.na(values)])
  numeric_values <- suppressWarnings(as.numeric(observed))
  levels <- if (length(observed) && all(is.finite(numeric_values))) {
    observed[order(numeric_values)]
  } else {
    sort(observed)
  }
  factor(values, levels = levels)
}

make_lightweight_seurat <- function(metadata, reductions) {
  cell_ids <- rownames(metadata)
  if (is.null(cell_ids) || anyDuplicated(cell_ids) || any(!nzchar(cell_ids))) {
    stop("Diagnostic metadata must have unique, nonblank cell IDs as row names.")
  }
  counts <- Matrix::sparseMatrix(
    i = rep(seq_len(2L), each = length(cell_ids)),
    j = rep.int(seq_along(cell_ids), times = 2L),
    x = rep.int(1, 2L * length(cell_ids)),
    dims = c(2L, length(cell_ids)),
    dimnames = list(c("diagnosticA", "diagnosticB"), cell_ids)
  )
  result <- CreateSeuratObject(
    counts = counts, assay = "Diagnostic", meta.data = metadata,
    min.cells = 0, min.features = 0
  )
  for (reduction_name in names(reductions)) {
    embeddings <- reductions[[reduction_name]]
    if (!identical(rownames(embeddings), cell_ids)) {
      stop("Reduction cell order differs from diagnostic metadata: ", reduction_name)
    }
    key <- paste0(toupper(gsub("[^A-Za-z0-9]", "", reduction_name)), "_")
    result[[reduction_name]] <- CreateDimReducObject(
      embeddings = embeddings, key = key, assay = "Diagnostic"
    )
  }
  result
}

save_tiff <- function(plot, path, width = 10, height = 8) {
  Cairo::CairoTIFF(
    filename = path, width = width, height = height, units = "in", res = 600
  )
  print(plot)
  grDevices::dev.off()
}

with_fork_workers <- function(expr, max_workers = workers) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  active_workers <- min(as.integer(max_workers), workers)
  if (active_workers > 1L && identical(.Platform$OS.type, "unix")) {
    plan(multicore, workers = active_workers)
  } else {
    plan(sequential)
  }
  force(expr)
}

staging_root <- file.path(
  diagnostic_root,
  paste0(".staging_", Sys.getenv("SLURM_JOB_ID", unset = "local"), "_", Sys.getpid())
)
staging_plot_dir <- file.path(staging_root, "plots")
staging_table_dir <- file.path(staging_root, "tables")
if (dir.exists(staging_root)) stop("Diagnostic staging directory already exists: ", staging_root)
dir.create(staging_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_table_dir, recursive = TRUE, showWarnings = FALSE)
completed <- FALSE
on.exit({
  if (!completed && dir.exists(staging_root)) {
    message("Diagnostic failed; partial files remain only in staging: ", staging_root)
  }
}, add = TRUE)

removal_decisions <- unique(as.character(removal_manifest$clusters_removed))
if (length(removal_decisions) != 1L || identical(removal_decisions, "none")) {
  stop("RL removal manifest does not contain one consistent nonzero decision.")
}
clusters_to_remove <- strsplit(removal_decisions, "|", fixed = TRUE)[[1]]
if (!all(focus_clusters %in% clusters_to_remove)) {
  stop(
    "The reviewed removal decision no longer includes focus cluster(s): ",
    paste(setdiff(focus_clusters, clusters_to_remove), collapse = ", ")
  )
}

message("Loading stage-3 RL object for diagnostics 1-3: ", stage3_path)
stage3 <- readRDS(stage3_path)
required_stage3_metadata <- c(
  "orig.ident", "consensus_label", "regional_consensus_label_pooled",
  "nCount_Xenium", "nFeature_Xenium", cluster_column
)
missing_stage3_metadata <- setdiff(required_stage3_metadata, colnames(stage3[[]]))
if (length(missing_stage3_metadata)) {
  stop("Stage-3 RL object lacks metadata: ", paste(missing_stage3_metadata, collapse = ", "))
}
if (!"umap_harmony" %in% Reductions(stage3)) {
  stop("Stage-3 RL object lacks the umap_harmony reduction.")
}
if (!setequal(unique(as.character(stage3$orig.ident)), sample_ids)) {
  stop("Stage-3 RL object does not contain exactly the configured sample IDs.")
}

stage3_clusters <- as.character(stage3[[cluster_column, drop = TRUE]])
if (anyNA(stage3_clusters) || any(!nzchar(stage3_clusters))) {
  stop("Stage-3 RL cluster assignments contain missing or blank values.")
}
unknown_removed_clusters <- setdiff(clusters_to_remove, unique(stage3_clusters))
if (length(unknown_removed_clusters)) {
  stop("Reviewed removal cluster(s) are absent from stage 3: ", paste(unknown_removed_clusters, collapse = ", "))
}
removed_index <- stage3_clusters %in% clusters_to_remove
retained_cell_ids <- colnames(stage3)[!removed_index]
removed_cell_ids <- colnames(stage3)[removed_index]
recorded_removed_cells <- sum(as.numeric(removal_manifest$removed_cells))
if (length(removed_cell_ids) != recorded_removed_cells) {
  stop(
    "Stage-3 removed-cell count (", length(removed_cell_ids),
    ") differs from the removal manifest (", recorded_removed_cells, ")."
  )
}

stage3_metadata <- stage3[[]][, required_stage3_metadata, drop = FALSE]
stage3_metadata[[cluster_column]] <- numeric_cluster_factor(stage3_metadata[[cluster_column]])
stage3_umap <- Embeddings(stage3, reduction = "umap_harmony")
if (!identical(rownames(stage3_umap), rownames(stage3_metadata))) {
  stop("Stage-3 UMAP cell order differs from metadata cell order.")
}
stage3_light <- make_lightweight_seurat(
  metadata = stage3_metadata,
  reductions = list(stage3_umap = stage3_umap)
)
rm(stage3, stage3_metadata, stage3_umap)
gc(full = TRUE)
check_mem("STAGE-3 LIGHTWEIGHT DIAGNOSTIC OBJECT READY")

retained_plot <- DimPlot(
  stage3_light,
  reduction = "stage3_umap",
  cells = retained_cell_ids,
  group.by = cluster_column,
  label = TRUE,
  label.size = 5,
  label.box = TRUE,
  raster = TRUE,
  raster.dpi = c(600, 600),
  pt.size = 0.45,
  alpha = 0.8
) +
  ggtitle("RL retained cells on the unchanged Stage-3 UMAP") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff(retained_plot, file.path(staging_plot_dir, retained_plot_name))
rm(retained_plot)

highlight_colors <- c("6" = "#D73027", "15" = "#4575B4")
for (focus_cluster in focus_clusters) {
  focus_cell_ids <- colnames(stage3_light)[
    as.character(stage3_light[[cluster_column, drop = TRUE]]) == focus_cluster
  ]
  if (!length(focus_cell_ids)) stop("Focus cluster has no cells: ", focus_cluster)
  highlight_plot <- DimPlot(
    stage3_light,
    reduction = "stage3_umap",
    cells.highlight = focus_cell_ids,
    cols = "grey85",
    cols.highlight = unname(highlight_colors[[focus_cluster]]),
    sizes.highlight = 0.5,
    raster = TRUE,
    raster.dpi = c(600, 600),
    pt.size = 0.25
  ) +
    ggtitle(paste0("Stage-3 UMAP: reviewed removal cluster ", focus_cluster)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  save_tiff(
    highlight_plot,
    file.path(
      staging_plot_dir,
      paste0("RL_Diagnostic_02_Stage3_RemovedCluster_", focus_cluster, "_OriginalUMAP.tif")
    )
  )
  rm(highlight_plot)
}

diagnostic_metadata <- stage3_light[[]] %>%
  mutate(
    cell_id = rownames(stage3_light[[]]),
    expression_cluster = as.character(.data[[cluster_column]])
  )

removed_overview <- diagnostic_metadata %>%
  filter(.data$expression_cluster %in% clusters_to_remove) %>%
  group_by(.data$expression_cluster) %>%
  summarise(
    cell_count = n(),
    reviewed_focus_cluster = first(.data$expression_cluster) %in% focus_clusters,
    median_nCount_Xenium = median(.data$nCount_Xenium),
    median_nFeature_Xenium = median(.data$nFeature_Xenium),
    samples_represented = n_distinct(.data$orig.ident),
    .groups = "drop"
  ) %>%
  mutate(fraction_of_stage3 = .data$cell_count / ncol(stage3_light)) %>%
  arrange(as.numeric(.data$expression_cluster))
write.csv(
  removed_overview,
  file.path(staging_table_dir, "RL_Diagnostic_03_RemovedCluster_Overview.csv"),
  row.names = FALSE
)

pooled_label_composition <- diagnostic_metadata %>%
  filter(.data$expression_cluster %in% focus_clusters) %>%
  count(
    .data$expression_cluster, .data$regional_consensus_label_pooled,
    name = "cell_count"
  ) %>%
  group_by(.data$expression_cluster) %>%
  mutate(fraction_of_cluster = .data$cell_count / sum(.data$cell_count)) %>%
  ungroup() %>%
  arrange(as.numeric(.data$expression_cluster), desc(.data$cell_count))
write.csv(
  pooled_label_composition,
  file.path(staging_table_dir, "RL_Diagnostic_03_PooledLabel_Composition.csv"),
  row.names = FALSE
)

sample_composition <- diagnostic_metadata %>%
  filter(.data$expression_cluster %in% focus_clusters) %>%
  count(.data$expression_cluster, sample_id = .data$orig.ident, name = "cell_count") %>%
  group_by(.data$expression_cluster) %>%
  mutate(fraction_of_cluster = .data$cell_count / sum(.data$cell_count)) %>%
  ungroup() %>%
  arrange(as.numeric(.data$expression_cluster), desc(.data$cell_count))
write.csv(
  sample_composition,
  file.path(staging_table_dir, "RL_Diagnostic_03_Sample_Composition.csv"),
  row.names = FALSE
)

qc_summary <- diagnostic_metadata %>%
  filter(.data$expression_cluster %in% focus_clusters) %>%
  group_by(.data$expression_cluster) %>%
  summarise(
    cell_count = n(),
    q25_nCount_Xenium = as.numeric(quantile(.data$nCount_Xenium, 0.25)),
    median_nCount_Xenium = median(.data$nCount_Xenium),
    mean_nCount_Xenium = mean(.data$nCount_Xenium),
    q75_nCount_Xenium = as.numeric(quantile(.data$nCount_Xenium, 0.75)),
    q25_nFeature_Xenium = as.numeric(quantile(.data$nFeature_Xenium, 0.25)),
    median_nFeature_Xenium = median(.data$nFeature_Xenium),
    mean_nFeature_Xenium = mean(.data$nFeature_Xenium),
    q75_nFeature_Xenium = as.numeric(quantile(.data$nFeature_Xenium, 0.75)),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(.data$expression_cluster))
write.csv(
  qc_summary,
  file.path(staging_table_dir, "RL_Diagnostic_03_QC_Summary.csv"),
  row.names = FALSE
)

stage4_markers <- read.csv(stage4_marker_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"cluster" %in% names(stage4_markers)) {
  stop("Stage-4 marker table lacks the required cluster column.")
}
focus_markers <- stage4_markers[
  as.character(stage4_markers$cluster) %in% focus_clusters, , drop = FALSE
]
if (!all(focus_clusters %in% unique(as.character(focus_markers$cluster)))) {
  warning("One or more focus clusters have no rows in the Stage-4 top-10 marker table.")
}
write.csv(
  focus_markers,
  file.path(staging_table_dir, "RL_Diagnostic_03_Stage4_Top10_Markers.csv"),
  row.names = FALSE
)

rm(stage3_light, diagnostic_metadata, stage4_markers, focus_markers)
gc(full = TRUE)
check_mem("DIAGNOSTICS 1-3 COMPLETE")

message("Loading stage-5 RL object for diagnostic 4: ", stage5_path)
stage5 <- readRDS(stage5_path)
stage5_manifest <- read.csv(stage5_manifest_path, stringsAsFactors = FALSE)
required_stage5_manifest_columns <- c(
  "cells_after_qc", "output_rds", "output_size", "output_mtime"
)
missing_manifest_columns <- setdiff(required_stage5_manifest_columns, names(stage5_manifest))
if (length(missing_manifest_columns) || nrow(stage5_manifest) != 1L) {
  stop(
    "Stage-5 manifest must contain exactly one row and columns: ",
    paste(required_stage5_manifest_columns, collapse = ", ")
  )
}
stage5_info <- file.info(stage5_path)
stage5_manifest_is_current <-
  isTRUE(all.equal(as.numeric(stage5_manifest$output_size), as.numeric(stage5_info$size))) &&
  isTRUE(all.equal(as.numeric(stage5_manifest$output_mtime), as.numeric(stage5_info$mtime))) &&
  identical(as.integer(stage5_manifest$cells_after_qc), as.integer(ncol(stage5)))
if (!stage5_manifest_is_current) {
  stop("Stage-5 RL object no longer matches its provenance manifest.")
}
if (!setequal(colnames(stage5), retained_cell_ids)) {
  stop("Stage-5 cell IDs are not exactly the retained Stage-3 cell IDs.")
}
if (!identical(colnames(stage5), retained_cell_ids)) {
  stop("Stage-5 cell order differs from the retained Stage-3 cell order.")
}
if (!"harmony_clean" %in% Reductions(stage5)) {
  stop("Stage-5 RL object lacks the harmony_clean reduction.")
}
stage5_harmony <- Embeddings(stage5, reduction = "harmony_clean")
if (ncol(stage5_harmony) < dimensions) {
  stop("Stage-5 harmony_clean has fewer than ", dimensions, " dimensions.")
}
stage5_harmony <- stage5_harmony[, seq_len(dimensions), drop = FALSE]
stage5_metadata <- data.frame(
  orig.ident = as.character(stage5$orig.ident),
  row.names = colnames(stage5),
  stringsAsFactors = FALSE
)
stage5_light <- make_lightweight_seurat(
  metadata = stage5_metadata,
  reductions = list(harmony_diag = stage5_harmony)
)
rm(stage5, stage5_harmony, stage5_metadata)
gc(full = TRUE)
check_mem("STAGE-5 LIGHTWEIGHT DIAGNOSTIC OBJECT READY")

message(
  "Calculating the controlled diagnostic UMAP with dimensions 1:", dimensions,
  " and n.neighbors=", umap_neighbors, "."
)
stage5_light <- with_fork_workers(
  RunUMAP(
    stage5_light,
    reduction = "harmony_diag",
    dims = seq_len(dimensions),
    n.neighbors = umap_neighbors,
    reduction.name = "umap_diag_d30",
    seed.use = random_seed,
    verbose = TRUE
  )
)

diagnostic_assignments <- data.frame(
  cell_id = colnames(stage5_light), stringsAsFactors = FALSE
)
for (k_value in neighbor_k_values) {
  nn_name <- paste0("RL_diag_nn_k", k_value)
  snn_name <- paste0("RL_diag_snn_k", k_value)
  message(
    Sys.time(), ": Building controlled dimensions-1:", dimensions,
    " graph with k.param=", k_value, "."
  )
  stage5_light <- FindNeighbors(
    stage5_light,
    reduction = "harmony_diag",
    dims = seq_len(dimensions),
    k.param = k_value,
    graph.name = c(nn_name, snn_name),
    verbose = TRUE
  )
  if (!all(c(nn_name, snn_name) %in% Graphs(stage5_light))) {
    stop("FindNeighbors did not create the expected k=", k_value, " graphs.")
  }

  stage5_light <- with_fork_workers(
    FindClusters(
      stage5_light,
      graph.name = snn_name,
      resolution = resolutions,
      algorithm = 1,
      n.start = 10,
      n.iter = 10,
      random.seed = random_seed,
      verbose = TRUE
    ),
    max_workers = length(resolutions)
  )
  source_columns <- paste0(snn_name, "_res.", resolutions)
  missing_source_columns <- setdiff(source_columns, colnames(stage5_light[[]]))
  if (length(missing_source_columns)) {
    stop(
      "Controlled clustering lacks expected result column(s): ",
      paste(missing_source_columns, collapse = ", ")
    )
  }

  for (i in seq_along(resolutions)) {
    resolution <- resolutions[[i]]
    result_column <- paste0("RL_diag_k", k_value, "_res.", resolution)
    cluster_values <- numeric_cluster_factor(
      stage5_light[[source_columns[[i]], drop = TRUE]]
    )
    stage5_light[[result_column]] <- cluster_values
    diagnostic_assignments[[result_column]] <- as.character(cluster_values)

    graph_plot <- DimPlot(
      stage5_light,
      reduction = "umap_diag_d30",
      group.by = result_column,
      label = TRUE,
      label.size = 5,
      label.box = TRUE,
      raster = TRUE,
      raster.dpi = c(600, 600),
      pt.size = 0.45,
      alpha = 0.8
    ) +
      ggtitle(
        paste0(
          "RL controlled comparison: 30 dimensions, k=", k_value,
          ", resolution ", resolution
        )
      ) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_tiff(
      graph_plot,
      file.path(
        staging_plot_dir,
        paste0(
          "RL_Diagnostic_04_Dims30_k", k_value, "_Res", resolution, ".tif"
        )
      )
    )
    rm(graph_plot)
  }

  stage5_light[[nn_name]] <- NULL
  stage5_light[[snn_name]] <- NULL
  gc(full = TRUE)
  check_mem(paste0("CONTROLLED k=", k_value, " COMPARISON COMPLETE"))
}

umap_coordinates <- Embeddings(stage5_light, reduction = "umap_diag_d30")
if (!identical(rownames(umap_coordinates), diagnostic_assignments$cell_id)) {
  stop("Diagnostic UMAP order differs from the assignment table.")
}
assignment_long <- bind_rows(lapply(neighbor_k_values, function(k_value) {
  bind_rows(lapply(resolutions, function(resolution) {
    column <- paste0("RL_diag_k", k_value, "_res.", resolution)
    data.frame(
      cell_id = diagnostic_assignments$cell_id,
      k = k_value,
      resolution = resolution,
      cluster = diagnostic_assignments[[column]],
      UMAP_1 = umap_coordinates[, 1],
      UMAP_2 = umap_coordinates[, 2],
      stringsAsFactors = FALSE
    )
  }))
}))
cluster_sizes <- assignment_long %>%
  group_by(.data$k, .data$resolution, .data$cluster) %>%
  summarise(
    cell_count = n(),
    median_UMAP_1 = median(.data$UMAP_1),
    median_UMAP_2 = median(.data$UMAP_2),
    .groups = "drop"
  ) %>%
  mutate(fraction_of_cells = .data$cell_count / ncol(stage5_light)) %>%
  arrange(.data$resolution, .data$k, as.numeric(.data$cluster))
write.csv(
  cluster_sizes,
  file.path(staging_table_dir, "RL_Diagnostic_04_Cluster_Sizes_And_Centroids.csv"),
  row.names = FALSE
)

contingency <- bind_rows(lapply(resolutions, function(resolution) {
  k30_column <- paste0("RL_diag_k30_res.", resolution)
  k50_column <- paste0("RL_diag_k50_res.", resolution)
  result <- as.data.frame(
    table(
      k30_cluster = diagnostic_assignments[[k30_column]],
      k50_cluster = diagnostic_assignments[[k50_column]]
    ),
    stringsAsFactors = FALSE
  )
  result$resolution <- resolution
  result
})) %>%
  select(.data$resolution, everything()) %>%
  filter(.data$Freq > 0) %>%
  arrange(.data$resolution, as.numeric(.data$k30_cluster), as.numeric(.data$k50_cluster))
write.csv(
  contingency,
  file.path(staging_table_dir, "RL_Diagnostic_04_k30_vs_k50_Contingency.csv"),
  row.names = FALSE
)

stage3_info <- file.info(stage3_path)
manifest <- data.frame(
  diagnostic = "RL stage-3 versus stage-5 post-QC comparison",
  stage3_rds = normalizePath(stage3_path),
  stage3_size = as.numeric(stage3_info$size),
  stage3_mtime = as.numeric(stage3_info$mtime),
  stage5_rds = normalizePath(stage5_path),
  stage5_size = as.numeric(stage5_info$size),
  stage5_mtime = as.numeric(stage5_info$mtime),
  removal_manifest = normalizePath(removal_path),
  clusters_removed = removal_decisions,
  focus_clusters = paste(focus_clusters, collapse = "|"),
  stage3_cells = length(retained_cell_ids) + length(removed_cell_ids),
  removed_cells = length(removed_cell_ids),
  retained_cells = length(retained_cell_ids),
  harmony_source = "stage5_harmony_clean",
  dimensions = dimensions,
  umap_neighbors = umap_neighbors,
  neighbor_k_values = paste(neighbor_k_values, collapse = "|"),
  resolutions = paste(resolutions, collapse = "|"),
  clustering_algorithm = 1L,
  clustering_n_start = 10L,
  clustering_n_iter = 10L,
  random_seed = random_seed,
  workers = workers,
  production_rds_written = FALSE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(staging_table_dir, "RL_Diagnostic_Manifest.csv"),
  row.names = FALSE
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
staged_outputs <- c(
  file.path(staging_plot_dir, c(retained_plot_name, highlight_plot_names, graph_plot_names)),
  file.path(staging_table_dir, table_names)
)
missing_staged_outputs <- staged_outputs[!file.exists(staged_outputs)]
if (length(missing_staged_outputs)) {
  stop("Diagnostic did not create expected staged output(s):\n- ", paste(missing_staged_outputs, collapse = "\n- "))
}
copy_ok <- file.copy(
  staged_outputs,
  expected_outputs,
  overwrite = overwrite,
  copy.mode = TRUE,
  copy.date = TRUE
)
if (!all(copy_ok)) {
  stop(
    "Failed to publish diagnostic output(s):\n- ",
    paste(expected_outputs[!copy_ok], collapse = "\n- ")
  )
}
unlink(staging_root, recursive = TRUE)
completed <- TRUE
plan(sequential)
check_mem("RL DIAGNOSTIC COMPLETE")
message(
  "Wrote ", length(expected_outputs), " protected diagnostic outputs under ",
  diagnostic_root, ". Production RDS files were not modified."
)
