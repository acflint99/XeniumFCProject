# Compare whole-tissue clustering at resolution 2.0, 3.0, 4.0, or 5.0 with the
# existing resolution-1.5 result. Pilot modes use three selected samples; the
# --all-samples-res4 mode uses every row in config/samples.csv. Use
# --resolution3, --resolution4, or --resolution5 for an isolated higher-resolution pilot;
# resolution 2.0 remains the default for compatibility with completed jobs.
#
# This script reclusters the saved SNN graph only. It does not repeat QC,
# normalization, PCA, UMAP, or neighbor finding, and it never writes into the
# production 03_clustered directory.

rm(list = ls())

library(here)

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--dry-run", "--overwrite", "--resolution3", "--resolution4",
  "--resolution5", "--all-samples-res4", "--list"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

resolution3_pilot <- "--resolution3" %in% args
resolution4_pilot <- "--resolution4" %in% args
resolution5_pilot <- "--resolution5" %in% args
all_samples_res4 <- "--all-samples-res4" %in% args
resolution_flags <- c(resolution3_pilot, resolution4_pilot, resolution5_pilot)
if (sum(resolution_flags) > 1L) {
  stop("Choose only one of --resolution3, --resolution4, or --resolution5.")
}
if (all_samples_res4 && any(resolution_flags)) {
  stop("--all-samples-res4 cannot be combined with a pilot-resolution option.")
}

baseline_resolution <- 1.5
candidate_resolution <- if (all_samples_res4) {
  4.0
} else if (resolution5_pilot) {
  5.0
} else if (resolution4_pilot) {
  4.0
} else if (resolution3_pilot) {
  3.0
} else {
  2.0
}

# PCW values were verified against metadata/samples_meta.xlsx when the pilot
# manifest was created. Every pilot sample must map to one authoritative
# config/samples.csv row.
pilot_samples <- load_resolution2_pilot_manifest(config)

manifest_counts <- vapply(
  pilot_samples$sample_id,
  function(sample_id) sum(samples$sample_id == sample_id),
  integer(1)
)
if (any(manifest_counts != 1L)) {
  stop(
    "Every pilot sample must map to exactly one config/samples.csv row; failures: ",
    paste(
      paste0(pilot_samples$sample_id[manifest_counts != 1L], "=", manifest_counts[manifest_counts != 1L]),
      collapse = ", "
    )
  )
}
manifest_matches <- match(pilot_samples$sample_id, samples$sample_id)
pilot_samples$input_layout <- samples$input_layout[manifest_matches]

analysis_samples <- if (all_samples_res4) {
  data.frame(
    task_id = seq_len(nrow(samples)),
    sample_id = samples$sample_id,
    PCW = "not_loaded_at_clustering",
    age_group = "not_loaded_at_clustering",
    input_layout = samples$input_layout,
    stringsAsFactors = FALSE
  )
} else {
  pilot_samples
}
if (anyDuplicated(analysis_samples$sample_id) || any(!nzchar(analysis_samples$sample_id))) {
  stop("Selected sample manifest contains duplicate or blank sample IDs.")
}

configured_resolution <- as.numeric(config$initial_clustering$resolution)
if (!isTRUE(all.equal(configured_resolution, baseline_resolution))) {
  stop(
    "This reclustering workflow expects initial_clustering.resolution = ", baseline_resolution,
    " in config/config.yml; found ", configured_resolution, "."
  )
}

format_resolution <- function(x) {
  format(x, scientific = FALSE, trim = TRUE)
}

baseline_tag <- format_resolution(baseline_resolution)
candidate_tag <- sprintf("%.1f", candidate_resolution)
candidate_identity_column <- paste0("whole_tissue_cluster_res", candidate_tag)
facet_point_size <- if (all_samples_res4) {
  0.02
} else if (any(resolution_flags)) {
  0.03
} else {
  0.2
}
assay_name <- as.character(config$initial_clustering$assay)
graph_name <- paste0(assay_name, "_snn")
baseline_column <- paste0(graph_name, "_res.", baseline_tag)
candidate_column <- paste0(graph_name, "_res.", format_resolution(candidate_resolution))

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "preprocess", "03_clustered", "rds")
pilot_stage <- if (all_samples_res4) {
  "03g_resolution4_all_samples"
} else if (resolution5_pilot) {
  "03e_resolution5_pilot"
} else if (resolution4_pilot) {
  "03d_resolution4_pilot"
} else if (resolution3_pilot) {
  "03c_resolution3_pilot"
} else {
  "03b_resolution2_pilot"
}
pilot_root <- file.path(output_root, "xenium", "preprocess", pilot_stage)
plot_dir <- file.path(pilot_root, "plots")
rds_dir <- file.path(pilot_root, "rds")
table_dir <- file.path(pilot_root, "tables")

analysis_samples$input_path <- file.path(
  input_dir,
  paste0(analysis_samples$sample_id, "_CB_QC_cluster.rds")
)

output_paths_for_sample <- function(sample_id) {
  stem <- paste0(
    sample_id, "_whole_tissue_Res", candidate_tag,
    if (all_samples_res4) "" else "_pilot"
  )
  c(
    rds = file.path(rds_dir, paste0(stem, ".rds")),
    umap = file.path(plot_dir, paste0(stem, "_UMAP_Res1.5_vs_Res", candidate_tag, ".tif")),
    spatial = file.path(plot_dir, paste0(stem, "_Spatial_Res1.5_vs_Res", candidate_tag, ".tif")),
    facet = file.path(plot_dir, paste0(stem, "_FacetSpatial_Res", candidate_tag, ".tif")),
    cluster_counts = file.path(table_dir, paste0(stem, "_cluster_counts.csv")),
    transitions = file.path(table_dir, paste0(stem, "_cluster_transitions.csv")),
    provenance = file.path(table_dir, paste0(stem, "_provenance.csv"))
  )
}

task_map <- analysis_samples[, c(
  "task_id", "sample_id", "PCW", "age_group", "input_layout", "input_path"
)]

if ("--list" %in% args) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept TASK_ID.")
  cat("Scope:", if (all_samples_res4) "all manifest samples" else "three-sample pilot", "\n")
  cat("Candidate resolution:", candidate_tag, "\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R ",
    "[--resolution3|--resolution4|--resolution5|--all-samples-res4] ",
    "[--list|--dry-run|--overwrite] [TASK_ID]"
  )
}

task_value <- if (length(task_args) == 1L) {
  task_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ". Use --list to inspect the mapping.")
}

task <- task_map[task_id, , drop = FALSE]
sample_id <- task$sample_id[[1]]
input_path <- task$input_path[[1]]
expected_outputs <- output_paths_for_sample(sample_id)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]

if (dry_run) {
  cat("Task:", task_id, "of", nrow(task_map), "\n")
  cat("Biological sample:", sample_id, "\n")
  cat("Developmental age:", task$PCW[[1]], "(", task$age_group[[1]], ")\n")
  cat("Baseline resolution:", baseline_tag, "\n")
  cat("Candidate resolution:", candidate_tag, "\n")
  cat("Input clustered RDS:", input_path, "\n")
  cat("Input exists:", file.exists(input_path), "\n")
  cat("Dry-run scope: path inspection only; no RDS is loaded.\n")
  write.table(
    data.frame(
      output_type = names(expected_outputs),
      output = unname(expected_outputs),
      exists = file.exists(expected_outputs),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  quit(save = "no", status = if (file.exists(input_path)) 0L else 1L)
}

if (!file.exists(input_path)) stop("Input clustered RDS not found: ", input_path)
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing resolution-", candidate_tag,
    if (all_samples_res4) " all-sample outputs for " else " pilot outputs for ",
    sample_id,
    ":\n- ", paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}
if (length(existing_outputs)) {
  warning("Overwriting ", length(existing_outputs), " reviewed output(s) for ", sample_id, ".")
}

required_packages <- c("Seurat", "ggplot2", "patchwork", "randomcoloR", "Cairo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Restore them explicitly through the project renv environment before running the pilot."
  )
}

for (path in c(plot_dir, rds_dir, table_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

set.seed(as.integer(config$runtime$random_seed))
message("Loading existing whole-tissue clustered object: ", input_path)
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", input_path)

cells_before <- SeuratObject::Cells(obj)
metadata_before <- colnames(obj[[]])
if (!length(cells_before) || anyDuplicated(cells_before)) {
  stop("Input object has no cells or contains duplicate cell IDs.")
}
if (!assay_name %in% SeuratObject::Assays(obj)) {
  stop("Configured assay '", assay_name, "' is absent from ", input_path, ".")
}
if (!graph_name %in% names(obj@graphs)) {
  stop(
    "Expected SNN graph '", graph_name, "' is absent. Available graphs: ",
    paste(names(obj@graphs), collapse = ", ")
  )
}
if (!"umap" %in% names(obj@reductions)) stop("Input object lacks the 'umap' reduction.")
if (!baseline_column %in% metadata_before) {
  stop("Input object lacks the baseline cluster column '", baseline_column, "'.")
}
if (!"seurat_clusters" %in% metadata_before) {
  stop("Input object lacks the 'seurat_clusters' metadata column.")
}
if (candidate_column %in% metadata_before) {
  stop(
    "Input production object already contains candidate column '", candidate_column,
    "'. Refusing to replace it during the pilot."
  )
}

baseline_clusters <- as.character(obj[[baseline_column]][, 1])
current_clusters <- as.character(obj$seurat_clusters)
if (anyNA(baseline_clusters) || !identical(baseline_clusters, current_clusters)) {
  stop(
    "The baseline column '", baseline_column,
    "' does not exactly match seurat_clusters for every cell."
  )
}
obj$whole_tissue_cluster_res1.5 <- factor(baseline_clusters)

SeuratObject::DefaultAssay(obj) <- assay_name
message(
  "Reclustering ", length(cells_before), " cells on ", graph_name,
  " at resolution ", candidate_tag, "."
)
obj <- Seurat::FindClusters(
  obj,
  resolution = candidate_resolution,
  graph.name = graph_name,
  algorithm = as.integer(config$initial_clustering$algorithm),
  n.iter = as.integer(config$initial_clustering$iterations),
  random.seed = as.integer(config$runtime$random_seed),
  verbose = TRUE
)

if (!candidate_column %in% colnames(obj[[]])) {
  stop("FindClusters did not create the expected column '", candidate_column, "'.")
}
candidate_clusters <- as.character(obj[[candidate_column]][, 1])
if (anyNA(candidate_clusters)) {
  stop("Resolution-", candidate_tag, " clustering produced missing cluster IDs.")
}
obj[[candidate_identity_column]] <- factor(candidate_clusters)
SeuratObject::Idents(obj) <- candidate_identity_column

cells_after <- SeuratObject::Cells(obj)
if (!identical(cells_before, cells_after)) {
  stop("Cell IDs or cell order changed during resolution-", candidate_tag, " clustering.")
}
missing_original_metadata <- setdiff(metadata_before, colnames(obj[[]]))
if (length(missing_original_metadata)) {
  stop(
    "Metadata columns were lost during clustering: ",
    paste(missing_original_metadata, collapse = ", ")
  )
}

sort_cluster_ids <- function(cluster_ids) {
  unique_ids <- unique(cluster_ids)
  numeric_ids <- suppressWarnings(as.numeric(unique_ids))
  if (!anyNA(numeric_ids)) unique_ids[order(numeric_ids)] else sort(unique_ids)
}

baseline_levels <- sort_cluster_ids(baseline_clusters)
candidate_levels <- sort_cluster_ids(candidate_clusters)

# Save the scientifically essential object and audit tables before plotting.
# This ensures that a graphics-device or coordinate-rendering error cannot
# discard a successfully completed clustering result.
cluster_counts <- rbind(
  data.frame(
    sample_id = sample_id,
    PCW = task$PCW[[1]],
    resolution = baseline_resolution,
    cluster = baseline_levels,
    cell_count = as.integer(table(factor(baseline_clusters, levels = baseline_levels))),
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = sample_id,
    PCW = task$PCW[[1]],
    resolution = candidate_resolution,
    cluster = candidate_levels,
    cell_count = as.integer(table(factor(candidate_clusters, levels = candidate_levels))),
    stringsAsFactors = FALSE
  )
)

transition_counts <- as.data.frame(
  table(
    baseline_clusters,
    candidate_clusters
  ),
  stringsAsFactors = FALSE,
  responseName = "cell_count"
)
names(transition_counts)[1:2] <- c(
  "cluster_res1.5",
  paste0("cluster_res", candidate_tag)
)
transition_counts <- transition_counts[transition_counts$cell_count > 0L, , drop = FALSE]
transition_counts$sample_id <- sample_id
transition_counts$PCW <- task$PCW[[1]]
transition_counts <- transition_counts[, c(
  "sample_id", "PCW", "cluster_res1.5", paste0("cluster_res", candidate_tag),
  "cell_count"
)]

provenance <- data.frame(
  field = c(
    "sample_id", "PCW", "age_group", "input_path", "output_rds",
    "baseline_resolution", "candidate_resolution", "graph_name", "algorithm",
    "iterations", "random_seed", "cells", "baseline_clusters", "candidate_clusters",
    "clustering_saved_at"
  ),
  value = c(
    sample_id, task$PCW[[1]], task$age_group[[1]], input_path,
    expected_outputs[["rds"]], baseline_tag, candidate_tag, graph_name,
    as.character(config$initial_clustering$algorithm),
    as.character(config$initial_clustering$iterations),
    as.character(config$runtime$random_seed), as.character(length(cells_after)),
    as.character(length(baseline_levels)), as.character(length(candidate_levels)),
    format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(cluster_counts, expected_outputs[["cluster_counts"]], row.names = FALSE)
utils::write.csv(transition_counts, expected_outputs[["transitions"]], row.names = FALSE)
utils::write.csv(provenance, expected_outputs[["provenance"]], row.names = FALSE)
saveRDS(obj, expected_outputs[["rds"]])
message("Saved resolution-", candidate_tag, " pilot RDS and audit tables before plotting.")

set.seed(as.integer(config$runtime$random_seed))
baseline_colors <- randomcoloR::distinctColorPalette(length(baseline_levels))
candidate_colors <- randomcoloR::distinctColorPalette(length(candidate_levels))
names(baseline_colors) <- baseline_levels
names(candidate_colors) <- candidate_levels

umap_baseline <- Seurat::DimPlot(
  obj,
  reduction = "umap",
  group.by = "whole_tissue_cluster_res1.5",
  label = TRUE,
  cols = baseline_colors
) + ggplot2::ggtitle(paste0(sample_id, " - whole tissue, resolution 1.5"))

umap_candidate <- Seurat::DimPlot(
  obj,
  reduction = "umap",
  group.by = candidate_identity_column,
  label = TRUE,
  cols = candidate_colors
) + ggplot2::ggtitle(paste0(sample_id, " - whole tissue, resolution ", candidate_tag))

Cairo::CairoTIFF(
  expected_outputs[["umap"]], width = 20, height = 10, units = "in",
  res = as.integer(config$plotting$tiff_dpi)
)
print(patchwork::wrap_plots(umap_baseline, umap_candidate, ncol = 2))
grDevices::dev.off()

if (!"fov" %in% SeuratObject::Images(obj)) {
  stop("Input object lacks the expected spatial FOV named 'fov'.")
}
spatial_baseline <- Seurat::ImageDimPlot(
  obj,
  fov = "fov",
  group.by = "whole_tissue_cluster_res1.5",
  cols = baseline_colors,
  size = 0.75
) + ggplot2::ggtitle(paste0(sample_id, " - spatial, resolution 1.5"))

spatial_candidate <- Seurat::ImageDimPlot(
  obj,
  fov = "fov",
  group.by = candidate_identity_column,
  cols = candidate_colors,
  size = 0.75
) + ggplot2::ggtitle(paste0(sample_id, " - spatial, resolution ", candidate_tag))

Cairo::CairoTIFF(
  expected_outputs[["spatial"]], width = 20, height = 10, units = "in",
  res = as.integer(config$plotting$tiff_dpi)
)
print(patchwork::wrap_plots(spatial_baseline, spatial_candidate, ncol = 2))
grDevices::dev.off()

coords <- Seurat::GetTissueCoordinates(obj, image = "fov", full = FALSE)
required_coordinate_columns <- c("x", "y", "cell")
missing_coordinate_columns <- setdiff(required_coordinate_columns, colnames(coords))
if (length(missing_coordinate_columns)) {
  stop(
    "Spatial coordinates lack required columns: ",
    paste(missing_coordinate_columns, collapse = ", ")
  )
}
coordinate_cells <- as.character(coords$cell)
if (anyNA(coordinate_cells) || any(!nzchar(coordinate_cells))) {
  stop("Spatial coordinates contain blank or missing cell IDs.")
}
if (anyDuplicated(coordinate_cells)) {
  stop("Spatial centroid coordinates contain duplicate cell IDs.")
}
missing_coordinate_cells <- setdiff(cells_after, coordinate_cells)
unexpected_coordinate_cells <- setdiff(coordinate_cells, cells_after)
if (length(missing_coordinate_cells) || length(unexpected_coordinate_cells)) {
  stop(
    "Spatial centroid cell-ID mismatch. Missing coordinates for ",
    length(missing_coordinate_cells), " object cells; unexpected coordinates for ",
    length(unexpected_coordinate_cells), " cells."
  )
}
coordinate_index <- match(cells_after, coordinate_cells)
if (anyNA(coordinate_index)) stop("Failed to construct the one-to-one spatial coordinate join.")
coords <- coords[coordinate_index, , drop = FALSE]
if (!identical(as.character(coords$cell), cells_after)) {
  stop("Spatial coordinates could not be reordered to match object cell IDs exactly.")
}
facet_data <- data.frame(
  x = coords$x,
  y = coords$y,
  cluster = factor(candidate_clusters, levels = candidate_levels),
  stringsAsFactors = FALSE
)
if (anyNA(facet_data$cluster)) {
  stop("Failed to map resolution-", candidate_tag, " clusters to spatial coordinates.")
}

facet_plot <- ggplot2::ggplot(
  facet_data,
  ggplot2::aes(x = y, y = x, color = cluster)
) +
  ggplot2::geom_point(size = facet_point_size) +
  ggplot2::facet_wrap(~cluster) +
  ggplot2::scale_color_manual(values = candidate_colors) +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "none",
    panel.background = ggplot2::element_rect(fill = "black", color = NA),
    plot.background = ggplot2::element_rect(fill = "black", color = NA),
    strip.text = ggplot2::element_text(color = "white", face = "bold"),
    plot.title = ggplot2::element_text(color = "white", hjust = 0.5)
  ) +
  ggplot2::ggtitle(
    paste0(sample_id, " - faceted spatial clusters, resolution ", candidate_tag)
  )

Cairo::CairoTIFF(
  expected_outputs[["facet"]], width = 14, height = 12, units = "in",
  res = as.integer(config$plotting$tiff_dpi)
)
print(facet_plot)
grDevices::dev.off()

message(
  "Resolution-", candidate_tag, " reclustering complete for ", sample_id, ": ",
  length(baseline_levels), " baseline clusters and ",
  length(candidate_levels), " candidate clusters."
)
message("Output root: ", pilot_root)
