#!/usr/bin/env Rscript

# Attach the cluster-level consensus labels produced by
# xenium_annotate_02_build_consensus.R to one individual Xenium object. The Aldinger object
# is used as the expression/spatial container; its original label columns are
# retained unchanged.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(readxl)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)

task_map <- data.frame(
  task_id = seq_len(nrow(samples)),
  sample_id = samples$sample_id,
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

valid_options <- c("--dry-run", "--overwrite", "--dotplot-only")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
dotplot_only <- "--dotplot-only" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (dotplot_only && overwrite) {
  stop("--dotplot-only already authorizes replacing the DotPlot files; do not combine it with --overwrite.")
}

task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_03_apply_consensus.R ",
    "[--dry-run|--overwrite|--dotplot-only] [TASK_ID]"
  )
}

task_value <- if (length(task_args) == 1L) {
  task_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(samples)) {
  stop("TASK_ID must be between 1 and ", nrow(samples), ". Use --list to inspect the mapping.")
}

sample_name <- samples$sample_id[[task_id]]
output_root <- here(config$project$outputs_dir)
input_path <- file.path(
  output_root,
  "xenium", "annotation", "01_label_transfer", "aldinger", "rds",
  paste0(sample_name, "_Aldinger_annotated.rds")
)
comparison_path <- file.path(
  output_root,
  "xenium", "annotation", "02_consensus", "tables",
  paste0(sample_name, "_comparison_merged.csv")
)
metadata_path <- resolve_config_path(config$manifests$sample_metadata, config)
output_dir <- file.path(
  output_root, "xenium", "annotation", "03_consensus_labels", "rds"
)
plot_dir <- file.path(
  output_root, "xenium", "annotation", "03_consensus_labels", "plots", "samples"
)
output_path <- file.path(output_dir, paste0(sample_name, "_Consensus_annotated.rds"))
plot_paths <- file.path(
  plot_dir,
  paste0(
    sample_name,
    c("_Consensus_UMAP.tif", "_Consensus_GlobalSpatial.tif", "_Consensus_FacetSpatial.tif")
  )
)
dotplot_paths <- file.path(
  plot_dir,
  paste0(sample_name, "_Consensus_Marker_DotPlot", c(".tif", ".pdf"))
)
expected_outputs <- c(output_path, plot_paths, dotplot_paths)

save_consensus_dotplot <- function(obj, sample_name, dotplot_paths) {
  if (!"consensus_label" %in% colnames(obj[[]])) {
    stop("Consensus object lacks the 'consensus_label' metadata column: ", sample_name)
  }

  consensus_values <- as.character(obj$consensus_label)
  if (anyNA(consensus_values) || any(!nzchar(consensus_values))) {
    stop("Consensus object contains blank or missing consensus labels: ", sample_name)
  }
  consensus_levels <- if (is.factor(obj$consensus_label)) {
    levels(obj$consensus_label)
  } else {
    unique(consensus_values)
  }
  consensus_levels <- consensus_levels[consensus_levels %in% unique(consensus_values)]
  obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
  Idents(obj) <- "consensus_label"

  existing_markers <- lapply(markers, function(features) intersect(features, rownames(obj)))
  existing_markers <- existing_markers[lengths(existing_markers) > 0L]
  if (!length(existing_markers)) {
    stop("None of the configured broad-cell markers are present in ", sample_name, ".")
  }

  p_dot <- DotPlot(
    obj,
    features = existing_markers,
    assay = "Xenium",
    cols = c("lightgrey", "red")
  ) +
    scale_y_discrete(limits = rev(consensus_levels)) +
    RotatedAxis() +
    ggtitle(paste(sample_name, "Consensus marker expression"))

  Cairo::CairoTIFF(dotplot_paths[[1]], width = 10, height = 6, units = "in", res = 600)
  print(p_dot)
  grDevices::dev.off()
  ggplot2::ggsave(
    filename = dotplot_paths[[2]],
    plot = p_dot,
    device = grDevices::cairo_pdf,
    width = 10,
    height = 6,
    units = "in"
  )
}

if (dotplot_only) {
  if (dry_run) {
    cat("Task:", task_id, "of", nrow(samples), "\n")
    cat("Biological sample:", sample_name, "\n")
    cat("Consensus input:", output_path, "\n")
    cat("Consensus input exists:", file.exists(output_path), "\n")
    write.table(
      data.frame(output = dotplot_paths, exists = file.exists(dotplot_paths)),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    quit(save = "no", status = 0L)
  }

  if (!file.exists(output_path)) stop("Consensus object not found: ", output_path)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(output_path)
  save_consensus_dotplot(obj, sample_name, dotplot_paths)
  message("Regenerated consensus DotPlots for ", sample_name, ".")
  quit(save = "no", status = 0L)
}

if (!file.exists(metadata_path)) stop("Sample metadata not found: ", metadata_path)
sample_map <- readxl::read_excel(metadata_path)
required_metadata_columns <- c("sample", "PCW")
if (!all(required_metadata_columns %in% names(sample_map))) {
  stop("Sample metadata must contain: ", paste(required_metadata_columns, collapse = ", "))
}
base_sample_name <- sub("_\\d+$", "", sample_name)
metadata_matches <- which(as.character(sample_map$sample) == base_sample_name)
if (length(metadata_matches) != 1L) {
  stop(
    "Expected one PCW metadata match for '", base_sample_name,
    "'; found ", length(metadata_matches), "."
  )
}
pcw_value <- as.character(sample_map$PCW[[metadata_matches]])
if (is.na(pcw_value) || !nzchar(pcw_value)) {
  stop("PCW metadata is blank for '", base_sample_name, "'.")
}

if (dry_run) {
  cat("Task:", task_id, "of", nrow(samples), "\n")
  cat("Biological sample:", sample_name, "\n")
  cat("Annotated input:", input_path, "\n")
  cat("Consensus table:", comparison_path, "\n")
  cat("Sample metadata:", metadata_path, "\n")
  cat("Matched PCW:", pcw_value, "\n")
  write.table(
    data.frame(output = expected_outputs, exists = file.exists(expected_outputs)),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("Annotated input not found: ", input_path)
if (!file.exists(comparison_path)) stop("Consensus table not found: ", comparison_path)

existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Refusing to overwrite existing consensus outputs for ", sample_name,
    ":\n- ", paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing files."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_path)
if (!"seurat_clusters" %in% colnames(obj[[]])) {
  stop("Annotated object lacks the 'seurat_clusters' metadata column: ", input_path)
}

obj$PCW <- pcw_value

comparison <- read.csv(comparison_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("seurat_clusters", "consensus_label")
if (!all(required_columns %in% names(comparison))) {
  stop("Consensus table must contain: ", paste(required_columns, collapse = ", "))
}
if (anyDuplicated(as.character(comparison$seurat_clusters))) {
  stop("Consensus table contains duplicate seurat_clusters values: ", comparison_path)
}

cell_clusters <- as.character(obj$seurat_clusters)
cluster_ids <- as.character(comparison$seurat_clusters)
cluster_lookup <- setNames(trimws(as.character(comparison$consensus_label)), cluster_ids)
unmapped_clusters <- setdiff(unique(cell_clusters), names(cluster_lookup))
if (length(unmapped_clusters) > 0L) {
  stop("No consensus label for cluster(s): ", paste(unmapped_clusters, collapse = ", "))
}

comparison_label_columns <- setdiff(names(comparison), "seurat_clusters")
for (column in comparison_label_columns) {
  column_lookup <- setNames(as.character(comparison[[column]]), cluster_ids)
  obj[[column]] <- unname(column_lookup[cell_clusters])
}

consensus_values <- as.character(obj$consensus_label)
consensus_values[is.na(consensus_values) | !nzchar(consensus_values)] <- "Unknown"
consensus_levels <- c(
  intersect(celltype_order, unique(consensus_values)),
  setdiff(unique(consensus_values), celltype_order)
)
validate_palette(setdiff(consensus_levels, "Unknown"))
obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
Idents(obj) <- "consensus_label"

if (!identical(as.character(Idents(obj)), as.character(obj$consensus_label))) {
  stop("Failed to set consensus_label as the active identity.")
}

plot_colors <- cluster_colors[consensus_levels]

p_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "consensus_label",
  label = TRUE,
  cols = plot_colors
) + ggtitle(paste(sample_name, "Consensus labels"))
Cairo::CairoTIFF(plot_paths[[1]], width = 8, height = 6, units = "in", res = 600)
print(p_umap)
grDevices::dev.off()

p_spatial <- ImageDimPlot(
  obj,
  group.by = "consensus_label",
  size = 0.75,
  cols = plot_colors
) + ggtitle(paste(sample_name, "Consensus labels"))
Cairo::CairoTIFF(
  plot_paths[[2]], width = 8, height = 8, units = "in", res = 600, bg = "black"
)
print(p_spatial)
grDevices::dev.off()

coords <- GetTissueCoordinates(obj)
if (!"cell" %in% colnames(coords)) {
  stop("Spatial coordinates lack the expected 'cell' identifier column.")
}
coordinate_cell_ids <- as.character(coords$cell)
if (anyDuplicated(coordinate_cell_ids)) {
  stop("Spatial coordinates contain duplicate cell identifiers.")
}
coordinate_cells <- match(coordinate_cell_ids, colnames(obj))
if (anyNA(coordinate_cells)) {
  missing_coordinate_cells <- unique(coordinate_cell_ids[is.na(coordinate_cells)])
  stop(
    "Spatial coordinates contain cell identifiers absent from the Seurat object: ",
    paste(head(missing_coordinate_cells, 10L), collapse = ", ")
  )
}
plot_data <- cbind(
  coords,
  consensus_label = factor(
    as.character(obj$consensus_label)[coordinate_cells],
    levels = consensus_levels
  )
)
p_facet <- ggplot(plot_data, aes(x = y, y = x, color = consensus_label)) +
  geom_point(size = 0.1) +
  facet_wrap(~consensus_label) +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  coord_fixed() +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black"),
    legend.position = "none",
    strip.text = element_text(color = "white")
  )
Cairo::CairoTIFF(
  plot_paths[[3]], width = 12, height = 8, units = "in", res = 600, bg = "black"
)
print(p_facet)
grDevices::dev.off()

save_consensus_dotplot(obj, sample_name, dotplot_paths)

saveRDS(obj, output_path, compress = FALSE)
message(
  "Saved consensus-labelled object with PCW metadata and plots for ",
  sample_name, "."
)
