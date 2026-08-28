#!/usr/bin/env Rscript

# Build one multi-page consensus-label report after all per-sample consensus
# objects have been generated. Each manifest-defined sample occupies one page
# containing its global spatial map, faceted spatial map, UMAP, and marker
# DotPlot. This script is deliberately separate from the per-sample Slurm array
# so only one process writes to the shared PDF.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)

metadata_path <- resolve_config_path(config$manifests$sample_metadata, config)
if (!file.exists(metadata_path)) stop("Sample metadata not found: ", metadata_path)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "The 'readxl' package is required to order report pages by PCW. ",
    "Restore it through the project renv environment before running this report."
  )
}
sample_metadata <- readxl::read_excel(metadata_path)
required_metadata_columns <- c("sample", "PCW")
if (!all(required_metadata_columns %in% names(sample_metadata))) {
  stop("Sample metadata must contain: ", paste(required_metadata_columns, collapse = ", "))
}

manifest_order <- seq_len(nrow(samples))
base_sample_ids <- sub("_\\d+$", "", samples$sample_id)
metadata_matches <- lapply(base_sample_ids, function(sample_id) {
  which(as.character(sample_metadata$sample) == sample_id)
})
bad_match_counts <- which(lengths(metadata_matches) != 1L)
if (length(bad_match_counts) > 0L) {
  stop(
    "Expected exactly one PCW metadata match for each biological sample; failures:\n- ",
    paste0(
      samples$sample_id[bad_match_counts], " (base ID ",
      base_sample_ids[bad_match_counts], "): ",
      lengths(metadata_matches)[bad_match_counts], " matches",
      collapse = "\n- "
    )
  )
}

pcw_values <- vapply(
  metadata_matches,
  function(i) trimws(as.character(sample_metadata$PCW[[i]])),
  character(1)
)
if (anyNA(pcw_values) || any(!nzchar(pcw_values))) {
  stop(
    "Blank or missing PCW metadata for: ",
    paste(samples$sample_id[is.na(pcw_values) | !nzchar(pcw_values)], collapse = ", ")
  )
}
pcw_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pcw_values)))
if (any(!is.finite(pcw_numeric))) {
  stop(
    "Could not parse numeric PCW values for: ",
    paste(samples$sample_id[!is.finite(pcw_numeric)], collapse = ", ")
  )
}
pcw_display <- format(pcw_numeric, trim = TRUE, scientific = FALSE)

page_order <- order(pcw_numeric, manifest_order)
page_map <- data.frame(
  page = seq_len(nrow(samples)),
  sample_id = samples$sample_id[page_order],
  PCW = pcw_values[page_order],
  PCW_numeric = pcw_numeric[page_order],
  PCW_display = pcw_display[page_order],
  manifest_order = manifest_order[page_order],
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
if (any(!startsWith(args, "--"))) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_03d_plot_report.R ",
    "[--list|--dry-run|--overwrite]"
  )
}
if ("--list" %in% args && length(args) > 1L) {
  stop("--list cannot be combined with another option.")
}

if ("--list" %in% args) {
  write.table(page_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(
  output_root, "xenium", "annotation", "03_consensus_labels", "rds"
)
plot_dir <- file.path(
  output_root, "xenium", "annotation", "03_consensus_labels", "plots"
)
input_paths <- file.path(
  input_dir, paste0(page_map$sample_id, "_Consensus_annotated.rds")
)
report_path <- file.path(plot_dir, "Xenium_Consensus_All_Samples_Report.pdf")

input_status <- data.frame(
  page = page_map$page,
  sample_id = page_map$sample_id,
  PCW = page_map$PCW,
  input = input_paths,
  exists = file.exists(input_paths),
  size_gb = round(file.info(input_paths)$size / 1024^3, 3),
  stringsAsFactors = FALSE
)

discovered_inputs <- if (dir.exists(input_dir)) {
  list.files(
    input_dir,
    pattern = "_Consensus_annotated\\.rds$",
    full.names = TRUE
  )
} else {
  character()
}
unexpected_inputs <- setdiff(basename(discovered_inputs), basename(input_paths))

if (dry_run) {
  cat("Manifest-defined report pages:", nrow(page_map), "\n")
  cat("Page order: numeric PCW, then manifest order within PCW\n")
  cat("Sample metadata:", metadata_path, "\n")
  cat("Complete consensus inputs:", sum(input_status$exists), "of", nrow(page_map), "\n")
  cat("Unexpected consensus inputs:", length(unexpected_inputs), "\n")
  cat("Report output:", report_path, "\n")
  cat("Report output exists:", file.exists(report_path), "\n")
  write.table(input_status, row.names = FALSE, quote = FALSE, sep = "\t")
  if (length(unexpected_inputs) > 0L) {
    cat("Unexpected input files:\n", paste(unexpected_inputs, collapse = "\n"), "\n")
  }
  quit(save = "no", status = 0L)
}

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing ", length(missing_inputs), " manifest-defined consensus object(s):\n- ",
    paste(missing_inputs, collapse = "\n- ")
  )
}
if (length(unexpected_inputs) > 0L) {
  stop(
    "Unexpected consensus object(s) are present in the input directory:\n- ",
    paste(unexpected_inputs, collapse = "\n- ")
  )
}
if (file.exists(report_path) && !overwrite) {
  stop(
    "Refusing to overwrite the existing report: ", report_path,
    "\nUse --overwrite only after reviewing the existing PDF."
  )
}

required_packages <- c("Seurat", "ggplot2", "patchwork", "Cairo", "png", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Restore them through the project renv environment before running this report."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
source(here("scripts", "color_palette.R"))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

partial_path <- file.path(
  plot_dir,
  paste0(".", basename(report_path), ".", Sys.getpid(), ".part")
)

build_page <- function(obj, sample_name, expected_pcw, pcw_title) {
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata) > 0L) {
    stop(
      sample_name, " is missing required metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  object_pcw <- unique(as.character(obj$PCW))
  if (
    length(object_pcw) != 1L || is.na(object_pcw) ||
      !nzchar(trimws(object_pcw)) || object_pcw != expected_pcw
  ) {
    stop(
      sample_name, " must contain exactly one PCW value matching sample metadata ('",
      expected_pcw, "'); found: ", paste(object_pcw, collapse = ", ")
    )
  }
  if (!"umap" %in% Reductions(obj)) {
    stop(sample_name, " lacks the expected 'umap' reduction.")
  }
  if (!"Xenium" %in% Assays(obj)) {
    stop(sample_name, " lacks the expected 'Xenium' assay.")
  }

  consensus_values <- as.character(obj$consensus_label)
  if (anyNA(consensus_values) || any(!nzchar(trimws(consensus_values)))) {
    stop(sample_name, " contains blank or missing consensus labels.")
  }
  consensus_levels <- if (is.factor(obj$consensus_label)) {
    levels(obj$consensus_label)
  } else {
    c(
      intersect(celltype_order, unique(consensus_values)),
      setdiff(unique(consensus_values), celltype_order)
    )
  }
  consensus_levels <- consensus_levels[consensus_levels %in% unique(consensus_values)]
  if (anyDuplicated(consensus_levels)) {
    stop(sample_name, " has duplicate consensus-label levels.")
  }
  validate_palette(consensus_levels)
  obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
  Idents(obj) <- "consensus_label"
  if (!identical(as.character(Idents(obj)), consensus_values)) {
    stop(sample_name, ": active identities do not match consensus_label.")
  }
  plot_colors <- cluster_colors[consensus_levels]

  p_spatial <- ImageDimPlot(
    obj,
    group.by = "consensus_label",
    size = 0.75,
    cols = plot_colors
  ) +
    ggtitle("Global spatial")

  coords <- GetTissueCoordinates(obj)
  required_coordinate_columns <- c("cell", "x", "y")
  missing_coordinate_columns <- setdiff(required_coordinate_columns, colnames(coords))
  if (length(missing_coordinate_columns) > 0L) {
    stop(
      sample_name, " spatial coordinates lack: ",
      paste(missing_coordinate_columns, collapse = ", ")
    )
  }
  coordinate_cell_ids <- as.character(coords$cell)
  if (anyDuplicated(coordinate_cell_ids)) {
    stop(sample_name, " spatial coordinates contain duplicate cell identifiers.")
  }
  coordinate_cells <- match(coordinate_cell_ids, colnames(obj))
  if (anyNA(coordinate_cells)) {
    missing_coordinate_cells <- unique(coordinate_cell_ids[is.na(coordinate_cells)])
    stop(
      sample_name, " spatial coordinates contain cell IDs absent from the object: ",
      paste(head(missing_coordinate_cells, 10L), collapse = ", ")
    )
  }
  plot_data <- cbind(
    coords,
    consensus_label = factor(
      consensus_values[coordinate_cells],
      levels = consensus_levels
    )
  )
  if (anyNA(plot_data$consensus_label)) {
    stop(sample_name, ": failed to join consensus labels to spatial coordinates.")
  }

  p_facet <- ggplot(plot_data, aes(x = y, y = x, color = consensus_label)) +
    geom_point(size = 0.1) +
    facet_wrap(~consensus_label) +
    scale_color_manual(values = plot_colors, drop = FALSE) +
    coord_fixed() +
    ggtitle("Faceted spatial") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white"),
      legend.position = "none",
      strip.text = element_text(color = "white")
    )

  p_umap <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "consensus_label",
    label = TRUE,
    cols = plot_colors
  ) +
    ggtitle("UMAP") +
    NoLegend()

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
    ggtitle("Consensus marker expression")

  ((p_spatial | p_facet) / (p_umap | p_dot)) +
    plot_layout(widths = c(1, 1.35), heights = c(1.15, 1)) +
    plot_annotation(
      title = paste0(sample_name, " (PCW ", pcw_title, ")"),
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    )
}

render_report <- function() {
  page_png <- NULL
  report_complete <- FALSE

  on.exit({
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
    if (!is.null(page_png) && file.exists(page_png)) unlink(page_png)
    if (!report_complete && file.exists(partial_path)) unlink(partial_path)
  }, add = TRUE)

  grDevices::cairo_pdf(
    filename = partial_path,
    width = 20,
    height = 14,
    onefile = TRUE,
    bg = "white"
  )

  for (i in seq_len(nrow(page_map))) {
    sample_name <- page_map$sample_id[[i]]
    expected_pcw <- page_map$PCW[[i]]
    pcw_title <- page_map$PCW_display[[i]]
    message(
      "Rendering page ", i, " of ", nrow(page_map), ": ",
      sample_name, " (PCW ", pcw_title, ")"
    )

    obj <- readRDS(input_paths[[i]])
    page_plot <- build_page(obj, sample_name, expected_pcw, pcw_title)

    # Rasterize the assembled page before embedding it in the PDF. Spatial
    # plots can contain hundreds of thousands of points, so vector pages would
    # be very large and slow to open.
    page_png <- tempfile(
      pattern = paste0("consensus_report_", i, "_"),
      tmpdir = plot_dir,
      fileext = ".png"
    )
    Cairo::CairoPNG(
      filename = page_png,
      width = 3000,
      height = 2100,
      units = "px",
      res = 150,
      bg = "white"
    )
    print(page_plot)
    grDevices::dev.off()

    page_image <- png::readPNG(page_png)
    grid::grid.newpage()
    grid::grid.raster(
      page_image,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    )

    unlink(page_png)
    page_png <- NULL
    rm(obj, page_plot, page_image)
    invisible(gc())
  }

  grDevices::dev.off()

  backup_path <- NULL
  if (file.exists(report_path)) {
    backup_path <- paste0(report_path, ".backup_", Sys.getpid())
    if (!file.rename(report_path, backup_path)) {
      stop("Could not preserve the existing report before overwrite: ", report_path)
    }
  }

  if (!file.rename(partial_path, report_path)) {
    if (!is.null(backup_path) && file.exists(backup_path)) {
      file.rename(backup_path, report_path)
    }
    stop("Could not move the completed report into place: ", report_path)
  }
  if (!is.null(backup_path) && file.exists(backup_path)) {
    unlink(backup_path)
  }

  report_complete <- TRUE
  message(
    "Saved ", nrow(page_map), "-page consensus report: ", report_path
  )
}

render_report()
