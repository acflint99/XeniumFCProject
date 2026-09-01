#!/usr/bin/env Rscript

# Build matched all-cell and known-only multi-page consensus-label reports after
# all per-sample consensus objects have been generated. Each manifest-defined
# sample occupies one page containing its global spatial map, faceted spatial
# map, UMAP, and marker DotPlot. This script renders pages independently before
# a single merge job writes each shared PDF.

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
valid_options <- c(
  "--list", "--dry-run", "--overwrite", "--all-samples-res4", "--all-samples-res5",
  "--weighted-2of3", "--render-page", "--merge-pages"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
if (any(!startsWith(args, "--"))) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_03d_plot_report.R ",
    "[--all-samples-res4|--all-samples-res5] [--weighted-2of3] ",
    "[--list|--dry-run|--render-page|--merge-pages] [--overwrite]"
  )
}
all_samples_res4 <- "--all-samples-res4" %in% args
all_samples_res5 <- "--all-samples-res5" %in% args
if (all_samples_res4 && all_samples_res5) {
  stop("Choose only one of --all-samples-res4 or --all-samples-res5.")
}
all_samples_mode <- all_samples_res4 || all_samples_res5
weighted_2of3 <- "--weighted-2of3" %in% args
if (xor(all_samples_mode, weighted_2of3)) {
  stop(
    "An all-samples resolution report requires its --all-samples-res* option and ",
    "--weighted-2of3 together; omit both for the original production report."
  )
}
all_samples_weighted <- all_samples_mode && weighted_2of3
report_resolution <- if (all_samples_res5) "5" else if (all_samples_res4) "4" else NA_character_
report_resolution_tag <- if (all_samples_res5) "5.0" else if (all_samples_res4) "4.0" else NA_character_
render_page <- "--render-page" %in% args
merge_pages <- "--merge-pages" %in% args
if (render_page && merge_pages) {
  stop("Choose only one of --render-page or --merge-pages.")
}
if ((render_page || merge_pages) && !all_samples_weighted) {
  stop("Parallel page rendering is supported only for all-sample weighted reports.")
}
report_mode <- if (all_samples_weighted) {
  paste0("resolution-", report_resolution_tag, " all-sample weighted 2-of-3 consensus")
} else {
  "production consensus"
}
facet_point_size <- if (all_samples_weighted) 0.01 else 0.1
facet_point_shape <- if (all_samples_weighted) "." else 19
facet_point_alpha <- if (all_samples_weighted) 0.55 else 1
if (
  "--list" %in% args &&
    any(c("--dry-run", "--overwrite", "--render-page", "--merge-pages") %in% args)
) {
  stop("--list cannot be combined with another report action.")
}

if ("--list" %in% args) {
  cat("Mode:", report_mode, "\n")
  write.table(page_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (dry_run && (render_page || merge_pages)) {
  stop("--dry-run cannot be combined with --render-page or --merge-pages.")
}

output_root <- here(config$project$outputs_dir)
annotation_root <- if (all_samples_weighted) {
  file.path(output_root, "xenium", "annotation", paste0("resolution", report_resolution, "_all_samples"))
} else {
  file.path(output_root, "xenium", "annotation")
}
consensus_stage <- if (all_samples_weighted) {
  "03_consensus_labels_weighted_2of3"
} else {
  "03_consensus_labels"
}
input_dir <- file.path(annotation_root, consensus_stage, "rds")
plot_dir <- file.path(annotation_root, consensus_stage, "plots")
input_paths <- file.path(
  input_dir, paste0(page_map$sample_id, "_Consensus_annotated.rds")
)
report_filename <- if (all_samples_weighted) {
  paste0("Xenium_Consensus_Res", report_resolution, "_Weighted2of3_All_Samples_Report.pdf")
} else {
  "Xenium_Consensus_All_Samples_Report.pdf"
}
known_only_report_filename <- sub(
  "_Report\\.pdf$", "_KnownOnly_Report.pdf", report_filename
)
report_variants <- c("with_unknown", "known_only")
report_paths <- stats::setNames(
  file.path(plot_dir, c(report_filename, known_only_report_filename)),
  report_variants
)
page_dir_base <- if (all_samples_weighted) {
  paste0("report_pages_res", report_resolution, "_weighted2of3")
} else {
  "report_pages"
}
page_dirs <- stats::setNames(
  file.path(plot_dir, c(page_dir_base, paste0(page_dir_base, "_known_only"))),
  report_variants
)
page_paths <- lapply(page_dirs, function(page_dir) {
  file.path(page_dir, sprintf("page_%02d_%s.png", page_map$page, page_map$sample_id))
})

input_status <- data.frame(
  page = page_map$page,
  sample_id = page_map$sample_id,
  PCW = page_map$PCW,
  input = input_paths,
  exists = file.exists(input_paths),
  size_gb = round(file.info(input_paths)$size / 1024^3, 3),
  page_output_with_unknown = page_paths$with_unknown,
  page_exists_with_unknown = file.exists(page_paths$with_unknown),
  page_output_known_only = page_paths$known_only,
  page_exists_known_only = file.exists(page_paths$known_only),
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
  compact_dry_run(
    paste0("Consensus report [", report_mode, "]"),
    inputs = c(metadata_path, input_paths),
    outputs = c(report_paths, unlist(page_paths, use.names = FALSE)),
    checks = c(no_unexpected_inputs = !length(unexpected_inputs))
  )
  if (length(unexpected_inputs) > 0L) {
    cat("  Unexpected inputs: ", paste(head(unexpected_inputs, 3L), collapse = "; "))
    if (length(unexpected_inputs) > 3L) cat("; +", length(unexpected_inputs) - 3L, " more")
    cat("\n")
  }
  inputs_ready <- all(input_status$exists) && !length(unexpected_inputs)
  quit(save = "no", status = if (inputs_ready) 0L else 1L)
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
existing_reports <- report_paths[file.exists(report_paths)]
if (length(existing_reports) && !overwrite) {
  stop(
    "Refusing to overwrite existing report(s):\n- ",
    paste(existing_reports, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing PDF."
  )
}

required_packages <- if (merge_pages) {
  c("Cairo", "png", "readxl")
} else {
  c("Seurat", "ggplot2", "patchwork", "Cairo", "png", "readxl")
}
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
  if (!merge_pages) {
    library(Seurat)
    library(ggplot2)
    library(patchwork)
  }
})
if (!merge_pages) source(here("scripts", "color_palette.R"))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
if (render_page || merge_pages) {
  invisible(lapply(page_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

partial_paths <- stats::setNames(
  file.path(
    plot_dir,
    paste0(".", basename(report_paths), ".", Sys.getpid(), ".part")
  ),
  report_variants
)

build_page <- function(
    obj, sample_name, expected_pcw, pcw_title, include_unknown = TRUE) {
  required_metadata <- c("consensus_label", "PCW")
  if (all_samples_weighted) {
    required_metadata <- c(
      required_metadata, "seurat_clusters",
      paste0("whole_tissue_cluster_res", report_resolution_tag),
      paste0("Xenium_snn_res.", report_resolution), "consensus_method"
    )
  }
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
  if (all_samples_weighted) {
    object_methods <- unique(trimws(as.character(obj$consensus_method)))
    if (!identical(object_methods, "weighted_2of3")) {
      stop(
        sample_name, " does not contain weighted_2of3 consensus provenance; found: ",
        paste(object_methods, collapse = ", ")
      )
    }
    if (!identical(
      as.character(obj$seurat_clusters),
      as.character(obj[[paste0("whole_tissue_cluster_res", report_resolution_tag)]][, 1])
    )) {
      stop(sample_name, ": seurat_clusters does not match resolution-", report_resolution, " identities.")
    }
  }

  consensus_values <- trimws(as.character(obj$consensus_label))
  if (anyNA(consensus_values) || any(!nzchar(consensus_values))) {
    stop(sample_name, " contains blank or missing consensus labels.")
  }
  if (!include_unknown) {
    known_cells <- Cells(obj)[!is_unknown_consensus_label(consensus_values)]
    if (!length(known_cells)) {
      stop(sample_name, " has no known cells for the known-only report.")
    }
    obj <- subset(obj, cells = known_cells)
    consensus_values <- trimws(as.character(obj$consensus_label))
  }
  consensus_levels <- order_consensus_levels(consensus_values, celltype_order)
  if (anyDuplicated(consensus_levels)) {
    stop(sample_name, " has duplicate consensus-label levels.")
  }
  validate_palette(consensus_levels)
  obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
  Idents(obj) <- "consensus_label"
  if (!identical(as.character(Idents(obj)), consensus_values)) {
    stop(sample_name, ": active identities do not match consensus_label.")
  }
  plot_colors <- resolve_cluster_colors(consensus_levels)

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
    geom_point(
      shape = facet_point_shape,
      size = facet_point_size,
      alpha = facet_point_alpha
    ) +
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
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot <- standardize_broad_dotplot(p_dot, consensus_levels) +
    ggtitle("Consensus marker expression")

  ((p_spatial | p_facet) / (p_umap | p_dot)) +
    plot_layout(widths = c(1, 1.35), heights = c(1.15, 1)) +
    plot_annotation(
      title = paste0(
        sample_name, " (PCW ", pcw_title, ")",
        if (all_samples_weighted) paste0(" — resolution ", report_resolution, ", weighted 2-of-3") else "",
        if (include_unknown) " — with Unknown clusters" else " — known cells only"
      ),
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    )
}

install_completed_file <- function(partial_file, destination) {
  backup_path <- NULL
  if (file.exists(destination)) {
    backup_path <- paste0(destination, ".backup_", Sys.getpid())
    if (!file.rename(destination, backup_path)) {
      stop("Could not preserve the existing file before overwrite: ", destination)
    }
  }

  if (!file.rename(partial_file, destination)) {
    if (!is.null(backup_path) && file.exists(backup_path)) {
      file.rename(backup_path, destination)
    }
    stop("Could not move the completed file into place: ", destination)
  }
  if (!is.null(backup_path) && file.exists(backup_path)) unlink(backup_path)
}

render_page_file <- function(obj, task_id, variant) {
  page_path <- page_paths[[variant]][[task_id]]
  if (file.exists(page_path) && !overwrite) {
    stop(
      "Refusing to overwrite the existing rendered page: ", page_path,
      "\nUse --overwrite only after reviewing it."
    )
  }
  sample_name <- page_map$sample_id[[task_id]]
  expected_pcw <- page_map$PCW[[task_id]]
  pcw_title <- page_map$PCW_display[[task_id]]
  include_unknown <- identical(variant, "with_unknown")
  partial_page <- paste0(page_path, ".", Sys.getpid(), ".part.png")
  page_complete <- FALSE
  device_open <- FALSE
  on.exit({
    if (device_open) grDevices::dev.off()
    if (!page_complete && file.exists(partial_page)) unlink(partial_page)
  }, add = TRUE)

  message(
    "Rendering ", gsub("_", " ", variant), " report page ", task_id,
    " of ", nrow(page_map), ": ", sample_name, " (PCW ", pcw_title, ")"
  )
  page_plot <- build_page(
    obj, sample_name, expected_pcw, pcw_title, include_unknown = include_unknown
  )
  Cairo::CairoPNG(
    filename = partial_page, width = 3000, height = 2100,
    units = "px", res = 150, bg = "white"
  )
  device_open <- TRUE
  print(page_plot)
  grDevices::dev.off()
  device_open <- FALSE

  page_info <- file.info(partial_page)
  if (!file.exists(partial_page) || is.na(page_info$size) || page_info$size <= 0) {
    stop("Rendered page is missing or empty: ", partial_page)
  }
  install_completed_file(partial_page, page_path)
  page_complete <- TRUE
  message("Saved rendered report page: ", page_path)
}

render_one_page <- function() {
  task_value <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  task_id <- suppressWarnings(as.integer(task_value))
  if (is.na(task_id) || task_id < 1L || task_id > nrow(page_map)) {
    stop("SLURM_ARRAY_TASK_ID must be between 1 and ", nrow(page_map), ".")
  }
  task_outputs <- vapply(page_paths, `[[`, character(1), task_id)
  existing_task_outputs <- task_outputs[file.exists(task_outputs)]
  if (length(existing_task_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite existing rendered page(s):\n- ",
      paste(existing_task_outputs, collapse = "\n- ")
    )
  }
  obj <- readRDS(input_paths[[task_id]])
  for (variant in names(page_paths)) render_page_file(obj, task_id, variant)
}

validate_page_set <- function(variant) {
  variant_paths <- page_paths[[variant]]
  missing_pages <- variant_paths[!file.exists(variant_paths)]
  empty_pages <- variant_paths[
    file.exists(variant_paths) &
      (is.na(file.info(variant_paths)$size) | file.info(variant_paths)$size <= 0)
  ]
  discovered_pages <- list.files(
    page_dirs[[variant]], pattern = "^page_[0-9]+_.*\\.png$", full.names = TRUE
  )
  unexpected_pages <- setdiff(basename(discovered_pages), basename(variant_paths))
  if (length(missing_pages)) {
    stop("Missing rendered report page(s):\n- ", paste(missing_pages, collapse = "\n- "))
  }
  if (length(empty_pages)) {
    stop("Empty rendered report page(s):\n- ", paste(empty_pages, collapse = "\n- "))
  }
  if (length(unexpected_pages)) {
    stop("Unexpected rendered report page(s):\n- ", paste(unexpected_pages, collapse = "\n- "))
  }
  invisible(TRUE)
}

merge_page_set <- function(variant) {
  variant_paths <- page_paths[[variant]]
  partial_path <- partial_paths[[variant]]
  report_path <- report_paths[[variant]]
  report_complete <- FALSE
  device_open <- FALSE
  on.exit({
    if (device_open) grDevices::dev.off()
    if (!report_complete && file.exists(partial_path)) unlink(partial_path)
  }, add = TRUE)

  grDevices::cairo_pdf(
    filename = partial_path,
    width = 20,
    height = 14,
    onefile = TRUE,
    bg = "white"
  )
  device_open <- TRUE
  for (i in seq_along(variant_paths)) {
    message(
      "Merging ", gsub("_", " ", variant), " report page ", i,
      " of ", length(variant_paths), ": ", basename(variant_paths[[i]])
    )
    page_image <- png::readPNG(variant_paths[[i]])
    grid::grid.newpage()
    grid::grid.raster(
      page_image,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    )
    rm(page_image)
    invisible(gc())
  }
  grDevices::dev.off()
  device_open <- FALSE

  install_completed_file(partial_path, report_path)
  report_complete <- TRUE
  message("Saved ", length(variant_paths), "-page consensus report: ", report_path)
}

merge_rendered_pages <- function() {
  for (variant in names(page_paths)) validate_page_set(variant)
  for (variant in names(page_paths)) merge_page_set(variant)
}

render_report <- function() {
  invisible(lapply(page_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  all_page_paths <- unlist(page_paths, use.names = FALSE)
  existing_pages <- all_page_paths[file.exists(all_page_paths)]
  if (length(existing_pages) && !overwrite) {
    stop(
      "Refusing to overwrite existing rendered page(s):\n- ",
      paste(existing_pages, collapse = "\n- "),
      "\nUse --overwrite only after reviewing them."
    )
  }
  for (i in seq_len(nrow(page_map))) {
    obj <- readRDS(input_paths[[i]])
    for (variant in names(page_paths)) render_page_file(obj, i, variant)
    rm(obj)
    invisible(gc())
  }
  merge_rendered_pages()
}

if (render_page) {
  render_one_page()
} else if (merge_pages) {
  merge_rendered_pages()
} else {
  render_report()
}
