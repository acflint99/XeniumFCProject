#!/usr/bin/env Rscript

# Assemble the saved resolution-5 Aldinger label-transfer TIFFs into one
# manifest-ordered PDF page per biological sample. This does not load Seurat
# objects or recompute any analysis or plots.

suppressPackageStartupMessages(library(here))

args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("--list", "--dry-run", "--overwrite")
unknown_args <- setdiff(args, allowed_args)
if (length(unknown_args)) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))
}

list_only <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (list_only && (dry_run || overwrite)) {
  stop("Use --list by itself; do not combine it with --dry-run or --overwrite.")
}

manifest_path <- here("config", "samples.csv")
if (!file.exists(manifest_path)) {
  stop("Missing authoritative sample manifest: ", manifest_path)
}
samples <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% names(samples)) {
  stop("Sample manifest lacks the required 'sample_id' column: ", manifest_path)
}
sample_ids <- trimws(as.character(samples$sample_id))
if (length(sample_ids) != 34L) {
  stop("Expected 34 biological samples; found ", length(sample_ids), " in ", manifest_path)
}
if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
  stop("Sample manifest contains missing, blank, or duplicate sample IDs.")
}

plot_dir <- here(
  "outputs", "xenium", "annotation",
  "01_label_transfer", "aldinger", "plots"
)
output_path <- file.path(
  plot_dir, "Xenium_Aldinger_Res5_All_Samples_Plot_Report.pdf"
)

plot_specs <- data.frame(
  key = c(
    "weighted_umap", "weighted_global", "weighted_facet",
    "majority_umap", "majority_global", "majority_facet",
    "weighted_dotplot", "majority_dotplot", "prediction_histogram"
  ),
  title = c(
    "Cluster-weighted UMAP", "Cluster-weighted global spatial",
    "Cluster-weighted faceted spatial", "Cluster-majority UMAP",
    "Cluster-majority global spatial", "Cluster-majority faceted spatial",
    "Cluster-weighted marker DotPlot", "Cluster-majority marker DotPlot",
    "Prediction-score histogram"
  ),
  suffix = c(
    "_Aldinger_Broad_ClusterWeighted_UMAP.tif",
    "_Broad_GlobalSpatial_ClusterWeighted.tif",
    "_Broad_FacetSpatial_ClusterWeighted.tif",
    "_Aldinger_Broad_ClusterMajority_UMAP.tif",
    "_Broad_GlobalSpatial_ClusterMajority.tif",
    "_Broad_FacetSpatial_ClusterMajority.tif",
    "_Broad_Marker_DotPlot_Weighted.tif",
    "_Broad_Marker_DotPlot_Majority.tif",
    "_Broad_PredictionScores_Hist.tif"
  ),
  stringsAsFactors = FALSE
)

page_map <- do.call(
  rbind,
  lapply(seq_along(sample_ids), function(page_index) {
    data.frame(
      page = page_index,
      sample_id = sample_ids[[page_index]],
      key = plot_specs$key,
      title = plot_specs$title,
      path = file.path(plot_dir, paste0(sample_ids[[page_index]], plot_specs$suffix)),
      stringsAsFactors = FALSE
    )
  })
)

if (list_only) {
  write.table(
    unique(page_map[c("page", "sample_id")]),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = 0L)
}

if (!dir.exists(plot_dir)) {
  stop("Missing Aldinger resolution-5 plot directory: ", plot_dir)
}
missing_inputs <- page_map$path[!file.exists(page_map$path)]
empty_inputs <- page_map$path[
  file.exists(page_map$path) & (is.na(file.info(page_map$path)$size) | file.info(page_map$path)$size <= 0)
]
discovered_tiffs <- list.files(
  plot_dir, pattern = "\\.(tif|tiff)$", full.names = TRUE, ignore.case = TRUE
)
unexpected_tiffs <- setdiff(
  normalizePath(discovered_tiffs, winslash = "/", mustWork = FALSE),
  normalizePath(page_map$path, winslash = "/", mustWork = FALSE)
)
if (length(missing_inputs)) {
  stop("Missing expected Aldinger TIFF(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
if (length(empty_inputs)) {
  stop("Empty expected Aldinger TIFF(s):\n- ", paste(empty_inputs, collapse = "\n- "))
}
if (length(unexpected_tiffs)) {
  stop("Unexpected TIFF(s) in the Aldinger plot directory:\n- ", paste(unexpected_tiffs, collapse = "\n- "))
}

if (dry_run) {
  cat("DRY-RUN PASS: Aldinger resolution-5 plot report\n")
  cat("  manifest: ", manifest_path, "\n", sep = "")
  cat("  samples/pages: ", length(sample_ids), "\n", sep = "")
  cat("  unique TIFF inputs: ", nrow(page_map), "\n", sep = "")
  cat("  output: ", output_path, "\n", sep = "")
  cat("  output exists: ", file.exists(output_path), "\n", sep = "")
  quit(save = "no", status = 0L)
}

required_packages <- c("Cairo", "tiff")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them explicitly in the project renv environment before submission; ",
    "the report script will not install packages."
  )
}
if (file.exists(output_path) && !overwrite) {
  stop(
    "Refusing to overwrite the existing report: ", output_path,
    "\nReview it first, then rerun with --overwrite if replacement is intentional."
  )
}

downsample_native_raster <- function(image, max_dimension = 1800L) {
  image_dimensions <- dim(image)
  if (length(image_dimensions) != 2L || any(image_dimensions < 1L)) {
    stop("TIFF reader returned an image with invalid dimensions.")
  }
  stride <- max(1L, ceiling(max(image_dimensions) / max_dimension))
  if (stride > 1L) {
    image <- image[
      seq.int(1L, image_dimensions[[1]], by = stride),
      seq.int(1L, image_dimensions[[2]], by = stride),
      drop = FALSE
    ]
    class(image) <- "nativeRaster"
  }
  image
}

draw_report_page <- function(sample_id, sample_plots) {
  grid::grid.newpage()
  page_layout <- grid::grid.layout(
    nrow = 2L, ncol = 1L,
    heights = grid::unit(c(0.055, 0.945), "npc")
  )
  grid::pushViewport(grid::viewport(layout = page_layout))
  grid::grid.text(
    paste0(sample_id, " — Aldinger label transfer, whole-tissue resolution 5"),
    vp = grid::viewport(layout.pos.row = 1L),
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2L,
      layout = grid::grid.layout(3L, 3L)
    )
  )

  cell_aspect <- (20 / 3) / ((14 * 0.945) / 3)
  for (plot_index in seq_len(nrow(sample_plots))) {
    row_index <- ((plot_index - 1L) %/% 3L) + 1L
    column_index <- ((plot_index - 1L) %% 3L) + 1L
    grid::pushViewport(
      grid::viewport(layout.pos.row = row_index, layout.pos.col = column_index)
    )
    grid::grid.text(
      sample_plots$title[[plot_index]],
      x = 0.5, y = 0.975,
      gp = grid::gpar(fontsize = 10, fontface = "bold")
    )

    image <- tiff::readTIFF(sample_plots$path[[plot_index]], native = TRUE)
    image <- downsample_native_raster(image)
    image_aspect <- ncol(image) / nrow(image)
    width_npc <- 0.96
    height_npc <- width_npc * cell_aspect / image_aspect
    if (height_npc > 0.88) {
      height_npc <- 0.88
      width_npc <- height_npc * image_aspect / cell_aspect
    }
    grid::grid.raster(
      image,
      x = 0.5, y = 0.455,
      width = grid::unit(width_npc, "npc"),
      height = grid::unit(height_npc, "npc"),
      interpolate = TRUE
    )
    rm(image)
    grid::popViewport()
  }
  grid::popViewport(2L)
  invisible(gc(verbose = FALSE))
}

partial_path <- paste0(output_path, ".", Sys.getpid(), ".part.pdf")
report_complete <- FALSE
device_open <- FALSE
on.exit({
  if (device_open) grDevices::dev.off()
  if (!report_complete && file.exists(partial_path)) unlink(partial_path)
}, add = TRUE)

Cairo::CairoPDF(
  file = partial_path, width = 20, height = 14,
  onefile = TRUE, bg = "white"
)
device_open <- TRUE
for (page_index in seq_along(sample_ids)) {
  message(
    "Rendering Aldinger resolution-5 report page ", page_index, "/",
    length(sample_ids), ": ", sample_ids[[page_index]]
  )
  draw_report_page(
    sample_ids[[page_index]],
    page_map[page_map$page == page_index, , drop = FALSE]
  )
}
grDevices::dev.off()
device_open <- FALSE

report_info <- file.info(partial_path)
if (!file.exists(partial_path) || is.na(report_info$size) || report_info$size <= 0) {
  stop("Rendered report is missing or empty: ", partial_path)
}

backup_path <- NULL
if (file.exists(output_path)) {
  backup_path <- paste0(output_path, ".backup_", Sys.getpid())
  if (!file.rename(output_path, backup_path)) {
    stop("Could not preserve the existing report before overwrite: ", output_path)
  }
}
if (!file.rename(partial_path, output_path)) {
  if (!is.null(backup_path) && file.exists(backup_path)) {
    file.rename(backup_path, output_path)
  }
  stop("Could not move the completed report into place: ", output_path)
}
if (!is.null(backup_path) && file.exists(backup_path)) unlink(backup_path)
report_complete <- TRUE

message("Saved 34-page Aldinger resolution-5 plot report: ", output_path)
