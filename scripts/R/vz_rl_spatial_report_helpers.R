# Build manifest-ordered reports from existing VZ/RL spatial TIFFs.

vzrl_report_specs <- function(output_root, page_map, labels = "both", point_size = "0.6") {
  definitions <- list(
    combined = list(stages = c("01_combined_labels", "01_combined_labels"),
                    plot_dirs = c("plots", "plots/combined"), column = "comb_subcluster",
                    title = "Combined subclusters", stem = "CombinedSubclusters",
                    suffix = "_combSubcluster_GlobalSpatial_pt"),
    `corrected-broad` = list(stages = c("02_corrected_broad", "01_corrected_broad"),
                           plot_dirs = c("plots", "plots"), column = "corrected_broad_label",
                           title = "Corrected broad labels", stem = "CorrectedBroad",
                           suffix = "_CorrectedBroad_GlobalSpatial_pt")
  )
  selected <- if (labels == "both") names(definitions) else labels
  if (!all(selected %in% names(definitions)) || !point_size %in% c("0.6", "0.85")) {
    stop("Choose labels=both|combined|corrected-broad and point_size=0.6|0.85.")
  }
  lapply(selected, function(label) {
    definition <- definitions[[label]]
    stages <- file.path(output_root, "xenium", "vz_rl", definition$stages)
    candidates <- file.path(stages, definition$plot_dirs)
    filenames <- paste0(page_map$sample_id, definition$suffix, point_size, ".tif")
    candidate_inputs <- lapply(candidates, function(directory) file.path(directory, filenames))
    found <- vapply(candidate_inputs, function(paths) sum(file.exists(paths)), integer(1))
    populated <- which(found > 0L)
    # Never silently combine files from different runs/layouts or choose between
    # duplicated plot sets. A partial set stays partial and fails the preflight.
    if (length(populated) > 1L) {
      stop(definition$title, " TIFFs occur in multiple layout locations; resolve duplicates before reporting:\n- ",
           paste0(candidates[populated], " (", found[populated], "/", nrow(page_map),
                  " inputs)", collapse = "\n- "))
    }
    chosen <- if (length(populated)) populated[[1]] else 1L
    stage <- stages[[chosen]]
    stem <- paste0("Xenium_VZRL_AllSamples_", definition$stem,
                   "_Spatial_Report_pt", point_size)
    list(
      label = label, column = definition$column, title = definition$title,
      point_size = point_size, input_dir = candidates[[chosen]],
      searched_dirs = candidates, inputs = candidate_inputs[[chosen]],
      pdf = file.path(stage, "reports", paste0(stem, ".pdf")),
      manifest = file.path(stage, "reports", paste0(stem, "_manifest.csv"))
    )
  })
}

write_vzrl_spatial_report <- function(
    spec, page_map, overwrite = FALSE,
    read_image = function(path) tiff::readTIFF(path, native = TRUE)) {
  if (!length(spec$inputs) || length(spec$inputs) != nrow(page_map) ||
      anyDuplicated(page_map$sample_id) ||
      !identical(page_map$page, seq_len(nrow(page_map))) ||
      any(!is.finite(page_map$PCW_numeric)) || is.unsorted(page_map$PCW_numeric)) {
    stop("Report inputs must match the unique, PCW-ordered sample page map.")
  }
  info <- file.info(spec$inputs)
  if (anyNA(info$size) || any(info$size <= 0)) stop("A report input is missing or empty.")
  outputs <- c(spec$pdf, spec$manifest)
  if (any(file.exists(outputs)) && !overwrite) {
    stop("Report outputs already exist; review before using --overwrite: ",
         paste(outputs[file.exists(outputs)], collapse = "; "))
  }
  if (!nzchar(Sys.which("pdfinfo"))) stop("Poppler pdfinfo is required to validate report page counts.")
  dir.create(dirname(spec$pdf), recursive = TRUE, showWarnings = FALSE)
  temporary_pdf <- tempfile(".vzrl_report_", dirname(spec$pdf), fileext = ".pdf")
  temporary_csv <- tempfile(".vzrl_report_", dirname(spec$manifest), fileext = ".csv")
  device_open <- FALSE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(c(temporary_pdf, temporary_csv))
  }, add = TRUE)
  grDevices::pdf(temporary_pdf, width = 14, height = 10.5, bg = "black",
                 onefile = TRUE, useDingbats = FALSE, compress = TRUE)
  device_open <- TRUE
  for (i in seq_len(nrow(page_map))) {
    message(spec$title, ": page ", i, "/", nrow(page_map), " [", page_map$sample_id[[i]], "]")
    image <- read_image(spec$inputs[[i]])
    dimensions <- dim(image)
    if (length(dimensions) != 2L || any(dimensions <= 0L)) {
      stop("Expected a single native raster image from: ", spec$inputs[[i]])
    }
    grid::grid.newpage()
    grid::grid.text(
      paste0(page_map$sample_id[[i]], " | PCW ", page_map$PCW_numeric[[i]], " | ",
             spec$title, " | ", i, "/", nrow(page_map)),
      x = 0.5, y = 0.975, gp = grid::gpar(col = "white", fontsize = 14)
    )
    # Fit the complete TIFF, including its legend, without cropping or resampling.
    image_aspect <- dimensions[[2]] / dimensions[[1]]
    width <- min(14, 9.9 * image_aspect)
    height <- width / image_aspect
    grid::grid.raster(image, x = 0.5, y = 0.48,
                      width = grid::unit(width, "in"), height = grid::unit(height, "in"),
                      interpolate = FALSE)
    rm(image)
    invisible(gc(verbose = FALSE))
  }
  grDevices::dev.off()
  device_open <- FALSE
  if (pdf_page_count(temporary_pdf) != nrow(page_map)) stop("Generated PDF page count is incorrect.")
  final_info <- file.info(spec$inputs)
  if (!identical(as.numeric(info$size), as.numeric(final_info$size)) ||
      !identical(as.numeric(info$mtime), as.numeric(final_info$mtime))) {
    stop("Input TIFFs changed while generating the report; wait for plotting jobs to finish.")
  }
  provenance <- data.frame(
    page_map, label_column = spec$column, point_size = spec$point_size,
    input_tiff = normalizePath(spec$inputs, mustWork = TRUE),
    input_bytes = as.numeric(info$size), input_mtime = as.numeric(info$mtime),
    output_pdf = normalizePath(spec$pdf, mustWork = FALSE),
    combined_page_count = nrow(page_map), stringsAsFactors = FALSE
  )
  write.csv(provenance, temporary_csv, row.names = FALSE)
  if (!file.rename(temporary_pdf, spec$pdf)) stop("Could not place report PDF: ", spec$pdf)
  if (!file.rename(temporary_csv, spec$manifest)) stop("Could not place report manifest: ", spec$manifest)
  message("Created ", spec$pdf, " and ", spec$manifest)
  invisible(provenance)
}
