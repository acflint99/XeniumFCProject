#!/usr/bin/env Rscript

# Compile the 34 production Unknown-QC audits into one evidence table and one
# blank review template, then merge the sample QC PDFs in manifest order.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) {
  stop("Usage: Rscript scripts/xenium_unknown_02_compile_qc_review.R [--dry-run|--overwrite]")
}
if (all(c("--dry-run", "--overwrite") %in% args)) {
  stop("--dry-run and --overwrite cannot be combined.")
}

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)
plot_dir <- file.path(
  production_paths$annotation_root, "04_unknown_qc_review", "plots"
)
sample_plot_paths <- file.path(
  plot_dir, paste0(sample_ids, "_Unknown_vs_Known_QC.pdf")
)
combined_plot_path <- file.path(plot_dir, "Xenium_Unknown_QC_All_Samples.pdf")
pdfunite <- Sys.which("pdfunite")
pdfinfo <- Sys.which("pdfinfo")

if ("--dry-run" %in% args) {
  compact_dry_run(
    "Production Unknown-QC PDF merge",
    inputs = sample_plot_paths,
    outputs = combined_plot_path,
    checks = c(
      pdfunite_available = nzchar(pdfunite),
      pdfinfo_available = nzchar(pdfinfo)
    )
  )
}

if (!"--dry-run" %in% args) {
  missing_plots <- sample_plot_paths[!file.exists(sample_plot_paths)]
  if (length(missing_plots)) {
    stop(
      "Cannot merge Unknown-QC plots; missing files:\n- ",
      paste(missing_plots, collapse = "\n- ")
    )
  }
  empty_plots <- sample_plot_paths[file.info(sample_plot_paths)$size <= 0]
  if (length(empty_plots)) {
    stop(
      "Cannot merge Unknown-QC plots; empty files:\n- ",
      paste(empty_plots, collapse = "\n- ")
    )
  }
  if (!nzchar(pdfunite) || !nzchar(pdfinfo)) {
    stop(
      "Merging Unknown-QC plots requires the Poppler commands pdfunite and ",
      "pdfinfo. Load the Cheaha module that provides Poppler before rerunning."
    )
  }
  if (file.exists(combined_plot_path) && !"--overwrite" %in% args) {
    stop(
      "Refusing to overwrite combined Unknown-QC plot: ", combined_plot_path,
      "\nUse --overwrite only after review."
    )
  }
}

driver <- here("scripts", "xenium_annotate_04_audit_unknown_qc.R")
driver_args <- c(
  shQuote(driver), "--all-samples-res5", "--weighted-2of3", "--combine", args
)
status <- system2(file.path(R.home("bin"), "Rscript"), driver_args)
if (!identical(status, 0L)) {
  stop("Unknown-QC compilation driver exited with status ", status, call. = FALSE)
}

if ("--dry-run" %in% args) quit(save = "no", status = 0L)

pdf_page_count <- function(path) {
  output <- suppressWarnings(system2(
    pdfinfo, shQuote(path), stdout = TRUE, stderr = TRUE
  ))
  command_status <- attr(output, "status")
  if (!is.null(command_status) && command_status != 0L) {
    stop("pdfinfo failed for ", path, ":\n", paste(output, collapse = "\n"))
  }
  page_line <- grep("^Pages:[[:space:]]+[0-9]+", output, value = TRUE)
  if (length(page_line) != 1L) stop("Could not determine PDF page count for ", path)
  as.integer(sub("^Pages:[[:space:]]+", "", page_line))
}

expected_pages <- sum(vapply(sample_plot_paths, pdf_page_count, integer(1)))
temporary_plot <- tempfile(
  pattern = "Xenium_Unknown_QC_All_Samples_", tmpdir = plot_dir, fileext = ".pdf"
)
on.exit(unlink(temporary_plot), add = TRUE)
merge_status <- system2(
  pdfunite,
  c(shQuote(sample_plot_paths), shQuote(temporary_plot))
)
if (!identical(merge_status, 0L) || !file.exists(temporary_plot) ||
    file.info(temporary_plot)$size <= 0) {
  stop("pdfunite failed to create the combined Unknown-QC plot.")
}
actual_pages <- pdf_page_count(temporary_plot)
if (actual_pages != expected_pages) {
  stop(
    "Combined Unknown-QC PDF page count mismatch: expected ", expected_pages,
    ", observed ", actual_pages, "."
  )
}
if (!file.rename(temporary_plot, combined_plot_path)) {
  stop("Could not move verified combined PDF to ", combined_plot_path)
}
message(
  "Merged ", length(sample_plot_paths), " sample QC PDFs (", actual_pages,
  " pages) in config/samples.csv order:\n- ", combined_plot_path
)
