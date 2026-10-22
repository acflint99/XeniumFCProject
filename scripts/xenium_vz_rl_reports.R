#!/usr/bin/env Rscript

# Report-only aggregation of completed sample TIFFs; never loads or writes RDS.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_vz_rl_reports.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))
source(here("scripts", "R", "vz_rl_spatial_report_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
valid <- args %in% c("--dry-run", "--overwrite") |
  grepl("^--labels=(both|combined|corrected-broad)$", args) |
  grepl("^--point-size=(0\\.6|0\\.85)$", args)
if (any(!valid) || anyDuplicated(sub("=.*", "", args))) {
  stop("Usage: Rscript scripts/xenium_vz_rl_reports.R [--labels=both|combined|corrected-broad] ",
       "[--point-size=0.6|0.85] [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
option_value <- function(prefix, default) {
  value <- args[startsWith(args, paste0(prefix, "="))]
  if (length(value)) sub("^[^=]*=", "", value) else default
}
labels <- option_value("--labels", "both")
point_size <- option_value("--point-size", "0.6")
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
page_map <- build_stage8_pcw_page_map(config, samples)
specs <- vzrl_report_specs(here(config$project$outputs_dir), page_map, labels, point_size)
inputs <- unlist(lapply(specs, `[[`, "inputs"), use.names = FALSE)
outputs <- unlist(lapply(specs, function(spec) c(spec$pdf, spec$manifest)), use.names = FALSE)
missing_inputs <- inputs[!file.exists(inputs)]
dependencies <- c(tiff = requireNamespace("tiff", quietly = TRUE),
                  pdfinfo = nzchar(Sys.which("pdfinfo")))
input_info <- file.info(inputs)
nonempty <- !length(missing_inputs) && all(!is.na(input_info$size) & input_info$size > 0)
if (dry_run) {
  compact_dry_run("VZ/RL merged spatial reports", inputs = inputs, outputs = outputs,
                 checks = c(dependencies, input_tiffs_nonempty = nonempty))
  cat("  Labels: ", labels, "; point size: ", point_size,
      "; pages per PDF: ", nrow(page_map), "; order: numeric PCW then sample manifest\n", sep = "")
  for (spec in specs) {
    cat("  ", spec$title, " inputs: ", spec$input_dir, "\n", sep = "")
    if (any(!file.exists(spec$inputs))) {
      cat("    Checked layout locations: ", paste(spec$searched_dirs, collapse = "; "), "\n", sep = "")
    }
  }
  cat("  PDF outputs:\n- ", paste(vapply(specs, `[[`, "", "pdf"), collapse = "\n- "), "\n", sep = "")
  quit(save = "no", status = if (!nonempty || !all(dependencies)) 1L else 0L)
}
if (length(missing_inputs)) stop("Missing spatial TIFFs:\n- ", paste(missing_inputs, collapse = "\n- "))
if (!nonempty) stop("A spatial TIFF is empty.")
if (!all(dependencies)) stop("Missing report dependencies: ",
                            paste(names(dependencies)[!dependencies], collapse = ", "),
                            ". Use the existing project renv environment and Poppler; no packages are installed by this script.")
if (any(file.exists(outputs)) && !overwrite) stop("Report outputs already exist; review before --overwrite.")
for (spec in specs) write_vzrl_spatial_report(spec, page_map, overwrite)
