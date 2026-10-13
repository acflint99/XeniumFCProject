#!/usr/bin/env Rscript

# Apply reviewed whole-population keep/remove decisions to one production
# sample. This entry point never combines tasks or selects an analysis mode.

rm(list = ls())

suppressPackageStartupMessages(library(here))

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
positional <- args[!args %in% valid_options]
if (length(positional) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_unknown_03_apply_qc_review.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]"
  )
}

driver <- here("scripts", "xenium_annotate_04b_apply_unknown_qc.R")
driver_args <- c(shQuote(driver), args)
status <- system2(file.path(R.home("bin"), "Rscript"), driver_args)
if (!identical(status, 0L)) {
  stop("Reviewed Unknown-QC application driver exited with status ", status, call. = FALSE)
}
