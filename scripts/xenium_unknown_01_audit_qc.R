#!/usr/bin/env Rscript

# Production entry point for one-sample-at-a-time QC auditing of all
# sample-qualified Unknown populations. This wrapper fixes the approved input
# to the resolution-5 weighted-two-of-three consensus; users do not select an
# analysis mode.

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
    "Usage: Rscript scripts/xenium_unknown_01_audit_qc.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]"
  )
}

driver <- here("scripts", "xenium_annotate_04_audit_unknown_qc.R")
driver_args <- c(
  shQuote(driver), "--all-samples-res5", "--weighted-2of3",
  args
)
status <- system2(file.path(R.home("bin"), "Rscript"), driver_args)
if (!identical(status, 0L)) {
  stop("Production Unknown-QC audit driver exited with status ", status, call. = FALSE)
}
