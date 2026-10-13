#!/usr/bin/env Rscript

rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_replicate_marker_evidence.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) stop("Usage: Rscript scripts/xenium_vz_06a_build_replicate_marker_evidence.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

input_path <- here(config$project$outputs_dir, "xenium", "vz", "05_post_qc", "rds", "Xenium_VZ_postQC.rds")
output_dir <- here(config$project$outputs_dir, "xenium", "vz", "06_subclusters", "tables", "replicate_marker_evidence")
cluster_column <- "Xenium_snn_res.0.8"
outputs <- file.path(output_dir, c(
  "Xenium_VZ_ReplicateAware_Marker_Evidence.csv",
  "Xenium_VZ_Positive_Marker_Candidates.csv",
  "Xenium_VZ_AntiMarker_Candidates.csv",
  "Xenium_VZ_Marker_Sample_Coverage.csv",
  "Xenium_VZ_ReplicateAware_Marker_Provenance.csv"
))
if (dry_run) {
  compact_dry_run("VZ replicate-aware positive and anti-marker evidence", inputs = input_path, outputs = outputs)
  quit(save = "no", status = if (file.exists(input_path)) 0L else 1L)
}
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
build_replicate_marker_evidence(
  "VZ", input_path, cluster_column, output_dir, overwrite = overwrite
)
