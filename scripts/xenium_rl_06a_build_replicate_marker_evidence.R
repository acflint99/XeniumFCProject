#!/usr/bin/env Rscript

# Optional standalone compatibility entry point. Production RL Stage-6 review
# now generates or reuses this evidence automatically; no separate job is needed.
rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_replicate_marker_evidence.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) stop("Usage: Rscript scripts/xenium_rl_06a_build_replicate_marker_evidence.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

input_path <- here(config$project$outputs_dir, "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC.rds")
output_dir <- here(config$project$outputs_dir, "xenium", "rl", "06_subclusters", "tables", "replicate_marker_evidence")
cluster_column <- "Xenium_snn_res.0.8"
outputs <- file.path(output_dir, c(
  "Xenium_RL_ReplicateAware_Marker_Evidence.csv",
  "Xenium_RL_Positive_Marker_Candidates.csv",
  "Xenium_RL_AntiMarker_Candidates.csv",
  "Xenium_RL_Marker_Sample_Coverage.csv",
  "Xenium_RL_ReplicateAware_Marker_Provenance.csv"
))
if (dry_run) {
  compact_dry_run("RL replicate-aware positive and anti-marker evidence", inputs = input_path, outputs = outputs)
  quit(save = "no", status = if (file.exists(input_path)) 0L else 1L)
}
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
build_replicate_marker_evidence(
  "RL", input_path, cluster_column, output_dir, overwrite = overwrite
)
