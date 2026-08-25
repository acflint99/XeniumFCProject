#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
task_map <- data.frame(task_id = seq_len(nrow(samples)), sample_id = samples$sample_id)

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop("Usage: Rscript scripts/xenium_vz_01_subset_samples.R [--dry-run|--overwrite] [TASK_ID]")
}
task_value <- if (length(task_args) == 1L) task_args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(samples)) {
  stop("TASK_ID must be between 1 and ", nrow(samples), ". Use --list to inspect the mapping.")
}

current_sample <- samples$sample_id[[task_id]]
output_root <- here(config$project$outputs_dir)
input_path <- file.path(
  output_root, "xenium", "annotation", "03_consensus_labels", "rds",
  paste0(current_sample, "_Consensus_annotated.rds")
)
output_dir <- file.path(output_root, "xenium", "vz", "01_subsets", "rds")
output_path <- file.path(output_dir, paste0(current_sample, "_VZsubset.rds"))

if (dry_run) {
  cat("Task:", task_id, "of", nrow(samples), "\n")
  cat("Biological sample:", current_sample, "\n")
  cat("Consensus input:", input_path, "\n")
  cat("Output:", output_path, "\n")
  cat("Output exists:", file.exists(output_path), "\n")
  quit(save = "no", status = 0L)
}
if (!file.exists(input_path)) stop("Consensus-labelled input not found: ", input_path)
if (file.exists(output_path) && !overwrite) {
  stop("Refusing to overwrite existing VZ subset: ", output_path, "\nUse --overwrite only after review.")
}

obj <- readRDS(input_path)
if (!"consensus_label" %in% colnames(obj[[]])) {
  stop("Input object lacks 'consensus_label': ", input_path)
}
Idents(obj) <- "consensus_label"
target_clusters <- unlist(config$regional_subsets$vz_broad_labels, use.names = FALSE)
existing_clusters <- intersect(target_clusters, levels(Idents(obj)))
if (!length(existing_clusters)) stop("No configured VZ consensus labels found in ", current_sample, ".")
missing_clusters <- setdiff(target_clusters, existing_clusters)
if (length(missing_clusters)) {
  warning("VZ labels absent from ", current_sample, ": ", paste(missing_clusters, collapse = ", "))
}

obj_subset <- subset(obj, idents = existing_clusters)
Idents(obj_subset) <- "consensus_label"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_subset, output_path, compress = FALSE)
message("Saved VZ consensus subset for ", current_sample, ": ", output_path)
