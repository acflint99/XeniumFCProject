#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
production_paths <- get_production_output_paths(config)
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
  production_paths$consensus_rds_dir,
  paste0(current_sample, "_Consensus_annotated.rds")
)
output_dir <- file.path(output_root, "xenium", "vz", "01_subsets", "rds")
output_path <- file.path(output_dir, paste0(current_sample, "_VZsubset.rds"))

if (dry_run) {
  compact_dry_run(
    paste0("VZ subset task ", task_id, "/", nrow(samples), " [", current_sample, "]"),
    inputs = input_path,
    outputs = output_path,
    checks = c(
      include_unknown_in_both_branches = identical(
        config$regional_subsets$include_unknown_in_both_branches, TRUE
      )
    )
  )
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

if (!identical(config$regional_subsets$include_unknown_in_both_branches, TRUE)) {
  stop("Configuration must include unresolved identities in both regional branches.")
}
unknown_clusters <- levels(Idents(obj))[is_unknown_consensus_label(levels(Idents(obj)))]
subset_clusters <- unique(c(existing_clusters, unknown_clusters))
obj_subset <- subset(obj, idents = subset_clusters)
unknown_cells <- is_unknown_consensus_label(obj_subset$consensus_label)
obj_subset$VZ_subset_reason <- ifelse(
  unknown_cells, "shared_unknown", "configured_broad_label"
)
Idents(obj_subset) <- "consensus_label"
input_info <- file.info(input_path)
obj_subset@misc$whole_tissue_production_source <- list(
  resolution = production_paths$resolution,
  consensus_stage = config$production$consensus_label_stage,
  input_path = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime)
)
obj_subset@misc$regional_subset_provenance <- list(
  branch = "VZ",
  configured_broad_labels = target_clusters,
  included_unknown_identities = unknown_clusters,
  configured_broad_cells = sum(!unknown_cells),
  shared_unknown_cells = sum(unknown_cells),
  include_unknown_in_both_branches = TRUE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_subset, output_path, compress = FALSE)
message(
  "Saved VZ consensus subset for ", current_sample, ": ", output_path,
  " (", sum(!unknown_cells), " configured-label cells; ",
  sum(unknown_cells), " shared Unknown cells)."
)
