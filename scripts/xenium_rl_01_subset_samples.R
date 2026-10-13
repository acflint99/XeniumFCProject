#!/usr/bin/env Rscript

# Build one RL regional subset from configured RL identities plus retained
# sample-qualified Unknown populations with two-of-three RL reference support.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
task_map <- data.frame(task_id = seq_along(sample_ids), sample_id = sample_ids)

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
  stop("Usage: Rscript scripts/xenium_rl_01_subset_samples.R [--dry-run|--overwrite] [TASK_ID]")
}
task_value <- if (length(task_args) == 1L) task_args[[1]] else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ". Use --list to inspect the mapping.")
}

current_sample <- sample_ids[[task_id]]
route_paths <- unknown_reference_route_paths(config)
output_root <- here(config$project$outputs_dir)
output_dir <- file.path(output_root, "xenium", "rl", "01_subsets", "rds")
output_path <- file.path(output_dir, paste0(current_sample, "_RLsubset.rds"))

manifest_state <- tryCatch({
  inputs <- load_unknown_reviewed_input_manifest(config)
  routes <- load_unknown_reference_routes(config)
  list(valid = TRUE, inputs = inputs, routes = routes)
}, error = function(error) list(valid = FALSE, error = conditionMessage(error)))

if (dry_run) {
  compact_dry_run(
    paste0("RL routed subset task ", task_id, "/", nrow(task_map), " [", current_sample, "]"),
    inputs = c(route_paths$input_manifest, route_paths$route_manifest),
    outputs = output_path,
    checks = c(route_manifests_valid = isTRUE(manifest_state$valid))
  )
  if (!isTRUE(manifest_state$valid)) cat("  Route error: ", manifest_state$error, "\n", sep = "")
  quit(save = "no", status = 0L)
}
if (!isTRUE(manifest_state$valid)) stop("Unknown routing is not ready: ", manifest_state$error)

input_row <- manifest_state$inputs[
  manifest_state$inputs$sample_id == current_sample, , drop = FALSE
]
if (nrow(input_row) != 1L) stop("Expected one reviewed input row for ", current_sample, ".")
input_path <- as.character(input_row$mapping_input_path[[1]])
if (!file.exists(input_path)) stop("Reviewed consensus input not found: ", input_path)
input_info <- file.info(input_path)
if (!isTRUE(all.equal(
      as.numeric(input_row$mapping_input_size[[1]]), as.numeric(input_info$size)
    )) || !isTRUE(all.equal(
      as.numeric(input_row$mapping_input_mtime[[1]]), as.numeric(input_info$mtime)
    ))) {
  stop("Reviewed input provenance changed for ", current_sample, ".")
}
if (file.exists(output_path) && !overwrite) {
  stop("Refusing to overwrite existing RL subset: ", output_path, "\nUse --overwrite only after review.")
}

obj <- readRDS(input_path)
if (anyDuplicated(colnames(obj))) stop("Duplicate cell IDs in ", current_sample, ".")
if (ncol(obj) != as.integer(input_row$remaining_cells[[1]])) {
  stop("Reviewed input cell count differs from its manifest for ", current_sample, ".")
}
if (!"consensus_label" %in% colnames(obj[[]])) {
  stop("Input object lacks consensus_label: ", input_path)
}
labels <- trimws(as.character(obj$consensus_label))
if (anyNA(labels) || any(!nzchar(labels))) stop("Blank consensus labels in ", current_sample, ".")

sample_routes <- manifest_state$routes[
  manifest_state$routes$sample_id == current_sample, , drop = FALSE
]
kept_routes <- sample_routes[tolower(sample_routes$qc_decision) == "keep", , drop = FALSE]
included_unknown_labels <- as.character(
  kept_routes$consensus_label[kept_routes$reference_route == "route_rl"]
)
missing_routed_labels <- setdiff(included_unknown_labels, unique(labels))
if (length(missing_routed_labels)) {
  stop("Kept RL-routed Unknown label(s) are absent from ", current_sample, ": ",
       paste(missing_routed_labels, collapse = ", "))
}
removed_labels_present <- intersect(
  as.character(sample_routes$consensus_label[tolower(sample_routes$qc_decision) == "remove"]),
  unique(labels)
)
if (length(removed_labels_present)) {
  stop("Reviewed-removed Unknown label(s) remain in the input: ",
       paste(removed_labels_present, collapse = ", "))
}

target_labels <- unlist(config$regional_subsets$rl_broad_labels, use.names = FALSE)
configured_present <- intersect(target_labels, unique(labels))
if (!length(configured_present)) stop("No configured RL consensus labels found in ", current_sample, ".")
missing_configured <- setdiff(target_labels, configured_present)
if (length(missing_configured)) {
  warning("RL labels absent from ", current_sample, ": ", paste(missing_configured, collapse = ", "))
}

include <- labels %in% c(configured_present, included_unknown_labels)
if (!any(include)) stop("RL routing retained no cells for ", current_sample, ".")
obj_subset <- subset(obj, cells = colnames(obj)[include])
subset_labels <- trimws(as.character(obj_subset$consensus_label))
unknown_cells <- is_unknown_consensus_label(subset_labels)
route_index <- match(subset_labels, sample_routes$consensus_label)
if (any(unknown_cells & is.na(route_index))) {
  stop("An RL-subset Unknown cell lacks a sample-qualified reference route.")
}
obj_subset$unknown_reference_route <- ifelse(
  unknown_cells, as.character(sample_routes$reference_route[route_index]), NA_character_
)
obj_subset$unknown_reference_route_support_n <- ifelse(
  unknown_cells,
  as.integer(sample_routes$reference_route_support_n[route_index]),
  NA_integer_
)
obj_subset$RL_subset_reason <- ifelse(
  unknown_cells, "reference_route_rl", "configured_broad_label"
)
Idents(obj_subset) <- "consensus_label"

production_paths <- get_production_output_paths(config)
obj_subset@misc$whole_tissue_production_source <- list(
  resolution = production_paths$resolution,
  consensus_stage = config$production$consensus_label_stage,
  input_path = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime)
)
obj_subset@misc$regional_subset_provenance <- list(
  branch = "RL",
  configured_broad_labels = target_labels,
  included_unknown_identities = sort(included_unknown_labels),
  configured_broad_cells = sum(!unknown_cells),
  reference_routed_unknown_cells = sum(unknown_cells),
  routing_method = config$regional_subsets$unknown_routing_method,
  minimum_reference_support = as.integer(
    config$regional_subsets$unknown_route_min_reference_support
  ),
  route_manifest = normalizePath(route_paths$route_manifest, mustWork = TRUE),
  input_manifest = normalizePath(route_paths$input_manifest, mustWork = TRUE)
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_subset, output_path, compress = FALSE)
message(
  "Saved RL subset for ", current_sample, ": ", output_path,
  " (", sum(!unknown_cells), " configured-label cells; ",
  sum(unknown_cells), " RL-routed Unknown cells)."
)
