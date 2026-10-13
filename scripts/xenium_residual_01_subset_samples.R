#!/usr/bin/env Rscript

# Extract retained Unknown populations routed to Other or unresolved for a
# separate residual analysis. Detailed consensus labels and cell IDs are kept.

rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
args <- commandArgs(trailingOnly = TRUE)
valid <- c("--list", "--dry-run", "--overwrite")
unknown <- args[startsWith(args, "--") & !args %in% valid]
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
positional <- args[!args %in% valid]
if (list_requested) {
  if (length(positional)) stop("--list does not accept TASK_ID.")
  write.table(
    data.frame(task_id = seq_along(sample_ids), sample_id = sample_ids),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = 0L)
}
if (length(positional) != 1L) {
  stop("Usage: Rscript scripts/xenium_residual_01_subset_samples.R [--list|--dry-run|--overwrite] TASK_ID")
}
task_id <- suppressWarnings(as.integer(positional[[1]]))
if (is.na(task_id) || task_id < 1L || task_id > length(sample_ids)) {
  stop("TASK_ID must be between 1 and ", length(sample_ids), ".")
}
sample_id <- sample_ids[[task_id]]
production_paths <- get_production_output_paths(config)
reviewed_inputs <- load_unknown_reviewed_input_manifest(config)
input_row <- reviewed_inputs[reviewed_inputs$sample_id == sample_id, , drop = FALSE]
if (nrow(input_row) != 1L) stop("Expected one reviewed input row for ", sample_id, ".")
required_input_columns <- c(
  "mapping_input_path", "mapping_input_size", "mapping_input_mtime",
  "remaining_cells"
)
missing_input_columns <- setdiff(required_input_columns, names(input_row))
if (length(missing_input_columns)) {
  stop(
    "Reviewed input manifest lacks required column(s): ",
    paste(missing_input_columns, collapse = ", ")
  )
}
input_path <- trimws(as.character(input_row$mapping_input_path[[1]]))
if (length(input_path) != 1L || is.na(input_path) || !nzchar(input_path)) {
  stop("Reviewed input path is blank or invalid for ", sample_id, ".")
}
route_path <- unknown_reference_route_paths(config)$route_manifest
output_root <- here(config$project$outputs_dir, "xenium", "residual", "01_subsets")
rds_path <- file.path(output_root, "rds", paste0(sample_id, "_ResidualSubset.rds"))
status_path <- file.path(output_root, "tables", paste0(sample_id, "_ResidualSubset_Status.csv"))

if (dry_run) {
  compact_dry_run(
    paste0("Residual subset task ", task_id, "/", length(sample_ids), " [", sample_id, "]"),
    inputs = c(input_path, route_path), outputs = c(rds_path, status_path)
  )
  quit(save = "no", status = if (all(file.exists(c(input_path, route_path)))) 0L else 1L)
}
missing <- c(input_path, route_path)[!file.exists(c(input_path, route_path))]
if (length(missing)) stop("Missing residual-subset input(s):\n- ", paste(missing, collapse = "\n- "))
input_info <- file.info(input_path)
if (!isTRUE(all.equal(
      as.numeric(input_row$mapping_input_size[[1]]), as.numeric(input_info$size)
    )) || !isTRUE(all.equal(
      as.numeric(input_row$mapping_input_mtime[[1]]), as.numeric(input_info$mtime)
    ))) {
  stop("Reviewed input provenance changed for ", sample_id, ".")
}
existing <- c(rds_path, status_path)[file.exists(c(rds_path, status_path))]
if (length(existing) && !overwrite) {
  stop("Refusing to overwrite residual-subset output(s):\n- ", paste(existing, collapse = "\n- "))
}
suppressPackageStartupMessages(library(Seurat))
routes <- load_unknown_reference_routes(config)
sample_routes <- routes[
  routes$sample_id == sample_id & tolower(routes$qc_decision) == "keep" &
    routes$reference_route %in% c("route_other", "route_unresolved"),
  , drop = FALSE
]
object <- readRDS(input_path)
if (anyDuplicated(Cells(object))) stop("Duplicate cell IDs in ", sample_id, ".")
if (ncol(object) != as.integer(input_row$remaining_cells[[1]])) {
  stop("Reviewed input cell count differs from its manifest for ", sample_id, ".")
}
if (!"consensus_label" %in% names(object[[]])) stop("Input lacks consensus_label.")
labels <- trimws(as.character(object$consensus_label))
target_labels <- as.character(sample_routes$consensus_label)
missing_labels <- setdiff(target_labels, unique(labels))
if (length(missing_labels)) {
  stop("Residual route label(s) are absent from ", sample_id, ": ", paste(missing_labels, collapse = ", "))
}
include <- labels %in% target_labels
status <- data.frame(
  task_id = task_id,
  sample_id = sample_id,
  input_rds = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(file.info(input_path)$size),
  input_mtime = as.numeric(file.info(input_path)$mtime),
  has_residual_cells = any(include),
  residual_cells = sum(include),
  residual_populations = length(target_labels),
  residual_labels = paste(sort(target_labels), collapse = "|"),
  route_manifest = normalizePath(route_path, mustWork = TRUE),
  route_manifest_size = as.numeric(file.info(route_path)$size),
  route_manifest_mtime = as.numeric(file.info(route_path)$mtime),
  stringsAsFactors = FALSE
)
dir.create(dirname(status_path), recursive = TRUE, showWarnings = FALSE)
if (any(include)) {
  subset_object <- subset(object, cells = colnames(object)[include])
  subset_labels <- as.character(subset_object$consensus_label)
  route_index <- match(subset_labels, sample_routes$consensus_label)
  if (anyNA(route_index)) stop("Residual subset contains a population without a route row.")
  subset_object$residual_reference_route <- sample_routes$reference_route[route_index]
  subset_object$residual_reference_support_n <- as.integer(
    sample_routes$reference_route_support_n[route_index]
  )
  subset_object$residual_reference_status <- sample_routes$reference_route_status[route_index]
  subset_object$residual_subset_reason <- "retained_high_quality_other_or_unresolved"
  subset_object@misc$residual_subset_provenance <- as.list(status[1, , drop = FALSE])
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(subset_object, rds_path, compress = FALSE)
} else if (file.exists(rds_path) && overwrite) {
  stop(
    "An obsolete residual RDS exists for a now-empty sample. Move or remove it ",
    "after review before rerunning: ", rds_path
  )
}
write.csv(status, status_path, row.names = FALSE)
message("Residual subset status saved for ", sample_id, ": ", sum(include), " cells.")
