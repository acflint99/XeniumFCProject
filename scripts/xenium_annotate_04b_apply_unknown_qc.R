#!/usr/bin/env Rscript

# Apply the reviewed post-consensus Unknown-QC decisions once, before either
# regional branch is subset or merged. Kept samples reuse the consensus RDS;
# samples with reviewed removals receive a protected filtered copy.
# Production users should invoke the single-purpose
# xenium_unknown_03_apply_qc_review.R entry point.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)
task_map <- data.frame(
  task_id = seq_along(sample_ids), sample_id = sample_ids,
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--combine", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
list_requested <- "--list" %in% args
combine <- "--combine" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

if (list_requested) {
  positional <- args[!args %in% valid_options]
  if (length(positional)) stop("--list does not accept TASK_ID.")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

output_root <- here(config$project$outputs_dir)
audit_table_dir <- file.path(
  output_root, "xenium", "annotation", "04_unknown_qc_review", "tables"
)
audit_template_path <- file.path(
  audit_table_dir, "Xenium_Unknown_QC_Review_Template.csv"
)
audit_manifest_path <- file.path(
  audit_table_dir, "Xenium_Unknown_QC_Audit_manifest.csv"
)
reviewed_path <- resolve_config_path(config$manifests$reviewed_unknown_qc, config)
stage_dir <- file.path(
  output_root, "xenium", "annotation", "04_unknown_qc_reviewed"
)
rds_dir <- file.path(stage_dir, "rds")
table_dir <- file.path(stage_dir, "tables")
per_sample_dir <- file.path(table_dir, "per_sample")
master_routing_path <- file.path(
  table_dir, "Xenium_Unknown_QC_routing_manifest.csv"
)

source_paths <- file.path(
  production_paths$consensus_rds_dir,
  paste0(sample_ids, "_Consensus_annotated.rds")
)
names(source_paths) <- sample_ids
filtered_paths <- file.path(
  rds_dir, paste0(sample_ids, "_Consensus_UnknownQC_reviewed.rds")
)
names(filtered_paths) <- sample_ids
per_sample_routing_paths <- file.path(
  per_sample_dir, paste0(sample_ids, "_Unknown_QC_routing.csv")
)
names(per_sample_routing_paths) <- sample_ids
removed_cell_paths <- file.path(
  per_sample_dir, paste0(sample_ids, "_Unknown_QC_removed_cells.csv")
)
names(removed_cell_paths) <- sample_ids

review_state <- tryCatch({
  if (!file.exists(audit_template_path)) stop("Missing audit template: ", audit_template_path)
  if (!file.exists(audit_manifest_path)) stop("Missing audit manifest: ", audit_manifest_path)
  if (!file.exists(reviewed_path)) stop("Missing reviewed manifest: ", reviewed_path)
  template <- read.csv(audit_template_path, stringsAsFactors = FALSE, na.strings = character())
  audit_manifest <- read.csv(audit_manifest_path, stringsAsFactors = FALSE, na.strings = character())
  review <- load_reviewed_unknown_qc(config)
  expected_keys <- unknown_qc_review_key(template)
  validate_unknown_qc_review(review, sample_ids, expected_keys, require_complete = TRUE)
  expected_provenance_columns <- c(
    "sample_id", "input_path", "input_size", "input_mtime", "input_cells"
  )
  missing_provenance <- setdiff(expected_provenance_columns, names(audit_manifest))
  if (length(missing_provenance)) {
    stop("Audit manifest lacks column(s): ", paste(missing_provenance, collapse = ", "))
  }
  if (!identical(as.character(audit_manifest$sample_id), sample_ids)) {
    stop("Audit manifest rows do not match config/samples.csv order.")
  }
  expected_sources <- normalizePath(source_paths, mustWork = FALSE)
  recorded_sources <- normalizePath(as.character(audit_manifest$input_path), mustWork = FALSE)
  if (!identical(recorded_sources, unname(expected_sources))) {
    stop("Audit manifest input paths do not match the production consensus inputs.")
  }
  list(valid = TRUE, template = template, audit = audit_manifest, review = review)
}, error = function(error) list(valid = FALSE, error = conditionMessage(error)))

positional <- args[!args %in% valid_options]
if (combine) {
  if (length(positional)) stop("--combine does not accept TASK_ID.")
  if (dry_run) {
    compact_dry_run(
      "Reviewed Unknown-QC routing combine",
      inputs = c(per_sample_routing_paths, reviewed_path, audit_manifest_path),
      outputs = master_routing_path,
      checks = c(review_complete_and_current = isTRUE(review_state$valid))
    )
    if (!isTRUE(review_state$valid)) cat("  Review error: ", review_state$error, "\n")
    quit(save = "no", status = 0L)
  }
  if (!isTRUE(review_state$valid)) stop("Unknown-QC review is not ready: ", review_state$error)
  missing <- per_sample_routing_paths[!file.exists(per_sample_routing_paths)]
  if (length(missing)) stop("Missing per-sample routing files:\n- ", paste(missing, collapse = "\n- "))
  if (file.exists(master_routing_path) && !overwrite) {
    stop("Refusing to overwrite routing manifest: ", master_routing_path, "\nUse --overwrite after review.")
  }
  routing <- do.call(rbind, lapply(per_sample_routing_paths, function(path) {
    read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  }))
  inspection <- inspect_unknown_qc_routing_manifest(routing, sample_ids)
  if (!inspection$valid) {
    stop("Combined Unknown-QC routing is invalid:\n- ", paste(inspection$errors, collapse = "\n- "))
  }
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(routing, master_routing_path, row.names = FALSE, na = "")
  message(
    "Wrote combined reviewed Unknown input manifest for ", nrow(routing), " samples: ",
    master_routing_path
  )
  quit(save = "no", status = 0L)
}

if (length(positional) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_04b_apply_unknown_qc.R ",
    "[--list|--combine|--dry-run|--overwrite] [TASK_ID]"
  )
}
task_value <- if (length(positional)) positional[[1]] else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ". Use --list to inspect mapping.")
}

sample_id <- task_map$sample_id[[task_id]]
source_path <- source_paths[[sample_id]]
filtered_path <- filtered_paths[[sample_id]]
routing_path <- per_sample_routing_paths[[sample_id]]
removed_path <- removed_cell_paths[[sample_id]]
if (dry_run) {
  compact_dry_run(
    paste0("Apply reviewed Unknown-QC task ", task_id, " [", sample_id, "]"),
    inputs = c(source_path, audit_template_path, audit_manifest_path, reviewed_path),
    outputs = c(routing_path, removed_path),
    checks = c(review_complete_and_current = isTRUE(review_state$valid))
  )
  if (!isTRUE(review_state$valid)) cat("  Review error: ", review_state$error, "\n")
  quit(save = "no", status = 0L)
}
if (!isTRUE(review_state$valid)) stop("Unknown-QC review is not ready: ", review_state$error)
if (!file.exists(source_path)) stop("Consensus input not found: ", source_path)

existing_tabular <- c(routing_path, removed_path)[file.exists(c(routing_path, removed_path))]
if (length(existing_tabular) && !overwrite) {
  stop("Refusing to overwrite Unknown-QC task outputs:\n- ", paste(existing_tabular, collapse = "\n- "))
}

audit_row <- review_state$audit[review_state$audit$sample_id == sample_id, , drop = FALSE]
source_info <- file.info(source_path)
source_current <- nrow(audit_row) == 1L &&
  isTRUE(all.equal(as.numeric(audit_row$input_size), as.numeric(source_info$size))) &&
  isTRUE(all.equal(as.numeric(audit_row$input_mtime), as.numeric(source_info$mtime)))
if (!source_current) stop(sample_id, " consensus object changed after the Unknown-QC audit.")

obj <- readRDS(source_path)
if (!"consensus_label" %in% colnames(obj[[]])) stop(sample_id, " lacks consensus_label metadata.")
labels <- trimws(as.character(obj$consensus_label))
if (anyNA(labels) || any(!nzchar(labels))) stop(sample_id, " contains blank consensus labels.")
sample_review <- review_state$review[review_state$review$sample_id == sample_id, , drop = FALSE]
sample_template <- review_state$template[
  review_state$template$sample_id == sample_id, , drop = FALSE
]
if (!setequal(unknown_qc_review_key(sample_review), unknown_qc_review_key(sample_template))) {
  stop(sample_id, " reviewed Unknown-QC keys differ from the audit template.")
}
remove_labels <- trimws(as.character(
  sample_review$consensus_label[tolower(trimws(sample_review$decision)) == "remove"]
))
missing_remove_labels <- setdiff(remove_labels, unique(labels))
if (length(missing_remove_labels)) {
  stop(sample_id, " reviewed removal label(s) are absent: ", paste(missing_remove_labels, collapse = ", "))
}
remove_cells <- colnames(obj)[labels %in% remove_labels]
removed_table <- data.frame(
  sample_id = rep(sample_id, length(remove_cells)),
  cell_id = remove_cells,
  consensus_label = labels[match(remove_cells, colnames(obj))],
  stringsAsFactors = FALSE
)
input_cells <- ncol(obj)
if (length(remove_cells) >= input_cells) stop(sample_id, " review would remove every cell.")

if (length(remove_cells)) {
  if (file.exists(filtered_path) && !overwrite) {
    stop("Refusing to overwrite filtered consensus object: ", filtered_path)
  }
  obj <- subset(obj, cells = setdiff(colnames(obj), remove_cells))
  obj@misc$unknown_qc_review <- list(
    reviewed_manifest = normalizePath(reviewed_path, mustWork = TRUE),
    source_path = normalizePath(source_path, mustWork = TRUE),
    removed_unknown_labels = remove_labels,
    removed_unknown_cells = length(remove_cells)
  )
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, filtered_path, compress = FALSE)
  mapping_path <- filtered_path
  mode <- "filtered_copy"
} else {
  if (file.exists(filtered_path)) {
    stop(
      "A stale filtered copy exists for a zero-removal decision: ", filtered_path,
      ". Archive or remove it explicitly before continuing."
    )
  }
  mapping_path <- source_path
  mode <- "reuse_consensus"
}

mapping_info <- file.info(mapping_path)
reviewed_info <- file.info(reviewed_path)
routing <- data.frame(
  sample_id = sample_id,
  source_path = normalizePath(source_path, mustWork = TRUE),
  mapping_input_path = normalizePath(mapping_path, mustWork = TRUE),
  materialization_mode = mode,
  reviewed_qc_path = normalizePath(reviewed_path, mustWork = TRUE),
  reviewed_qc_size = as.numeric(reviewed_info$size),
  reviewed_qc_mtime = as.numeric(reviewed_info$mtime),
  source_size = as.numeric(source_info$size),
  source_mtime = as.numeric(source_info$mtime),
  mapping_input_size = as.numeric(mapping_info$size),
  mapping_input_mtime = as.numeric(mapping_info$mtime),
  input_cells = input_cells,
  removed_unknown_cells = length(remove_cells),
  remaining_cells = input_cells - length(remove_cells),
  removed_unknown_labels = if (length(remove_labels)) {
    paste(sort(remove_labels), collapse = "|")
  } else "none",
  stringsAsFactors = FALSE
)
dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(removed_table, removed_path, row.names = FALSE, na = "")
write.csv(routing, routing_path, row.names = FALSE, na = "")
message(
  sample_id, ": removed ", length(remove_cells), " reviewed low-QC Unknown cells; ",
  "routing mode = ", mode, "."
)
