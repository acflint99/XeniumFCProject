#!/usr/bin/env Rscript

# Create a compact, non-mutating score audit for user-reviewed Xenium clusters
# across the existing Aldinger, Sepp, and Science resolution-5 transfer RDS
# files. This script never writes an RDS or changes a label or threshold.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
references <- c("Aldinger", "Sepp", "Science")

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
positional_args <- args[!args %in% valid_options]
if (length(positional_args)) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_00c_audit_reviewed_clusters.R ",
    "[--list|--dry-run|--overwrite]"
  )
}

list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

review_manifest_path <- here("scripts", "config", "reviewed_annotation_clusters.csv")
if (!file.exists(review_manifest_path)) {
  review_manifest_path <- here("config", "reviewed_annotation_clusters.csv")
}
if (!file.exists(review_manifest_path)) {
  stop("Reviewed-cluster manifest not found: config/reviewed_annotation_clusters.csv")
}
reviewed <- read.csv(
  review_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = character()
)
required_manifest_columns <- c(
  "sample_id", "seurat_cluster", "current_consensus_label",
  "expected_label", "review_source"
)
missing_manifest_columns <- setdiff(required_manifest_columns, names(reviewed))
if (length(missing_manifest_columns)) {
  stop(
    "Reviewed-cluster manifest lacks columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}
for (column in required_manifest_columns) {
  reviewed[[column]] <- trimws(as.character(reviewed[[column]]))
  if (anyNA(reviewed[[column]]) || any(!nzchar(reviewed[[column]]))) {
    stop("Reviewed-cluster manifest has a blank or missing ", column, ".")
  }
}
review_keys <- paste(reviewed$sample_id, reviewed$seurat_cluster, sep = "|")
if (anyDuplicated(review_keys)) {
  stop("Reviewed-cluster manifest has duplicate sample/cluster rows.")
}
sample_matches <- vapply(
  reviewed$sample_id,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(sample_matches != 1L)) {
  stop(
    "Every reviewed sample must match exactly one config/samples.csv row. Invalid: ",
    paste(unique(reviewed$sample_id[sample_matches != 1L]), collapse = ", ")
  )
}

output_root <- here(config$project$outputs_dir)
annotation_root <- file.path(
  output_root, "xenium", "annotation"
)
input_grid <- do.call(rbind, lapply(unique(reviewed$sample_id), function(sample_id) {
  data.frame(
    sample_id = sample_id,
    reference = references,
    input_path = file.path(
      annotation_root,
      "01_label_transfer",
      tolower(references),
      "rds",
      paste0(sample_id, "_", references, "_annotated.rds")
    ),
    stringsAsFactors = FALSE
  )
}))

output_dir <- file.path(
  output_root,
  "validation",
  "resolution5_method_development",
  "04_transfer_threshold_calibration_audit",
  "reviewed_clusters",
  "tables"
)
output_paths <- c(
  compact_table = file.path(
    output_dir, "Reviewed_cluster_reference_score_audit.csv"
  ),
  provenance = file.path(
    output_dir, "Reviewed_cluster_reference_score_audit_provenance.csv"
  )
)

if (list_requested) {
  cat("Reviewed clusters:\n")
  write.table(
    reviewed[required_manifest_columns],
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  cat("\nReference-specific inputs:\n")
  write.table(input_grid, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

if (dry_run) {
  ok <- compact_dry_run(
    "Reviewed-cluster reference-score audit [5 clusters x 3 references]",
    inputs = c(review_manifest_path, input_grid$input_path),
    outputs = output_paths,
    checks = c(
      reviewed_rows_equal_5 = nrow(reviewed) == 5L,
      unique_reviewed_clusters = !anyDuplicated(review_keys),
      samples_match_manifest = all(sample_matches == 1L)
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- input_grid$input_path[!file.exists(input_grid$input_path)]
if (length(missing_inputs)) {
  stop("Missing reference-specific transfer RDS files:\n- ",
       paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite reviewed-cluster audit outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing the existing audit."
  )
}

legacy_threshold <- as.numeric(config$label_transfer$prediction_score_threshold)
focus_classes <- c("RL", "Granule", "Glia", "VZ")
audit_rows <- list()
row_index <- 1L

single_cluster_label <- function(metadata, column, cell_index, sample_id,
                                 cluster_id, reference) {
  if (!column %in% names(metadata)) return(NA_character_)
  labels <- unique(trimws(as.character(metadata[[column]][cell_index])))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) != 1L) {
    stop(
      sample_id, " cluster ", cluster_id, " has ", length(labels),
      " distinct ", column, " values in the ", reference, " object."
    )
  }
  labels[[1]]
}

for (grid_index in seq_len(nrow(input_grid))) {
  sample_id <- input_grid$sample_id[[grid_index]]
  reference <- input_grid$reference[[grid_index]]
  input_path <- input_grid$input_path[[grid_index]]
  message("Loading ", sample_id, " / ", reference, ": ", input_path)
  object <- readRDS(input_path)
  cell_ids <- Cells(object)
  metadata <- object[[]]
  if (anyDuplicated(cell_ids)) {
    stop("Duplicate cell IDs in ", sample_id, " / ", reference, ".")
  }
  if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Missing or duplicate metadata row names in ", sample_id, " / ", reference, ".")
  }
  if (!setequal(cell_ids, rownames(metadata))) {
    stop("Cell IDs and metadata rows do not match in ", sample_id, " / ", reference, ".")
  }
  metadata <- metadata[cell_ids, , drop = FALSE]
  required_metadata <- c("seurat_clusters", "predicted.id", "prediction.score.max")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata)) {
    stop(
      sample_id, " / ", reference, " lacks metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  score_columns <- grep("^prediction\\.score\\.", names(metadata), value = TRUE)
  score_columns <- setdiff(
    score_columns, c("prediction.score.max", "prediction.score.id")
  )
  if (length(score_columns) < 2L) {
    stop("Fewer than two class-score columns in ", sample_id, " / ", reference, ".")
  }
  class_labels <- sub("^prediction\\.score\\.", "", score_columns)
  if (anyDuplicated(class_labels)) {
    stop("Duplicate class-score labels in ", sample_id, " / ", reference, ".")
  }
  score_matrix <- as.matrix(metadata[, score_columns, drop = FALSE])
  storage.mode(score_matrix) <- "double"
  colnames(score_matrix) <- class_labels
  if (any(!is.finite(score_matrix))) {
    stop("Non-finite prediction scores in ", sample_id, " / ", reference, ".")
  }
  max_scores <- suppressWarnings(as.numeric(metadata$prediction.score.max))
  if (any(!is.finite(max_scores))) {
    stop("Non-finite prediction.score.max in ", sample_id, " / ", reference, ".")
  }
  predicted_ids <- trimws(as.character(metadata$predicted.id))
  if (anyNA(predicted_ids) || any(!nzchar(predicted_ids))) {
    stop("Missing or blank predicted.id in ", sample_id, " / ", reference, ".")
  }
  cluster_ids <- as.character(metadata$seurat_clusters)

  sample_targets <- reviewed[reviewed$sample_id == sample_id, , drop = FALSE]
  for (target_index in seq_len(nrow(sample_targets))) {
    target <- sample_targets[target_index, , drop = FALSE]
    cluster_id <- target$seurat_cluster[[1]]
    cell_index <- which(cluster_ids == cluster_id)
    if (!length(cell_index)) {
      stop(
        "Reviewed cluster ", cluster_id, " is absent from ", sample_id,
        " / ", reference, "."
      )
    }
    cluster_scores <- score_matrix[cell_index, , drop = FALSE]
    mean_scores <- colMeans(cluster_scores)
    score_order <- order(mean_scores, decreasing = TRUE)
    top_label <- names(mean_scores)[score_order[[1]]]
    runner_up_label <- names(mean_scores)[score_order[[2]]]
    current_label <- target$current_consensus_label[[1]]
    expected_label <- target$expected_label[[1]]
    score_or_na <- function(label) {
      if (label %in% names(mean_scores)) unname(mean_scores[[label]]) else NA_real_
    }

    output_row <- data.frame(
      sample_id = sample_id,
      seurat_cluster = cluster_id,
      current_consensus_label = current_label,
      expected_label = expected_label,
      reference = reference,
      n_cells = length(cell_index),
      reference_cluster_majority = single_cluster_label(
        metadata, "cluster_majority", cell_index, sample_id, cluster_id, reference
      ),
      reference_cluster_weighted = single_cluster_label(
        metadata, "cluster_weighted", cell_index, sample_id, cluster_id, reference
      ),
      top_mean_score_label = top_label,
      top_mean_score = unname(mean_scores[[top_label]]),
      runner_up_mean_score_label = runner_up_label,
      runner_up_mean_score = unname(mean_scores[[runner_up_label]]),
      top2_mean_score_margin = unname(
        mean_scores[[top_label]] - mean_scores[[runner_up_label]]
      ),
      expected_label_mean_score = score_or_na(expected_label),
      current_label_mean_score = score_or_na(current_label),
      expected_minus_current_mean_score = (
        score_or_na(expected_label) - score_or_na(current_label)
      ),
      fraction_cells_predicted_expected = mean(predicted_ids[cell_index] == expected_label),
      fraction_cells_predicted_current = mean(predicted_ids[cell_index] == current_label),
      fraction_cells_max_score_above_legacy_threshold = mean(
        max_scores[cell_index] > legacy_threshold
      ),
      mean_prediction_score_max = mean(max_scores[cell_index]),
      median_prediction_score_max = median(max_scores[cell_index]),
      legacy_threshold = legacy_threshold,
      threshold_alone_can_select_expected = identical(top_label, expected_label),
      review_source = target$review_source[[1]],
      stringsAsFactors = FALSE
    )
    for (class_label in focus_classes) {
      output_row[[paste0("mean_score_", class_label)]] <- score_or_na(class_label)
      output_row[[paste0("fraction_predicted_", class_label)]] <- mean(
        predicted_ids[cell_index] == class_label
      )
    }
    audit_rows[[row_index]] <- output_row
    row_index <- row_index + 1L
  }
  rm(object, metadata, score_matrix)
  invisible(gc(verbose = FALSE))
}

audit_table <- do.call(rbind, audit_rows)
expected_rows <- nrow(reviewed) * length(references)
if (nrow(audit_table) != expected_rows) {
  stop("Expected ", expected_rows, " audit rows; created ", nrow(audit_table), ".")
}
audit_key <- paste(
  audit_table$sample_id, audit_table$seurat_cluster, audit_table$reference,
  sep = "|"
)
if (anyDuplicated(audit_key)) stop("Audit unexpectedly duplicated a sample/cluster/reference row.")

sample_order <- match(audit_table$sample_id, reviewed$sample_id)
cluster_order <- match(
  paste(audit_table$sample_id, audit_table$seurat_cluster, sep = "|"),
  review_keys
)
reference_order <- match(audit_table$reference, references)
audit_table <- audit_table[
  order(sample_order, cluster_order, reference_order), , drop = FALSE
]

provenance <- data.frame(
  field = c(
    "purpose", "review_manifest", "resolution_mode", "references",
    "reviewed_clusters", "expected_rows", "legacy_threshold",
    "score_summary", "production_effect"
  ),
  value = c(
    "compact audit of user-reviewed mislabels",
    review_manifest_path,
    "resolution5_all_samples",
    paste(references, collapse = "|"),
    as.character(nrow(reviewed)),
    as.character(expected_rows),
    as.character(legacy_threshold),
    "unweighted cell mean within each existing Seurat cluster",
    "none; no object, label, threshold, or consensus table is modified"
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(audit_table, output_paths[["compact_table"]], row.names = FALSE, na = "")
write.csv(provenance, output_paths[["provenance"]], row.names = FALSE, na = "")
message("Wrote reviewed-cluster compact audit: ", output_paths[["compact_table"]])
