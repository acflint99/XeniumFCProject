#!/usr/bin/env Rscript

# Focused resolution-5 RPCA sensitivity analysis for k.anchor = 5 and 10.
# Runs three requested Xenium samples against Aldinger, Sepp, and Science while
# preserving every other production transfer setting. Outputs are compact
# tables only; no Seurat object or production annotation output is written.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)

requested_samples <- c(
  "GZFB_9_X_G_1",
  "GZFB_20_X_G_1",
  "GZFB_22_X_G_3"
)
references <- c("Aldinger", "Sepp", "Science")
k_anchor_values <- c(5L, 10L)

sample_counts <- vapply(
  requested_samples,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(sample_counts != 1L)) {
  stop(
    "Every k.anchor sensitivity sample must match exactly one config/samples.csv row. Invalid: ",
    paste(requested_samples[sample_counts != 1L], collapse = ", ")
  )
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

task_map <- expand.grid(
  sample_id = requested_samples,
  reference = references,
  k_anchor = k_anchor_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_map <- task_map[order(
  match(task_map$sample_id, requested_samples),
  match(task_map$reference, references),
  match(task_map$k_anchor, k_anchor_values)
), , drop = FALSE]
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c("task_id", "sample_id", "reference", "k_anchor")]
if (nrow(task_map) != 18L) stop("Expected exactly 18 k.anchor sensitivity tasks.")

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite", "--combine")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
combine <- "--combine" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

positional_args <- args[!args %in% valid_options]
if (list_requested) {
  if (length(positional_args)) stop("--list does not accept TASK_ID.")
  cat("Mode: resolution-5 all-sample k.anchor sensitivity\n")
  cat("Production settings retained: RPCA, k.score=30, k.weight=50, 30 dimensions, 1000 reference cells/identity, threshold=0.4\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

preprocess_root <- file.path(
  here(config$project$outputs_dir), "xenium", "preprocess",
  "04_resolution5_clustered", "rds"
)
annotation_root <- file.path(
  here(config$project$outputs_dir), "validation",
  "resolution5_method_development", "k_anchor_sensitivity"
)

task_input_path <- function(sample_id) {
  file.path(preprocess_root, paste0(sample_id, "_whole_tissue_Res5.0.rds"))
}

task_output_paths <- function(sample_id, reference, k_anchor) {
  output_dir <- file.path(
    annotation_root,
    paste0("k_anchor_", k_anchor),
    tolower(reference),
    "tables"
  )
  prefix <- paste0(sample_id, "_", reference, "_kanchor", k_anchor)
  c(
    cluster_scores = file.path(
      output_dir, paste0(prefix, "_cluster_score_summary.csv")
    ),
    prediction_counts = file.path(
      output_dir, paste0(prefix, "_prediction_cellcounts.csv")
    ),
    provenance = file.path(
      output_dir, paste0(prefix, "_provenance.csv")
    )
  )
}

all_task_outputs <- do.call(rbind, lapply(seq_len(nrow(task_map)), function(i) {
  paths <- task_output_paths(
    task_map$sample_id[[i]], task_map$reference[[i]], task_map$k_anchor[[i]]
  )
  data.frame(
    task_id = task_map$task_id[[i]],
    sample_id = task_map$sample_id[[i]],
    reference = task_map$reference[[i]],
    k_anchor = task_map$k_anchor[[i]],
    output_type = names(paths),
    output_path = unname(paths),
    stringsAsFactors = FALSE
  )
}))

combined_output_dir <- file.path(annotation_root, "combined_tables")
combined_outputs <- c(
  cluster_scores = file.path(
    combined_output_dir, "k_anchor_5_10_combined_cluster_score_summary.csv"
  ),
  provenance = file.path(
    combined_output_dir, "k_anchor_5_10_combined_provenance.csv"
  )
)

if (combine) {
  if (length(positional_args)) stop("--combine does not accept TASK_ID.")
  cluster_inputs <- all_task_outputs$output_path[
    all_task_outputs$output_type == "cluster_scores"
  ]
  if (dry_run) {
    ok <- compact_dry_run(
      "Combine 18 k.anchor sensitivity cluster-score tables",
      inputs = cluster_inputs,
      outputs = combined_outputs,
      checks = c(expected_task_tables = length(cluster_inputs) == 18L)
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }
  missing_inputs <- cluster_inputs[!file.exists(cluster_inputs)]
  if (length(missing_inputs)) {
    stop("Missing k.anchor sensitivity tables:\n- ",
         paste(missing_inputs, collapse = "\n- "))
  }
  existing_outputs <- combined_outputs[file.exists(combined_outputs)]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite combined k.anchor sensitivity outputs:\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing the existing outputs."
    )
  }
  tables <- lapply(cluster_inputs, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  combined <- do.call(rbind, tables)
  combined_key <- paste(
    combined$sample_id, combined$reference, combined$k_anchor,
    combined$seurat_cluster, sep = "|"
  )
  if (anyDuplicated(combined_key)) {
    stop("Combined k.anchor table has duplicate sample/reference/k/cluster rows.")
  }
  provenance <- data.frame(
    field = c(
      "purpose", "samples", "references", "k_anchor_values",
      "task_tables", "combined_rows", "production_effect"
    ),
    value = c(
      "focused resolution-5 RPCA k.anchor sensitivity",
      paste(requested_samples, collapse = "|"),
      paste(references, collapse = "|"),
      paste(k_anchor_values, collapse = "|"),
      as.character(length(cluster_inputs)),
      as.character(nrow(combined)),
      "none; no RDS, production transfer, or consensus output is modified"
    ),
    stringsAsFactors = FALSE
  )
  dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(combined, combined_outputs[["cluster_scores"]], row.names = FALSE, na = "")
  write.csv(provenance, combined_outputs[["provenance"]], row.names = FALSE, na = "")
  message("Wrote combined k.anchor sensitivity table: ", combined_outputs[["cluster_scores"]])
  quit(save = "no", status = 0L)
}

if (length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]\n",
    "   or: Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R ",
    "--combine [--dry-run|--overwrite]"
  )
}
task_value <- if (length(positional_args)) {
  positional_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ". Use --list.")
}

sample_id <- task_map$sample_id[[task_id]]
reference <- task_map$reference[[task_id]]
k_anchor <- as.integer(task_map$k_anchor[[task_id]])
input_path <- task_input_path(sample_id)
reference_path <- unname(reference_paths[[reference]])
output_paths <- task_output_paths(sample_id, reference, k_anchor)

if (dry_run) {
  ok <- compact_dry_run(
    paste0(
      "k.anchor sensitivity task ", task_id, "/", nrow(task_map), " [",
      sample_id, " | ", reference, " | k.anchor=", k_anchor, "]"
    ),
    inputs = c(input_path, reference_path),
    outputs = output_paths,
    checks = c(
      requested_k_anchor = k_anchor %in% k_anchor_values,
      sample_matches_manifest = sum(sample_manifest$sample_id == sample_id) == 1L
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- c(input_path, reference_path)[!file.exists(c(input_path, reference_path))]
if (length(missing_inputs)) {
  stop("Missing k.anchor sensitivity input(s):\n- ",
       paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite k.anchor sensitivity outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing tables."
  )
}

set.seed(as.integer(config$runtime$random_seed))
query <- readRDS(input_path)
reference_object <- readRDS(reference_path)

if (anyDuplicated(Cells(query))) stop("Query contains duplicate cell IDs: ", sample_id)
if (anyDuplicated(Cells(reference_object))) {
  stop("Reference contains duplicate cell IDs: ", reference)
}
if (!"clusters_refined" %in% colnames(reference_object[[]])) {
  stop(reference, " reference lacks clusters_refined metadata.")
}
reference_labels <- trimws(as.character(reference_object$clusters_refined))
if (anyNA(reference_labels) || any(!nzchar(reference_labels))) {
  stop(reference, " reference has blank or missing clusters_refined labels.")
}
reference_object$clusters_refined <- reference_labels

required_query_metadata <- c(
  "whole_tissue_cluster_res1.5", "whole_tissue_cluster_res5.0",
  "Xenium_snn_res.5", "seurat_clusters"
)
missing_query_metadata <- setdiff(required_query_metadata, colnames(query[[]]))
if (length(missing_query_metadata)) {
  stop(
    sample_id, " resolution-5 input lacks metadata: ",
    paste(missing_query_metadata, collapse = ", ")
  )
}
resolution_clusters <- as.character(query$whole_tissue_cluster_res5.0)
active_clusters <- as.character(query$seurat_clusters)
if (anyNA(resolution_clusters) || !identical(resolution_clusters, active_clusters)) {
  stop("Query seurat_clusters does not exactly match whole_tissue_cluster_res5.0.")
}

if ("RNA" %in% Assays(reference_object)) {
  DefaultAssay(reference_object) <- "RNA"
} else {
  DefaultAssay(reference_object) <- Assays(reference_object)[[1]]
}
if ("Xenium" %in% Assays(query)) {
  DefaultAssay(query) <- "Xenium"
} else {
  DefaultAssay(query) <- Assays(query)[[1]]
}

variable_features <- as.integer(config$label_transfer$variable_features)
reference_object <- FindVariableFeatures(
  reference_object,
  selection.method = "vst",
  nfeatures = variable_features,
  verbose = FALSE
)
shared_genes <- intersect(rownames(reference_object), rownames(query))
Idents(reference_object) <- "clusters_refined"
reference_balanced <- subset(
  reference_object,
  downsample = as.integer(config$label_transfer$reference_cells_per_identity)
)
transfer_features <- intersect(VariableFeatures(reference_balanced), shared_genes)
transfer_dimensions <- seq_len(as.integer(config$label_transfer$dimensions))
if (length(transfer_features) <= max(transfer_dimensions)) {
  stop(
    "Only ", length(transfer_features), " transfer features are available; more than ",
    max(transfer_dimensions), " are required."
  )
}

reference_balanced <- ScaleData(
  reference_balanced, features = transfer_features, verbose = FALSE
)
query <- ScaleData(query, features = transfer_features, verbose = FALSE)
reference_balanced <- RunPCA(
  reference_balanced, features = transfer_features, verbose = FALSE
)
query <- RunPCA(query, features = transfer_features, verbose = FALSE)

k_score <- as.integer(config$label_transfer$k_score)
k_weight <- 50L
anchors <- FindTransferAnchors(
  reference = reference_balanced,
  query = query,
  normalization.method = "LogNormalize",
  reduction = config$label_transfer$method,
  features = transfer_features,
  dims = transfer_dimensions,
  k.anchor = k_anchor,
  k.score = k_score,
  approx.pca = TRUE
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference_balanced$clusters_refined,
  dims = transfer_dimensions,
  k.weight = k_weight,
  store.weights = FALSE
)
predictions <- as.data.frame(predictions, check.names = FALSE)
if (nrow(predictions) != ncol(query)) {
  stop("Transfer predictions do not align one-to-one with query cells.")
}

score_columns <- grep("^prediction\\.score\\.", names(predictions), value = TRUE)
score_columns <- setdiff(score_columns, "prediction.score.max")
if (length(score_columns) < 2L) stop("TransferData returned fewer than two class scores.")
score_matrix <- as.matrix(predictions[, score_columns, drop = FALSE])
storage.mode(score_matrix) <- "double"
class_labels <- sub("^prediction\\.score\\.", "", score_columns)
colnames(score_matrix) <- class_labels
if (any(!is.finite(score_matrix))) stop("TransferData returned non-finite class scores.")
predicted_ids <- trimws(as.character(predictions$predicted.id))
max_scores <- suppressWarnings(as.numeric(predictions$prediction.score.max))
if (anyNA(predicted_ids) || any(!nzchar(predicted_ids)) || any(!is.finite(max_scores))) {
  stop("TransferData returned invalid predicted labels or maximum scores.")
}

cluster_ids <- as.character(query$seurat_clusters)
cluster_levels <- unique(cluster_ids)
numeric_clusters <- suppressWarnings(as.numeric(cluster_levels))
cluster_levels <- cluster_levels[
  order(is.na(numeric_clusters), numeric_clusters, cluster_levels)
]
legacy_threshold <- as.numeric(config$label_transfer$prediction_score_threshold)
focus_classes <- c("RL", "Granule", "Glia", "VZ")
cluster_rows <- vector("list", length(cluster_levels))

for (cluster_index in seq_along(cluster_levels)) {
  cluster_id <- cluster_levels[[cluster_index]]
  cell_index <- which(cluster_ids == cluster_id)
  cluster_scores <- score_matrix[cell_index, , drop = FALSE]
  mean_scores <- colMeans(cluster_scores)
  score_order <- order(mean_scores, decreasing = TRUE)
  weighted_winner <- names(mean_scores)[score_order[[1]]]
  runner_up <- names(mean_scores)[score_order[[2]]]
  majority_counts <- sort(table(predicted_ids[cell_index]), decreasing = TRUE)
  majority_winner <- names(majority_counts)[[1]]
  mean_max_score <- mean(max_scores[cell_index])
  weighted_score <- unname(mean_scores[[weighted_winner]])
  score_or_na <- function(label) {
    if (label %in% names(mean_scores)) unname(mean_scores[[label]]) else NA_real_
  }

  row <- data.frame(
    sample_id = sample_id,
    reference = reference,
    k_anchor = k_anchor,
    k_score = k_score,
    k_weight = k_weight,
    seurat_cluster = cluster_id,
    n_cells = length(cell_index),
    cluster_majority = if (mean_max_score >= legacy_threshold) majority_winner else "Unknown",
    cluster_weighted = if (weighted_score >= legacy_threshold) weighted_winner else "Unknown",
    weighted_winner_before_threshold = weighted_winner,
    weighted_winner_mean_score = weighted_score,
    runner_up_mean_score_label = runner_up,
    runner_up_mean_score = unname(mean_scores[[runner_up]]),
    top2_mean_score_margin = weighted_score - unname(mean_scores[[runner_up]]),
    mean_prediction_score_max = mean_max_score,
    median_prediction_score_max = median(max_scores[cell_index]),
    fraction_cells_max_score_above_legacy_threshold = mean(
      max_scores[cell_index] > legacy_threshold
    ),
    legacy_threshold = legacy_threshold,
    stringsAsFactors = FALSE
  )
  for (class_label in focus_classes) {
    row[[paste0("mean_score_", class_label)]] <- score_or_na(class_label)
    row[[paste0("fraction_predicted_", class_label)]] <- mean(
      predicted_ids[cell_index] == class_label
    )
  }
  cluster_rows[[cluster_index]] <- row
}
cluster_summary <- do.call(rbind, cluster_rows)
if (anyDuplicated(cluster_summary$seurat_cluster)) {
  stop("Cluster summary unexpectedly duplicated seurat_clusters.")
}
if (sum(cluster_summary$n_cells) != ncol(query)) {
  stop("Cluster summary counts do not sum to the query cell count.")
}

prediction_count_matrix <- as.data.frame.matrix(table(cluster_ids, predicted_ids))
prediction_counts <- data.frame(
  seurat_cluster = rownames(prediction_count_matrix),
  prediction_count_matrix,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
count_cluster_number <- suppressWarnings(as.numeric(prediction_counts$seurat_cluster))
prediction_counts <- prediction_counts[
  order(is.na(count_cluster_number), count_cluster_number, prediction_counts$seurat_cluster),
  , drop = FALSE
]

anchor_count <- tryCatch(
  nrow(methods::slot(anchors, "anchors")),
  error = function(e) NA_integer_
)
provenance <- data.frame(
  field = c(
    "sample_id", "reference", "input_path", "reference_path", "random_seed",
    "query_cells", "balanced_reference_cells", "reference_cells_per_identity_max",
    "shared_genes", "transfer_features", "dimensions", "reduction",
    "k_anchor", "k_score", "k_weight", "prediction_score_threshold",
    "anchor_count", "outputs", "production_effect"
  ),
  value = c(
    sample_id, reference, input_path, reference_path,
    as.character(config$runtime$random_seed), as.character(ncol(query)),
    as.character(ncol(reference_balanced)),
    as.character(config$label_transfer$reference_cells_per_identity),
    as.character(length(shared_genes)), as.character(length(transfer_features)),
    paste(transfer_dimensions, collapse = ","), config$label_transfer$method,
    as.character(k_anchor), as.character(k_score), as.character(k_weight),
    as.character(legacy_threshold), as.character(anchor_count),
    "cluster score summary, prediction cell counts, provenance; no RDS",
    "none; isolated sensitivity tables only"
  ),
  stringsAsFactors = FALSE
)

output_dir <- dirname(output_paths[["cluster_scores"]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_summary, output_paths[["cluster_scores"]], row.names = FALSE, na = "")
write.csv(prediction_counts, output_paths[["prediction_counts"]], row.names = FALSE, na = "")
write.csv(provenance, output_paths[["provenance"]], row.names = FALSE, na = "")
message(
  "Wrote isolated k.anchor sensitivity tables for ", sample_id, " / ",
  reference, " / k.anchor=", k_anchor, "."
)
