#!/usr/bin/env Rscript

# Focused resolution-5 sensitivity analysis for two reference-sampling modes
# and four cluster-level abstention rules. The three query samples are derived
# from config/reviewed_annotation_clusters.csv and are mapped independently to
# Aldinger, Sepp, and Science. Outputs are compact tables only; no RDS object,
# production transfer output, consensus label, or threshold is modified.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
references <- c("Aldinger", "Sepp", "Science")
sampling_modes <- c("cap1000", "full")
legacy_threshold <- as.numeric(config$label_transfer$prediction_score_threshold)
reference_cap <- as.integer(config$label_transfer$reference_cells_per_identity)
k_anchor <- as.integer(config$label_transfer$k_anchor)
k_score <- as.integer(config$label_transfer$k_score)
k_weight <- 50L
transfer_dimensions <- seq_len(as.integer(config$label_transfer$dimensions))

if (!is.finite(legacy_threshold) || legacy_threshold <= 0 || legacy_threshold >= 1) {
  stop("label_transfer.prediction_score_threshold must be between 0 and 1.")
}
if (is.na(reference_cap) || reference_cap < 1L) {
  stop("label_transfer.reference_cells_per_identity must be a positive integer.")
}

# These are prespecified sensitivity rules, not production recommendations.
# production_score reproduces the current weighted per-reference vote.
abstention_rules <- data.frame(
  rule_id = c(
    "production_score",
    "score_margin_0.05",
    "score_margin_0.10",
    "score_margin_0.05_cell_support_0.60"
  ),
  min_winner_mean_score = rep(legacy_threshold, 4L),
  min_top2_mean_score_margin = c(0, 0.05, 0.10, 0.05),
  min_winner_cell_fraction = c(0, 0, 0, 0.60),
  stringsAsFactors = FALSE
)
if (anyDuplicated(abstention_rules$rule_id)) {
  stop("Abstention rule IDs must be unique.")
}

review_manifest_candidates <- c(
  here("config", "reviewed_annotation_clusters.csv"),
  here("scripts", "config", "reviewed_annotation_clusters.csv")
)
review_manifest_existing <- review_manifest_candidates[
  file.exists(review_manifest_candidates)
]
if (!length(review_manifest_existing)) {
  stop(
    "Missing reviewed-cluster manifest. Checked:\n- ",
    paste(review_manifest_candidates, collapse = "\n- ")
  )
}
review_manifest_path <- normalizePath(review_manifest_existing[[1]], mustWork = TRUE)
reviewed <- read.csv(
  review_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = character()
)
required_review_columns <- c(
  "sample_id", "seurat_cluster", "current_consensus_label",
  "expected_label", "review_source"
)
missing_review_columns <- setdiff(required_review_columns, names(reviewed))
if (length(missing_review_columns)) {
  stop(
    "Reviewed-cluster manifest lacks columns: ",
    paste(missing_review_columns, collapse = ", ")
  )
}
for (column in required_review_columns) {
  reviewed[[column]] <- trimws(as.character(reviewed[[column]]))
  if (anyNA(reviewed[[column]]) || any(!nzchar(reviewed[[column]]))) {
    stop("Reviewed-cluster manifest has blank or missing values in ", column, ".")
  }
}
review_keys <- paste(reviewed$sample_id, reviewed$seurat_cluster, sep = "|")
if (anyDuplicated(review_keys)) {
  stop("Reviewed-cluster manifest has duplicate sample/cluster rows.")
}

requested_samples <- unique(reviewed$sample_id)
if (length(requested_samples) != 3L) {
  stop(
    "This focused test expects exactly three unique reviewed samples; found ",
    length(requested_samples), "."
  )
}
sample_counts <- vapply(
  requested_samples,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(sample_counts != 1L)) {
  stop(
    "Every sensitivity sample must match exactly one config/samples.csv row. Invalid: ",
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
  sampling_mode = sampling_modes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_map <- task_map[order(
  match(task_map$sample_id, requested_samples),
  match(task_map$reference, references),
  match(task_map$sampling_mode, sampling_modes)
), , drop = FALSE]
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c("task_id", "sample_id", "reference", "sampling_mode")]
if (nrow(task_map) != 18L) stop("Expected exactly 18 sensitivity tasks.")

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
  cat("Scope: three reviewed resolution-5 samples x three references x two sampling modes\n")
  cat(
    "Fixed transfer settings: RPCA, dimensions=1:", max(transfer_dimensions),
    ", k.anchor=", k_anchor, ", k.score=", k_score,
    ", k.weight=", k_weight, ", seed=", config$runtime$random_seed, "\n",
    sep = ""
  )
  cat("Sampling modes: cap", reference_cap, " per identity; full reference\n", sep = "")
  cat("Abstention rules:\n")
  write.table(abstention_rules, row.names = FALSE, quote = FALSE, sep = "\t")
  cat("\nTask map:\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

output_root <- here(config$project$outputs_dir)
query_root <- file.path(
  output_root, "xenium", "preprocess", "04_resolution5_clustered", "rds"
)
sensitivity_root <- file.path(
  output_root, "validation", "resolution5_method_development",
  "reference_sampling_abstention_sensitivity"
)

query_path_for <- function(sample_id) {
  file.path(query_root, paste0(sample_id, "_whole_tissue_Res5.0.rds"))
}

task_output_paths <- function(sample_id, reference, sampling_mode) {
  output_dir <- file.path(
    sensitivity_root, sampling_mode, tolower(reference), "tables"
  )
  prefix <- paste(sample_id, reference, sampling_mode, sep = "_")
  c(
    cluster_metrics = file.path(
      output_dir, paste0(prefix, "_cluster_metrics.csv")
    ),
    reference_votes = file.path(
      output_dir, paste0(prefix, "_reference_votes.csv")
    ),
    reference_balance = file.path(
      output_dir, paste0(prefix, "_reference_class_balance.csv")
    ),
    provenance = file.path(
      output_dir, paste0(prefix, "_provenance.csv")
    )
  )
}

all_task_outputs <- do.call(rbind, lapply(seq_len(nrow(task_map)), function(i) {
  paths <- task_output_paths(
    task_map$sample_id[[i]],
    task_map$reference[[i]],
    task_map$sampling_mode[[i]]
  )
  data.frame(
    task_id = task_map$task_id[[i]],
    sample_id = task_map$sample_id[[i]],
    reference = task_map$reference[[i]],
    sampling_mode = task_map$sampling_mode[[i]],
    output_type = names(paths),
    output_path = unname(paths),
    stringsAsFactors = FALSE
  )
}))

combined_dir <- file.path(sensitivity_root, "combined_tables")
combined_outputs <- c(
  cluster_metrics = file.path(
    combined_dir, "sampling_abstention_cluster_metrics.csv"
  ),
  reference_votes = file.path(
    combined_dir, "sampling_abstention_reference_votes.csv"
  ),
  consensus = file.path(
    combined_dir, "sampling_abstention_consensus_results.csv"
  ),
  unknown_summary = file.path(
    combined_dir, "sampling_abstention_unknown_summary.csv"
  ),
  reviewed_clusters = file.path(
    combined_dir, "sampling_abstention_reviewed_cluster_results.csv"
  ),
  sampling_comparison = file.path(
    combined_dir, "sampling_abstention_sampling_comparison.csv"
  ),
  reference_balance = file.path(
    combined_dir, "sampling_abstention_reference_class_balance.csv"
  ),
  provenance = file.path(
    combined_dir, "sampling_abstention_combined_provenance.csv"
  )
)

read_task_outputs <- function(output_type) {
  paths <- all_task_outputs$output_path[
    all_task_outputs$output_type == output_type
  ]
  tables <- lapply(paths, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  do.call(rbind, tables)
}

consensus_from_votes <- function(group) {
  if (nrow(group) != length(references) ||
      !setequal(as.character(group$reference), references)) {
    stop(
      "Consensus group does not contain exactly one vote from each reference: ",
      paste(unique(group$sample_id), unique(group$sampling_mode),
            unique(group$rule_id), unique(group$seurat_cluster), collapse = " | ")
    )
  }
  group <- group[match(references, group$reference), , drop = FALSE]
  votes <- trimws(as.character(group$reference_vote))
  valid_votes <- votes[
    !is.na(votes) & nzchar(votes) & tolower(votes) != "unknown"
  ]
  if (!length(valid_votes)) {
    label <- "Unknown"
    support <- 0L
    status <- "insufficient_support"
  } else {
    counts <- table(valid_votes)
    support <- as.integer(max(counts))
    winners <- names(counts)[counts == support]
    if (support < 2L || length(winners) != 1L) {
      label <- "Unknown"
      status <- "insufficient_support"
    } else {
      label <- winners[[1]]
      status <- if (support == 3L) "unanimous_3of3" else "majority_2of3"
    }
  }
  n_cells <- unique(as.integer(group$n_cells))
  if (length(n_cells) != 1L || is.na(n_cells)) {
    stop("References disagree on query cluster cell count.")
  }
  data.frame(
    sample_id = group$sample_id[[1]],
    sampling_mode = group$sampling_mode[[1]],
    rule_id = group$rule_id[[1]],
    seurat_cluster = as.character(group$seurat_cluster[[1]]),
    n_cells = n_cells,
    aldinger_vote = votes[[1]],
    sepp_vote = votes[[2]],
    science_vote = votes[[3]],
    consensus_label = label,
    consensus_support_n = support,
    consensus_nonunknown_n = length(valid_votes),
    consensus_status = status,
    stringsAsFactors = FALSE
  )
}

if (combine) {
  if (length(positional_args)) stop("--combine does not accept TASK_ID.")
  required_inputs <- all_task_outputs$output_path
  if (dry_run) {
    ok <- compact_dry_run(
      "Combine reference-sampling and abstention sensitivity outputs",
      inputs = required_inputs,
      outputs = combined_outputs,
      checks = c(
        expected_task_outputs = length(required_inputs) == 18L * 4L,
        expected_rules = nrow(abstention_rules) == 4L,
        expected_samples = length(requested_samples) == 3L
      )
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }

  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop(
      "Missing sensitivity outputs:\n- ",
      paste(missing_inputs, collapse = "\n- ")
    )
  }
  existing_outputs <- combined_outputs[file.exists(combined_outputs)]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite combined sensitivity outputs:\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing the existing tables."
    )
  }

  cluster_metrics <- read_task_outputs("cluster_metrics")
  reference_votes <- read_task_outputs("reference_votes")
  reference_balance <- read_task_outputs("reference_balance")
  cluster_metrics$seurat_cluster <- as.character(cluster_metrics$seurat_cluster)
  reference_votes$seurat_cluster <- as.character(reference_votes$seurat_cluster)

  metric_key <- paste(
    cluster_metrics$sample_id, cluster_metrics$reference,
    cluster_metrics$sampling_mode, cluster_metrics$seurat_cluster,
    sep = "|"
  )
  if (anyDuplicated(metric_key)) {
    stop("Combined cluster metrics contain duplicate task/cluster rows.")
  }
  vote_key <- paste(
    reference_votes$sample_id, reference_votes$reference,
    reference_votes$sampling_mode, reference_votes$rule_id,
    reference_votes$seurat_cluster, sep = "|"
  )
  if (anyDuplicated(vote_key)) {
    stop("Combined reference votes contain duplicate task/rule/cluster rows.")
  }

  consensus_group_key <- paste(
    reference_votes$sample_id, reference_votes$sampling_mode,
    reference_votes$rule_id, reference_votes$seurat_cluster,
    sep = "|"
  )
  consensus_groups <- split(reference_votes, consensus_group_key, drop = TRUE)
  consensus <- do.call(rbind, lapply(consensus_groups, consensus_from_votes))
  rownames(consensus) <- NULL

  sample_order <- match(consensus$sample_id, requested_samples)
  sampling_order <- match(consensus$sampling_mode, sampling_modes)
  rule_order <- match(consensus$rule_id, abstention_rules$rule_id)
  numeric_cluster <- suppressWarnings(as.numeric(consensus$seurat_cluster))
  consensus <- consensus[order(
    sample_order, sampling_order, rule_order,
    is.na(numeric_cluster), numeric_cluster, consensus$seurat_cluster
  ), , drop = FALSE]

  summary_group_key <- paste(
    consensus$sample_id, consensus$sampling_mode, consensus$rule_id,
    sep = "|"
  )
  unknown_summary <- do.call(rbind, lapply(
    split(consensus, summary_group_key, drop = TRUE),
    function(group) {
      unknown <- tolower(group$consensus_label) == "unknown"
      data.frame(
        sample_id = group$sample_id[[1]],
        sampling_mode = group$sampling_mode[[1]],
        rule_id = group$rule_id[[1]],
        n_clusters = nrow(group),
        n_known_clusters = sum(!unknown),
        n_unknown_clusters = sum(unknown),
        fraction_unknown_clusters = mean(unknown),
        n_cells = sum(group$n_cells),
        n_unknown_cells = sum(group$n_cells[unknown]),
        fraction_cells_in_unknown_clusters = sum(group$n_cells[unknown]) / sum(group$n_cells),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(unknown_summary) <- NULL
  unknown_summary <- unknown_summary[order(
    match(unknown_summary$sample_id, requested_samples),
    match(unknown_summary$sampling_mode, sampling_modes),
    match(unknown_summary$rule_id, abstention_rules$rule_id)
  ), , drop = FALSE]

  reviewed_clusters <- merge(
    consensus,
    reviewed[required_review_columns],
    by = c("sample_id", "seurat_cluster"),
    all = FALSE,
    sort = FALSE
  )
  reviewed_clusters$matches_expected_label <- (
    reviewed_clusters$consensus_label == reviewed_clusters$expected_label
  )
  reviewed_clusters$matches_current_label <- (
    reviewed_clusters$consensus_label == reviewed_clusters$current_consensus_label
  )
  reviewed_clusters <- reviewed_clusters[order(
    match(reviewed_clusters$sample_id, requested_samples),
    match(
      paste(reviewed_clusters$sample_id, reviewed_clusters$seurat_cluster, sep = "|"),
      review_keys
    ),
    match(reviewed_clusters$sampling_mode, sampling_modes),
    match(reviewed_clusters$rule_id, abstention_rules$rule_id)
  ), , drop = FALSE]

  cap_results <- consensus[consensus$sampling_mode == "cap1000", , drop = FALSE]
  full_results <- consensus[consensus$sampling_mode == "full", , drop = FALSE]
  cap_results <- cap_results[c(
    "sample_id", "rule_id", "seurat_cluster", "consensus_label",
    "consensus_support_n", "consensus_nonunknown_n"
  )]
  full_results <- full_results[c(
    "sample_id", "rule_id", "seurat_cluster", "consensus_label",
    "consensus_support_n", "consensus_nonunknown_n"
  )]
  names(cap_results)[4:6] <- paste0("cap1000_", names(cap_results)[4:6])
  names(full_results)[4:6] <- paste0("full_", names(full_results)[4:6])
  sampling_comparison <- merge(
    cap_results, full_results,
    by = c("sample_id", "rule_id", "seurat_cluster"),
    all = TRUE, sort = FALSE
  )
  if (anyNA(sampling_comparison$cap1000_consensus_label) ||
      anyNA(sampling_comparison$full_consensus_label)) {
    stop("Cap1000 and full-reference consensus cluster sets do not match.")
  }
  sampling_comparison$consensus_label_changed <- (
    sampling_comparison$cap1000_consensus_label !=
      sampling_comparison$full_consensus_label
  )
  sampling_comparison <- sampling_comparison[order(
    match(sampling_comparison$sample_id, requested_samples),
    match(sampling_comparison$rule_id, abstention_rules$rule_id),
    suppressWarnings(as.numeric(sampling_comparison$seurat_cluster)),
    sampling_comparison$seurat_cluster
  ), , drop = FALSE]

  combined_provenance <- data.frame(
    field = c(
      "purpose", "review_manifest", "samples", "references",
      "sampling_modes", "abstention_rules", "transfer_settings",
      "task_count", "production_effect"
    ),
    value = c(
      "focused resolution-5 reference-sampling and abstention sensitivity",
      review_manifest_path,
      paste(requested_samples, collapse = "|"),
      paste(references, collapse = "|"),
      paste(sampling_modes, collapse = "|"),
      paste(abstention_rules$rule_id, collapse = "|"),
      paste0(
        "RPCA;dims=1:", max(transfer_dimensions),
        ";k.anchor=", k_anchor, ";k.score=", k_score,
        ";k.weight=", k_weight, ";seed=", config$runtime$random_seed
      ),
      as.character(nrow(task_map)),
      "none; compact isolated tables only"
    ),
    stringsAsFactors = FALSE
  )

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(cluster_metrics, combined_outputs[["cluster_metrics"]], row.names = FALSE, na = "")
  write.csv(reference_votes, combined_outputs[["reference_votes"]], row.names = FALSE, na = "")
  write.csv(consensus, combined_outputs[["consensus"]], row.names = FALSE, na = "")
  write.csv(unknown_summary, combined_outputs[["unknown_summary"]], row.names = FALSE, na = "")
  write.csv(reviewed_clusters, combined_outputs[["reviewed_clusters"]], row.names = FALSE, na = "")
  write.csv(sampling_comparison, combined_outputs[["sampling_comparison"]], row.names = FALSE, na = "")
  write.csv(reference_balance, combined_outputs[["reference_balance"]], row.names = FALSE, na = "")
  write.csv(combined_provenance, combined_outputs[["provenance"]], row.names = FALSE, na = "")
  message("Wrote combined reference-sampling and abstention sensitivity tables: ", combined_dir)
  quit(save = "no", status = 0L)
}

if (length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]\n",
    "   or: Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R ",
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
reference_name <- task_map$reference[[task_id]]
sampling_mode <- task_map$sampling_mode[[task_id]]
query_path <- query_path_for(sample_id)
reference_path <- unname(reference_paths[[reference_name]])
output_paths <- task_output_paths(sample_id, reference_name, sampling_mode)

if (dry_run) {
  ok <- compact_dry_run(
    paste0(
      "Reference sampling/abstention task ", task_id, "/", nrow(task_map),
      " [", sample_id, " | ", reference_name, " | ", sampling_mode, "]"
    ),
    inputs = c(query_path, reference_path, review_manifest_path),
    outputs = output_paths,
    checks = c(
      sample_matches_manifest = sum(sample_manifest$sample_id == sample_id) == 1L,
      valid_reference = reference_name %in% references,
      valid_sampling_mode = sampling_mode %in% sampling_modes,
      expected_rules = nrow(abstention_rules) == 4L
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- c(query_path, reference_path)[
  !file.exists(c(query_path, reference_path))
]
if (length(missing_inputs)) {
  stop("Missing sensitivity input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite sensitivity outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing tables."
  )
}

set.seed(as.integer(config$runtime$random_seed))
query <- readRDS(query_path)
reference <- readRDS(reference_path)

if (anyDuplicated(Cells(query))) stop("Query contains duplicate cell IDs: ", sample_id)
if (anyDuplicated(Cells(reference))) {
  stop("Reference contains duplicate cell IDs: ", reference_name)
}
if (!"clusters_refined" %in% colnames(reference[[]])) {
  stop(reference_name, " reference lacks clusters_refined metadata.")
}
reference_labels <- trimws(as.character(reference$clusters_refined))
if (anyNA(reference_labels) || any(!nzchar(reference_labels))) {
  stop(reference_name, " reference has blank or missing clusters_refined labels.")
}
reference$clusters_refined <- reference_labels

required_query_metadata <- c(
  "whole_tissue_cluster_res5.0", "Xenium_snn_res.5", "seurat_clusters"
)
missing_query_metadata <- setdiff(required_query_metadata, colnames(query[[]]))
if (length(missing_query_metadata)) {
  stop(
    sample_id, " resolution-5 input lacks metadata: ",
    paste(missing_query_metadata, collapse = ", ")
  )
}
cluster_ids <- as.character(query$seurat_clusters)
resolution_clusters <- as.character(query$whole_tissue_cluster_res5.0)
if (anyNA(cluster_ids) || !identical(cluster_ids, resolution_clusters)) {
  stop("Query seurat_clusters does not exactly match whole_tissue_cluster_res5.0.")
}
sample_reviewed_clusters <- reviewed$seurat_cluster[reviewed$sample_id == sample_id]
missing_reviewed_clusters <- setdiff(sample_reviewed_clusters, unique(cluster_ids))
if (length(missing_reviewed_clusters)) {
  stop(
    sample_id, " lacks reviewed resolution-5 cluster(s): ",
    paste(missing_reviewed_clusters, collapse = ", ")
  )
}

DefaultAssay(reference) <- if ("RNA" %in% Assays(reference)) {
  "RNA"
} else {
  Assays(reference)[[1]]
}
DefaultAssay(query) <- if ("Xenium" %in% Assays(query)) {
  "Xenium"
} else {
  Assays(query)[[1]]
}

reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures = as.integer(config$label_transfer$variable_features),
  verbose = FALSE
)
shared_genes <- intersect(rownames(reference), rownames(query))
Idents(reference) <- "clusters_refined"
set.seed(as.integer(config$runtime$random_seed))
reference_for_transfer <- if (sampling_mode == "cap1000") {
  subset(reference, downsample = reference_cap)
} else if (sampling_mode == "full") {
  reference
} else {
  stop("Unsupported sampling mode: ", sampling_mode)
}

class_counts <- function(object, count_name) {
  counts <- as.data.frame(
    table(reference_class = as.character(object$clusters_refined)),
    stringsAsFactors = FALSE
  )
  names(counts)[[2]] <- count_name
  counts
}
reference_balance <- merge(
  class_counts(reference, "full_reference_cells"),
  class_counts(reference_for_transfer, "used_reference_cells"),
  by = "reference_class", all = TRUE, sort = FALSE
)
reference_balance$full_reference_cells[is.na(reference_balance$full_reference_cells)] <- 0L
reference_balance$used_reference_cells[is.na(reference_balance$used_reference_cells)] <- 0L
reference_balance$retained_fraction <- (
  reference_balance$used_reference_cells / reference_balance$full_reference_cells
)
reference_balance$sample_id <- sample_id
reference_balance$reference <- reference_name
reference_balance$sampling_mode <- sampling_mode
reference_balance$reference_cap_per_identity <- if (sampling_mode == "cap1000") {
  reference_cap
} else {
  NA_integer_
}
reference_balance <- reference_balance[c(
  "sample_id", "reference", "sampling_mode", "reference_class",
  "full_reference_cells", "used_reference_cells", "retained_fraction",
  "reference_cap_per_identity"
)]
full_reference_cells_total <- ncol(reference)
if (sampling_mode == "cap1000" &&
    any(reference_balance$used_reference_cells > reference_cap)) {
  stop("Capped reference retained more than ", reference_cap, " cells in an identity.")
}
if (sampling_mode == "full" &&
    any(reference_balance$used_reference_cells != reference_balance$full_reference_cells)) {
  stop("Full-reference mode did not retain every reference cell.")
}

# Release the second reference binding before ScaleData/RunPCA. This avoids
# retaining both the complete input and its scaled copy in full-reference mode.
rm(reference)
invisible(gc(verbose = FALSE))

transfer_features <- intersect(VariableFeatures(reference_for_transfer), shared_genes)
if (length(transfer_features) <= max(transfer_dimensions)) {
  stop(
    "Only ", length(transfer_features), " transfer features are available; more than ",
    max(transfer_dimensions), " are required."
  )
}

reference_for_transfer <- ScaleData(
  reference_for_transfer, features = transfer_features, verbose = FALSE
)
query <- ScaleData(query, features = transfer_features, verbose = FALSE)
set.seed(as.integer(config$runtime$random_seed))
reference_for_transfer <- RunPCA(
  reference_for_transfer, features = transfer_features, verbose = FALSE
)
set.seed(as.integer(config$runtime$random_seed))
query <- RunPCA(query, features = transfer_features, verbose = FALSE)

set.seed(as.integer(config$runtime$random_seed))
anchors <- FindTransferAnchors(
  reference = reference_for_transfer,
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
  refdata = reference_for_transfer$clusters_refined,
  dims = transfer_dimensions,
  k.weight = k_weight,
  store.weights = FALSE
)
predictions <- as.data.frame(predictions, check.names = FALSE)
if (!identical(rownames(predictions), Cells(query))) {
  if (!setequal(rownames(predictions), Cells(query))) {
    stop("Transfer predictions do not align one-to-one with query cells.")
  }
  predictions <- predictions[Cells(query), , drop = FALSE]
}

score_columns <- grep("^prediction\\.score\\.", names(predictions), value = TRUE)
score_columns <- setdiff(
  score_columns, c("prediction.score.max", "prediction.score.id")
)
if (length(score_columns) < 2L) {
  stop("TransferData returned fewer than two class-score columns.")
}
score_matrix <- as.matrix(predictions[, score_columns, drop = FALSE])
storage.mode(score_matrix) <- "double"
class_labels <- sub("^prediction\\.score\\.", "", score_columns)
if (anyDuplicated(class_labels)) stop("TransferData returned duplicate class labels.")
colnames(score_matrix) <- class_labels
if (any(!is.finite(score_matrix))) stop("TransferData returned non-finite class scores.")
predicted_ids <- trimws(as.character(predictions$predicted.id))
max_scores <- suppressWarnings(as.numeric(predictions$prediction.score.max))
if (anyNA(predicted_ids) || any(!nzchar(predicted_ids)) ||
    any(!is.finite(max_scores))) {
  stop("TransferData returned invalid predicted labels or maximum scores.")
}

cluster_levels <- unique(cluster_ids)
numeric_clusters <- suppressWarnings(as.numeric(cluster_levels))
cluster_levels <- cluster_levels[order(
  is.na(numeric_clusters), numeric_clusters, cluster_levels
)]
focus_classes <- c("RL", "Granule", "Glia", "VZ")
cluster_rows <- vector("list", length(cluster_levels))
vote_rows <- vector("list", length(cluster_levels) * nrow(abstention_rules))
vote_index <- 1L

for (cluster_index in seq_along(cluster_levels)) {
  cluster_id <- cluster_levels[[cluster_index]]
  cell_index <- which(cluster_ids == cluster_id)
  cluster_scores <- score_matrix[cell_index, , drop = FALSE]
  mean_scores <- colMeans(cluster_scores)
  score_order <- order(mean_scores, decreasing = TRUE)
  winner <- names(mean_scores)[score_order[[1]]]
  runner_up <- names(mean_scores)[score_order[[2]]]
  winner_score <- unname(mean_scores[[winner]])
  runner_up_score <- unname(mean_scores[[runner_up]])
  top2_margin <- winner_score - runner_up_score
  winner_cell_fraction <- mean(predicted_ids[cell_index] == winner)
  majority_counts <- sort(table(predicted_ids[cell_index]), decreasing = TRUE)
  majority_winner <- names(majority_counts)[[1]]
  score_or_na <- function(label) {
    if (label %in% names(mean_scores)) unname(mean_scores[[label]]) else NA_real_
  }

  metric_row <- data.frame(
    sample_id = sample_id,
    reference = reference_name,
    sampling_mode = sampling_mode,
    seurat_cluster = cluster_id,
    n_cells = length(cell_index),
    weighted_winner_before_abstention = winner,
    weighted_winner_mean_score = winner_score,
    runner_up_mean_score_label = runner_up,
    runner_up_mean_score = runner_up_score,
    top2_mean_score_margin = top2_margin,
    winner_cell_fraction = winner_cell_fraction,
    majority_winner_before_abstention = majority_winner,
    mean_prediction_score_max = mean(max_scores[cell_index]),
    median_prediction_score_max = median(max_scores[cell_index]),
    fraction_cells_max_score_at_least_production_threshold = mean(
      max_scores[cell_index] >= legacy_threshold
    ),
    production_threshold = legacy_threshold,
    stringsAsFactors = FALSE
  )
  for (class_label in focus_classes) {
    metric_row[[paste0("mean_score_", class_label)]] <- score_or_na(class_label)
    metric_row[[paste0("fraction_predicted_", class_label)]] <- mean(
      predicted_ids[cell_index] == class_label
    )
  }
  cluster_rows[[cluster_index]] <- metric_row

  for (rule_index in seq_len(nrow(abstention_rules))) {
    rule <- abstention_rules[rule_index, , drop = FALSE]
    failed_score <- winner_score < rule$min_winner_mean_score[[1]]
    failed_margin <- top2_margin < rule$min_top2_mean_score_margin[[1]]
    failed_cell_support <- (
      winner_cell_fraction < rule$min_winner_cell_fraction[[1]]
    )
    abstained <- failed_score || failed_margin || failed_cell_support
    vote_rows[[vote_index]] <- data.frame(
      sample_id = sample_id,
      reference = reference_name,
      sampling_mode = sampling_mode,
      rule_id = rule$rule_id[[1]],
      seurat_cluster = cluster_id,
      n_cells = length(cell_index),
      weighted_winner_before_abstention = winner,
      reference_vote = if (abstained) "Unknown" else winner,
      abstained = abstained,
      failed_score = failed_score,
      failed_margin = failed_margin,
      failed_cell_support = failed_cell_support,
      winner_mean_score = winner_score,
      top2_mean_score_margin = top2_margin,
      winner_cell_fraction = winner_cell_fraction,
      min_winner_mean_score = rule$min_winner_mean_score[[1]],
      min_top2_mean_score_margin = rule$min_top2_mean_score_margin[[1]],
      min_winner_cell_fraction = rule$min_winner_cell_fraction[[1]],
      stringsAsFactors = FALSE
    )
    vote_index <- vote_index + 1L
  }
}

cluster_metrics <- do.call(rbind, cluster_rows)
reference_votes <- do.call(rbind, vote_rows)
if (anyDuplicated(cluster_metrics$seurat_cluster)) {
  stop("Cluster metrics unexpectedly duplicated seurat_clusters.")
}
if (sum(cluster_metrics$n_cells) != ncol(query)) {
  stop("Cluster metric cell counts do not sum to the query cell count.")
}
if (nrow(reference_votes) != nrow(cluster_metrics) * nrow(abstention_rules)) {
  stop("Reference-vote table has an unexpected number of rows.")
}

production_votes <- reference_votes[
  reference_votes$rule_id == "production_score", , drop = FALSE
]
expected_production_votes <- ifelse(
  cluster_metrics$weighted_winner_mean_score >= legacy_threshold,
  cluster_metrics$weighted_winner_before_abstention,
  "Unknown"
)
if (!identical(
  as.character(production_votes$reference_vote),
  as.character(expected_production_votes)
)) {
  stop("production_score rule does not reproduce the current weighted threshold rule.")
}

anchor_count <- tryCatch(
  nrow(methods::slot(anchors, "anchors")),
  error = function(e) NA_integer_
)
provenance <- data.frame(
  field = c(
    "sample_id", "reference", "sampling_mode", "query_path",
    "reference_path", "review_manifest", "random_seed", "query_cells",
    "full_reference_cells", "used_reference_cells", "reference_cap",
    "shared_genes", "transfer_features", "dimensions", "reduction",
    "k_anchor", "k_score", "k_weight", "production_threshold",
    "abstention_rules", "seed_control", "anchor_count", "outputs",
    "production_effect"
  ),
  value = c(
    sample_id, reference_name, sampling_mode, query_path,
    reference_path, review_manifest_path, as.character(config$runtime$random_seed),
    as.character(ncol(query)), as.character(full_reference_cells_total),
    as.character(ncol(reference_for_transfer)),
    if (sampling_mode == "cap1000") as.character(reference_cap) else "none",
    as.character(length(shared_genes)), as.character(length(transfer_features)),
    paste(transfer_dimensions, collapse = ","), config$label_transfer$method,
    as.character(k_anchor), as.character(k_score), as.character(k_weight),
    as.character(legacy_threshold),
    paste(abstention_rules$rule_id, collapse = "|"),
    "seed reset before sampling, each PCA, and anchor finding to isolate sampling mode",
    as.character(anchor_count),
    "cluster metrics, reference votes, reference balance, provenance; no RDS",
    "none; isolated sensitivity tables only"
  ),
  stringsAsFactors = FALSE
)

output_dir <- dirname(output_paths[["cluster_metrics"]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_metrics, output_paths[["cluster_metrics"]], row.names = FALSE, na = "")
write.csv(reference_votes, output_paths[["reference_votes"]], row.names = FALSE, na = "")
write.csv(reference_balance, output_paths[["reference_balance"]], row.names = FALSE, na = "")
write.csv(provenance, output_paths[["provenance"]], row.names = FALSE, na = "")
message(
  "Wrote isolated sampling/abstention tables for ", sample_id, " / ",
  reference_name, " / ", sampling_mode, "."
)
