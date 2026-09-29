#!/usr/bin/env Rscript

# Focused resolution-5 cap-convergence screen. Test intermediate per-identity
# caps of 2,500, 5,000, and 10,000 cells on the three reviewed query samples
# and three published references. Reuse the completed cap1000 and full results
# when combining. Outputs are compact tables only; production is untouched.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
references <- c("Aldinger", "Sepp", "Science")
tested_caps <- c(2500L, 5000L, 10000L)
sampling_modes <- paste0("cap", tested_caps)
cap_by_mode <- setNames(tested_caps, sampling_modes)
comparison_modes <- c("cap1000", sampling_modes, "full")
random_seed <- as.integer(config$runtime$random_seed)
legacy_threshold <- as.numeric(
  config$label_transfer$prediction_score_threshold
)
k_anchor <- as.integer(config$label_transfer$k_anchor)
k_score <- as.integer(config$label_transfer$k_score)
k_weight <- 50L
transfer_dimensions <- seq_len(
  as.integer(config$label_transfer$dimensions)
)

if (!is.finite(legacy_threshold) ||
    legacy_threshold <= 0 || legacy_threshold >= 1) {
  stop("label_transfer.prediction_score_threshold must be between 0 and 1.")
}
if (anyNA(tested_caps) || any(tested_caps < 1L) ||
    anyDuplicated(tested_caps)) {
  stop("Tested reference caps must be unique positive integers.")
}
if (is.na(random_seed)) stop("runtime.random_seed must be an integer.")

# Preserve all four prior rules so cap selection and abstention can be reviewed
# independently. score_margin_0.05 remains the leading candidate, not a fixed
# production decision.
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
review_manifest_path <- normalizePath(
  review_manifest_existing[[1]], mustWork = TRUE
)
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
    "This focused test expects exactly three reviewed samples; found ",
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
    "Every cap-convergence sample must match exactly one samples.csv row: ",
    paste(requested_samples[sample_counts != 1L], collapse = ", ")
  )
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

# Ordering makes each group of nine tasks one cap, which is easy to audit from
# Slurm task IDs.
task_map <- expand.grid(
  sample_id = requested_samples,
  reference = references,
  sampling_mode = sampling_modes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_map <- task_map[order(
  match(task_map$sampling_mode, sampling_modes),
  match(task_map$sample_id, requested_samples),
  match(task_map$reference, references)
), , drop = FALSE]
task_map$reference_cap_per_identity <- unname(
  cap_by_mode[task_map$sampling_mode]
)
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c(
  "task_id", "sample_id", "reference", "sampling_mode",
  "reference_cap_per_identity"
)]
if (nrow(task_map) != 27L ||
    !identical(task_map$task_id, seq_len(27L))) {
  stop("Expected exactly 27 consecutive cap-convergence tasks.")
}

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
  cat("Scope: three reviewed samples x three references x three intermediate caps\n")
  cat("New caps: ", paste(tested_caps, collapse = ", "), " cells per identity\n", sep = "")
  cat("Combine modes: ", paste(comparison_modes, collapse = ", "), "\n", sep = "")
  cat(
    "Fixed transfer settings: RPCA, dimensions=1:",
    max(transfer_dimensions), ", k.anchor=", k_anchor,
    ", k.score=", k_score, ", k.weight=", k_weight,
    ", seed=", random_seed, "\n", sep = ""
  )
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

task_output_paths <- function(sample_id, reference_name, sampling_mode) {
  output_dir <- file.path(
    sensitivity_root, sampling_mode, tolower(reference_name), "tables"
  )
  prefix <- paste(sample_id, reference_name, sampling_mode, sep = "_")
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

existing_combined_dir <- file.path(sensitivity_root, "combined_tables")
existing_combined_inputs <- c(
  cluster_metrics = file.path(
    existing_combined_dir, "sampling_abstention_cluster_metrics.csv"
  ),
  reference_votes = file.path(
    existing_combined_dir, "sampling_abstention_reference_votes.csv"
  ),
  reference_balance = file.path(
    existing_combined_dir, "sampling_abstention_reference_class_balance.csv"
  )
)

cap_combined_dir <- file.path(sensitivity_root, "cap_convergence_screen_tables")
combined_outputs <- c(
  cluster_metrics = file.path(cap_combined_dir, "cap_convergence_cluster_metrics.csv"),
  reference_votes = file.path(cap_combined_dir, "cap_convergence_reference_votes.csv"),
  consensus = file.path(cap_combined_dir, "cap_convergence_consensus_results.csv"),
  unknown_summary = file.path(cap_combined_dir, "cap_convergence_unknown_summary.csv"),
  reviewed_clusters = file.path(cap_combined_dir, "cap_convergence_reviewed_cluster_results.csv"),
  cluster_comparison = file.path(cap_combined_dir, "cap_convergence_cluster_comparison.csv"),
  convergence_summary = file.path(cap_combined_dir, "cap_convergence_to_full_summary.csv"),
  reference_vote_summary = file.path(cap_combined_dir, "cap_convergence_reference_vote_summary.csv"),
  reference_balance = file.path(cap_combined_dir, "cap_convergence_reference_class_balance.csv"),
  provenance = file.path(cap_combined_dir, "cap_convergence_combined_provenance.csv")
)

read_new_task_outputs <- function(output_type) {
  paths <- all_task_outputs$output_path[
    all_task_outputs$output_type == output_type
  ]
  do.call(rbind, lapply(paths, function(path) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

bind_same_columns <- function(first, second, first_name, second_name) {
  if (!setequal(names(first), names(second))) {
    stop(
      "Cannot combine ", first_name, " and ", second_name,
      "; column sets differ. Only in ", first_name, ": ",
      paste(setdiff(names(first), names(second)), collapse = ", "),
      "; only in ", second_name, ": ",
      paste(setdiff(names(second), names(first)), collapse = ", ")
    )
  }
  second <- second[names(first)]
  rbind(first, second)
}

consensus_from_votes <- function(group) {
  if (nrow(group) != length(references) ||
      !setequal(as.character(group$reference), references)) {
    stop("Consensus group must contain exactly one vote from each reference.")
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
  required_inputs <- c(
    existing_combined_inputs,
    all_task_outputs$output_path
  )
  if (dry_run) {
    ok <- compact_dry_run(
      "Combine five-mode reference-cap convergence screen",
      inputs = required_inputs,
      outputs = combined_outputs,
      checks = c(
        expected_new_task_outputs = nrow(all_task_outputs) == 27L * 4L,
        expected_comparison_modes = length(comparison_modes) == 5L,
        expected_rules = nrow(abstention_rules) == 4L,
        expected_samples = length(requested_samples) == 3L
      )
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }

  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop(
      "Missing cap-convergence combine input(s):\n- ",
      paste(missing_inputs, collapse = "\n- ")
    )
  }
  existing_outputs <- combined_outputs[file.exists(combined_outputs)]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite cap-convergence combined outputs:\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing the existing tables."
    )
  }

  old_metrics <- read.csv(
    existing_combined_inputs[["cluster_metrics"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  old_votes <- read.csv(
    existing_combined_inputs[["reference_votes"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  old_balance <- read.csv(
    existing_combined_inputs[["reference_balance"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  old_metrics <- old_metrics[
    old_metrics$sampling_mode %in% c("cap1000", "full"), , drop = FALSE
  ]
  old_votes <- old_votes[
    old_votes$sampling_mode %in% c("cap1000", "full"), , drop = FALSE
  ]
  old_balance <- old_balance[
    old_balance$sampling_mode %in% c("cap1000", "full"), , drop = FALSE
  ]

  new_metrics <- read_new_task_outputs("cluster_metrics")
  new_votes <- read_new_task_outputs("reference_votes")
  new_balance <- read_new_task_outputs("reference_balance")
  cluster_metrics <- bind_same_columns(
    old_metrics, new_metrics,
    "completed cap1000/full cluster metrics", "intermediate-cap cluster metrics"
  )
  reference_votes <- bind_same_columns(
    old_votes, new_votes,
    "completed cap1000/full votes", "intermediate-cap votes"
  )
  reference_balance <- bind_same_columns(
    old_balance, new_balance,
    "completed cap1000/full balance", "intermediate-cap balance"
  )
  cluster_metrics$seurat_cluster <- as.character(
    cluster_metrics$seurat_cluster
  )
  reference_votes$seurat_cluster <- as.character(
    reference_votes$seurat_cluster
  )

  observed_modes <- sort(unique(as.character(cluster_metrics$sampling_mode)))
  if (!setequal(observed_modes, comparison_modes)) {
    stop(
      "Combined metrics do not contain exactly the five expected modes. Found: ",
      paste(observed_modes, collapse = ", ")
    )
  }
  metric_key <- paste(
    cluster_metrics$sample_id, cluster_metrics$reference,
    cluster_metrics$sampling_mode, cluster_metrics$seurat_cluster,
    sep = "|"
  )
  if (anyDuplicated(metric_key)) {
    stop("Combined cap-convergence metrics contain duplicate task/cluster rows.")
  }
  vote_key <- paste(
    reference_votes$sample_id, reference_votes$reference,
    reference_votes$sampling_mode, reference_votes$rule_id,
    reference_votes$seurat_cluster, sep = "|"
  )
  if (anyDuplicated(vote_key)) {
    stop("Combined cap-convergence votes contain duplicate task/rule/cluster rows.")
  }

  consensus_group_key <- paste(
    reference_votes$sample_id, reference_votes$sampling_mode,
    reference_votes$rule_id, reference_votes$seurat_cluster,
    sep = "|"
  )
  consensus <- do.call(rbind, lapply(
    split(reference_votes, consensus_group_key, drop = TRUE),
    consensus_from_votes
  ))
  rownames(consensus) <- NULL
  consensus$seurat_cluster <- as.character(consensus$seurat_cluster)
  numeric_cluster <- suppressWarnings(as.numeric(consensus$seurat_cluster))
  consensus <- consensus[order(
    match(consensus$sample_id, requested_samples),
    match(consensus$sampling_mode, comparison_modes),
    match(consensus$rule_id, abstention_rules$rule_id),
    is.na(numeric_cluster), numeric_cluster, consensus$seurat_cluster
  ), , drop = FALSE]

  full_cluster_keys <- paste(
    consensus$sample_id[consensus$sampling_mode == "full"],
    consensus$rule_id[consensus$sampling_mode == "full"],
    consensus$seurat_cluster[consensus$sampling_mode == "full"],
    sep = "|"
  )
  mode_cluster_sets <- split(
    paste(consensus$sample_id, consensus$rule_id,
          consensus$seurat_cluster, sep = "|"),
    consensus$sampling_mode
  )
  if (!all(vapply(
    mode_cluster_sets,
    function(keys) setequal(keys, full_cluster_keys),
    logical(1)
  ))) {
    stop("The five sampling modes do not contain identical cluster/rule sets.")
  }

  unknown_group_key <- paste(
    consensus$sample_id, consensus$sampling_mode, consensus$rule_id,
    sep = "|"
  )
  unknown_summary <- do.call(rbind, lapply(
    split(consensus, unknown_group_key, drop = TRUE),
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
        fraction_cells_in_unknown_clusters = (
          sum(group$n_cells[unknown]) / sum(group$n_cells)
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(unknown_summary) <- NULL
  unknown_summary <- unknown_summary[order(
    match(unknown_summary$sample_id, requested_samples),
    match(unknown_summary$sampling_mode, comparison_modes),
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
    reviewed_clusters$consensus_label ==
      reviewed_clusters$current_consensus_label
  )
  expected_reviewed_rows <- (
    nrow(reviewed) * length(comparison_modes) * nrow(abstention_rules)
  )
  if (nrow(reviewed_clusters) != expected_reviewed_rows) {
    stop(
      "Reviewed-cluster join returned ", nrow(reviewed_clusters),
      " rows; expected ", expected_reviewed_rows, "."
    )
  }
  reviewed_clusters <- reviewed_clusters[order(
    match(reviewed_clusters$sample_id, requested_samples),
    match(
      paste(reviewed_clusters$sample_id,
            reviewed_clusters$seurat_cluster, sep = "|"),
      review_keys
    ),
    match(reviewed_clusters$sampling_mode, comparison_modes),
    match(reviewed_clusters$rule_id, abstention_rules$rule_id)
  ), , drop = FALSE]

  full_consensus <- consensus[
    consensus$sampling_mode == "full",
    c(
      "sample_id", "rule_id", "seurat_cluster", "consensus_label",
      "consensus_support_n", "consensus_nonunknown_n"
    ),
    drop = FALSE
  ]
  names(full_consensus)[4:6] <- paste0(
    "full_", names(full_consensus)[4:6]
  )
  cluster_comparison <- merge(
    consensus,
    full_consensus,
    by = c("sample_id", "rule_id", "seurat_cluster"),
    all = FALSE,
    sort = FALSE
  )
  cluster_comparison$matches_full_consensus <- (
    cluster_comparison$consensus_label ==
      cluster_comparison$full_consensus_label
  )
  cluster_comparison$support_difference_from_full <- (
    cluster_comparison$consensus_support_n -
      cluster_comparison$full_consensus_support_n
  )
  cluster_comparison$nonunknown_vote_difference_from_full <- (
    cluster_comparison$consensus_nonunknown_n -
      cluster_comparison$full_consensus_nonunknown_n
  )
  numeric_comparison_cluster <- suppressWarnings(
    as.numeric(cluster_comparison$seurat_cluster)
  )
  cluster_comparison <- cluster_comparison[order(
    match(cluster_comparison$sample_id, requested_samples),
    match(cluster_comparison$sampling_mode, comparison_modes),
    match(cluster_comparison$rule_id, abstention_rules$rule_id),
    is.na(numeric_comparison_cluster), numeric_comparison_cluster,
    cluster_comparison$seurat_cluster
  ), , drop = FALSE]

  convergence_group_key <- paste(
    cluster_comparison$sampling_mode,
    cluster_comparison$rule_id,
    sep = "|"
  )
  convergence_summary <- do.call(rbind, lapply(
    split(cluster_comparison, convergence_group_key, drop = TRUE),
    function(group) {
      reviewed_group <- reviewed_clusters[
        reviewed_clusters$sampling_mode == group$sampling_mode[[1]] &
          reviewed_clusters$rule_id == group$rule_id[[1]],
        , drop = FALSE
      ]
      matches_full <- group$matches_full_consensus
      unknown <- tolower(group$consensus_label) == "unknown"
      data.frame(
        sampling_mode = group$sampling_mode[[1]],
        rule_id = group$rule_id[[1]],
        n_clusters = nrow(group),
        n_clusters_matching_full = sum(matches_full),
        fraction_clusters_matching_full = mean(matches_full),
        n_cells = sum(group$n_cells),
        n_cells_in_clusters_matching_full = sum(group$n_cells[matches_full]),
        fraction_cells_in_clusters_matching_full = (
          sum(group$n_cells[matches_full]) / sum(group$n_cells)
        ),
        n_unknown_clusters = sum(unknown),
        fraction_unknown_clusters = mean(unknown),
        n_unknown_cells = sum(group$n_cells[unknown]),
        fraction_cells_in_unknown_clusters = (
          sum(group$n_cells[unknown]) / sum(group$n_cells)
        ),
        n_unanimous_clusters = sum(group$consensus_support_n == 3L),
        n_majority_clusters = sum(group$consensus_support_n == 2L),
        mean_consensus_support_n = mean(group$consensus_support_n),
        mean_nonunknown_reference_votes = mean(
          group$consensus_nonunknown_n
        ),
        reviewed_cluster_count = nrow(reviewed_group),
        reviewed_matches_expected = sum(
          reviewed_group$matches_expected_label
        ),
        reviewed_matches_current = sum(
          reviewed_group$matches_current_label
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(convergence_summary) <- NULL
  convergence_summary <- convergence_summary[order(
    match(convergence_summary$sampling_mode, comparison_modes),
    match(convergence_summary$rule_id, abstention_rules$rule_id)
  ), , drop = FALSE]

  full_votes <- reference_votes[
    reference_votes$sampling_mode == "full",
    c(
      "sample_id", "reference", "rule_id", "seurat_cluster",
      "reference_vote"
    ),
    drop = FALSE
  ]
  names(full_votes)[[5]] <- "full_reference_vote"
  vote_comparison <- merge(
    reference_votes,
    full_votes,
    by = c("sample_id", "reference", "rule_id", "seurat_cluster"),
    all = FALSE,
    sort = FALSE
  )
  vote_comparison$matches_full_reference_vote <- (
    vote_comparison$reference_vote == vote_comparison$full_reference_vote
  )
  reference_vote_group_key <- paste(
    vote_comparison$sampling_mode,
    vote_comparison$rule_id,
    vote_comparison$reference,
    sep = "|"
  )
  reference_vote_summary <- do.call(rbind, lapply(
    split(vote_comparison, reference_vote_group_key, drop = TRUE),
    function(group) {
      data.frame(
        sampling_mode = group$sampling_mode[[1]],
        rule_id = group$rule_id[[1]],
        reference = group$reference[[1]],
        n_reference_votes = nrow(group),
        n_votes_matching_full = sum(group$matches_full_reference_vote),
        fraction_votes_matching_full = mean(
          group$matches_full_reference_vote
        ),
        n_abstentions = sum(group$abstained),
        fraction_abstained = mean(group$abstained),
        median_winner_mean_score = median(group$winner_mean_score),
        median_top2_mean_score_margin = median(
          group$top2_mean_score_margin
        ),
        median_winner_cell_fraction = median(
          group$winner_cell_fraction
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(reference_vote_summary) <- NULL
  reference_vote_summary <- reference_vote_summary[order(
    match(reference_vote_summary$sampling_mode, comparison_modes),
    match(reference_vote_summary$rule_id, abstention_rules$rule_id),
    match(reference_vote_summary$reference, references)
  ), , drop = FALSE]

  combined_provenance <- data.frame(
    field = c(
      "purpose", "review_manifest", "samples", "references",
      "comparison_modes", "new_caps", "random_seed",
      "abstention_rules", "transfer_settings", "new_task_count",
      "suggested_screen_interpretation", "production_effect"
    ),
    value = c(
      "screen intermediate per-identity caps for convergence toward full-reference labels",
      review_manifest_path,
      paste(requested_samples, collapse = "|"),
      paste(references, collapse = "|"),
      paste(comparison_modes, collapse = "|"),
      paste(tested_caps, collapse = "|"),
      as.character(random_seed),
      paste(abstention_rules$rule_id, collapse = "|"),
      paste0(
        "RPCA;dims=1:", max(transfer_dimensions),
        ";k.anchor=", k_anchor, ";k.score=", k_score,
        ";k.weight=", k_weight
      ),
      as.character(nrow(task_map)),
      "identify the smallest cap with stable agreement to full; confirm the selected cap with additional seeds before production",
      "none; new isolated tables only"
    ),
    stringsAsFactors = FALSE
  )

  dir.create(cap_combined_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(cluster_metrics, combined_outputs[["cluster_metrics"]], row.names = FALSE, na = "")
  write.csv(reference_votes, combined_outputs[["reference_votes"]], row.names = FALSE, na = "")
  write.csv(consensus, combined_outputs[["consensus"]], row.names = FALSE, na = "")
  write.csv(unknown_summary, combined_outputs[["unknown_summary"]], row.names = FALSE, na = "")
  write.csv(reviewed_clusters, combined_outputs[["reviewed_clusters"]], row.names = FALSE, na = "")
  write.csv(cluster_comparison, combined_outputs[["cluster_comparison"]], row.names = FALSE, na = "")
  write.csv(convergence_summary, combined_outputs[["convergence_summary"]], row.names = FALSE, na = "")
  write.csv(reference_vote_summary, combined_outputs[["reference_vote_summary"]], row.names = FALSE, na = "")
  write.csv(reference_balance, combined_outputs[["reference_balance"]], row.names = FALSE, na = "")
  write.csv(combined_provenance, combined_outputs[["provenance"]], row.names = FALSE, na = "")
  message("Wrote five-mode cap-convergence summaries: ", cap_combined_dir)
  quit(save = "no", status = 0L)
}
if (length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]\n",
    "   or: Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R ",
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
reference_cap <- task_map$reference_cap_per_identity[[task_id]]
query_path <- query_path_for(sample_id)
reference_path <- unname(reference_paths[[reference_name]])
output_paths <- task_output_paths(sample_id, reference_name, sampling_mode)

if (dry_run) {
  ok <- compact_dry_run(
    paste0(
      "Reference-cap convergence task ", task_id, "/", nrow(task_map),
      " [", sample_id, " | ", reference_name, " | ", sampling_mode, "]"
    ),
    inputs = c(query_path, reference_path, review_manifest_path),
    outputs = output_paths,
    checks = c(
      sample_matches_manifest = sum(sample_manifest$sample_id == sample_id) == 1L,
      valid_reference = reference_name %in% references,
      valid_sampling_mode = sampling_mode %in% sampling_modes,
      valid_reference_cap = is.finite(reference_cap) && reference_cap >= 1L,
      expected_rules = nrow(abstention_rules) == 4L
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- c(query_path, reference_path)[
  !file.exists(c(query_path, reference_path))
]
if (length(missing_inputs)) {
  stop("Missing cap-convergence input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite cap-convergence outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing tables."
  )
}

set.seed(random_seed)
query <- readRDS(query_path)
reference <- readRDS(reference_path)

if (anyDuplicated(Cells(query))) stop("Query contains duplicate cell IDs: ", sample_id)
if (anyDuplicated(Cells(reference))) {
  stop("Reference contains duplicate cell IDs: ", reference_name)
}
if (!identical(rownames(query[[]]), Cells(query))) {
  stop("Query metadata rows do not exactly align with query cell IDs.")
}
if (!identical(rownames(reference[[]]), Cells(reference))) {
  stop("Reference metadata rows do not exactly align with reference cell IDs.")
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
set.seed(random_seed)
reference_for_transfer <- subset(reference, downsample = reference_cap)

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
reference_balance$reference_cap_per_identity <- reference_cap
reference_balance <- reference_balance[c(
  "sample_id", "reference", "sampling_mode", "reference_class",
  "full_reference_cells", "used_reference_cells", "retained_fraction",
  "reference_cap_per_identity"
)]
full_reference_cells_total <- ncol(reference)
expected_used_cells <- pmin(
  reference_balance$full_reference_cells,
  reference_cap
)
if (!identical(
  as.integer(reference_balance$used_reference_cells),
  as.integer(expected_used_cells)
)) {
  stop(
    "Downsampling did not retain exactly min(class size, cap) cells for every identity."
  )
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
set.seed(random_seed)
reference_for_transfer <- RunPCA(
  reference_for_transfer, features = transfer_features, verbose = FALSE
)
set.seed(random_seed)
query <- RunPCA(query, features = transfer_features, verbose = FALSE)

set.seed(random_seed)
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
    reference_path, review_manifest_path, as.character(random_seed),
    as.character(ncol(query)), as.character(full_reference_cells_total),
    as.character(ncol(reference_for_transfer)), as.character(reference_cap),
    as.character(length(shared_genes)), as.character(length(transfer_features)),
    paste(transfer_dimensions, collapse = ","), config$label_transfer$method,
    as.character(k_anchor), as.character(k_score), as.character(k_weight),
    as.character(legacy_threshold),
    paste(abstention_rules$rule_id, collapse = "|"),
    "seed reset before sampling, each PCA, and anchor finding to isolate sampling mode",
    as.character(anchor_count),
    "cluster metrics, reference votes, reference balance, provenance; no RDS",
    "none; isolated cap-convergence tables only"
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
  "Wrote isolated cap-convergence tables for ", sample_id, " / ",
  reference_name, " / ", sampling_mode, "."
)
