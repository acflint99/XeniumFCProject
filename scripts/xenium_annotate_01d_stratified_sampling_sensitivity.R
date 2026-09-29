#!/usr/bin/env Rscript

# Add an age/stage-stratified cap of 1,000 cells per clusters_refined identity
# to the completed cap1000/full reference-sampling sensitivity analysis. The
# three reviewed resolution-5 query samples are mapped to Aldinger, Sepp, and
# Science. Outputs are compact tables only; production objects are untouched.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
references <- c("Aldinger", "Sepp", "Science")
sampling_mode <- "stratified_cap1000"
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
    "Every stratified-sampling query must match exactly one samples.csv row: ",
    paste(requested_samples[sample_counts != 1L], collapse = ", ")
  )
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

stage_candidates <- list(
  Aldinger = c(
    "age", "Age", "PCW", "pcw", "stage", "Stage",
    "gestational_age", "gestational.age", "week", "Week"
  ),
  Sepp = c(
    "Stage", "stage", "age", "Age", "PCW", "pcw",
    "gestational_age", "gestational.age", "week", "Week"
  ),
  Science = c(
    "PCW", "pcw", "age", "Age", "stage", "Stage",
    "gestational_age", "gestational.age", "week", "Week"
  )
)

clean_stratum_values <- function(values, field, reference_name) {
  values <- trimws(as.character(values))
  if (anyNA(values) || any(!nzchar(values))) {
    stop(
      reference_name, " reference stratum field '", field,
      "' contains ", sum(is.na(values) | !nzchar(values)),
      " missing or blank value(s)."
    )
  }
  values
}

select_stage_field <- function(object, reference_name) {
  metadata_names <- colnames(object[[]])
  explicit <- intersect(stage_candidates[[reference_name]], metadata_names)
  regex_matches <- metadata_names[grepl(
    "age|stage|pcw|gest|week|wpc",
    metadata_names,
    ignore.case = TRUE
  )]
  candidates <- unique(c(explicit, regex_matches))
  if (!length(candidates)) {
    stop(
      reference_name,
      " reference lacks an age/stage metadata field. Available metadata: ",
      paste(metadata_names, collapse = ", ")
    )
  }

  candidate_summary <- lapply(candidates, function(field) {
    values <- trimws(as.character(object[[field]][, 1]))
    valid <- !is.na(values) & nzchar(values)
    data.frame(
      field = field,
      missing_n = sum(!valid),
      unique_n = length(unique(values[valid])),
      stringsAsFactors = FALSE
    )
  })
  candidate_summary <- do.call(rbind, candidate_summary)
  usable <- candidate_summary[
    candidate_summary$missing_n == 0L &
      candidate_summary$unique_n >= 2L &
      candidate_summary$unique_n <= 50L,
    , drop = FALSE
  ]
  if (!nrow(usable)) {
    stop(
      reference_name,
      " reference has no complete, informative age/stage field with 2-50 levels. Candidates:\n",
      paste(
        paste0(
          candidate_summary$field,
          " [levels=", candidate_summary$unique_n,
          ", missing=", candidate_summary$missing_n, "]"
        ),
        collapse = "\n"
      )
    )
  }
  usable$field[[1]]
}

allocate_evenly <- function(stratum_counts, target_n) {
  stratum_names <- names(stratum_counts)
  stratum_counts <- as.integer(stratum_counts)
  names(stratum_counts) <- stratum_names
  if (!length(stratum_counts) || any(stratum_counts <= 0L)) {
    stop("Stratum counts must be positive.")
  }
  if (target_n < 1L || target_n > sum(stratum_counts)) {
    stop("Invalid stratified-sampling target: ", target_n)
  }
  allocation <- setNames(integer(length(stratum_counts)), names(stratum_counts))
  for (step in seq_len(target_n)) {
    eligible <- which(allocation < stratum_counts)
    if (!length(eligible)) stop("Stratified allocation exhausted cells too early.")
    minimum_allocated <- min(allocation[eligible])
    tied <- eligible[allocation[eligible] == minimum_allocated]
    selected_index <- tied[[sample.int(length(tied), size = 1L)]]
    allocation[[selected_index]] <- allocation[[selected_index]] + 1L
  }
  if (sum(allocation) != target_n || any(allocation > stratum_counts)) {
    stop("Stratified allocation failed its quota checks.")
  }
  allocation
}

sample_stratified_reference <- function(reference, stratum_field, cap, seed) {
  cells <- Cells(reference)
  labels <- trimws(as.character(reference$clusters_refined))
  strata <- clean_stratum_values(
    reference[[stratum_field]][, 1], stratum_field, "Reference"
  )
  if (length(cells) != length(labels) || length(cells) != length(strata)) {
    stop("Reference cells, labels, and strata do not align one-to-one.")
  }
  if (anyDuplicated(cells)) stop("Reference contains duplicate cell IDs.")

  set.seed(seed)
  selected_cells <- character()
  allocation_rows <- list()
  allocation_index <- 1L
  for (label in sort(unique(labels))) {
    label_index <- which(labels == label)
    label_strata <- strata[label_index]
    stratum_levels <- sort(unique(label_strata))
    stratum_counts <- table(factor(label_strata, levels = stratum_levels))
    target_n <- min(cap, length(label_index))
    allocation <- allocate_evenly(stratum_counts, target_n)

    for (stratum_value in names(allocation)) {
      stratum_index <- label_index[label_strata == stratum_value]
      take_n <- allocation[[stratum_value]]
      chosen <- if (take_n == length(stratum_index)) {
        cells[stratum_index]
      } else {
        sample(cells[stratum_index], size = take_n, replace = FALSE)
      }
      selected_cells <- c(selected_cells, chosen)
      allocation_rows[[allocation_index]] <- data.frame(
        reference_class = label,
        stratum_field = stratum_field,
        stratum_value = stratum_value,
        full_reference_cells = length(stratum_index),
        used_reference_cells = take_n,
        retained_fraction = take_n / length(stratum_index),
        stringsAsFactors = FALSE
      )
      allocation_index <- allocation_index + 1L
    }
  }

  if (anyDuplicated(selected_cells)) {
    stop("Stratified sampling selected duplicate cell IDs.")
  }
  expected_n <- sum(pmin(as.integer(table(labels)), cap))
  if (length(selected_cells) != expected_n) {
    stop(
      "Stratified sampling selected ", length(selected_cells),
      " cells; expected ", expected_n, "."
    )
  }
  selected_cells <- cells[cells %in% selected_cells]
  list(
    object = subset(reference, cells = selected_cells),
    strata_balance = do.call(rbind, allocation_rows)
  )
}

task_map <- expand.grid(
  sample_id = requested_samples,
  reference = references,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_map <- task_map[order(
  match(task_map$sample_id, requested_samples),
  match(task_map$reference, references)
), , drop = FALSE]
task_map$sampling_mode <- sampling_mode
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c("task_id", "sample_id", "reference", "sampling_mode")]
if (nrow(task_map) != 9L) stop("Expected exactly nine stratified-sampling tasks.")

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--list", "--dry-run", "--overwrite", "--combine", "--inspect-strata"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
combine <- "--combine" %in% args
inspect_strata <- "--inspect-strata" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (combine && inspect_strata) stop("Choose only one of --combine or --inspect-strata.")
positional_args <- args[!args %in% valid_options]

if (list_requested) {
  if (length(positional_args)) stop("--list does not accept TASK_ID.")
  cat("Scope: three reviewed resolution-5 samples x three references\n")
  cat("Sampling: maximum ", reference_cap,
      " cells per identity, allocated evenly across developmental age/stage strata\n",
      sep = "")
  cat(
    "Fixed transfer settings: RPCA, dimensions=1:", max(transfer_dimensions),
    ", k.anchor=", k_anchor, ", k.score=", k_score,
    ", k.weight=", k_weight, ", seed=", config$runtime$random_seed, "\n",
    sep = ""
  )
  cat("Candidate age/stage fields by reference:\n")
  for (reference_name in references) {
    cat(reference_name, ": ", paste(stage_candidates[[reference_name]], collapse = ", "), "\n", sep = "")
  }
  cat("\nAbstention rules:\n")
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
preflight_dir <- file.path(sensitivity_root, "stratified_preflight")
preflight_output <- file.path(preflight_dir, "reference_age_stage_strata.csv")

query_path_for <- function(sample_id) {
  file.path(query_root, paste0(sample_id, "_whole_tissue_Res5.0.rds"))
}

task_output_paths <- function(sample_id, reference_name) {
  output_dir <- file.path(
    sensitivity_root, sampling_mode, tolower(reference_name), "tables"
  )
  prefix <- paste(sample_id, reference_name, sampling_mode, sep = "_")
  c(
    cluster_metrics = file.path(output_dir, paste0(prefix, "_cluster_metrics.csv")),
    reference_votes = file.path(output_dir, paste0(prefix, "_reference_votes.csv")),
    reference_balance = file.path(output_dir, paste0(prefix, "_reference_class_balance.csv")),
    strata_balance = file.path(output_dir, paste0(prefix, "_reference_strata_balance.csv")),
    provenance = file.path(output_dir, paste0(prefix, "_provenance.csv"))
  )
}

all_task_outputs <- do.call(rbind, lapply(seq_len(nrow(task_map)), function(i) {
  paths <- task_output_paths(task_map$sample_id[[i]], task_map$reference[[i]])
  data.frame(
    task_id = task_map$task_id[[i]],
    sample_id = task_map$sample_id[[i]],
    reference = task_map$reference[[i]],
    output_type = names(paths),
    output_path = unname(paths),
    stringsAsFactors = FALSE
  )
}))

existing_combined_dir <- file.path(sensitivity_root, "combined_tables")
existing_combined_inputs <- c(
  cluster_metrics = file.path(existing_combined_dir, "sampling_abstention_cluster_metrics.csv"),
  reference_votes = file.path(existing_combined_dir, "sampling_abstention_reference_votes.csv"),
  consensus = file.path(existing_combined_dir, "sampling_abstention_consensus_results.csv"),
  reference_balance = file.path(existing_combined_dir, "sampling_abstention_reference_class_balance.csv")
)
stratified_combined_dir <- file.path(sensitivity_root, "stratified_combined_tables")
combined_outputs <- c(
  cluster_metrics = file.path(stratified_combined_dir, "three_mode_cluster_metrics.csv"),
  reference_votes = file.path(stratified_combined_dir, "three_mode_reference_votes.csv"),
  consensus = file.path(stratified_combined_dir, "three_mode_consensus_results.csv"),
  unknown_summary = file.path(stratified_combined_dir, "three_mode_unknown_summary.csv"),
  reviewed_clusters = file.path(stratified_combined_dir, "three_mode_reviewed_cluster_results.csv"),
  sampling_comparison = file.path(stratified_combined_dir, "three_mode_sampling_comparison.csv"),
  reference_balance = file.path(stratified_combined_dir, "three_mode_reference_class_balance.csv"),
  strata_balance = file.path(stratified_combined_dir, "stratified_reference_strata_balance.csv"),
  provenance = file.path(stratified_combined_dir, "three_mode_combined_provenance.csv")
)

if (inspect_strata) {
  if (length(positional_args)) stop("--inspect-strata does not accept TASK_ID.")
  if (dry_run) {
    ok <- compact_dry_run(
      "Inspect reference age/stage strata",
      inputs = unname(reference_paths),
      outputs = preflight_output,
      checks = c(three_references = length(reference_paths) == 3L)
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }
  missing_references <- reference_paths[!file.exists(reference_paths)]
  if (length(missing_references)) {
    stop("Missing reference RDS file(s):\n- ", paste(missing_references, collapse = "\n- "))
  }
  if (file.exists(preflight_output) && !overwrite) {
    stop(
      "Refusing to overwrite stratification preflight: ", preflight_output,
      "\nUse --overwrite only after reviewing the existing file."
    )
  }
  preflight_rows <- list()
  for (reference_name in references) {
    reference <- readRDS(reference_paths[[reference_name]])
    if (!"clusters_refined" %in% colnames(reference[[]])) {
      stop(reference_name, " reference lacks clusters_refined metadata.")
    }
    field <- select_stage_field(reference, reference_name)
    values <- clean_stratum_values(
      reference[[field]][, 1], field, reference_name
    )
    counts <- as.data.frame(table(stratum_value = values), stringsAsFactors = FALSE)
    names(counts)[[2]] <- "n_cells"
    counts$reference <- reference_name
    counts$stratum_field <- field
    counts$n_strata <- nrow(counts)
    counts$n_reference_cells <- ncol(reference)
    preflight_rows[[reference_name]] <- counts[c(
      "reference", "stratum_field", "stratum_value", "n_cells",
      "n_strata", "n_reference_cells"
    )]
    rm(reference)
    invisible(gc(verbose = FALSE))
  }
  preflight <- do.call(rbind, preflight_rows)
  dir.create(preflight_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(preflight, preflight_output, row.names = FALSE, na = "")
  message("Wrote validated reference age/stage strata: ", preflight_output)
  quit(save = "no", status = 0L)
}

read_task_output <- function(output_type) {
  paths <- all_task_outputs$output_path[all_task_outputs$output_type == output_type]
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
    stop("Consensus group does not contain exactly one vote from each reference.")
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
    sampling_mode = sampling_mode,
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
  required_inputs <- c(existing_combined_inputs, all_task_outputs$output_path, preflight_output)
  if (dry_run) {
    ok <- compact_dry_run(
      "Combine capped, stratified, and full-reference sensitivity results",
      inputs = required_inputs,
      outputs = combined_outputs,
      checks = c(
        expected_new_task_outputs = nrow(all_task_outputs) == 9L * 5L,
        expected_rules = nrow(abstention_rules) == 4L
      )
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Missing combine input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
  }
  existing_outputs <- combined_outputs[file.exists(combined_outputs)]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite stratified combined outputs:\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing the existing tables."
    )
  }

  old_metrics <- read.csv(existing_combined_inputs[["cluster_metrics"]], stringsAsFactors = FALSE, check.names = FALSE)
  old_votes <- read.csv(existing_combined_inputs[["reference_votes"]], stringsAsFactors = FALSE, check.names = FALSE)
  old_consensus <- read.csv(existing_combined_inputs[["consensus"]], stringsAsFactors = FALSE, check.names = FALSE)
  old_balance <- read.csv(existing_combined_inputs[["reference_balance"]], stringsAsFactors = FALSE, check.names = FALSE)
  new_metrics <- read_task_output("cluster_metrics")
  new_votes <- read_task_output("reference_votes")
  new_balance <- read_task_output("reference_balance")
  strata_balance <- read_task_output("strata_balance")

  old_metrics <- old_metrics[old_metrics$sampling_mode %in% c("cap1000", "full"), , drop = FALSE]
  old_votes <- old_votes[old_votes$sampling_mode %in% c("cap1000", "full"), , drop = FALSE]
  old_consensus <- old_consensus[old_consensus$sampling_mode %in% c("cap1000", "full"), , drop = FALSE]
  old_balance <- old_balance[old_balance$sampling_mode %in% c("cap1000", "full"), , drop = FALSE]

  all_metrics <- bind_same_columns(
    old_metrics, new_metrics, "existing cluster metrics", "stratified cluster metrics"
  )
  all_votes <- bind_same_columns(
    old_votes, new_votes, "existing reference votes", "stratified reference votes"
  )
  all_balance <- bind_same_columns(
    old_balance, new_balance, "existing class balance", "stratified class balance"
  )
  all_metrics$seurat_cluster <- as.character(all_metrics$seurat_cluster)
  all_votes$seurat_cluster <- as.character(all_votes$seurat_cluster)
  old_consensus$seurat_cluster <- as.character(old_consensus$seurat_cluster)

  new_group_key <- paste(
    new_votes$sample_id, new_votes$rule_id, new_votes$seurat_cluster, sep = "|"
  )
  new_consensus <- do.call(rbind, lapply(
    split(new_votes, new_group_key, drop = TRUE), consensus_from_votes
  ))
  rownames(new_consensus) <- NULL
  all_consensus <- bind_same_columns(
    old_consensus, new_consensus, "existing consensus", "stratified consensus"
  )
  sampling_order <- c("cap1000", sampling_mode, "full")
  numeric_cluster <- suppressWarnings(as.numeric(all_consensus$seurat_cluster))
  all_consensus <- all_consensus[order(
    match(all_consensus$sample_id, requested_samples),
    match(all_consensus$sampling_mode, sampling_order),
    match(all_consensus$rule_id, abstention_rules$rule_id),
    is.na(numeric_cluster), numeric_cluster, all_consensus$seurat_cluster
  ), , drop = FALSE]

  unknown_group_key <- paste(
    all_consensus$sample_id, all_consensus$sampling_mode,
    all_consensus$rule_id, sep = "|"
  )
  unknown_summary <- do.call(rbind, lapply(
    split(all_consensus, unknown_group_key, drop = TRUE),
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

  reviewed_results <- merge(
    all_consensus,
    reviewed[required_review_columns],
    by = c("sample_id", "seurat_cluster"),
    all = FALSE,
    sort = FALSE
  )
  reviewed_results$matches_expected_label <- (
    reviewed_results$consensus_label == reviewed_results$expected_label
  )
  reviewed_results$matches_current_label <- (
    reviewed_results$consensus_label == reviewed_results$current_consensus_label
  )

  comparison_columns <- c(
    "sample_id", "rule_id", "seurat_cluster", "consensus_label",
    "consensus_support_n", "consensus_nonunknown_n"
  )
  mode_table <- function(mode, prefix) {
    tab <- all_consensus[all_consensus$sampling_mode == mode, comparison_columns, drop = FALSE]
    names(tab)[4:6] <- paste0(prefix, "_", names(tab)[4:6])
    tab
  }
  comparison <- Reduce(
    function(x, y) merge(
      x, y, by = c("sample_id", "rule_id", "seurat_cluster"),
      all = TRUE, sort = FALSE
    ),
    list(
      mode_table("cap1000", "cap1000"),
      mode_table(sampling_mode, "stratified"),
      mode_table("full", "full")
    )
  )
  label_columns <- c(
    "cap1000_consensus_label", "stratified_consensus_label",
    "full_consensus_label"
  )
  if (anyNA(comparison[label_columns])) {
    stop("The three sampling modes do not contain identical cluster sets.")
  }
  comparison$stratified_vs_cap1000_changed <- (
    comparison$stratified_consensus_label != comparison$cap1000_consensus_label
  )
  comparison$stratified_vs_full_changed <- (
    comparison$stratified_consensus_label != comparison$full_consensus_label
  )

  combined_provenance <- data.frame(
    field = c(
      "purpose", "review_manifest", "samples", "references",
      "sampling_modes", "stratification", "abstention_rules",
      "transfer_settings", "new_task_count", "production_effect"
    ),
    value = c(
      "add developmental age/stage-stratified cap1000 to completed sampling sensitivity",
      review_manifest_path,
      paste(requested_samples, collapse = "|"),
      paste(references, collapse = "|"),
      paste(c("cap1000", sampling_mode, "full"), collapse = "|"),
      "within each clusters_refined identity, allocate at most 1000 cells as evenly as possible across validated age/stage levels",
      paste(abstention_rules$rule_id, collapse = "|"),
      paste0(
        "RPCA;dims=1:", max(transfer_dimensions),
        ";k.anchor=", k_anchor, ";k.score=", k_score,
        ";k.weight=", k_weight, ";seed=", config$runtime$random_seed
      ),
      as.character(nrow(task_map)),
      "none; new isolated tables only"
    ),
    stringsAsFactors = FALSE
  )

  dir.create(stratified_combined_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_metrics, combined_outputs[["cluster_metrics"]], row.names = FALSE, na = "")
  write.csv(all_votes, combined_outputs[["reference_votes"]], row.names = FALSE, na = "")
  write.csv(all_consensus, combined_outputs[["consensus"]], row.names = FALSE, na = "")
  write.csv(unknown_summary, combined_outputs[["unknown_summary"]], row.names = FALSE, na = "")
  write.csv(reviewed_results, combined_outputs[["reviewed_clusters"]], row.names = FALSE, na = "")
  write.csv(comparison, combined_outputs[["sampling_comparison"]], row.names = FALSE, na = "")
  write.csv(all_balance, combined_outputs[["reference_balance"]], row.names = FALSE, na = "")
  write.csv(strata_balance, combined_outputs[["strata_balance"]], row.names = FALSE, na = "")
  write.csv(combined_provenance, combined_outputs[["provenance"]], row.names = FALSE, na = "")
  message("Wrote three-mode stratified sensitivity summaries: ", stratified_combined_dir)
  quit(save = "no", status = 0L)
}

if (length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R ",
    "[--list|--dry-run|--overwrite] [TASK_ID]\n",
    "   or: Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R ",
    "--inspect-strata [--dry-run|--overwrite]\n",
    "   or: Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R ",
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
query_path <- query_path_for(sample_id)
reference_path <- unname(reference_paths[[reference_name]])
output_paths <- task_output_paths(sample_id, reference_name)

if (dry_run) {
  ok <- compact_dry_run(
    paste0(
      "Stratified sampling task ", task_id, "/", nrow(task_map),
      " [", sample_id, " | ", reference_name, "]"
    ),
    inputs = c(query_path, reference_path, review_manifest_path, preflight_output),
    outputs = output_paths,
    checks = c(
      sample_matches_manifest = sum(sample_manifest$sample_id == sample_id) == 1L,
      valid_reference = reference_name %in% references,
      expected_rules = nrow(abstention_rules) == 4L
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- c(query_path, reference_path, preflight_output)[
  !file.exists(c(query_path, reference_path, preflight_output))
]
if (length(missing_inputs)) {
  stop("Missing stratified-sampling input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite stratified-sampling outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing tables."
  )
}

preflight <- read.csv(preflight_output, stringsAsFactors = FALSE, check.names = FALSE)
required_preflight_columns <- c(
  "reference", "stratum_field", "stratum_value", "n_cells",
  "n_strata", "n_reference_cells"
)
missing_preflight_columns <- setdiff(required_preflight_columns, names(preflight))
if (length(missing_preflight_columns)) {
  stop(
    "Stratification preflight lacks columns: ",
    paste(missing_preflight_columns, collapse = ", ")
  )
}
preflight_fields <- unique(preflight$stratum_field[preflight$reference == reference_name])
if (length(preflight_fields) != 1L || is.na(preflight_fields)) {
  stop("Preflight does not define exactly one stratum field for ", reference_name, ".")
}
stratum_field <- preflight_fields[[1]]

set.seed(as.integer(config$runtime$random_seed))
query <- readRDS(query_path)
reference <- readRDS(reference_path)
if (anyDuplicated(Cells(query))) stop("Query contains duplicate cell IDs: ", sample_id)
if (anyDuplicated(Cells(reference))) stop("Reference contains duplicate cell IDs: ", reference_name)
if (!identical(rownames(query[[]]), Cells(query))) {
  stop("Query metadata rows do not exactly align with cell IDs: ", sample_id)
}
if (!identical(rownames(reference[[]]), Cells(reference))) {
  stop("Reference metadata rows do not exactly align with cell IDs: ", reference_name)
}
if (!"clusters_refined" %in% colnames(reference[[]])) {
  stop(reference_name, " reference lacks clusters_refined metadata.")
}
reference_labels <- trimws(as.character(reference$clusters_refined))
if (anyNA(reference_labels) || any(!nzchar(reference_labels))) {
  stop(reference_name, " reference has blank or missing clusters_refined labels.")
}
reference$clusters_refined <- reference_labels
if (!stratum_field %in% colnames(reference[[]])) {
  stop(reference_name, " reference no longer contains preflight field '", stratum_field, "'.")
}
detected_field <- select_stage_field(reference, reference_name)
if (!identical(detected_field, stratum_field)) {
  stop(
    reference_name, " preflight selected '", stratum_field,
    "' but task detection selected '", detected_field, "'."
  )
}
current_strata <- clean_stratum_values(
  reference[[stratum_field]][, 1], stratum_field, reference_name
)
current_stratum_counts <- as.data.frame(
  table(stratum_value = current_strata), stringsAsFactors = FALSE
)
names(current_stratum_counts)[[2]] <- "n_cells"
preflight_reference <- preflight[
  preflight$reference == reference_name,
  c("stratum_value", "n_cells", "n_strata", "n_reference_cells"),
  drop = FALSE
]
if (anyDuplicated(preflight_reference$stratum_value)) {
  stop("Preflight contains duplicate strata for ", reference_name, ".")
}
preflight_reference <- preflight_reference[
  order(preflight_reference$stratum_value), , drop = FALSE
]
current_stratum_counts <- current_stratum_counts[
  order(current_stratum_counts$stratum_value), , drop = FALSE
]
counts_match <- identical(
  as.character(preflight_reference$stratum_value),
  as.character(current_stratum_counts$stratum_value)
) && identical(
  as.integer(preflight_reference$n_cells),
  as.integer(current_stratum_counts$n_cells)
)
summary_match <- nrow(preflight_reference) > 0L &&
  !anyNA(preflight_reference) &&
  length(unique(preflight_reference$n_strata)) == 1L &&
  unique(as.integer(preflight_reference$n_strata)) == nrow(current_stratum_counts) &&
  length(unique(preflight_reference$n_reference_cells)) == 1L &&
  unique(as.integer(preflight_reference$n_reference_cells)) == ncol(reference)
if (!counts_match || !summary_match) {
  stop(
    reference_name,
    " age/stage counts no longer match the preflight; rerun the preflight."
  )
}

required_query_metadata <- c(
  "whole_tissue_cluster_res5.0", "Xenium_snn_res.5", "seurat_clusters"
)
missing_query_metadata <- setdiff(required_query_metadata, colnames(query[[]]))
if (length(missing_query_metadata)) {
  stop(sample_id, " resolution-5 input lacks metadata: ", paste(missing_query_metadata, collapse = ", "))
}
cluster_ids <- as.character(query$seurat_clusters)
if (anyNA(cluster_ids) ||
    !identical(cluster_ids, as.character(query$whole_tissue_cluster_res5.0))) {
  stop("Query seurat_clusters does not exactly match whole_tissue_cluster_res5.0.")
}
missing_reviewed_clusters <- setdiff(
  reviewed$seurat_cluster[reviewed$sample_id == sample_id], unique(cluster_ids)
)
if (length(missing_reviewed_clusters)) {
  stop(sample_id, " lacks reviewed cluster(s): ", paste(missing_reviewed_clusters, collapse = ", "))
}

DefaultAssay(reference) <- if ("RNA" %in% Assays(reference)) "RNA" else Assays(reference)[[1]]
DefaultAssay(query) <- if ("Xenium" %in% Assays(query)) "Xenium" else Assays(query)[[1]]
reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures = as.integer(config$label_transfer$variable_features),
  verbose = FALSE
)
shared_genes <- intersect(rownames(reference), rownames(query))

stratified <- sample_stratified_reference(
  reference,
  stratum_field = stratum_field,
  cap = reference_cap,
  seed = as.integer(config$runtime$random_seed)
)
reference_for_transfer <- stratified$object
strata_balance <- stratified$strata_balance
strata_balance$sample_id <- sample_id
strata_balance$reference <- reference_name
strata_balance$sampling_mode <- sampling_mode
strata_balance <- strata_balance[c(
  "sample_id", "reference", "sampling_mode", "reference_class",
  "stratum_field", "stratum_value", "full_reference_cells",
  "used_reference_cells", "retained_fraction"
)]

class_count <- function(object, count_name) {
  counts <- as.data.frame(
    table(reference_class = as.character(object$clusters_refined)),
    stringsAsFactors = FALSE
  )
  names(counts)[[2]] <- count_name
  counts
}
reference_balance <- merge(
  class_count(reference, "full_reference_cells"),
  class_count(reference_for_transfer, "used_reference_cells"),
  by = "reference_class", all = TRUE, sort = FALSE
)
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
expected_class_counts <- pmin(
  reference_balance$full_reference_cells, reference_cap
)
if (anyNA(reference_balance) ||
    !identical(
      as.integer(reference_balance$used_reference_cells),
      as.integer(expected_class_counts)
    )) {
  stop("Stratified reference did not retain the expected capped class counts.")
}
full_reference_cells_total <- ncol(reference)
rm(reference, stratified)
invisible(gc(verbose = FALSE))

transfer_features <- intersect(VariableFeatures(reference_for_transfer), shared_genes)
if (length(transfer_features) <= max(transfer_dimensions)) {
  stop("Too few transfer features: ", length(transfer_features))
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
score_columns <- setdiff(score_columns, c("prediction.score.max", "prediction.score.id"))
if (length(score_columns) < 2L) stop("TransferData returned fewer than two class scores.")
score_matrix <- as.matrix(predictions[, score_columns, drop = FALSE])
storage.mode(score_matrix) <- "double"
class_labels <- sub("^prediction\\.score\\.", "", score_columns)
if (anyDuplicated(class_labels)) stop("TransferData returned duplicate class labels.")
colnames(score_matrix) <- class_labels
if (any(!is.finite(score_matrix))) stop("TransferData returned non-finite class scores.")
predicted_ids <- trimws(as.character(predictions$predicted.id))
max_scores <- suppressWarnings(as.numeric(predictions$prediction.score.max))
if (anyNA(predicted_ids) || any(!nzchar(predicted_ids)) || any(!is.finite(max_scores))) {
  stop("TransferData returned invalid labels or scores.")
}

cluster_levels <- unique(cluster_ids)
numeric_clusters <- suppressWarnings(as.numeric(cluster_levels))
cluster_levels <- cluster_levels[order(is.na(numeric_clusters), numeric_clusters, cluster_levels)]
focus_classes <- c("RL", "Granule", "Glia", "VZ")
cluster_rows <- vector("list", length(cluster_levels))
vote_rows <- vector("list", length(cluster_levels) * nrow(abstention_rules))
vote_index <- 1L

for (cluster_index in seq_along(cluster_levels)) {
  cluster_id <- cluster_levels[[cluster_index]]
  cell_index <- which(cluster_ids == cluster_id)
  mean_scores <- colMeans(score_matrix[cell_index, , drop = FALSE])
  score_order <- order(mean_scores, decreasing = TRUE)
  winner <- names(mean_scores)[score_order[[1]]]
  runner_up <- names(mean_scores)[score_order[[2]]]
  winner_score <- unname(mean_scores[[winner]])
  runner_up_score <- unname(mean_scores[[runner_up]])
  margin <- winner_score - runner_up_score
  winner_fraction <- mean(predicted_ids[cell_index] == winner)
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
    top2_mean_score_margin = margin,
    winner_cell_fraction = winner_fraction,
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
    failed_margin <- margin < rule$min_top2_mean_score_margin[[1]]
    failed_cell_support <- winner_fraction < rule$min_winner_cell_fraction[[1]]
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
      top2_mean_score_margin = margin,
      winner_cell_fraction = winner_fraction,
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
if (sum(cluster_metrics$n_cells) != ncol(query)) {
  stop("Cluster metric cell counts do not sum to the query cell count.")
}
if (anyDuplicated(cluster_metrics$seurat_cluster)) {
  stop("Cluster metrics contain duplicate cluster IDs.")
}
production_votes <- reference_votes[reference_votes$rule_id == "production_score", , drop = FALSE]
expected_votes <- ifelse(
  cluster_metrics$weighted_winner_mean_score >= legacy_threshold,
  cluster_metrics$weighted_winner_before_abstention,
  "Unknown"
)
if (!identical(as.character(production_votes$reference_vote), as.character(expected_votes))) {
  stop("production_score does not reproduce the current weighted threshold rule.")
}

anchor_count <- tryCatch(
  nrow(methods::slot(anchors, "anchors")),
  error = function(e) NA_integer_
)
provenance <- data.frame(
  field = c(
    "sample_id", "reference", "sampling_mode", "stratum_field",
    "sampling_rule", "query_path", "reference_path", "review_manifest",
    "random_seed", "query_cells", "full_reference_cells",
    "used_reference_cells", "reference_cap", "shared_genes",
    "transfer_features", "dimensions", "reduction", "k_anchor",
    "k_score", "k_weight", "production_threshold", "abstention_rules",
    "seed_control", "anchor_count", "outputs", "production_effect"
  ),
  value = c(
    sample_id, reference_name, sampling_mode, stratum_field,
    "within each clusters_refined identity, allocate cap evenly across age/stage strata; redistribute when a stratum has fewer cells",
    query_path, reference_path, review_manifest_path,
    as.character(config$runtime$random_seed), as.character(ncol(query)),
    as.character(full_reference_cells_total), as.character(ncol(reference_for_transfer)),
    as.character(reference_cap), as.character(length(shared_genes)),
    as.character(length(transfer_features)), paste(transfer_dimensions, collapse = ","),
    config$label_transfer$method, as.character(k_anchor), as.character(k_score),
    as.character(k_weight), as.character(legacy_threshold),
    paste(abstention_rules$rule_id, collapse = "|"),
    "seed reset before sampling, each PCA, and anchor finding",
    as.character(anchor_count),
    "cluster metrics, reference votes, class balance, stratum balance, provenance; no RDS",
    "none; isolated sensitivity tables only"
  ),
  stringsAsFactors = FALSE
)

output_dir <- dirname(output_paths[["cluster_metrics"]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_metrics, output_paths[["cluster_metrics"]], row.names = FALSE, na = "")
write.csv(reference_votes, output_paths[["reference_votes"]], row.names = FALSE, na = "")
write.csv(reference_balance, output_paths[["reference_balance"]], row.names = FALSE, na = "")
write.csv(strata_balance, output_paths[["strata_balance"]], row.names = FALSE, na = "")
write.csv(provenance, output_paths[["provenance"]], row.names = FALSE, na = "")
message(
  "Wrote stratified-sampling sensitivity tables for ", sample_id, " / ",
  reference_name, " using '", stratum_field, "'."
)
