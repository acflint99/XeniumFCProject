#!/usr/bin/env Rscript

# Combine all broad cluster scores and evaluate prespecified per-reference
# score, margin, and cell-support grids. Agreement of the other two references
# is used as an internal concordance target; it is not claimed as ground truth.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_broad_reference.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
references <- broad_reference_names()

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) {
  stop("Usage: Rscript scripts/xenium_broad_03_compile_calibration.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

validation_root <- file.path(
  here(config$project$outputs_dir), "validation", "broad_reference_workflow"
)
score_dir <- file.path(validation_root, "02_cluster_scores", "per_task")
input_grid <- expand.grid(
  sample_id = sample_ids,
  reference = references,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
input_grid <- input_grid[order(
  match(input_grid$sample_id, sample_ids),
  match(input_grid$reference, references)
), , drop = FALSE]
input_grid$path <- file.path(
  score_dir,
  paste0(input_grid$sample_id, "_", input_grid$reference, "_Broad_Cluster_Scores.csv")
)
reference_validation_path <- file.path(
  validation_root, "01_reference_validation", "tables",
  "Xenium_Broad_Reference_Validation.csv"
)
output_dir <- file.path(validation_root, "03_calibration", "tables")
outputs <- c(
  combined_scores = file.path(output_dir, "Xenium_Broad_All_Cluster_Scores.csv"),
  threshold_grid = file.path(output_dir, "Xenium_Broad_Reference_Threshold_Sensitivity.csv"),
  current_rule = file.path(output_dir, "Xenium_Broad_Current_Numeric_Threshold_Consensus.csv"),
  current_summary = file.path(output_dir, "Xenium_Broad_Current_Numeric_Threshold_Summary.csv"),
  review_template = file.path(output_dir, "Xenium_Broad_Reference_Threshold_Review_Template.csv"),
  provenance = file.path(output_dir, "Xenium_Broad_Calibration_Provenance.csv")
)

if (dry_run) {
  compact_dry_run(
    "Compile broad-reference calibration",
    inputs = c(reference_validation_path, input_grid$path),
    outputs = outputs,
    checks = c(
      expected_samples = length(sample_ids) == nrow(samples),
      expected_references = identical(references, c("Aldinger", "Sepp", "Science")),
      expected_tasks = nrow(input_grid) == length(sample_ids) * 3L
    )
  )
  quit(
    save = "no",
    status = if (all(file.exists(c(reference_validation_path, input_grid$path)))) 0L else 1L
  )
}

missing_inputs <- c(reference_validation_path, input_grid$path)[
  !file.exists(c(reference_validation_path, input_grid$path))
]
if (length(missing_inputs)) {
  stop("Missing calibration input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- outputs[file.exists(outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite broad calibration output(s):\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing the existing files."
  )
}

required_score_columns <- c(
  "sample_id", "reference", "seurat_clusters", "n_cells",
  "consensus_label", "existing_consensus_broad", "broad_winner",
  "broad_winner_route", "broad_winner_mean_score",
  "broad_top2_mean_score_margin", "broad_winner_cell_fraction",
  "broad_score_entropy", "broad_mean_score_vz", "broad_mean_score_rl",
  "broad_mean_score_other"
)
score_tables <- vector("list", nrow(input_grid))
for (i in seq_len(nrow(input_grid))) {
  table <- read.csv(
    input_grid$path[[i]], stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  missing <- setdiff(required_score_columns, names(table))
  if (length(missing)) {
    stop(input_grid$path[[i]], " lacks column(s): ", paste(missing, collapse = ", "))
  }
  if (any(as.character(table$sample_id) != input_grid$sample_id[[i]]) ||
      any(as.character(table$reference) != input_grid$reference[[i]])) {
    stop("Broad-score file identity does not match its expected task: ", input_grid$path[[i]])
  }
  key <- paste(table$sample_id, table$reference, table$seurat_clusters, sep = "\u001f")
  if (!nrow(table) || anyDuplicated(key)) {
    stop("Broad-score file is empty or contains duplicate keys: ", input_grid$path[[i]])
  }
  score_tables[[i]] <- table
}
all_scores <- do.call(rbind, score_tables)
rownames(all_scores) <- NULL
all_key <- paste(
  all_scores$sample_id, all_scores$reference, all_scores$seurat_clusters,
  sep = "\u001f"
)
if (anyDuplicated(all_key)) stop("Combined broad scores contain duplicate keys.")

population_key <- paste(all_scores$sample_id, all_scores$seurat_clusters, sep = "\u001f")
population_counts <- table(population_key)
if (any(population_counts != length(references))) {
  stop("At least one sample/cluster does not have exactly three reference score rows.")
}

wide_rows <- lapply(unique(population_key), function(key) {
  rows <- all_scores[population_key == key, , drop = FALSE]
  rows <- rows[match(references, rows$reference), , drop = FALSE]
  if (anyNA(rows$reference) || !identical(as.character(rows$reference), references)) {
    stop("A sample/cluster cannot be ordered across all three references.")
  }
  stable_columns <- c("sample_id", "seurat_clusters", "n_cells", "consensus_label", "existing_consensus_broad")
  for (column in stable_columns) {
    if (length(unique(as.character(rows[[column]]))) != 1L) {
      stop("References disagree on ", column, " for key ", key, ".")
    }
  }
  result <- rows[1, stable_columns, drop = FALSE]
  for (i in seq_along(references)) {
    prefix <- tolower(references[[i]])
    result[[paste0(prefix, "_winner")]] <- rows$broad_winner[[i]]
    result[[paste0(prefix, "_score")]] <- rows$broad_winner_mean_score[[i]]
    result[[paste0(prefix, "_margin")]] <- rows$broad_top2_mean_score_margin[[i]]
    result[[paste0(prefix, "_cell_fraction")]] <- rows$broad_winner_cell_fraction[[i]]
    result[[paste0(prefix, "_entropy")]] <- rows$broad_score_entropy[[i]]
  }
  result
})
wide <- do.call(rbind, wide_rows)
rownames(wide) <- NULL

score_grid <- seq(0.30, 0.80, by = 0.05)
margin_grid <- seq(0, 0.20, by = 0.025)
cell_fraction_grid <- c(0, 0.50, 0.60, 0.70)
threshold_grid <- expand.grid(
  min_broad_winner_score = score_grid,
  min_broad_top2_margin = margin_grid,
  min_broad_winner_cell_fraction = cell_fraction_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

calibration_rows <- vector("list", length(references) * nrow(threshold_grid))
row_index <- 0L
for (reference in references) {
  target <- tolower(reference)
  other <- tolower(setdiff(references, reference))
  target_winner <- as.character(wide[[paste0(target, "_winner")]])
  target_score <- as.numeric(wide[[paste0(target, "_score")]])
  target_margin <- as.numeric(wide[[paste0(target, "_margin")]])
  target_fraction <- as.numeric(wide[[paste0(target, "_cell_fraction")]])
  other_agree <- wide[[paste0(other[[1]], "_winner")]] ==
    wide[[paste0(other[[2]], "_winner")]]
  internal_target <- as.character(wide[[paste0(other[[1]], "_winner")]])
  established_known <- as.character(wide$existing_consensus_broad) != "Unknown"
  original_unknown <- !established_known

  for (j in seq_len(nrow(threshold_grid))) {
    rule <- threshold_grid[j, , drop = FALSE]
    pass <- target_score >= rule$min_broad_winner_score &
      target_margin >= rule$min_broad_top2_margin &
      target_fraction >= rule$min_broad_winner_cell_fraction
    eligible <- other_agree
    evaluated <- eligible & pass
    row_index <- row_index + 1L
    calibration_rows[[row_index]] <- data.frame(
      reference = reference,
      min_broad_winner_score = rule$min_broad_winner_score,
      min_broad_top2_margin = rule$min_broad_top2_margin,
      min_broad_winner_cell_fraction = rule$min_broad_winner_cell_fraction,
      n_all_clusters = nrow(wide),
      n_internal_target_clusters = sum(eligible),
      n_votes_retained = sum(pass),
      overall_vote_coverage = mean(pass),
      n_internal_target_votes_retained = sum(evaluated),
      internal_target_coverage = if (sum(eligible)) mean(pass[eligible]) else NA_real_,
      internal_target_concordance = if (sum(evaluated)) {
        mean(target_winner[evaluated] == internal_target[evaluated])
      } else NA_real_,
      internal_target_conflict_rate = if (sum(evaluated)) {
        mean(target_winner[evaluated] != internal_target[evaluated])
      } else NA_real_,
      known_vote_coverage = if (sum(established_known)) mean(pass[established_known]) else NA_real_,
      known_consensus_concordance = if (sum(pass & established_known)) {
        mean(target_winner[pass & established_known] ==
          wide$existing_consensus_broad[pass & established_known])
      } else NA_real_,
      unknown_vote_coverage = if (sum(original_unknown)) mean(pass[original_unknown]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
calibration <- do.call(rbind, calibration_rows)

# The current numeric thresholds are shown only as a diagnostic at the broad
# level. They have not been calibrated for summed broad probabilities.
current_score <- as.numeric(config$label_transfer$prediction_score_threshold)
current_margin <- as.numeric(
  config$label_transfer$resolution5_all_samples$min_top2_mean_score_margin
)
vote_matrix <- matrix(
  "abstain", nrow = nrow(wide), ncol = length(references),
  dimnames = list(NULL, tolower(references))
)
for (reference in references) {
  key <- tolower(reference)
  pass <- wide[[paste0(key, "_score")]] >= current_score &
    wide[[paste0(key, "_margin")]] >= current_margin
  vote_matrix[pass, key] <- as.character(wide[[paste0(key, "_winner")]])[pass]
  wide[[paste0(key, "_diagnostic_vote")]] <- vote_matrix[, key]
}

consensus_result <- lapply(seq_len(nrow(wide)), function(i) {
  valid <- vote_matrix[i, vote_matrix[i, ] != "abstain"]
  if (!length(valid)) {
    return(c(label = "Unknown", support = "0", nonabstain = "0", status = "all_references_abstained"))
  }
  counts <- table(valid)
  support <- max(counts)
  winners <- names(counts)[counts == support]
  if (support < 2L || length(winners) != 1L) {
    return(c(
      label = "Unknown", support = as.character(support),
      nonabstain = as.character(length(valid)), status = "insufficient_broad_support"
    ))
  }
  c(
    label = winners[[1]], support = as.character(support),
    nonabstain = as.character(length(valid)),
    status = if (support == 3L) "unanimous_3of3" else "majority_2of3"
  )
})
consensus_matrix <- do.call(rbind, consensus_result)
wide$diagnostic_broad_consensus <- consensus_matrix[, "label"]
wide$diagnostic_broad_support_n <- as.integer(consensus_matrix[, "support"])
wide$diagnostic_broad_nonabstain_n <- as.integer(consensus_matrix[, "nonabstain"])
wide$diagnostic_broad_status <- consensus_matrix[, "status"]
wide$diagnostic_matches_existing_known <- ifelse(
  wide$existing_consensus_broad == "Unknown", NA,
  wide$diagnostic_broad_consensus == wide$existing_consensus_broad
)

current_summary <- as.data.frame(
  table(
    existing_consensus_broad = wide$existing_consensus_broad,
    diagnostic_broad_consensus = wide$diagnostic_broad_consensus,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
names(current_summary)[[3]] <- "n_populations"
current_summary$current_score_threshold <- current_score
current_summary$current_margin_threshold <- current_margin
current_summary$interpretation <- paste(
  "diagnostic only; current fine-label numeric thresholds have not been",
  "approved for summed broad probabilities"
)

review_template <- data.frame(
  reference = references,
  min_broad_winner_score = NA_real_,
  min_broad_top2_margin = NA_real_,
  min_broad_winner_cell_fraction = NA_real_,
  approved = FALSE,
  review_notes = "",
  stringsAsFactors = FALSE
)
provenance <- data.frame(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  n_samples = length(sample_ids),
  n_references = length(references),
  n_cluster_reference_rows = nrow(all_scores),
  n_unique_sample_clusters = nrow(wide),
  calibration_target = "agreement of the other two raw broad winners; internal concordance target, not ground truth",
  score_grid = paste(score_grid, collapse = "|"),
  margin_grid = paste(margin_grid, collapse = "|"),
  cell_fraction_grid = paste(cell_fraction_grid, collapse = "|"),
  production_labels_modified = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(all_scores, outputs[["combined_scores"]], row.names = FALSE)
write.csv(calibration, outputs[["threshold_grid"]], row.names = FALSE)
write.csv(wide, outputs[["current_rule"]], row.names = FALSE)
write.csv(current_summary, outputs[["current_summary"]], row.names = FALSE)
write.csv(review_template, outputs[["review_template"]], row.names = FALSE, na = "")
write.csv(provenance, outputs[["provenance"]], row.names = FALSE)
message("Compiled broad-reference calibration evidence for ", nrow(wide), " sample/cluster populations.")
