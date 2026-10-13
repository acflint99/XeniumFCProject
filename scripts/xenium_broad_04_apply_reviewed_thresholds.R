#!/usr/bin/env Rscript

# Apply explicitly reviewed per-reference broad thresholds to the compiled
# VZ/RL/Other scores. This writes versioned tables only and preserves the
# established detailed consensus labels for auditing and rollback.

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
  stop("Usage: Rscript scripts/xenium_broad_04_apply_reviewed_thresholds.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
thresholds <- load_reviewed_broad_thresholds(
  config, require_approved = !dry_run
)

score_path <- file.path(
  here(config$project$outputs_dir), "validation", "broad_reference_workflow",
  "03_calibration", "tables", "Xenium_Broad_All_Cluster_Scores.csv"
)
threshold_path <- attr(thresholds, "path")
output_dir <- file.path(
  here(config$project$outputs_dir), "xenium", "annotation",
  "05_broad_reference_consensus", "tables"
)
per_sample_dir <- file.path(output_dir, "per_sample")
outputs <- c(
  combined = file.path(output_dir, "Xenium_Broad_Reference_Consensus.csv"),
  summary = file.path(output_dir, "Xenium_Broad_Reference_Consensus_Summary.csv"),
  conflict = file.path(output_dir, "Xenium_Broad_Known_Consensus_Conflicts.csv"),
  provenance = file.path(output_dir, "Xenium_Broad_Reference_Consensus_Provenance.csv"),
  file.path(per_sample_dir, paste0(sample_ids, "_Broad_Reference_Consensus.csv"))
)

if (dry_run) {
  compact_dry_run(
    "Apply reviewed broad-reference thresholds",
    inputs = c(score_path, threshold_path),
    outputs = outputs,
    checks = c(
      all_thresholds_approved = all(thresholds$approved),
      thresholds_complete = all(is.finite(as.matrix(thresholds[c(
        "min_broad_winner_score", "min_broad_top2_margin",
        "min_broad_winner_cell_fraction"
      )]))),
      three_reference_rows = nrow(thresholds) == 3L
    )
  )
  quit(
    save = "no",
    status = if (
      file.exists(score_path) && all(thresholds$approved) &&
        all(is.finite(as.matrix(thresholds[c(
          "min_broad_winner_score", "min_broad_top2_margin",
          "min_broad_winner_cell_fraction"
        )])))
    ) 0L else 1L
  )
}

if (!file.exists(score_path)) stop("Combined broad-score table not found: ", score_path)
existing_outputs <- outputs[file.exists(outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite reviewed broad-consensus output(s):\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing and archiving the existing outputs."
  )
}

scores <- read.csv(
  score_path, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = character()
)
required <- c(
  "sample_id", "reference", "seurat_clusters", "n_cells",
  "consensus_label", "existing_consensus_broad", "broad_winner",
  "broad_winner_mean_score", "broad_top2_mean_score_margin",
  "broad_winner_cell_fraction"
)
missing <- setdiff(required, names(scores))
if (length(missing)) stop("Combined broad scores lack: ", paste(missing, collapse = ", "))
if (any(!scores$sample_id %in% sample_ids) ||
    any(!scores$reference %in% references)) {
  stop("Combined broad scores contain unexpected sample or reference names.")
}
key <- paste(scores$sample_id, scores$reference, scores$seurat_clusters, sep = "\u001f")
if (anyDuplicated(key)) stop("Combined broad scores contain duplicate keys.")

thresholds <- thresholds[match(references, thresholds$reference), , drop = FALSE]
rows <- vector("list", length(unique(paste(scores$sample_id, scores$seurat_clusters))))
population_keys <- unique(paste(scores$sample_id, scores$seurat_clusters, sep = "\u001f"))
for (i in seq_along(population_keys)) {
  population <- scores[
    paste(scores$sample_id, scores$seurat_clusters, sep = "\u001f") == population_keys[[i]],
    , drop = FALSE
  ]
  population <- population[match(references, population$reference), , drop = FALSE]
  if (nrow(population) != 3L || anyNA(population$reference)) {
    stop("Population does not contain exactly three ordered references: ", population_keys[[i]])
  }
  stable <- c("sample_id", "seurat_clusters", "n_cells", "consensus_label", "existing_consensus_broad")
  if (any(vapply(stable, function(column) {
    length(unique(as.character(population[[column]]))) != 1L
  }, logical(1)))) {
    stop("References disagree on stable population metadata: ", population_keys[[i]])
  }
  passed <- population$broad_winner_mean_score >= thresholds$min_broad_winner_score &
    population$broad_top2_mean_score_margin >= thresholds$min_broad_top2_margin &
    population$broad_winner_cell_fraction >= thresholds$min_broad_winner_cell_fraction
  votes <- ifelse(passed, as.character(population$broad_winner), "Abstain")
  valid <- votes[votes != "Abstain"]
  if (!length(valid)) {
    broad_label <- "Unknown"
    support <- 0L
    status <- "all_references_abstained"
  } else {
    counts <- table(valid)
    support <- as.integer(max(counts))
    winners <- names(counts)[counts == support]
    if (support < 2L || length(winners) != 1L) {
      broad_label <- "Unknown"
      status <- "insufficient_broad_support"
    } else {
      broad_label <- winners[[1]]
      status <- if (support == 3L) "unanimous_3of3" else "majority_2of3"
    }
  }
  row <- population[1, stable, drop = FALSE]
  for (j in seq_along(references)) {
    prefix <- tolower(references[[j]])
    row[[paste0(prefix, "_broad_winner")]] <- population$broad_winner[[j]]
    row[[paste0(prefix, "_broad_score")]] <- population$broad_winner_mean_score[[j]]
    row[[paste0(prefix, "_broad_margin")]] <- population$broad_top2_mean_score_margin[[j]]
    row[[paste0(prefix, "_broad_cell_fraction")]] <- population$broad_winner_cell_fraction[[j]]
    row[[paste0(prefix, "_broad_vote")]] <- votes[[j]]
    row[[paste0(prefix, "_broad_abstained")]] <- !passed[[j]]
  }
  row$broad_reference_consensus <- broad_label
  row$broad_reference_route <- if (broad_label == "Unknown") {
    "route_unresolved"
  } else {
    broad_category_to_route(broad_label)
  }
  row$broad_reference_support_n <- support
  row$broad_reference_nonabstain_n <- length(valid)
  row$broad_reference_status <- status
  row$broad_matches_existing_known <- if (
    row$existing_consensus_broad[[1]] == "Unknown"
  ) NA else broad_label == row$existing_consensus_broad[[1]]
  rows[[i]] <- row
}
result <- do.call(rbind, rows)
rownames(result) <- NULL
result$sample_id <- factor(result$sample_id, levels = sample_ids)
cluster_number <- suppressWarnings(as.numeric(as.character(result$seurat_clusters)))
result <- result[order(result$sample_id, is.na(cluster_number), cluster_number, result$seurat_clusters), ]
result$sample_id <- as.character(result$sample_id)

summary_table <- as.data.frame(
  table(
    existing_consensus_broad = result$existing_consensus_broad,
    broad_reference_consensus = result$broad_reference_consensus,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
names(summary_table)[[3]] <- "n_populations"
conflicts <- result[
  !is.na(result$broad_matches_existing_known) & !result$broad_matches_existing_known,
  , drop = FALSE
]
score_info <- file.info(score_path)
threshold_info <- file.info(threshold_path)
provenance <- data.frame(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  combined_score_path = normalizePath(score_path, mustWork = TRUE),
  combined_score_size = as.numeric(score_info$size),
  combined_score_mtime = as.numeric(score_info$mtime),
  reviewed_threshold_path = threshold_path,
  reviewed_threshold_size = as.numeric(threshold_info$size),
  reviewed_threshold_mtime = as.numeric(threshold_info$mtime),
  consensus_method = "per-reference reviewed broad abstention followed by broad 2-of-3",
  original_consensus_label_preserved = TRUE,
  stringsAsFactors = FALSE
)

dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, outputs[["combined"]], row.names = FALSE)
write.csv(summary_table, outputs[["summary"]], row.names = FALSE)
write.csv(conflicts, outputs[["conflict"]], row.names = FALSE)
write.csv(provenance, outputs[["provenance"]], row.names = FALSE)
for (sample_id in sample_ids) {
  sample_table <- result[result$sample_id == sample_id, , drop = FALSE]
  if (!nrow(sample_table)) stop("Broad consensus has no rows for ", sample_id, ".")
  write.csv(
    sample_table,
    file.path(per_sample_dir, paste0(sample_id, "_Broad_Reference_Consensus.csv")),
    row.names = FALSE
  )
}
message("Applied reviewed broad thresholds without replacing detailed consensus labels.")
