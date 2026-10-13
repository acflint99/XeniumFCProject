#!/usr/bin/env Rscript

# Re-express existing per-class transfer probabilities at the VZ/RL/Other
# level for every whole-tissue cluster. This diagnostic does not rerun label
# transfer and does not change consensus labels, reference votes, or routes.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_broad_reference.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
references <- broad_reference_names()
taxonomy <- load_reference_broad_taxonomy(config)
validate_taxonomy_against_config(taxonomy, config)
production_paths <- get_production_output_paths(config)

task_map <- expand.grid(
  sample_id = sample_ids,
  reference = references,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
task_map <- task_map[order(
  match(task_map$sample_id, sample_ids),
  match(task_map$reference, references)
), , drop = FALSE]
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c("task_id", "sample_id", "reference")]
if (nrow(task_map) != length(sample_ids) * length(references)) {
  stop("Broad-score task map has an unexpected size.")
}

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
positional <- args[!args %in% valid_options]

if (list_requested) {
  if (length(positional)) stop("--list does not accept TASK_ID.")
  cat(
    "Scope:", length(sample_ids), "samples x", length(references),
    "references =", nrow(task_map), "tasks\n"
  )
  cat("Broad categories: VZ, RL, Other\n")
  cat("This stage reads existing annotated RDS objects; it does not rerun TransferData.\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

if (length(positional) != 1L) {
  stop(
    "Usage: Rscript scripts/xenium_broad_02_score_clusters.R ",
    "[--list|--dry-run|--overwrite] TASK_ID"
  )
}
task_id <- suppressWarnings(as.integer(positional[[1]]))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ".")
}

task <- task_map[task_id, , drop = FALSE]
sample_id <- task$sample_id[[1]]
reference <- task$reference[[1]]
reference_key <- tolower(reference)
annotated_path <- file.path(
  production_paths$label_transfer_dir, reference_key, "rds",
  paste0(sample_id, "_", reference, "_annotated.rds")
)
consensus_path <- file.path(
  production_paths$consensus_table_dir,
  paste0(sample_id, "_comparison_merged.csv")
)
validation_path <- file.path(
  here(config$project$outputs_dir), "validation", "broad_reference_workflow",
  "01_reference_validation", "tables", "Xenium_Broad_Reference_Validation.csv"
)
output_dir <- file.path(
  here(config$project$outputs_dir), "validation", "broad_reference_workflow",
  "02_cluster_scores", "per_task"
)
output_path <- file.path(
  output_dir, paste0(sample_id, "_", reference, "_Broad_Cluster_Scores.csv")
)
provenance_path <- file.path(
  output_dir, paste0(sample_id, "_", reference, "_Broad_Cluster_Scores_Provenance.csv")
)

if (dry_run) {
  compact_dry_run(
    paste0(
      "Broad cluster scores task ", task_id, "/", nrow(task_map),
      " [", sample_id, " / ", reference, "]"
    ),
    inputs = c(annotated_path, consensus_path, validation_path, attr(taxonomy, "path")),
    outputs = c(output_path, provenance_path)
  )
  quit(
    save = "no",
    status = if (all(file.exists(c(
      annotated_path, consensus_path, validation_path, attr(taxonomy, "path")
    )))) 0L else 1L
  )
}

inputs <- c(annotated_path, consensus_path, validation_path, attr(taxonomy, "path"))
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("Missing broad-score input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- c(output_path, provenance_path)[
  file.exists(c(output_path, provenance_path))
]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite broad-score output(s):\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing the existing files."
  )
}

suppressPackageStartupMessages(library(Seurat))

object <- readRDS(annotated_path)
if (anyDuplicated(Cells(object))) stop("Annotated object contains duplicate cell IDs.")
if (!identical(rownames(object[[]]), Cells(object))) {
  stop("Annotated object metadata does not align exactly with cell IDs.")
}
score_summary <- summarize_cluster_broad_scores(object, taxonomy)

consensus <- read.csv(
  consensus_path, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = character()
)
required_consensus <- c(
  "seurat_clusters", "consensus_label", paste0(reference_key, "_reference_vote")
)
missing_consensus <- setdiff(required_consensus, names(consensus))
if (length(missing_consensus)) {
  stop("Consensus table lacks column(s): ", paste(missing_consensus, collapse = ", "))
}
consensus$seurat_clusters <- trimws(as.character(consensus$seurat_clusters))
if (anyNA(consensus$seurat_clusters) || any(!nzchar(consensus$seurat_clusters)) ||
    anyDuplicated(consensus$seurat_clusters)) {
  stop("Consensus table contains invalid or duplicate cluster IDs.")
}
if (!setequal(score_summary$seurat_clusters, consensus$seurat_clusters)) {
  stop("Annotated object and consensus table contain different cluster sets.")
}
index <- match(score_summary$seurat_clusters, consensus$seurat_clusters)
score_summary$sample_id <- sample_id
score_summary$reference <- reference
score_summary$consensus_label <- trimws(as.character(consensus$consensus_label[index]))
score_summary$existing_reference_vote <- trimws(as.character(
  consensus[[paste0(reference_key, "_reference_vote")]][index]
))
if (anyNA(score_summary$consensus_label) || any(!nzchar(score_summary$consensus_label))) {
  stop("Consensus table contains blank or missing consensus labels.")
}

score_summary$existing_consensus_broad <- "Unknown"
known <- !is_unknown_consensus_label(score_summary$consensus_label)
known_index <- match(score_summary$consensus_label[known], taxonomy$fine_label)
if (anyNA(known_index)) {
  stop(
    "Known consensus label(s) are absent from the broad taxonomy: ",
    paste(sort(unique(score_summary$consensus_label[known][is.na(known_index)])), collapse = ", ")
  )
}
score_summary$existing_consensus_broad[known] <- taxonomy$broad_label[known_index]

score_summary$existing_reference_broad <- "Abstain"
existing_known <- !is.na(score_summary$existing_reference_vote) &
  nzchar(score_summary$existing_reference_vote) &
  tolower(score_summary$existing_reference_vote) != "unknown"
vote_index <- match(score_summary$existing_reference_vote[existing_known], taxonomy$fine_label)
if (anyNA(vote_index)) {
  stop("Existing reference vote contains a label absent from the broad taxonomy.")
}
score_summary$existing_reference_broad[existing_known] <- taxonomy$broad_label[vote_index]
score_summary$broad_winner_matches_existing_consensus <- ifelse(
  known, score_summary$broad_winner == score_summary$existing_consensus_broad, NA
)
score_summary$broad_winner_matches_existing_reference <- ifelse(
  existing_known, score_summary$broad_winner == score_summary$existing_reference_broad, NA
)

score_summary <- score_summary[c(
  "sample_id", "reference", "seurat_clusters", "n_cells",
  "consensus_label", "existing_consensus_broad", "existing_reference_vote",
  "existing_reference_broad", "broad_winner", "broad_winner_route",
  "broad_winner_mean_score", "broad_runner_up",
  "broad_runner_up_mean_score", "broad_top2_mean_score_margin",
  "broad_winner_cell_fraction", "broad_score_entropy",
  "broad_mean_score_vz", "broad_mean_score_rl", "broad_mean_score_other",
  "broad_winner_matches_existing_consensus",
  "broad_winner_matches_existing_reference"
)]
if (anyDuplicated(paste(score_summary$sample_id, score_summary$seurat_clusters))) {
  stop("Broad-score output contains duplicate sample/cluster keys.")
}

annotated_info <- file.info(annotated_path)
consensus_info <- file.info(consensus_path)
provenance <- data.frame(
  sample_id = sample_id,
  reference = reference,
  task_id = task_id,
  total_tasks = nrow(task_map),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  annotated_rds = normalizePath(annotated_path, mustWork = TRUE),
  annotated_rds_size = as.numeric(annotated_info$size),
  annotated_rds_mtime = as.numeric(annotated_info$mtime),
  consensus_table = normalizePath(consensus_path, mustWork = TRUE),
  consensus_table_size = as.numeric(consensus_info$size),
  consensus_table_mtime = as.numeric(consensus_info$mtime),
  taxonomy_path = attr(taxonomy, "path"),
  score_method = "sum cell-level fine-label probabilities within VZ/RL/Other then average by cluster",
  production_labels_modified = FALSE,
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(score_summary, output_path, row.names = FALSE)
write.csv(provenance, provenance_path, row.names = FALSE)
message("Saved broad cluster-score audit: ", output_path)
