#!/usr/bin/env Rscript
# Apply the explicitly approved pre-QC cluster decision; retain cluster IDs/embedding.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_04b_apply_qc.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite"))) stop("Usage: xenium_residual_04b_apply_qc.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
p <- residual_qc_removal_paths(config)
inputs <- unlist(p[c("source", "decision", "qc_cells", "qc_provenance")])
outputs <- unlist(p[c("filtered", "removed", "manifest", "provenance")])
if (dry_run) {
  compact_dry_run("Approved residual QC removal (paths/decision CSV only)", inputs, outputs)
  if (file.exists(p$decision)) {
    d <- read_residual_qc_decision(p$decision)
    cat("  Remove cluster(s): ", paste(d$cluster, collapse = ", "),
        "; expected cells: ", d$expected_input_cells[[1]], " -> ",
        d$expected_input_cells[[1]] - sum(d$expected_removed_cells), "\n", sep = "")
  }
  quit(save = "no", status = if (all(file.exists(inputs))) 0L else 1L)
}
if (any(!file.exists(inputs))) stop("Missing QC decision inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
if (any(file.exists(outputs)) && !overwrite) stop("Refusing to overwrite QC-removal outputs; review before --overwrite.")
d <- read_residual_qc_decision(p$decision)
decision_digest <- unname(tools::md5sum(p$decision))
if (!residual_qc_file_matches(p$source, d$source_size[[1]], d$source_mtime[[1]])) {
  stop("Integrated residual object differs from the approved source snapshot. Do not apply the old cluster number.")
}
ids <- as.character(load_sample_manifest(config)$sample_id)
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid manifest sample IDs.")
provenance <- read.csv(p$qc_provenance, stringsAsFactors = FALSE)
source_row <- provenance[provenance$input_path == normalizePath(p$source, mustWork = TRUE), , drop = FALSE]
if (nrow(source_row) != 1L || !identical(as.character(source_row$cluster_column), d$cluster_column[[1]]) ||
    !residual_qc_file_matches(p$source, source_row$input_size, source_row$input_mtime)) {
  stop("QC evidence is not current for the approved source and clustering.")
}
suppressPackageStartupMessages(library(SeuratObject))
object <- readRDS(p$source)
metadata <- object[[]]
column <- d$cluster_column[[1]]
if (!all(c("orig.ident", "consensus_label", "residual_reference_route", column) %in% names(metadata)) ||
    !identical(colnames(object), rownames(metadata)) || anyDuplicated(colnames(object)) ||
    ncol(object) != d$expected_input_cells[[1]]) stop("Invalid pre-QC metadata or cell count.")
qc <- read.csv(p$qc_cells, stringsAsFactors = FALSE)
if (!all(c("cell_id", "sample_id", "cluster") %in% names(qc)) || anyDuplicated(qc$cell_id) ||
    !setequal(qc$cell_id, colnames(object))) stop("QC per-cell IDs differ from the approved source.")
qc <- qc[match(colnames(object), qc$cell_id), , drop = FALSE]
sample_ids <- as.character(metadata$orig.ident)
clusters <- as.character(metadata[[column]])
if (anyNA(sample_ids) || any(!sample_ids %in% ids) || anyNA(clusters) || any(!nzchar(clusters)) ||
    !identical(as.character(qc$sample_id), sample_ids) || !identical(as.character(qc$cluster), clusters) ||
    any(!is_unknown_consensus_label(metadata$consensus_label)) ||
    anyNA(metadata$residual_reference_route) ||
    any(!metadata$residual_reference_route %in% c("route_other", "route_unresolved"))) {
  stop("QC sample/cluster identities or residual routes do not match the source.")
}
remove <- clusters %in% d$cluster
for (i in seq_len(nrow(d))) {
  if (sum(clusters == d$cluster[[i]]) != d$expected_removed_cells[[i]]) stop("Approved removal count differs for cluster ", d$cluster[[i]])
}
removed <- data.frame(cell_id = colnames(object)[remove], sample_id = sample_ids[remove],
                       cluster = clusters[remove], consensus_label = as.character(metadata$consensus_label[remove]),
                       review_notes = d$review_notes[match(clusters[remove], d$cluster)], stringsAsFactors = FALSE)
for (i in seq_len(nrow(removed))) {
  if (!startsWith(removed$cell_id[[i]], paste0(removed$sample_id[[i]], "_"))) stop("Invalid sample-qualified cell ID.")
}
manifest <- data.frame(sample_id = ids,
  input_cells = as.integer(table(factor(sample_ids, levels = ids))),
  removed_cells = as.integer(table(factor(removed$sample_id, levels = ids))))
manifest$remaining_cells <- manifest$input_cells - manifest$removed_cells
keep_ids <- colnames(object)[!remove]
object <- subset(object, cells = keep_ids)
if (!identical(colnames(object), keep_ids) ||
    !identical(as.character(object[[column, drop = TRUE]]), clusters[!remove])) {
  stop("Filtering changed surviving IDs/order or cluster labels.")
}
digest <- unname(tools::md5sum(p$decision))
if (!identical(digest, decision_digest) ||
    !residual_qc_file_matches(p$source, d$source_size[[1]], d$source_mtime[[1]])) {
  stop("Approved decision or source changed during removal.")
}
object@misc$residual_qc_removal <- list(decision_md5 = digest,
  source_size = d$source_size[[1]], source_mtime = d$source_mtime[[1]],
  removed_clusters = d$cluster, cells_before = nrow(metadata), removed_cells = nrow(removed),
  cells_after = ncol(object), reintegrated = FALSE, cluster_column = column)
for (path in outputs) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
saveRDS(object, p$filtered, compress = FALSE)
write.csv(removed, p$removed, row.names = FALSE)
write.csv(manifest, p$manifest, row.names = FALSE)
filtered_info <- file.info(p$filtered)
write.csv(data.frame(created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  decision_md5 = digest, source_size = d$source_size[[1]], source_mtime = d$source_mtime[[1]],
  filtered_size = as.numeric(filtered_info$size), filtered_mtime = as.numeric(filtered_info$mtime),
  cells_before = nrow(metadata), removed_cells = nrow(removed), cells_after = ncol(object),
  cluster_column = column, removed_clusters = paste(d$cluster, collapse = "|"),
  reintegrated = FALSE, embedding = "pre-QC coordinates retained; no reclustering"),
  p$provenance, row.names = FALSE)
state <- load_residual_qc_state(config)
validate_residual_qc_object(object, state)
message("Approved residual QC removal complete: ", nrow(removed), " removed; ", ncol(object),
        " retained. Original integrated RDS preserved. Regenerate post-QC marker review before labels/mapping.")
