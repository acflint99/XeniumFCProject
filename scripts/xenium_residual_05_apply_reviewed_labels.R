#!/usr/bin/env Rscript

# Apply a complete, explicitly approved residual-cluster label manifest. Labels
# may remain Unknown when marker evidence does not support a defensible winner.

rm(list = ls())
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_05_apply_reviewed_labels.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly=TRUE)
valid <- c("--dry-run", "--overwrite")
if (any(!args %in% valid)) stop("Usage: Rscript scripts/xenium_residual_05_apply_reviewed_labels.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
root <- here(config$project$outputs_dir, "xenium", "residual")
qc_state <- load_residual_analysis_state(config)
input <- qc_state$input
evidence <- file.path(qc_state$markers, "Xenium_Residual_ReplicateAware_Marker_Evidence.csv")
positive <- file.path(qc_state$markers, "Xenium_Residual_Positive_Marker_Candidates.csv")
anti <- file.path(qc_state$markers, "Xenium_Residual_AntiMarker_Candidates.csv")
marker_provenance <- file.path(qc_state$markers, "Xenium_Residual_ReplicateAware_Marker_Provenance.csv")
review_path <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels, config)
output <- file.path(root, "05_subclusters", "rds", "Xenium_Residual_Subclusters.rds")
manifest <- file.path(root, "05_subclusters", "tables", "Xenium_Residual_Subcluster_Labels.csv")
inputs <- c(input, evidence, positive, anti, marker_provenance, review_path)
if (dry_run) { compact_dry_run("Apply reviewed residual subcluster labels", inputs=inputs, outputs=c(output, manifest)); quit(save="no", status=if(all(file.exists(inputs))) 0L else 1L) }
missing <- inputs[!file.exists(inputs)]
if (length(missing)) stop("Missing residual label input(s):\n- ", paste(missing, collapse="\n- "))
existing <- c(output, manifest)[file.exists(c(output, manifest))]
if (length(existing) && !overwrite) stop("Refusing to overwrite residual label output(s):\n- ", paste(existing, collapse="\n- "))
suppressPackageStartupMessages(library(Seurat))
object <- readRDS(input)
validate_residual_qc_object(object, qc_state)
marker_source <- read.csv(marker_provenance, stringsAsFactors = FALSE)
if (nrow(marker_source) != 1L ||
    !all(c("input_rds", "input_size", "input_mtime") %in% names(marker_source)) ||
    !identical(normalizePath(marker_source$input_rds, mustWork = TRUE), normalizePath(input, mustWork = TRUE)) ||
    !residual_qc_file_matches(input, marker_source$input_size, marker_source$input_mtime)) {
  stop("Residual marker evidence does not match the current analysis input.")
}
cluster_column <- "residual_postQC_Xenium_snn_res.0.8"
if (!cluster_column %in% names(object[[]])) stop("Residual object lacks ", cluster_column, ".")
clusters <- unique(trimws(as.character(object[[cluster_column, drop=TRUE]])))
review <- load_reviewed_residual_subcluster_labels(config)
required <- c("cluster", "final_label", "annotation_confidence", "positive_marker_support", "anti_marker_support", "cross_sample_support", "review_notes", "approved")
if (!all(required %in% names(review))) stop("Reviewed residual-label manifest is incomplete.")
if (!nrow(review) || !all(c("source_input_md5", "source_cluster_column") %in% names(review)) ||
    anyNA(review$source_input_md5) || anyNA(review$source_cluster_column) ||
    any(review$source_input_md5 != unname(tools::md5sum(input))) ||
    any(review$source_cluster_column != cluster_column)) {
  stop("Generate and approve a fresh residual label template bound to the current post-QC input; old cluster IDs cannot be reused.")
}
review$cluster <- trimws(as.character(review$cluster)); review$final_label <- trimws(as.character(review$final_label))
review$approved <- suppressWarnings(as.logical(review$approved))
if (!nrow(review) || anyNA(review$approved) || any(!review$approved)) stop("Every residual cluster requires an explicitly approved review row.")
if (anyDuplicated(review$cluster) || !setequal(review$cluster, clusters)) stop("Reviewed residual-label clusters do not exactly match the integrated object.")
if (any(vapply(review[setdiff(required, "approved")], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1)))) stop("Approved residual-label rows contain blank evidence or notes.")
index <- match(as.character(object[[cluster_column, drop=TRUE]]), review$cluster)
if (anyNA(index)) stop("A residual cell cannot be matched to its reviewed cluster.")
object$residual_subcluster <- review$final_label[index]
object$residual_annotation_confidence <- review$annotation_confidence[index]
object$residual_annotation_status <- ifelse(tolower(review$final_label[index]) == "unknown", "retained_unknown", "reviewed_label")
Idents(object) <- "residual_subcluster"
input_info <- file.info(input); review_info <- file.info(review_path)
review$source_cluster_column <- cluster_column
review$input_rds <- normalizePath(input, mustWork=TRUE)
review$input_size <- as.numeric(input_info$size); review$input_mtime <- as.numeric(input_info$mtime)
review$review_manifest <- normalizePath(review_path, mustWork=TRUE)
review$review_manifest_size <- as.numeric(review_info$size); review$review_manifest_mtime <- as.numeric(review_info$mtime)
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE); dir.create(dirname(manifest), recursive=TRUE, showWarnings=FALSE)
saveRDS(object, output, compress=FALSE); write.csv(review, manifest, row.names=FALSE)
message("Applied reviewed residual labels; retained Unknown is permitted when documented.")
