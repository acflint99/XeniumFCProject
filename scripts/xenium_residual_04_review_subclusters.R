#!/usr/bin/env Rscript

# Generate non-spatial residual-subcluster review evidence and a blank label
# review template. Spatial plots remain a later validation step.

rm(list = ls())
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_04_review_subclusters.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_replicate_marker_evidence.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly=TRUE)
valid <- c("--dry-run", "--overwrite")
if (any(!args %in% valid)) stop("Usage: Rscript scripts/xenium_residual_04_review_subclusters.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
root <- here(config$project$outputs_dir, "xenium", "residual")
qc_state <- load_residual_qc_state(config)
input <- qc_state$input
out <- qc_state$markers
cluster_column <- "Xenium_snn_res.0.8"
outputs <- c(
  file.path(out, "Xenium_Residual_ReplicateAware_Marker_Evidence.csv"),
  file.path(out, "Xenium_Residual_Positive_Marker_Candidates.csv"),
  file.path(out, "Xenium_Residual_AntiMarker_Candidates.csv"),
  file.path(out, "Xenium_Residual_Marker_Sample_Coverage.csv"),
  file.path(out, "Xenium_Residual_ReplicateAware_Marker_Provenance.csv"),
  file.path(out, "Xenium_Residual_Subcluster_Label_Review_Template.csv")
)
if (dry_run) { compact_dry_run("Residual non-spatial subcluster review", inputs=input, outputs=outputs); quit(save="no", status=if(file.exists(input)) 0L else 1L) }
existing <- outputs[file.exists(outputs)]
if (length(existing) && !overwrite) stop("Refusing to overwrite residual review output(s):\n- ", paste(existing, collapse="\n- "))
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
object <- readRDS(input)
validate_residual_qc_object(object, qc_state)
build_replicate_marker_evidence("Residual", input, cluster_column, out, overwrite=overwrite)
if (!cluster_column %in% names(object[[]])) stop("Residual object lacks ", cluster_column, ".")
clusters <- unique(trimws(as.character(object[[cluster_column, drop=TRUE]])))
numeric_clusters <- suppressWarnings(as.numeric(clusters))
clusters <- clusters[order(is.na(numeric_clusters), numeric_clusters, clusters)]
template <- data.frame(
  cluster=clusters, final_label="", annotation_confidence="",
  positive_marker_support="", anti_marker_support="",
  cross_sample_support="", review_notes="", approved=FALSE,
  stringsAsFactors=FALSE
)
write.csv(template, outputs[[6]], row.names=FALSE)
message("Residual marker evidence and label-review template are ready. No labels were applied.")
