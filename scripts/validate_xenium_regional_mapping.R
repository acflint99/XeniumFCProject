#!/usr/bin/env Rscript
# Lightweight audit of regional mapping manifests; no large RDS is loaded.
# Expected runtime: under one minute, depending on filesystem metadata latency.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/validate_xenium_regional_mapping.R")
source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()
sample_ids <- as.character(load_sample_manifest(config)$sample_id)
args <- commandArgs(trailingOnly = TRUE)
branches <- if (length(args)) toupper(args) else c("VZ", "RL", "RESIDUAL")
if (any(!branches %in% c("VZ", "RL", "RESIDUAL")) || anyDuplicated(branches)) {
  stop("Usage: Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R [VZ RL RESIDUAL]")
}
root <- here(config$project$outputs_dir, "xenium")
for (branch in branches) {
  region <- tolower(branch)
  if (branch == "RESIDUAL") {
    manifest_path <- file.path(root, region, "06_mapped", "tables", "Xenium_Residual_Mapping_Manifest.csv")
    master_path <- file.path(root, region, "05_subclusters", "rds", "Xenium_Residual_Subclusters.rds")
    label_path <- file.path(root, region, "05_subclusters", "tables", "Xenium_Residual_Subcluster_Labels.csv")
    count_column <- "mapped_residual_cells"
  } else {
    manifest_path <- file.path(root, region, "07_mapped", "rds", paste0("Xenium_", branch, "_Mapping_manifest.csv"))
    master_path <- file.path(root, region, "06_subclusters", "rds", paste0("Xenium_", branch, "_subclusters.rds"))
    label_path <- file.path(root, region, "06_subclusters", "tables", paste0("Xenium_", branch, "_subcluster_labels.csv"))
    count_column <- paste0("mapped_", region, "_cells")
  }
  required_files <- c(manifest_path, master_path, label_path)
  if (any(!file.exists(required_files))) stop(branch, " is missing:\n- ", paste(required_files[!file.exists(required_files)], collapse = "\n- "))
  mapping <- read.csv(manifest_path, stringsAsFactors = FALSE)
  required_columns <- c("sample_id", "output_path", count_column)
  if (!all(required_columns %in% names(mapping))) stop(branch, " mapping manifest lacks required columns.")
  if (!identical(as.character(mapping$sample_id), sample_ids) || anyDuplicated(mapping$sample_id)) {
    stop(branch, " mapping sample IDs/order do not match the active manifest.")
  }
  paths <- as.character(mapping$output_path)
  if (anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) || any(!file.exists(paths))) {
    stop(branch, " mapped output paths are missing, duplicated, or invalid.")
  }
  counts <- suppressWarnings(as.numeric(mapping[[count_column]]))
  if (any(!is.finite(counts)) || any(counts < 0 | counts != floor(counts))) stop(branch, " has invalid mapped-cell counts.")
  output_info <- file.info(paths)
  source_time <- max(file.info(c(master_path, label_path))$mtime)
  if (anyNA(output_info$size) || any(output_info$size <= 0) ||
      any(output_info$mtime < source_time) || file.info(manifest_path)$mtime < source_time) {
    stop(branch, " mapped outputs or manifest predate their label source or are empty; mapping needs review/regeneration.")
  }
  if (branch == "VZ") {
    approved <- c("Patterning PCs", "BG", "Maturing PCs", "MLIs", "RG Progenitors", "Early-born PCs", "Late-born PCs", "eCN", "Pre-OPCs", "GABA Progenitors", "Cycling Glia", "iCN", "Astrocytes/Ependyma", "Cycling Glia", "RL VZ", "Maturing PCs", "Golgi Cells", "Differentiating GCs", "OPCs", "VZPs", "GABA Interneurons", "iCN", "Unresolved", "Patterning PCs", "Maturing PCs")
    labels <- read.csv(label_path, stringsAsFactors = FALSE)
    if (!all(c("cluster_id", "VZ_subcluster") %in% names(labels)) ||
        anyDuplicated(labels$cluster_id) || !setequal(as.character(labels$cluster_id), as.character(0:24))) {
      stop("VZ label manifest does not contain the approved 0-24 cluster map.")
    }
    if (!identical(as.character(labels$VZ_subcluster[match(as.character(0:24), as.character(labels$cluster_id))]), approved) ||
        sum(counts) != 1447477L) {
      stop("VZ mapping does not match the approved VZ labels/cell total (cluster 17 = Differentiating GCs).")
    }
  }
  cat("PASS | ", branch, " | samples=", length(sample_ids), "/", length(sample_ids),
      " | mapped cells=", format(sum(counts), scientific = FALSE),
      " | output paths and timestamps checked; RDS contents not loaded\n", sep = "")
}
