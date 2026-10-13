#!/usr/bin/env Rscript

# Memory-conscious spatial merge of all 34 manifest-defined combined objects.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_spatial_01_merge.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
rds_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "01_merged", "rds")
table_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "01_merged", "tables")
input_names <- paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds")
input_paths <- file.path(input_dir, input_names)
merged_path <- file.path(rds_dir, "Xenium_VZRL_spatial_merged.rds")
manifest_path <- file.path(table_dir, "Xenium_VZRL_spatial_merged_manifest.csv")

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  compact_dry_run(
    "Spatial regional merge",
    inputs = input_paths,
    outputs = c(merged_path, manifest_path),
    checks = c(no_unexpected_inputs = !length(unexpected_names))
  )
  if (length(unexpected_names)) cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  quit(save = "no", status = 0L)
}

if (length(missing_names)) {
  stop("Cannot spatially merge; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot spatially merge; unexpected RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
existing_outputs <- c(merged_path, manifest_path)[file.exists(c(merged_path, manifest_path))]
if (length(existing_outputs) && !overwrite) {
  stop("Refusing to overwrite spatial merge outputs:\n- ",
       paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
}

options(future.globals.maxSize = Inf)
plan(sequential)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
merge_manifest <- data.frame(
  sample_id = sample_ids, input_path = input_paths,
  cells = integer(length(sample_ids)), images = integer(length(sample_ids)),
  PCW = character(length(sample_ids)),
  vz_cells = integer(length(sample_ids)), rl_cells = integer(length(sample_ids)),
  residual_cells = integer(length(sample_ids)),
  residual_retained_unknown_cells = integer(length(sample_ids)),
  unassigned_consensus_cells = integer(length(sample_ids)),
  stringsAsFactors = FALSE
)
merged_obj <- NULL

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Reading spatial object ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])
  required_metadata <- c(
    "comb_subcluster", "consensus_label", "VZ_subcluster", "RL_subcluster",
    "VZ_subset_reason", "RL_subset_reason", "residual_subcluster",
    "residual_reference_route", "residual_subset_reason",
    "residual_annotation_confidence", "residual_annotation_status",
    "regional_assignment", "PCW"
  )
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (ncol(obj) == 0L) stop(sample_id, " contains zero cells.")
  if (!length(Images(obj))) stop(sample_id, " has no spatial image/FOV data.")
  pcw_values <- unique(as.character(obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")
  branch_count <- as.integer(!is.na(obj$VZ_subcluster)) +
    as.integer(!is.na(obj$RL_subcluster)) +
    as.integer(!is.na(obj$residual_subcluster))
  if (any(branch_count > 1L)) {
    stop(sample_id, " contains an overlap among VZ, RL, and residual branches.")
  }
  assignment <- as.character(obj$regional_assignment)
  allowed_assignments <- c("VZ", "RL", "Residual", "unassigned_consensus")
  if (anyNA(assignment) || any(!assignment %in% allowed_assignments)) {
    stop(sample_id, " has missing or invalid regional_assignment values.")
  }
  expected_assignment <- rep("unassigned_consensus", ncol(obj))
  expected_assignment[!is.na(obj$VZ_subcluster)] <- "VZ"
  expected_assignment[!is.na(obj$RL_subcluster)] <- "RL"
  expected_assignment[!is.na(obj$residual_subcluster)] <- "Residual"
  if (!identical(assignment, expected_assignment)) {
    stop(sample_id, " regional_assignment does not match branch metadata.")
  }

  retained_layers <- intersect(c("counts", "data"), Layers(obj[["Xenium"]]))
  if (!"counts" %in% retained_layers) stop(sample_id, " Xenium assay lacks a counts layer.")
  obj <- DietSeurat(obj, assays = "Xenium", layers = retained_layers,
                    dimreducs = NULL, graphs = NULL)
  obj$sample_id <- sample_id
  obj$orig.ident <- sample_id
  obj <- RenameCells(obj, add.cell.id = sample_id)
  merge_manifest$cells[[i]] <- ncol(obj)
  merge_manifest$images[[i]] <- length(Images(obj))
  merge_manifest$PCW[[i]] <- pcw_values[[1]]
  merge_manifest$vz_cells[[i]] <- sum(assignment == "VZ")
  merge_manifest$rl_cells[[i]] <- sum(assignment == "RL")
  merge_manifest$residual_cells[[i]] <- sum(assignment == "Residual")
  merge_manifest$residual_retained_unknown_cells[[i]] <- sum(
    assignment == "Residual" & obj$residual_annotation_status == "retained_unknown",
    na.rm = TRUE
  )
  merge_manifest$unassigned_consensus_cells[[i]] <- sum(assignment == "unassigned_consensus")

  if (is.null(merged_obj)) {
    merged_obj <- obj
  } else {
    merged_obj <- merge(merged_obj, obj, project = "Xenium_VZRL_Spatial")
  }
  rm(obj)
  gc()
}

if (anyDuplicated(colnames(merged_obj))) stop("Spatial merged object contains duplicate cell names.")
if (ncol(merged_obj) != sum(merge_manifest$cells)) {
  stop("Spatial merged cell count does not equal the per-sample manifest total.")
}
assignment_totals <- with(
  merge_manifest,
  vz_cells + rl_cells + residual_cells + unassigned_consensus_cells
)
if (!identical(as.integer(assignment_totals), as.integer(merge_manifest$cells))) {
  stop("Spatial regional-assignment counts do not equal per-sample cell counts.")
}
if (!setequal(unique(as.character(merged_obj$sample_id)), sample_ids)) {
  stop("Spatial merged object does not contain exactly the configured sample IDs.")
}
saveRDS(merged_obj, merged_path, compress = FALSE)
write.csv(merge_manifest, manifest_path, row.names = FALSE)
message("Saved complete spatial merge: ", merged_path)
