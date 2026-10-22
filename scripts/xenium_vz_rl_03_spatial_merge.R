#!/usr/bin/env Rscript

# Memory-conscious spatial merge of all manifest-defined corrected combined objects.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "downstream_resource_logging.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))
source(here("scripts", "R", "vz_rl_batch_merge.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
batch_option <- grep("^--batch-size=", args, value = TRUE)
if (length(batch_option) > 1L) stop("Specify --batch-size only once.")
other_args <- args[!startsWith(args, "--batch-size=")]
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- other_args[startsWith(other_args, "--") & !other_args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!other_args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_03_spatial_merge.R [--dry-run|--overwrite] [--batch-size=2..10]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
batch_size <- if (length(batch_option)) {
  suppressWarnings(as.integer(sub("^--batch-size=", "", batch_option)))
} else 6L
if (length(batch_size) != 1L || is.na(batch_size) ||
    batch_size < 2L || batch_size > 10L) {
  stop("--batch-size must be an integer from 2 through 10.")
}

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "vz_rl", "02_corrected_broad", "rds")
source_input_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
corrected_manifest_path <- file.path(
  output_root, "xenium", "vz_rl", "02_corrected_broad", "tables",
  "Xenium_VZRL_corrected_broad_manifest.csv"
)
rds_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "03_merged", "rds")
table_dir <- file.path(output_root, "xenium", "vz_rl", "spatial", "03_merged", "tables")
input_names <- paste0(sample_ids, "_VZ_RL_QC_CorrectedBroad.rds")
input_paths <- file.path(input_dir, input_names)
source_input_paths <- file.path(source_input_dir, paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds"))
merged_path <- file.path(rds_dir, "Xenium_VZRL_spatial_merged.rds")
manifest_path <- file.path(table_dir, "Xenium_VZRL_spatial_merged_manifest.csv")

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)
corrected_manifest <- if (file.exists(corrected_manifest_path)) {
  tryCatch(read.csv(corrected_manifest_path, stringsAsFactors = FALSE),
           error = function(error) NULL)
} else NULL
required_corrected_manifest <- c(
  "sample_id", "input_path", "output_path", "input_size", "input_mtime",
  "output_size", "cells", "mapping_md5"
)
corrected_manifest_valid <- !is.null(corrected_manifest) &&
  all(required_corrected_manifest %in% names(corrected_manifest)) &&
  nrow(corrected_manifest) == length(sample_ids) &&
  identical(as.character(corrected_manifest$sample_id), sample_ids) &&
  identical(as.character(corrected_manifest$input_path), source_input_paths) &&
  identical(as.character(corrected_manifest$output_path), input_paths) &&
  identical(as.character(corrected_manifest$mapping_md5), rep(
    unname(tools::md5sum(here("scripts", "R", "corrected_broad_labels.R"))),
    length(sample_ids)
  )) &&
  all(file.exists(c(source_input_paths, input_paths)))
if (corrected_manifest_valid) {
  original_info <- file.info(source_input_paths)
  corrected_info <- file.info(input_paths)
  corrected_manifest_valid <-
    isTRUE(all.equal(as.numeric(corrected_manifest$input_size), as.numeric(original_info$size))) &&
    isTRUE(all.equal(as.numeric(corrected_manifest$input_mtime), as.numeric(original_info$mtime))) &&
    isTRUE(all.equal(as.numeric(corrected_manifest$output_size), as.numeric(corrected_info$size))) &&
    all(!is.na(corrected_manifest$cells) & corrected_manifest$cells > 0L)
}

if (dry_run) {
  compact_dry_run(
    "Spatial regional merge",
    inputs = c(input_paths, corrected_manifest_path),
    outputs = c(merged_path, manifest_path),
    checks = c(no_unexpected_inputs = !length(unexpected_names),
               corrected_broad_provenance_current = corrected_manifest_valid)
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
if (!corrected_manifest_valid) {
  stop("Corrected-broad manifest is missing, incomplete, or stale. Regenerate Stage 02 before spatial merging.")
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
  corrected_broad_mapping_md5 = as.character(corrected_manifest$mapping_md5),
  merge_strategy = rep("ordered_batch_list", length(sample_ids)),
  merge_batch_size = rep(batch_size, length(sample_ids)),
  merge_batch_index = as.integer(ceiling(seq_along(sample_ids) / batch_size)),
  input_image_names = character(length(sample_ids)),
  merged_image_names = character(length(sample_ids)),
  cells = integer(length(sample_ids)), images = integer(length(sample_ids)),
  PCW = character(length(sample_ids)),
  vz_cells = integer(length(sample_ids)), rl_cells = integer(length(sample_ids)),
  residual_cells = integer(length(sample_ids)),
  residual_retained_unknown_cells = integer(length(sample_ids)),
  unassigned_consensus_cells = integer(length(sample_ids)),
  stringsAsFactors = FALSE
)
message("Merging ", length(sample_ids), " spatial samples in ",
        ceiling(length(sample_ids) / batch_size), " ordered batches of at most ",
        batch_size, ".")
batch_members <- list()
merged_batches <- list()
expected_cell_ids <- vector("list", length(sample_ids))
expected_image_names <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Reading spatial object ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])
  log_downstream_memory("OBJECT LOADED", obj)
  required_metadata <- c(
    "comb_subcluster", "consensus_label", "corrected_broad_label",
    "VZ_subcluster", "RL_subcluster",
    "VZ_subset_reason", "RL_subset_reason", "residual_subcluster",
    "residual_reference_route", "residual_subset_reason",
    "residual_annotation_confidence", "residual_annotation_status",
    "regional_assignment", "PCW"
  )
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (ncol(obj) != corrected_manifest$cells[[i]]) {
    stop(sample_id, " cell count differs from the corrected-broad manifest.")
  }
  if (!identical(as.character(obj$corrected_broad_label),
                 map_corrected_broad_labels(obj$comb_subcluster))) {
    stop(sample_id, " corrected broad labels differ from the approved fine-label mapping.")
  }
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
  merge_manifest$input_image_names[[i]] <- paste(Images(obj), collapse = "|")
  obj <- namespace_vz_rl_images(obj, sample_id)
  expected_image_names[[i]] <- Images(obj)
  merge_manifest$merged_image_names[[i]] <- paste(Images(obj), collapse = "|")
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
  expected_cell_ids[[i]] <- colnames(obj)
  batch_members[[length(batch_members) + 1L]] <- obj
  rm(obj)
  if (length(batch_members) == batch_size || i == length(sample_ids)) {
    batch_index <- length(merged_batches) + 1L
    member_indices <- which(merge_manifest$merge_batch_index == batch_index)
    batch_obj <- merge_vz_rl_objects(batch_members)
    if (ncol(batch_obj) != sum(merge_manifest$cells[member_indices]) ||
        length(Images(batch_obj)) != sum(merge_manifest$images[member_indices]) ||
        !identical(colnames(batch_obj), unlist(expected_cell_ids[member_indices], use.names = FALSE))) {
      stop("Batch ", batch_index, " lost or reordered cells or spatial images.")
    }
    merged_batches[[batch_index]] <- batch_obj
    rm(batch_obj, batch_members)
    batch_members <- list()
    gc()
    log_downstream_memory(paste("MERGED BATCH", batch_index), merged_batches[[batch_index]])
  }
}

merged_obj <- merge_vz_rl_objects(merged_batches)
rm(merged_batches, batch_members)
gc()
log_downstream_memory("MERGED ALL BATCHES", merged_obj)

if (anyDuplicated(colnames(merged_obj))) stop("Spatial merged object contains duplicate cell names.")
if (ncol(merged_obj) != sum(merge_manifest$cells)) {
  stop("Spatial merged cell count does not equal the per-sample manifest total.")
}
if (!identical(colnames(merged_obj), unlist(expected_cell_ids, use.names = FALSE))) {
  stop("Spatial merged cell IDs or order differ from the validated inputs.")
}
if (length(Images(merged_obj)) != sum(merge_manifest$images)) {
  stop("Spatial merged image/FOV count differs from the validated inputs.")
}
if (!identical(Images(merged_obj), unlist(expected_image_names, use.names = FALSE))) {
  stop("Spatial merged image/FOV names or order differ from the validated inputs.")
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
