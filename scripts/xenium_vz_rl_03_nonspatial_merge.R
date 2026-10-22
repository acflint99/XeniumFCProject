#!/usr/bin/env Rscript

# Remove spatial/QC overhead and merge the complete manifest-defined set of
# combined VZ/RL/residual whole-tissue objects.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "downstream_resource_logging.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_03_nonspatial_merge.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "vz_rl", "02_corrected_broad", "rds")
source_input_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
corrected_manifest_path <- file.path(
  output_root, "xenium", "vz_rl", "02_corrected_broad", "tables",
  "Xenium_VZRL_corrected_broad_manifest.csv"
)
rds_dir <- file.path(output_root, "xenium", "vz_rl", "nonspatial", "03_merged", "rds")
table_dir <- file.path(output_root, "xenium", "vz_rl", "nonspatial", "03_merged", "tables")
input_names <- paste0(sample_ids, "_VZ_RL_QC_CorrectedBroad.rds")
input_paths <- file.path(input_dir, input_names)
source_input_paths <- file.path(source_input_dir, paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds"))
clean_names <- paste0("clean_", sample_ids, ".rds")
clean_paths <- file.path(rds_dir, clean_names)
merged_path <- file.path(rds_dir, "Xenium_VZRL_merged.rds")
manifest_path <- file.path(table_dir, "Xenium_VZRL_merge_manifest.csv")

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
    "Combined regional merge",
    inputs = c(input_paths, corrected_manifest_path),
    outputs = c(clean_paths, merged_path, manifest_path),
    checks = c(no_unexpected_inputs = !length(unexpected_names),
               corrected_broad_provenance_current = corrected_manifest_valid)
  )
  if (length(unexpected_names)) cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  quit(save = "no", status = 0L)
}

if (length(missing_names)) {
  stop("Cannot clean/merge combined objects; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot clean/merge combined objects; unexpected RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
if (!corrected_manifest_valid) {
  stop("Corrected-broad manifest is missing, incomplete, or stale. Regenerate Stage 02 before merging.")
}
expected_outputs <- c(clean_paths, merged_path, manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop("Refusing to overwrite existing combined clean/merge outputs:\n- ",
       paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
}

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
merge_manifest <- data.frame(
  sample_id = sample_ids, input_path = input_paths, clean_path = clean_paths,
  corrected_broad_mapping_md5 = as.character(corrected_manifest$mapping_md5),
  cells = integer(length(sample_ids)), PCW = character(length(sample_ids)),
  vz_cells = integer(length(sample_ids)), rl_cells = integer(length(sample_ids)),
  residual_cells = integer(length(sample_ids)),
  residual_retained_unknown_cells = integer(length(sample_ids)),
  unassigned_consensus_cells = integer(length(sample_ids)),
  stringsAsFactors = FALSE
)
merged_obj <- NULL

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Cleaning and merging sample ", i, " of ", length(sample_ids), ": ", sample_id)
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
  expected_broad <- map_corrected_broad_labels(obj$comb_subcluster)
  if (!identical(as.character(obj$corrected_broad_label), expected_broad) ||
      ncol(obj) != corrected_manifest$cells[[i]]) {
    stop(sample_id, " corrected broad labels or cell count differ from the approved mapping manifest.")
  }
  if (ncol(obj) == 0L) stop(sample_id, " contains zero cells.")
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

  if ("scale.data" %in% Layers(obj[["Xenium"]])) obj[["Xenium"]]@layers$scale.data <- NULL
  for (assay_name in c("BlankCodeword", "ControlCodeword", "ControlProbe", "GenomicControl")) {
    if (assay_name %in% Assays(obj)) obj[[assay_name]] <- NULL
  }
  obj@images <- list()
  obj@project.name <- sample_id
  obj$orig.ident <- sample_id
  Idents(obj) <- "comb_subcluster"
  saveRDS(obj, clean_paths[[i]], compress = FALSE)

  merge_manifest$cells[[i]] <- ncol(obj)
  merge_manifest$PCW[[i]] <- pcw_values[[1]]
  merge_manifest$vz_cells[[i]] <- sum(assignment == "VZ")
  merge_manifest$rl_cells[[i]] <- sum(assignment == "RL")
  merge_manifest$residual_cells[[i]] <- sum(assignment == "Residual")
  merge_manifest$residual_retained_unknown_cells[[i]] <- sum(
    assignment == "Residual" & obj$residual_annotation_status == "retained_unknown",
    na.rm = TRUE
  )
  merge_manifest$unassigned_consensus_cells[[i]] <- sum(assignment == "unassigned_consensus")
  obj <- RenameCells(obj, add.cell.id = sample_id)
  if (is.null(merged_obj)) {
    merged_obj <- obj
  } else {
    merged_obj <- merge(merged_obj, obj, project = "Xenium_VZRL_Combined")
  }
  rm(obj)
  gc()
  log_downstream_memory(paste("MERGED SAMPLE", sample_id), merged_obj)
}

if (anyDuplicated(colnames(merged_obj))) stop("Combined merged object contains duplicate cell names.")
if (ncol(merged_obj) != sum(merge_manifest$cells)) {
  stop("Combined merged cell count does not equal the per-sample manifest total.")
}
assignment_totals <- with(
  merge_manifest,
  vz_cells + rl_cells + residual_cells + unassigned_consensus_cells
)
if (!identical(as.integer(assignment_totals), as.integer(merge_manifest$cells))) {
  stop("Combined regional-assignment counts do not equal per-sample cell counts.")
}
if (!setequal(unique(as.character(merged_obj$orig.ident)), sample_ids)) {
  stop("Combined merged object does not contain exactly the configured sample IDs.")
}
Idents(merged_obj) <- "comb_subcluster"
saveRDS(merged_obj, merged_path, compress = FALSE)
write.csv(merge_manifest, manifest_path, row.names = FALSE)
message("Saved complete ", length(sample_ids), "-sample combined merge: ", merged_path)
