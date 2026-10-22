#!/usr/bin/env Rscript

# Lightweight post-job audit of the 41 corrected-broad sample outputs.
# Reads CSVs and file metadata only; no RDS object or TIFF pixels are loaded.

suppressPackageStartupMessages(library(here))
here::i_am("scripts/validate_xenium_vz_rl_02_correct_broad.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))
config <- load_pipeline_config()
sample_ids <- as.character(load_sample_manifest(config)$sample_id)
stage_root <- here(config$project$outputs_dir, "xenium", "vz_rl", "02_corrected_broad")
source_root <- here(config$project$outputs_dir, "xenium", "vz_rl", "01_combined_labels", "rds")
manifest_path <- file.path(stage_root, "tables", "Xenium_VZRL_corrected_broad_manifest.csv")
if (!file.exists(manifest_path)) stop("Corrected-broad manifest is missing: ", manifest_path)
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required <- c(
  "sample_id", "input_path", "output_path", "input_size", "input_mtime",
  "output_size", "cells", "mapping_md5", "corrected_unknown_cells",
  "changed_broad_cells"
)
if (!all(required %in% names(manifest)) ||
    !identical(as.character(manifest$sample_id), sample_ids) ||
    anyDuplicated(manifest$sample_id)) {
  stop("Corrected-broad manifest columns or sample IDs/order are invalid.")
}
expected_input <- file.path(source_root, paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds"))
expected_output <- file.path(stage_root, "rds", paste0(sample_ids, "_VZ_RL_QC_CorrectedBroad.rds"))
if (!identical(as.character(manifest$input_path), expected_input) ||
    !identical(as.character(manifest$output_path), expected_output) ||
    !identical(as.character(manifest$mapping_md5), rep(
      unname(tools::md5sum(here("scripts", "R", "corrected_broad_labels.R"))),
      length(sample_ids)
    ))) {
  stop("Corrected-broad paths or approved mapping signature are stale.")
}
all_rds <- c(expected_input, expected_output)
if (any(!file.exists(all_rds))) stop("Combined or corrected-broad RDS files are missing.")
source_info <- file.info(expected_input)
output_info <- file.info(expected_output)
if (!isTRUE(all.equal(as.numeric(manifest$input_size), as.numeric(source_info$size))) ||
    !isTRUE(all.equal(as.numeric(manifest$input_mtime), as.numeric(source_info$mtime))) ||
    !isTRUE(all.equal(as.numeric(manifest$output_size), as.numeric(output_info$size))) ||
    anyNA(output_info$size) || any(output_info$size <= 0) ||
    any(output_info$mtime < source_info$mtime)) {
  stop("Corrected-broad RDS size or timestamp differs from the provenance manifest.")
}
if (anyNA(manifest$cells) || any(manifest$cells <= 0L) ||
    anyNA(manifest$corrected_unknown_cells) || any(manifest$corrected_unknown_cells < 0L) ||
    anyNA(manifest$changed_broad_cells) || any(manifest$changed_broad_cells < 0L)) {
  stop("Corrected-broad manifest has invalid cell counts.")
}

plot_suffixes <- c("_CorrectedBroad_GlobalSpatial_pt0.6.tif",
                   "_CorrectedBroad_GlobalSpatial_pt0.85.tif")
plot_paths <- file.path(
  stage_root, "plots",
  paste0(rep(sample_ids, each = length(plot_suffixes)),
         rep(plot_suffixes, times = length(sample_ids)))
)
count_paths <- file.path(stage_root, "tables", paste0(sample_ids, "_counts_corrected_broad_label.csv"))
transition_paths <- file.path(stage_root, "tables", paste0(sample_ids, "_consensus_to_corrected_broad.csv"))
all_small_outputs <- c(plot_paths, count_paths, transition_paths)
if (any(!file.exists(all_small_outputs)) ||
    anyNA(file.info(all_small_outputs)$size) ||
    any(file.info(all_small_outputs)$size <= 0)) {
  stop("A corrected-broad plot or audit table is missing or empty.")
}
for (i in seq_along(sample_ids)) {
  counts <- read.csv(count_paths[[i]], stringsAsFactors = FALSE)
  transitions <- read.csv(transition_paths[[i]], stringsAsFactors = FALSE)
  transition_by_label <- if (all(c("Corrected_Broad_Label", "Cell_Count") %in% names(transitions))) {
    tapply(transitions$Cell_Count, transitions$Corrected_Broad_Label, sum)
  } else numeric()
  if (!all(c("Corrected_Broad_Label", "Cell_Count") %in% names(counts)) ||
      !all(c("Original_Consensus_Label", "Combined_Subcluster",
             "Corrected_Broad_Label", "Cell_Count") %in% names(transitions)) ||
      anyDuplicated(counts$Corrected_Broad_Label) ||
      anyNA(counts$Cell_Count) || any(counts$Cell_Count <= 0L) ||
      anyNA(transitions$Cell_Count) || any(transitions$Cell_Count <= 0L) ||
      sum(counts$Cell_Count) != manifest$cells[[i]] ||
      sum(transitions$Cell_Count) != manifest$cells[[i]] ||
      !setequal(names(transition_by_label), as.character(counts$Corrected_Broad_Label)) ||
      !identical(as.numeric(transition_by_label[as.character(counts$Corrected_Broad_Label)]),
                 as.numeric(counts$Cell_Count)) ||
      !identical(map_corrected_broad_labels(transitions$Combined_Subcluster),
                 as.character(transitions$Corrected_Broad_Label)) ||
      sum(counts$Cell_Count[counts$Corrected_Broad_Label == "Unknown"]) !=
        manifest$corrected_unknown_cells[[i]] ||
      sum(transitions$Cell_Count[
        transitions$Original_Consensus_Label != transitions$Corrected_Broad_Label
      ]) != manifest$changed_broad_cells[[i]]) {
    stop("Corrected-broad tables do not reconcile for sample ", sample_ids[[i]], ".")
  }
}
cat("PASS | corrected broad | samples=", length(sample_ids), "/", length(sample_ids),
    " | cells=", format(sum(manifest$cells), scientific = FALSE),
    " | plots=", length(plot_paths),
    " | original combined files retained; RDS contents not loaded\n", sep = "")
