#!/usr/bin/env Rscript

# Merge the complete manifest-defined set of VZ consensus subsets.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- samples$sample_id

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_02_merge_samples.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
subset_dir <- file.path(output_root, "xenium", "vz", "01_subsets", "rds")
merged_rds_dir <- file.path(output_root, "xenium", "vz", "02_merged", "rds")
merged_table_dir <- file.path(output_root, "xenium", "vz", "02_merged", "tables")
subset_names <- paste0(sample_ids, "_VZsubset.rds")
subset_paths <- file.path(subset_dir, subset_names)
output_path <- file.path(merged_rds_dir, "Xenium_Merged_VZSubsets.rds")
manifest_path <- file.path(merged_table_dir, "Xenium_Merged_VZSubsets_manifest.csv")

observed_names <- if (dir.exists(subset_dir)) {
  list.files(subset_dir, pattern = "\\.rds$", full.names = FALSE)
} else {
  character()
}
missing_names <- subset_names[!file.exists(subset_paths)]
unexpected_names <- setdiff(observed_names, subset_names)

if (dry_run) {
  compact_dry_run(
    "VZ subset merge",
    inputs = subset_paths,
    outputs = c(output_path, manifest_path),
    checks = c(no_unexpected_inputs = !length(unexpected_names))
  )
  if (length(unexpected_names)) cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  quit(save = "no", status = 0L)
}

if (length(missing_names)) {
  stop("Cannot merge VZ subsets; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot merge VZ subsets; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}

existing_outputs <- c(output_path, manifest_path)[file.exists(c(output_path, manifest_path))]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ merge outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages(library(Seurat))

merge_manifest <- data.frame(
  sample_id = sample_ids,
  input_path = subset_paths,
  cells = integer(length(sample_ids)),
  PCW = character(length(sample_ids)),
  stringsAsFactors = FALSE
)
merged_obj <- NULL

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Reading VZ subset ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(subset_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  }
  pcw_values <- unique(as.character(obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")
  if (ncol(obj) == 0L) stop(sample_id, " VZ subset contains zero cells.")

  obj$orig.ident <- sample_id
  Idents(obj) <- "consensus_label"
  obj@images <- list()
  obj <- RenameCells(obj, add.cell.id = sample_id)
  merge_manifest$cells[[i]] <- ncol(obj)
  merge_manifest$PCW[[i]] <- pcw_values[[1]]

  if (is.null(merged_obj)) {
    merged_obj <- obj
  } else {
    merged_obj <- merge(
      x = merged_obj,
      y = obj,
      project = "Xenium_VZ_Refinement"
    )
  }
  rm(obj)
  gc()
}

if (anyDuplicated(colnames(merged_obj))) stop("Merged VZ object contains duplicate cell names.")
if (ncol(merged_obj) != sum(merge_manifest$cells)) {
  stop("Merged VZ cell count does not equal the sum of per-sample cell counts.")
}
Idents(merged_obj) <- "consensus_label"

dir.create(merged_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_table_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(merged_obj, output_path, compress = FALSE)
write.csv(merge_manifest, manifest_path, row.names = FALSE)
message("Saved complete 34-sample VZ merge: ", output_path)
message("Saved VZ merge manifest: ", manifest_path)
