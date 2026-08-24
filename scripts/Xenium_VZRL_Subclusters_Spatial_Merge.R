#!/usr/bin/env Rscript

# Memory-conscious spatial merge of all 34 manifest-defined combined objects.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
})
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/Xenium_VZRL_Subclusters_Spatial_Merge.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
output_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Merged_RDS")
input_names <- paste0(sample_ids, "_Ald_VZ_RL_QC_Subclusters.rds")
input_paths <- file.path(input_dir, input_names)
merged_path <- file.path(output_dir, "XenAld_VZRL_spatial_merged.rds")
manifest_path <- file.path(output_dir, "XenAld_VZRL_spatial_merged_manifest.csv")

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  cat("Expected spatial-merge inputs:", length(input_paths), "\n")
  cat("Existing spatial-merge inputs:", sum(file.exists(input_paths)), "\n")
  cat("Unexpected top-level RDS files:", length(unexpected_names), "\n")
  cat("Spatial merged output:", merged_path, "\n")
  cat("Spatial merged output exists:", file.exists(merged_path), "\n")
  cat("Spatial merge manifest exists:", file.exists(manifest_path), "\n")
  write.table(data.frame(sample_id = sample_ids, input_exists = file.exists(input_paths)),
              row.names = FALSE, quote = FALSE, sep = "\t")
  if (length(unexpected_names)) cat("Unexpected files:\n- ", paste(unexpected_names, collapse = "\n- "), "\n")
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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
merge_manifest <- data.frame(
  sample_id = sample_ids, input_path = input_paths,
  cells = integer(length(sample_ids)), images = integer(length(sample_ids)),
  PCW = character(length(sample_ids)), stringsAsFactors = FALSE
)
merged_obj <- NULL

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Reading spatial object ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])
  required_metadata <- c("comb_subcluster", "consensus_label", "VZ_subcluster", "RL_subcluster", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (ncol(obj) == 0L) stop(sample_id, " contains zero cells.")
  if (!length(Images(obj))) stop(sample_id, " has no spatial image/FOV data.")
  pcw_values <- unique(as.character(obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

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
if (!setequal(unique(as.character(merged_obj$sample_id)), sample_ids)) {
  stop("Spatial merged object does not contain exactly the configured sample IDs.")
}
saveRDS(merged_obj, merged_path, compress = FALSE)
write.csv(merge_manifest, manifest_path, row.names = FALSE)
message("Saved complete spatial merge: ", merged_path)
