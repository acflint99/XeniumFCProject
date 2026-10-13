#!/usr/bin/env Rscript

rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
args <- commandArgs(trailingOnly = TRUE)
valid <- c("--dry-run", "--overwrite")
if (any(!args %in% valid)) stop("Usage: Rscript scripts/xenium_residual_02_merge_samples.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
root <- here(config$project$outputs_dir, "xenium", "residual")
status_paths <- file.path(root, "01_subsets", "tables", paste0(sample_ids, "_ResidualSubset_Status.csv"))
output_path <- file.path(root, "02_merged", "rds", "Xenium_Residual_Merged.rds")
manifest_path <- file.path(root, "02_merged", "tables", "Xenium_Residual_Merge_Manifest.csv")
if (dry_run) {
  compact_dry_run("Residual subset merge", inputs = status_paths, outputs = c(output_path, manifest_path))
  quit(save="no", status=if(all(file.exists(status_paths))) 0L else 1L)
}
missing <- status_paths[!file.exists(status_paths)]
if (length(missing)) stop("Missing residual status file(s):\n- ", paste(missing, collapse="\n- "))
existing <- c(output_path, manifest_path)[file.exists(c(output_path, manifest_path))]
if (length(existing) && !overwrite) stop("Refusing to overwrite residual merge output(s):\n- ", paste(existing, collapse="\n- "))
status <- do.call(rbind, lapply(status_paths, function(path) read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)))
if (nrow(status) != length(sample_ids) || anyDuplicated(status$sample_id) ||
    !setequal(status$sample_id, sample_ids)) stop("Residual status files do not map one-to-one to config/samples.csv.")
status <- status[match(sample_ids, status$sample_id), , drop=FALSE]
has_cells <- as.logical(status$has_residual_cells)
if (anyNA(has_cells)) stop("Residual status has invalid has_residual_cells values.")
subset_paths <- file.path(root, "01_subsets", "rds", paste0(sample_ids, "_ResidualSubset.rds"))
if (any(!file.exists(subset_paths[has_cells]))) stop("A residual status row declares cells but its RDS is missing.")
if (!any(has_cells)) stop("No retained Other or unresolved cells are available for residual analysis.")
suppressPackageStartupMessages(library(Seurat))
merged <- NULL
manifest <- status
manifest$subset_rds <- ifelse(has_cells, subset_paths, "")
for (i in which(has_cells)) {
  sample_id <- sample_ids[[i]]
  object <- readRDS(subset_paths[[i]])
  if (ncol(object) != as.integer(status$residual_cells[[i]])) stop(sample_id, " residual cell count differs from status.")
  object$orig.ident <- sample_id
  object <- add_regional_unknown_labels(object, sample_id)
  object@images <- list()
  object <- RenameCells(object, add.cell.id = sample_id)
  if (is.null(merged)) merged <- object else merged <- merge(merged, object, project="Xenium_Residual")
  rm(object); invisible(gc(verbose=FALSE))
}
if (anyDuplicated(colnames(merged))) stop("Residual merge contains duplicate cell IDs.")
if (ncol(merged) != sum(as.integer(status$residual_cells[has_cells]))) stop("Residual merge cell count differs from status totals.")
dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(manifest_path), recursive=TRUE, showWarnings=FALSE)
saveRDS(merged, output_path, compress=FALSE)
write.csv(manifest, manifest_path, row.names=FALSE)
message("Saved residual merge with ", ncol(merged), " cells from ", sum(has_cells), " samples.")
