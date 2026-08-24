#!/usr/bin/env Rscript

# Remove spatial/QC overhead and merge the complete manifest-defined set of
# combined VZ/RL whole-tissue objects.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/Xenium_CombSubclusters_Clean&Merge.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
output_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS")
input_names <- paste0(sample_ids, "_Ald_VZ_RL_QC_Subclusters.rds")
input_paths <- file.path(input_dir, input_names)
clean_names <- paste0("clean_", sample_ids, ".rds")
clean_paths <- file.path(output_dir, clean_names)
merged_path <- file.path(output_dir, "XenAld_VZRL_clean_merge.rds")
manifest_path <- file.path(output_dir, "XenAld_VZRL_clean_merge_manifest.csv")

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  cat("Expected combined inputs:", length(input_paths), "\n")
  cat("Existing combined inputs:", sum(file.exists(input_paths)), "\n")
  cat("Unexpected top-level RDS files:", length(unexpected_names), "\n")
  cat("Clean sample outputs existing:", sum(file.exists(clean_paths)), "of", length(clean_paths), "\n")
  cat("Merged output:", merged_path, "\n")
  cat("Merged output exists:", file.exists(merged_path), "\n")
  write.table(
    data.frame(sample_id = sample_ids, input_exists = file.exists(input_paths),
               clean_output_exists = file.exists(clean_paths)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  if (length(unexpected_names)) cat("Unexpected files:\n- ", paste(unexpected_names, collapse = "\n- "), "\n")
  quit(save = "no", status = 0L)
}

if (length(missing_names)) {
  stop("Cannot clean/merge combined objects; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot clean/merge combined objects; unexpected RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
expected_outputs <- c(clean_paths, merged_path, manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop("Refusing to overwrite existing combined clean/merge outputs:\n- ",
       paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
merge_manifest <- data.frame(
  sample_id = sample_ids, input_path = input_paths, clean_path = clean_paths,
  cells = integer(length(sample_ids)), PCW = character(length(sample_ids)),
  stringsAsFactors = FALSE
)
merged_obj <- NULL

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Cleaning and merging sample ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])
  required_metadata <- c("comb_subcluster", "consensus_label", "VZ_subcluster", "RL_subcluster", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (ncol(obj) == 0L) stop(sample_id, " contains zero cells.")
  pcw_values <- unique(as.character(obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

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
  obj <- RenameCells(obj, add.cell.id = sample_id)
  if (is.null(merged_obj)) {
    merged_obj <- obj
  } else {
    merged_obj <- merge(merged_obj, obj, project = "Xenium_VZRL_Combined")
  }
  rm(obj)
  gc()
}

if (anyDuplicated(colnames(merged_obj))) stop("Combined merged object contains duplicate cell names.")
if (ncol(merged_obj) != sum(merge_manifest$cells)) {
  stop("Combined merged cell count does not equal the per-sample manifest total.")
}
if (!setequal(unique(as.character(merged_obj$orig.ident)), sample_ids)) {
  stop("Combined merged object does not contain exactly the configured sample IDs.")
}
Idents(merged_obj) <- "comb_subcluster"
saveRDS(merged_obj, merged_path, compress = FALSE)
write.csv(merge_manifest, manifest_path, row.names = FALSE)
message("Saved complete ", length(sample_ids), "-sample combined merge: ", merged_path)
