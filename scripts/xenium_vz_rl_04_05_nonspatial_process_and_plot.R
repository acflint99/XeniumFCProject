#!/usr/bin/env Rscript

# Process the complete combined VZ/RL/residual merge, or create plots from the saved
# processed object without repeating integration.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(future)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "downstream_resource_logging.R"))
source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))
source(here("scripts", "R", "nonspatial_merge_plots.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")

# Keep future evaluation sequential even for UMAP. Allocated CPUs may still
# serve compiled routines, but do not create full-object R worker sessions.
run_umap_with_workers <- function(object, ...) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  Seurat::RunUMAP(object, ...)
}

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--process-only", "--plots-only", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) stop("All arguments must be named options.")
process_only <- "--process-only" %in% args
plots_only <- "--plots-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (process_only && plots_only) stop("Choose --process-only or --plots-only, not both.")
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
nonspatial_root <- file.path(output_root, "xenium", "vz_rl", "nonspatial")
# Project-relative overrides support isolated reruns from an archived merge.
# The backup remains read-only; both jobs must use the same overrides.
rerun_directory <- function(variable, default) {
  value <- Sys.getenv(variable, unset = "")
  if (!nzchar(value)) return(default)
  if (!nzchar(trimws(value))) stop(variable, " must not contain only whitespace.")
  normalizePath(here(value), mustWork = FALSE)
}
merge_dir <- rerun_directory("COMBINED_MERGE_INPUT_DIR", file.path(nonspatial_root, "03_merged"))
nonspatial_root <- rerun_directory("COMBINED_NONSPATIAL_OUTPUT_DIR", nonspatial_root)
input_path <- file.path(merge_dir, "rds", "Xenium_VZRL_merged.rds")
input_manifest_path <- file.path(merge_dir, "tables", "Xenium_VZRL_merge_manifest.csv")
rds_dir <- file.path(nonspatial_root, "04_processed", "rds")
plot_dir <- rerun_directory("COMBINED_PLOT_OUTPUT_DIR", file.path(nonspatial_root, "05_plots"))
table_dir <- file.path(nonspatial_root, "04_processed", "tables")
processed_path <- file.path(rds_dir, "Xenium_VZRL_processed.rds")
processed_manifest_path <- file.path(table_dir, "Xenium_VZRL_processed_manifest.csv")
message("Combined merge input: ", input_path)
message("Combined processed object: ", processed_path)
message("Combined plot directory: ", plot_dir)

plot_names <- c(nonspatial_plot_spec()$filename, "Xenium_Merged_plot_manifest.csv",
                "Xenium_Merged_Heatmap_marker_summary.csv")
plot_paths <- file.path(plot_dir, plot_names)

if (dry_run) {
  compact_dry_run(
    paste0(
      "Combined regional analysis [",
      if (process_only) "processing" else if (plots_only) "plotting" else "inspection",
      "]"
    ),
    inputs = if (plots_only) c(processed_path, processed_manifest_path, input_path) else c(input_path, input_manifest_path),
    outputs = c(processed_path, processed_manifest_path, plot_paths)
  )
  quit(save = "no", status = 0L)
}
if (!process_only && !plots_only) {
  stop("Choose --process-only for integration or --plots-only after processing.")
}

if (process_only) {
  if (!file.exists(input_path)) stop("Combined merge not found: ", input_path)
  if (!file.exists(input_manifest_path)) stop("Combined merge manifest not found: ", input_manifest_path)
  existing_outputs <- c(processed_path, processed_manifest_path)[
    file.exists(c(processed_path, processed_manifest_path))
  ]
  if (length(existing_outputs) && !overwrite) {
    stop("Refusing to overwrite combined processed outputs:\n- ",
         paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
  }

  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  input_manifest <- read.csv(input_manifest_path, stringsAsFactors = FALSE)
  if (nrow(input_manifest) != length(sample_ids) ||
      !setequal(input_manifest$sample_id, sample_ids)) {
    stop("Combined input manifest must contain exactly the samples in config/samples.csv.")
  }
  merged_obj <- readRDS(input_path)
  log_downstream_memory("OBJECT LOADED", merged_obj)
  required_metadata <- c(
    "orig.ident", "comb_subcluster", "consensus_label", "corrected_broad_label",
    "VZ_subcluster",
    "RL_subcluster", "VZ_subset_reason", "RL_subset_reason",
    "residual_subcluster", "residual_reference_route",
    "residual_subset_reason", "residual_annotation_confidence",
    "residual_annotation_status",
    "regional_assignment", "PCW"
  )
  missing_metadata <- setdiff(required_metadata, colnames(merged_obj[[]]))
  if (length(missing_metadata)) stop("Combined merge lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if (!identical(as.character(merged_obj$corrected_broad_label),
                 map_corrected_broad_labels(merged_obj$comb_subcluster))) {
    stop("Combined merge corrected broad labels differ from the approved fine-label mapping.")
  }
  if (ncol(merged_obj) != sum(input_manifest$cells)) {
    stop("Combined merge cell count does not match its manifest.")
  }
  if (!setequal(unique(as.character(merged_obj$orig.ident)), sample_ids)) {
    stop("Combined merge does not contain exactly the configured sample IDs.")
  }

  set.seed(config$runtime$random_seed)
  plan(sequential)
  DefaultAssay(merged_obj) <- "Xenium"
  merged_obj <- NormalizeData(merged_obj)
  log_downstream_memory("NormalizeData COMPLETE", merged_obj)
  merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000)
  merged_obj <- JoinLayers(merged_obj)
  log_downstream_memory("JoinLayers COMPLETE", merged_obj)
  merged_obj <- ScaleData(merged_obj)
  log_downstream_memory("ScaleData COMPLETE", merged_obj)
  merged_obj <- RunPCA(merged_obj, npcs = 30, seed.use = config$runtime$random_seed)
  log_downstream_memory("RunPCA COMPLETE", merged_obj)
  merged_obj <- run_umap_with_workers(
    merged_obj, reduction = "pca", dims = 1:30,
    reduction.name = "umap.unintegrated", seed.use = config$runtime$random_seed
  )
  merged_obj <- RunHarmony(
    merged_obj, group.by.vars = "orig.ident", reduction.use = "pca",
    reduction.save = "harmony"
  )
  merged_obj <- run_umap_with_workers(
    merged_obj, reduction = "harmony", dims = 1:30,
    reduction.name = "umap.harmony", seed.use = config$runtime$random_seed
  )
  saveRDS(merged_obj, processed_path, compress = FALSE)
  input_info <- file.info(input_path)
  processed_manifest <- data.frame(
    input_path = input_path, input_size = as.numeric(input_info$size),
    input_mtime = as.numeric(input_info$mtime), cells = ncol(merged_obj),
    samples = length(unique(merged_obj$orig.ident)),
    random_seed = config$runtime$random_seed, stringsAsFactors = FALSE
  )
  write.csv(processed_manifest, processed_manifest_path, row.names = FALSE)
  message("Saved processed combined object: ", processed_path)
  quit(save = "no", status = 0L)
}

if (!file.exists(processed_path)) stop("Processed combined object not found: ", processed_path)
if (!file.exists(processed_manifest_path)) stop("Processed manifest not found: ", processed_manifest_path)
existing_plots <- plot_paths[file.exists(plot_paths)]
if (length(existing_plots) && !overwrite) {
  stop("Refusing to overwrite combined plot outputs:\n- ",
       paste(existing_plots, collapse = "\n- "), "\nUse --overwrite only after review.")
}
processed_manifest <- read.csv(processed_manifest_path, stringsAsFactors = FALSE)
input_info <- file.info(input_path)
input_is_current <- file.exists(input_path) && nrow(processed_manifest) == 1L &&
  isTRUE(all.equal(as.numeric(processed_manifest$input_size), as.numeric(input_info$size))) &&
  isTRUE(all.equal(as.numeric(processed_manifest$input_mtime), as.numeric(input_info$mtime)))
if (!input_is_current) stop("Combined merge changed after processing. Rerun --process-only.")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
merged_obj <- readRDS(processed_path)
log_downstream_memory("OBJECT LOADED", merged_obj)
if (ncol(merged_obj) != processed_manifest$cells[[1]]) stop("Processed object cell count does not match its manifest.")
if (!"corrected_broad_label" %in% colnames(merged_obj[[]]) ||
    !identical(as.character(merged_obj$corrected_broad_label),
               map_corrected_broad_labels(merged_obj$comb_subcluster))) {
  stop("Processed object lacks current corrected broad labels.")
}
required_reductions <- c("umap.unintegrated", "umap.harmony")
if (!all(required_reductions %in% Reductions(merged_obj))) {
  stop("Processed object lacks reduction(s): ", paste(setdiff(required_reductions, Reductions(merged_obj)), collapse = ", "))
}
set.seed(config$runtime$random_seed)

render_nonspatial_merge_plots(
  merged_obj, plot_dir = plot_dir, sample_ids = sample_ids,
  seed = config$runtime$random_seed, processed_path = processed_path
)
message("Saved combined plots without rerunning integration.")
