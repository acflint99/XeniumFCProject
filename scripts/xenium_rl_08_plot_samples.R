#!/usr/bin/env Rscript

# Slurm-array driver for one manifest-defined RL spatial report per task.
rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
task_map <- data.frame(task_id = seq_along(sample_ids), sample_id = sample_ids)
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop("Usage: Rscript scripts/xenium_rl_08_plot_samples.R [--list|--dry-run|--overwrite] [TASK_ID]")
}
task_value <- if (length(task_args)) task_args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > length(sample_ids)) {
  stop("TASK_ID must be between 1 and ", length(sample_ids), ". Use --list to inspect the mapping.")
}

sample_id <- sample_ids[[task_id]]
output_root <- here(config$project$outputs_dir)
input_path <- file.path(
  output_root, "xenium", "rl", "07_mapped", "rds",
  paste0(sample_id, "_Ald_VZ_RL_QC_Subclusters.rds")
)
output_base <- file.path(output_root, "xenium", "rl", "08_sample_reports")
sample_dir <- file.path(output_base, sample_id)
output_suffixes <- c(
  "_cluster_counts.csv", "_Global_Spatial2.tif", "_RL_Spatial.tif",
  "_Granule_Spatial2.tif", "_Granule&RL_Spatial2.tif", "_UBC_Spatial2.tif",
  "_Faceted_Clusters0.3.tif", "_Markers_DotPlot.tif", "_Markers_DotPlot.pdf"
)
output_paths <- file.path(sample_dir, paste0(sample_id, output_suffixes))

if (dry_run) {
  compact_dry_run(
    paste0("RL report task ", task_id, "/", length(sample_ids), " [", sample_id, "]"),
    inputs = input_path,
    outputs = output_paths
  )
  quit(save = "no", status = 0L)
}
if (!file.exists(input_path)) stop("RL mapped input not found: ", input_path)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop("Refusing to overwrite existing RL report outputs:\n- ",
       paste(existing_outputs, collapse = "\n- "), "\nUse --overwrite only after review.")
}

source(here("scripts", "xenium_rl_plot_sample.R"))
generate_spatial_reports(sample_id, input_path, output_base)
missing_outputs <- output_paths[!file.exists(output_paths)]
if (length(missing_outputs)) stop("RL report did not create:\n- ", paste(missing_outputs, collapse = "\n- "))
message("Completed RL spatial report for ", sample_id, ".")
