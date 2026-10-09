#!/usr/bin/env Rscript

# Slurm-array driver for one manifest-defined RL spatial report per task.
rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- samples$sample_id
task_map <- data.frame(task_id = seq_along(sample_ids), sample_id = sample_ids)
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
valid_options <- c("--dry-run", "--overwrite", "--combine")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
combine_reports <- "--combine" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L || (combine_reports && length(task_args))) {
  stop(
    "Usage: Rscript scripts/xenium_rl_08_plot_samples.R ",
    "[--list | --combine [--dry-run|--overwrite] | [--dry-run|--overwrite] TASK_ID]"
  )
}
output_root <- here(config$project$outputs_dir)
output_base <- file.path(output_root, "xenium", "rl", "08_sample_reports")
if (combine_reports) {
  source(here("scripts", "R", "spatial_report_helpers.R"))
  page_map <- build_stage8_pcw_page_map(config, samples)
  combine_stage8_sample_reports(
    page_map = page_map,
    report_root = output_base,
    branch = "RL",
    dry_run = dry_run,
    overwrite = overwrite
  )
  quit(save = "no", status = 0L)
}
task_value <- if (length(task_args)) task_args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > length(sample_ids)) {
  stop("TASK_ID must be between 1 and ", length(sample_ids), ". Use --list to inspect the mapping.")
}

sample_id <- sample_ids[[task_id]]
input_path <- file.path(
  output_root, "xenium", "rl", "07_mapped", "rds",
  paste0(sample_id, "_RL_QC_Subclusters.rds")
)
sample_dir <- file.path(output_base, sample_id)
output_stem <- paste0(sample_id, "_RL")
output_suffixes <- c(
  "_cluster_counts.csv", "_Global_Spatial_pt0.75.tif", "_Global_Spatial_NoOther_pt0.75.tif",
  "_Lineage_Spatial_pt0.75.tif",
  "_Granule_Spatial_pt0.75.tif", "_Granule&RL_Spatial_pt0.75.tif", "_UBC_Spatial_pt0.75.tif",
  "_Faceted_Clusters_pt0.6.tif", "_Faceted_Other_pt0.6.tif",
  "_Markers_DotPlot.tif", "_Markers_DotPlot.pdf", "_Spatial_Report.pdf"
)
output_paths <- file.path(sample_dir, paste0(output_stem, output_suffixes))

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
tiff_outputs <- output_paths[grepl("\\.tif$", output_paths)]
tiff_sizes <- file.info(tiff_outputs)$size
undersized_tiffs <- tiff_outputs[is.na(tiff_sizes) | tiff_sizes < 10000L]
if (length(undersized_tiffs)) {
  stop(
    "RL report created implausibly small TIFF output(s); check device dimensions:\n- ",
    paste(undersized_tiffs, collapse = "\n- ")
  )
}
message("Completed RL spatial report for ", sample_id, ".")
