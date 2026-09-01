# Preprocess one biological sample from a manually separated multi-sample slide.
# A Slurm array task ID selects one of the manifest rows. Use --list to print
# the complete task-to-sample mapping without loading Xenium data.

rm(list = ls())

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
split_samples <- samples[
  samples$input_layout == "multi_sample_slide",
  ,
  drop = FALSE
]

if (nrow(split_samples) == 0L) {
  stop("No multi-sample-slide rows were found in config/samples.csv.")
}

task_map <- data.frame(
  task_id = seq_len(nrow(split_samples)),
  sample_id = split_samples$sample_id,
  input_directory = split_samples$input_directory,
  cell_stats_file = split_samples$cell_stats_file,
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)

if (identical(args, "--list")) {
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) {
  stop("--dry-run and --overwrite cannot be used together.")
}

task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_preprocess_split_slides.R ",
    "[--dry-run|--overwrite] [TASK_ID]"
  )
}

task_value <- if (length(task_args) == 1L) {
  task_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}

task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(split_samples)) {
  stop(
    "TASK_ID must be between 1 and ", nrow(split_samples),
    ". Run with --list to see the mapping."
  )
}

sample_record <- split_samples[task_id, , drop = FALSE]
sample_paths <- resolve_sample_paths(sample_record, config)
sample_name <- sample_paths$sample_id

if (!dir.exists(sample_paths$input_dir)) {
  stop("Input slide directory does not exist: ", sample_paths$input_dir)
}
if (!file.exists(sample_paths$cell_stats_path)) {
  stop("Cell-stat CSV does not exist: ", sample_paths$cell_stats_path)
}

output_root <- here(config$project$outputs_dir)
expected_outputs <- c(
  file.path(output_root, "xenium", "preprocess", "01_cropped", "plots", paste0(sample_name, "_CB_nCount_FeatPlot.tif")),
  file.path(output_root, "xenium", "preprocess", "01_cropped", "rds", paste0(sample_name, "_CB.rds")),
  file.path(output_root, "xenium", "preprocess", "02_qc", "reports", paste0(sample_name, "_QCplots.pdf")),
  file.path(output_root, "xenium", "preprocess", "02_qc", "reports", paste0(sample_name, "_QC_thresholds.txt")),
  file.path(output_root, "xenium", "preprocess", "02_qc", "rds", paste0(sample_name, "_CB_QC.rds")),
  file.path(output_root, "xenium", "preprocess", "03_clustered", "plots", paste0(sample_name, "_UMAP.tif")),
  file.path(output_root, "xenium", "preprocess", "03_clustered", "plots", paste0(sample_name, "_RawCluster_UMAP.tif")),
  file.path(output_root, "xenium", "preprocess", "03_clustered", "plots", paste0(sample_name, "_GlobalRawClustersSpatialPlot.tif")),
  file.path(output_root, "xenium", "preprocess", "03_clustered", "plots", paste0(sample_name, "_FacetRawClustersSpatialPlot.tif")),
  file.path(output_root, "xenium", "preprocess", "03_clustered", "rds", paste0(sample_name, "_CB_QC_cluster.rds"))
)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]

if (dry_run) {
  compact_dry_run(
    paste0(
      "Split-slide preprocessing task ", task_id, "/", nrow(split_samples),
      " [", sample_name, "]"
    ),
    inputs = c(sample_paths$input_dir, sample_paths$cell_stats_path),
    outputs = expected_outputs
  )
  quit(save = "no", status = 0L)
}

if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Refusing to overwrite existing outputs for ", sample_name, ":\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}
if (length(existing_outputs) > 0L) {
  warning("Overwriting ", length(existing_outputs), " existing outputs for ", sample_name, ".")
}

worker_cores <- suppressWarnings(as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")
))
if (is.na(worker_cores) || worker_cores < 1L) worker_cores <- 8L

Sys.setenv(
  OMP_NUM_THREADS = worker_cores,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1
)

library(Seurat)
library(future)

source(here("scripts", "xenium_preprocess_01_crop_cerebellum.R"))
source(here("scripts", "xenium_preprocess_02_qc_cells.R"))
source(here("scripts", "xenium_preprocess_03_normalize_cluster.R"))

future::plan(future::sequential)
options(
  future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3,
  future.fork.enable = FALSE
)

message("Task ", task_id, " of ", nrow(split_samples))
message("Biological sample: ", sample_name)
message("Input slide: ", sample_record$input_directory[[1]])
message("Cell-stat CSV: ", sample_paths$cell_stats_file)
message("Started: ", Sys.time())

xenium_cereb <- XeniumCropCerebellum(
  sample_name = sample_name,
  cell_stat_file = sample_paths$cell_stats_file,
  input_directory = sample_record$input_directory[[1]],
  sample_path = sample_paths$input_dir
)
xenium_cereb_qc <- qc_xenium(xenium_cereb, sample_name = sample_name)
process_xenium_clusters(xenium_cereb_qc, sample_name = sample_name)

rm(xenium_cereb, xenium_cereb_qc)
gc()

message("Finished: ", sample_name, " at ", Sys.time())
