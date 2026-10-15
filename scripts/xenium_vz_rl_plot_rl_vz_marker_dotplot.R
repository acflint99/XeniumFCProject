#!/usr/bin/env Rscript

# Plot the same RL-VZ-focused marker panel across either the RL or VZ
# resolution-0.8 post-QC clusters. This is a read-only diagnostic and does not
# recalculate clusters or write a Seurat object.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "subcluster_review_markers.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)

branch_map <- data.frame(
  task_id = 1:2,
  branch = c("RL", "VZ"),
  stringsAsFactors = FALSE
)
if (identical(args, "--list")) {
  write.table(branch_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

valid_flags <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") &
  !args %in% valid_flags &
  !grepl("^--branch=", args)]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
branch_args <- args[grepl("^--branch=", args)]
if (length(branch_args) > 1L) stop("Specify --branch only once.")
branch <- if (length(branch_args)) {
  toupper(sub("^--branch=", "", branch_args))
} else {
  task_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")))
  if (is.na(task_id) || !task_id %in% branch_map$task_id) {
    stop("Specify --branch=RL|VZ or set SLURM_ARRAY_TASK_ID to 1 or 2.")
  }
  branch_map$branch[match(task_id, branch_map$task_id)]
}
if (!branch %in% branch_map$branch) stop("--branch must be RL or VZ.")
args <- args[!grepl("^--branch=", args)]
if (any(!args %in% valid_flags)) {
  stop(
    "Usage: Rscript scripts/xenium_vz_rl_plot_rl_vz_marker_dotplot.R ",
    "--branch=RL|VZ [--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

resolution_label <- "0.8"
cluster_column <- paste0("Xenium_snn_res.", resolution_label)
branch_lower <- tolower(branch)
output_root <- here(config$project$outputs_dir)
input_path <- file.path(
  output_root, "xenium", branch_lower, "05_post_qc", "rds",
  paste0("Xenium_", branch, "_postQC.rds")
)
diagnostic_root <- file.path(
  output_root, "xenium", "vz_rl", "diagnostics", "rl_vz_marker_dotplots"
)
plot_dir <- file.path(diagnostic_root, "plots")
table_dir <- file.path(diagnostic_root, "tables")
file_prefix <- paste0("Xenium_", branch, "_Res0.8_RLVZMarker")
plot_tiff_path <- file.path(plot_dir, paste0(file_prefix, "_DotPlot.tif"))
plot_pdf_path <- file.path(plot_dir, paste0(file_prefix, "_DotPlot.pdf"))
availability_path <- file.path(
  table_dir, paste0(file_prefix, "_Availability.csv")
)
cluster_summary_path <- file.path(
  table_dir, paste0(file_prefix, "_ClusterSummary.csv")
)
manifest_path <- file.path(
  table_dir, paste0(file_prefix, "_Manifest.csv")
)
expected_outputs <- c(
  plot_tiff_path, plot_pdf_path, availability_path,
  cluster_summary_path, manifest_path
)

rl_vz_diagnostic_markers <- list(
  "RL VZ / roof plate core" = c(
    "LMX1A", "LMX1B", "WLS", "WNT1", "PAX6", "SOX2",
    "CALCB", "CRYAB", "COL2A1", "MSX1", "MSX2",
    "GDF7", "BMP4", "BMP7"
  ),
  "RL VZ proliferation" = c(
    "HES1", "NES", "MKI67", "TOP2A", "PCNA", "CCND1", "CCND2"
  ),
  "RL-SVZ / early RL output" = c(
    "ATOH1", "EOMES", "OTX2", "RSPO3", "LHX9", "BARHL1",
    "ZIC1", "ZIC4", "TRPC3", "BCAN", "ZFHX4"
  ),
  "BG / radial glia contrast" = c(
    "SLC1A3", "GDF10", "TNC", "PTPRZ1", "FABP7", "SOX9",
    "EGFR", "PAX3", "HOPX", "ATP1A2", "MEGF10", "NDRG2"
  ),
  "Cerebellar VZ / GABA contrast" = c(
    "PTF1A", "PRDM13", "ASCL1", "PAX2", "GAD1", "GAD2"
  ),
  "Ependymal / choroid contrast" = c(
    "FOXJ1", "PIFO", "TTR", "AQP1", "KCNJ13", "KRT8", "KRT18"
  )
)
validate_focused_marker_sets(rl_vz_diagnostic_markers)

if (dry_run) {
  compact_dry_run(
    paste0(branch, " resolution-0.8 RL-VZ marker DotPlot"),
    inputs = input_path,
    outputs = expected_outputs,
    checks = c(
      branch_is_supported = branch %in% c("RL", "VZ"),
      resolution_is_0.8 = identical(resolution_label, "0.8"),
      marker_panel_is_valid = TRUE
    )
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop(branch, " post-QC input not found: ", input_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing ", branch, " RL-VZ marker outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(config$runtime$random_seed)

message(Sys.time(), ": loading ", branch, " post-QC object...")
obj <- readRDS(input_path)
if (anyDuplicated(colnames(obj))) stop(branch, " object contains duplicate cell IDs.")
if (!"Xenium" %in% Assays(obj)) stop(branch, " object lacks the Xenium assay.")
required_metadata <- c("orig.ident", cluster_column)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) {
  stop(branch, " object lacks metadata: ", paste(missing_metadata, collapse = ", "))
}
cluster_values <- as.character(obj[[cluster_column, drop = TRUE]])
if (anyNA(cluster_values) || any(!nzchar(cluster_values))) {
  stop(branch, " resolution-0.8 cluster assignments contain missing values.")
}
numeric_clusters <- suppressWarnings(as.numeric(cluster_values))
if (anyNA(numeric_clusters)) {
  stop(branch, " resolution-0.8 cluster assignments are not numeric.")
}
cluster_levels <- as.character(sort(unique(numeric_clusters)))
obj[[cluster_column]] <- factor(cluster_values, levels = cluster_levels)
DefaultAssay(obj) <- "Xenium"

cluster_summary <- data.frame(
  branch = branch,
  resolution = resolution_label,
  cluster = cluster_levels,
  cell_count = vapply(
    cluster_levels,
    function(cluster) sum(cluster_values == cluster),
    integer(1)
  ),
  sample_count = vapply(
    cluster_levels,
    function(cluster) {
      length(unique(as.character(obj$orig.ident)[cluster_values == cluster]))
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(cluster_summary, cluster_summary_path, row.names = FALSE)

write_focused_marker_dotplot(
  object = obj,
  marker_sets = rl_vz_diagnostic_markers,
  assay = "Xenium",
  group_by = cluster_column,
  branch = paste0(branch, " RL-VZ diagnostic"),
  resolution_label = resolution_label,
  availability_path = availability_path,
  tiff_path = plot_tiff_path,
  pdf_path = plot_pdf_path
)

input_info <- file.info(input_path)
manifest <- data.frame(
  branch = branch,
  input_rds = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime),
  cells = ncol(obj),
  samples = length(unique(as.character(obj$orig.ident))),
  cluster_column = cluster_column,
  clusters = paste(cluster_levels, collapse = "|"),
  marker_panel = "rl_vz_diagnostic_markers",
  marker_signature = focused_marker_set_signature(rl_vz_diagnostic_markers),
  random_seed = config$runtime$random_seed,
  writes_rds = FALSE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
write.csv(manifest, manifest_path, row.names = FALSE)

missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing_outputs)) {
  stop("RL-VZ marker diagnostic did not create output(s):\n- ",
       paste(missing_outputs, collapse = "\n- "))
}
message(branch, " resolution-0.8 RL-VZ marker DotPlot is ready: ", plot_tiff_path)
