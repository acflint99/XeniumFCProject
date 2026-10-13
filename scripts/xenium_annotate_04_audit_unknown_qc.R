#!/usr/bin/env Rscript

# Descriptive QC audit for cells assigned to an Unknown-## consensus label.
# This script does not filter, relabel, or save a modified Seurat object.
# Production users should invoke xenium_unknown_01_audit_qc.R and
# xenium_unknown_02_compile_qc_review.R. This multi-mode driver remains shared
# with the established resolution-pilot launchers.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))

config <- load_pipeline_config()
pilot_manifest <- load_resolution2_pilot_manifest(config)
sample_manifest <- load_sample_manifest(config)
production_paths <- get_production_output_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--pilot-res3", "--pilot-res4", "--pilot-res5", "--all-samples-res5",
  "--weighted-2of3", "--combine", "--list", "--dry-run", "--overwrite"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

pilot_res3 <- "--pilot-res3" %in% args
pilot_res4 <- "--pilot-res4" %in% args
pilot_res5 <- "--pilot-res5" %in% args
all_samples_res5 <- "--all-samples-res5" %in% args
mode_flags <- c(pilot_res3, pilot_res4, pilot_res5, all_samples_res5)
if (sum(mode_flags) != 1L) {
  stop(
    "Choose exactly one of --pilot-res3, --pilot-res4, --pilot-res5, ",
    "or --all-samples-res5."
  )
}
weighted_2of3 <- "--weighted-2of3" %in% args
combine <- "--combine" %in% args
if (all_samples_res5 && !weighted_2of3) {
  stop("--all-samples-res5 requires --weighted-2of3 for the production consensus.")
}
resolution_tag <- if (pilot_res5) "5.0" else if (pilot_res4) "4.0" else "3.0"
pilot_stage <- if (pilot_res5) {
  "03e_resolution5_pilot"
} else if (pilot_res4) {
  "03d_resolution4_pilot"
} else {
  "03c_resolution3_pilot"
}
consensus_stage <- if (weighted_2of3) {
  "03_consensus_labels_weighted_2of3"
} else {
  "03_consensus_labels"
}

task_samples <- if (all_samples_res5) sample_manifest$sample_id else pilot_manifest$sample_id
task_map <- data.frame(
  task_id = seq_along(task_samples),
  sample_id = task_samples,
  stringsAsFactors = FALSE
)

output_root <- here(config$project$outputs_dir)
if (all_samples_res5) {
  annotation_root <- production_paths$annotation_root
  consensus_stage <- config$production$consensus_label_stage
  input_dir <- production_paths$consensus_rds_dir
  audit_root <- file.path(annotation_root, "04_unknown_qc_review")
} else {
  annotation_root <- file.path(
    output_root, "xenium", "preprocess", pilot_stage, "annotation"
  )
  input_dir <- file.path(annotation_root, consensus_stage, "rds")
  audit_root <- file.path(annotation_root, consensus_stage, "unknown_qc_audit")
}
table_dir <- file.path(audit_root, "tables")
plot_dir <- file.path(audit_root, "plots")

output_paths_for <- function(sample_name) {
  c(
    group_summary = file.path(
      table_dir, paste0(sample_name, "_Unknown_vs_Known_QC_summary.csv")
    ),
    statistical_summary = file.path(
      table_dir, paste0(sample_name, "_Unknown_vs_Known_QC_effects.csv")
    ),
    cluster_summary = file.path(
      table_dir, paste0(sample_name, "_Unknown_cluster_QC_summary.csv")
    ),
    thresholds = file.path(
      table_dir, paste0(sample_name, "_QC_thresholds_and_known_tails.csv")
    ),
    provenance = file.path(
      table_dir, paste0(sample_name, "_Unknown_QC_audit_provenance.csv")
    ),
    plot = file.path(
      plot_dir, paste0(sample_name, "_Unknown_vs_Known_QC.pdf")
    )
  )
}

resolve_xenium_cells_path <- function(input_dir, must_exist = TRUE) {
  candidates <- file.path(input_dir, c("cells.csv.gz", "cells.csv"))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 1L) {
    stop(
      "Expected one Xenium cell table but found both compressed and ",
      "uncompressed files:\n- ", paste(existing, collapse = "\n- ")
    )
  }
  if (!length(existing)) {
    if (must_exist) {
      stop(
        "Xenium cell table not found; expected one of:\n- ",
        paste(candidates, collapse = "\n- ")
      )
    }
    return(candidates[[1]])
  }
  existing[[1]]
}

read_xenium_cell_morphology <- function(path) {
  required <- c("cell_id", "cell_area", "nucleus_area", "nucleus_count")
  header <- names(utils::read.csv(
    path, nrows = 0L, check.names = FALSE, stringsAsFactors = FALSE
  ))
  missing <- setdiff(required, header)
  if (length(missing)) {
    stop(
      "Xenium cell table lacks required morphology column(s): ",
      paste(missing, collapse = ", "), ". File: ", path
    )
  }
  column_classes <- rep("NULL", length(header))
  column_classes[match(required, header)] <- c(
    "character", "numeric", "numeric", "numeric"
  )
  morphology <- utils::read.csv(
    path,
    colClasses = column_classes,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (anyNA(morphology$cell_id) || any(!nzchar(morphology$cell_id))) {
    stop("Xenium cell table contains missing or blank cell_id values: ", path)
  }
  if (anyDuplicated(morphology$cell_id)) {
    stop("Xenium cell table contains duplicate cell_id values: ", path)
  }
  if (any(!is.finite(morphology$cell_area)) ||
      any(!is.finite(morphology$nucleus_count))) {
    stop(
      "Xenium cell table contains non-finite cell_area or nucleus_count ",
      "values: ", path
    )
  }
  if (any(is.infinite(morphology$nucleus_area))) {
    stop("Xenium cell table contains infinite nucleus_area values: ", path)
  }
  finite_nucleus_area <- is.finite(morphology$nucleus_area)
  if (any(morphology$cell_area <= 0) ||
      any(morphology$nucleus_area[finite_nucleus_area] < 0) ||
      any(morphology$nucleus_count < 0) ||
      any(abs(morphology$nucleus_count - round(morphology$nucleus_count)) > 1e-8)) {
    stop(
      "Xenium morphology values must have positive cell_area, nonnegative ",
      "finite nucleus_area when present, and nonnegative integer ",
      "nucleus_count: ", path
    )
  }
  morphology
}

combined_paths <- c(
  evidence = file.path(table_dir, "Xenium_Unknown_QC_Evidence.csv"),
  review_template = file.path(table_dir, "Xenium_Unknown_QC_Review_Template.csv"),
  provenance = file.path(table_dir, "Xenium_Unknown_QC_Audit_manifest.csv")
)

if ("--list" %in% args) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept TASK_ID.")
  cat(
    "Mode: ",
    if (all_samples_res5) "resolution-5 production" else paste0("resolution-", resolution_tag, " pilot"),
    "\n", sep = ""
  )
  cat(
    "Consensus method: ",
    if (weighted_2of3) "weighted_2of3" else "legacy_six_vote",
    "\n",
    sep = ""
  )
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

task_args <- args[!args %in% valid_options]
if (combine) {
  if (length(task_args)) stop("--combine does not accept TASK_ID.")
  per_sample_paths <- lapply(task_samples, output_paths_for)
  cluster_paths <- vapply(per_sample_paths, `[[`, character(1), "cluster_summary")
  provenance_paths <- vapply(per_sample_paths, `[[`, character(1), "provenance")
  combine_inputs <- c(cluster_paths, provenance_paths)
  if (dry_run) {
    compact_dry_run(
      "Production Unknown-QC audit combine",
      inputs = combine_inputs,
      outputs = combined_paths,
      checks = c(all_samples_res5 = all_samples_res5)
    )
    quit(save = "no", status = 0L)
  }
  if (!all_samples_res5) stop("--combine is supported only with --all-samples-res5.")
  missing_inputs <- combine_inputs[!file.exists(combine_inputs)]
  if (length(missing_inputs)) {
    stop("Cannot combine Unknown-QC audit; missing files:\n- ", paste(missing_inputs, collapse = "\n- "))
  }
  existing <- combined_paths[file.exists(combined_paths)]
  if (length(existing) && !overwrite) {
    stop(
      "Refusing to overwrite combined Unknown-QC audit outputs:\n- ",
      paste(existing, collapse = "\n- "), "\nUse --overwrite only after review."
    )
  }
  evidence <- do.call(rbind, lapply(cluster_paths, function(path) {
    read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  }))
  provenance <- do.call(rbind, lapply(provenance_paths, function(path) {
    read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  }))
  if (!identical(unique(as.character(provenance$sample_id)), task_samples)) {
    stop("Unknown-QC audit provenance is not in config/samples.csv order.")
  }
  evidence <- evidence[order(
    match(evidence$sample_id, task_samples), evidence$consensus_label
  ), , drop = FALSE]
  keys <- unknown_qc_review_key(evidence)
  if (anyDuplicated(keys)) stop("Combined Unknown-QC evidence contains duplicate keys.")
  review_template <- evidence
  review_template$decision <- ""
  review_template$review_notes <- ""
  review_template <- review_template[c(
    unknown_qc_review_key_columns(), "decision", "review_notes",
    setdiff(names(review_template), c(
      unknown_qc_review_key_columns(), "decision", "review_notes"
    ))
  )]
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(evidence, combined_paths[["evidence"]], row.names = FALSE, na = "")
  write.csv(review_template, combined_paths[["review_template"]], row.names = FALSE, na = "")
  write.csv(provenance, combined_paths[["provenance"]], row.names = FALSE, na = "")
  message(
    "Combined ", nrow(evidence), " sample-qualified Unknown populations. Copy ",
    combined_paths[["review_template"]], " to ",
    resolve_config_path(config$manifests$reviewed_unknown_qc, config),
    " and fill decision plus review_notes."
  )
  quit(save = "no", status = 0L)
}

if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_04_audit_unknown_qc.R ",
    "--pilot-res3|--pilot-res4|--pilot-res5|--all-samples-res5 ",
    "[--weighted-2of3] [--combine|--list|--dry-run|--overwrite] [TASK_ID]"
  )
}
task_value <- if (length(task_args) == 1L) {
  task_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop(
    "TASK_ID must be between 1 and ", nrow(task_map),
    ". Use --list to inspect the mapping."
  )
}

sample_name <- task_map$sample_id[[task_id]]
sample_matches <- which(sample_manifest$sample_id == sample_name)
if (length(sample_matches) != 1L) {
  stop(
    "Expected exactly one config/samples.csv row for ", sample_name,
    "; found ", length(sample_matches), "."
  )
}
sample_record <- sample_manifest[sample_matches, , drop = FALSE]
sample_paths <- resolve_sample_paths(sample_record, config)
xenium_cells_path <- resolve_xenium_cells_path(
  sample_paths$input_dir, must_exist = !dry_run
)
input_path <- file.path(
  input_dir, paste0(sample_name, "_Consensus_annotated.rds")
)
output_paths <- output_paths_for(sample_name)

if (dry_run) {
  compact_dry_run(
    paste0("Unknown-QC audit task ", task_id, " [", sample_name, "]"),
    inputs = c(input_path, xenium_cells_path),
    outputs = output_paths
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("Missing consensus object: ", input_path)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing Unknown-QC audit outputs for ", sample_name,
    ":\n- ", paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}

message("Loading consensus object: ", input_path)
obj <- readRDS(input_path)
cell_ids <- Cells(obj)
metadata <- obj[[]]
if (anyDuplicated(cell_ids)) stop("Duplicate cell IDs in object: ", sample_name)
if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
  stop("Metadata row names are missing or duplicated: ", sample_name)
}
if (!setequal(cell_ids, rownames(metadata))) {
  stop("Object cell IDs and metadata row names do not match: ", sample_name)
}
metadata <- metadata[cell_ids, , drop = FALSE]

base_required_columns <- c(
  "seurat_clusters", "consensus_label", "nCount_Xenium",
  "nFeature_Xenium", "cell_area", "percent_control"
)
missing_base_columns <- setdiff(base_required_columns, colnames(metadata))
if (length(missing_base_columns)) {
  stop(
    "Consensus object is missing required metadata columns: ",
    paste(missing_base_columns, collapse = ", ")
  )
}

message("Loading Xenium morphology: ", xenium_cells_path)
morphology <- read_xenium_cell_morphology(xenium_cells_path)
morphology_index <- match(cell_ids, morphology$cell_id)
if (anyNA(morphology_index)) {
  missing_cell_ids <- cell_ids[is.na(morphology_index)]
  stop(
    "Consensus object contains ", length(missing_cell_ids),
    " cell ID(s) absent from the Xenium cell table for ", sample_name,
    ". First missing IDs: ",
    paste(utils::head(missing_cell_ids, 10L), collapse = ", ")
  )
}
morphology <- morphology[morphology_index, , drop = FALSE]
if (!identical(cell_ids, morphology$cell_id)) {
  stop("Xenium morphology join did not preserve cell order for ", sample_name, ".")
}

object_cell_area <- suppressWarnings(as.numeric(metadata$cell_area))
raw_cell_area <- morphology$cell_area
cell_area_absolute_difference <- abs(object_cell_area - raw_cell_area)
# cerebellum_cells_stats.csv stores Area (µm^2) rounded to two decimals,
# whereas the raw Xenium cells table retains additional decimal places. The
# largest legitimate rounding difference is therefore 0.005 µm^2; the small
# numerical margin avoids rejecting a boundary value because of binary
# floating-point representation.
cell_area_rounding_tolerance <- 0.00501
cell_area_matches_raw <- is.finite(object_cell_area) &
  cell_area_absolute_difference <= cell_area_rounding_tolerance
if (any(!cell_area_matches_raw)) {
  mismatch <- which(!cell_area_matches_raw)
  stop(
    "Existing cell_area does not match raw Xenium cell_area for ",
    length(mismatch), " cell(s) in ", sample_name,
    " within the allowed two-decimal rounding tolerance of ",
    cell_area_rounding_tolerance, " square micrometers. Refusing to mix area ",
    "sources. First mismatched cell IDs: ",
    paste(utils::head(cell_ids[mismatch], 10L), collapse = ", ")
  )
}

metadata$nucleus_area <- morphology$nucleus_area
metadata$nucleus_count <- morphology$nucleus_count
metadata$cell_to_nucleus_area_ratio <- ifelse(
  morphology$nucleus_area > 0,
  raw_cell_area / morphology$nucleus_area,
  NA_real_
)
metadata$transcripts_per_cell_area <- suppressWarnings(
  as.numeric(metadata$nCount_Xenium) / raw_cell_area
)

qc_metrics <- c(
  "nCount_Xenium", "nFeature_Xenium", "cell_area", "nucleus_area",
  "cell_to_nucleus_area_ratio", "transcripts_per_cell_area",
  "percent_control"
)
required_columns <- c(
  "seurat_clusters", "consensus_label", qc_metrics, "nucleus_count"
)
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns)) {
  stop(
    "Consensus object is missing required metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

labels <- trimws(as.character(metadata$consensus_label))
if (anyNA(labels) || any(!nzchar(labels))) {
  stop("Blank or missing consensus labels in ", sample_name, ".")
}
unknown_cells <- is_unknown_consensus_label(labels)
has_unknown <- any(unknown_cells)
if (!any(!unknown_cells)) {
  stop("All cells have an Unknown-## consensus label in ", sample_name, ".")
}

qc_data <- data.frame(
  cell_id = cell_ids,
  seurat_clusters = as.character(metadata$seurat_clusters),
  consensus_label = labels,
  qc_group = ifelse(unknown_cells, "Unknown", "Known"),
  stringsAsFactors = FALSE
)
for (metric in qc_metrics) {
  values <- suppressWarnings(as.numeric(metadata[[metric]]))
  invalid <- !is.finite(values) & !is.na(values)
  if (any(invalid)) {
    stop("Non-finite values in ", metric, " for ", sample_name, ".")
  }
  metrics_allowing_missing <- c(
    "nucleus_area", "cell_to_nucleus_area_ratio"
  )
  if (anyNA(values) && !metric %in% metrics_allowing_missing) {
    stop("Missing values in ", metric, " for ", sample_name, ".")
  }
  qc_data[[metric]] <- values
}
qc_data$nucleus_count <- suppressWarnings(as.numeric(metadata$nucleus_count))
if (any(!is.finite(qc_data$nucleus_count))) {
  stop("Non-finite values in nucleus_count for ", sample_name, ".")
}

summarize_values <- function(values) {
  finite_values <- values[is.finite(values)]
  if (!length(finite_values)) {
    return(c(
      n = length(values), mean = NA_real_, sd = NA_real_, q10 = NA_real_,
      q25 = NA_real_, median = NA_real_, q75 = NA_real_, q90 = NA_real_
    ))
  }
  c(
    n = length(values), mean = mean(finite_values), sd = stats::sd(finite_values),
    q10 = unname(stats::quantile(finite_values, 0.10)),
    q25 = unname(stats::quantile(finite_values, 0.25)),
    median = stats::median(finite_values),
    q75 = unname(stats::quantile(finite_values, 0.75)),
    q90 = unname(stats::quantile(finite_values, 0.90))
  )
}

metric_is_informative <- function(values) {
  finite_values <- values[is.finite(values)]
  length(finite_values) > 0L && length(unique(finite_values)) > 1L
}

finite_median <- function(values) {
  finite_values <- values[is.finite(values)]
  if (length(finite_values)) stats::median(finite_values) else NA_real_
}

tail_fraction <- function(flags) {
  if (all(is.na(flags))) NA_real_ else mean(flags, na.rm = TRUE)
}

if (!has_unknown) {
  group_summary <- do.call(rbind, lapply(qc_metrics, function(metric) {
    data.frame(
      sample = sample_name, qc_group = "Known", metric = metric,
      as.list(summarize_values(qc_data[[metric]])),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
  effect_summary <- data.frame(
    sample = character(), metric = character(), median_known = numeric(),
    median_unknown = numeric(), median_difference_unknown_minus_known = numeric(),
    rank_biserial_unknown_vs_known = numeric(), wilcoxon_p = numeric(),
    known_tail_direction = character(), known_tail_threshold = numeric(),
    known_tail_fraction = numeric(), unknown_tail_fraction = numeric(),
    metric_is_informative = logical(), wilcoxon_p_bh = numeric(),
    stringsAsFactors = FALSE
  )
  cluster_summary <- data.frame(
    sample_id = character(), seurat_clusters = character(),
    consensus_label = character(), n_cells = integer(),
    fraction_of_sample = numeric(), median_nCount_Xenium = numeric(),
    median_nFeature_Xenium = numeric(), median_cell_area = numeric(),
    median_nucleus_area = numeric(),
    median_cell_to_nucleus_area_ratio = numeric(),
    median_transcripts_per_cell_area = numeric(),
    median_nucleus_count = numeric(), median_percent_control = numeric(),
    fraction_low_nCount_vs_known_q10 = numeric(),
    fraction_low_nFeature_vs_known_q10 = numeric(),
    fraction_low_cell_area_vs_known_q10 = numeric(),
    fraction_low_nucleus_area_vs_known_q10 = numeric(),
    fraction_low_transcript_density_vs_known_q10 = numeric(),
    fraction_low_cell_to_nucleus_ratio_vs_known_q10 = numeric(),
    fraction_high_cell_to_nucleus_ratio_vs_known_q90 = numeric(),
    fraction_nucleus_area_missing = numeric(),
    fraction_nucleus_area_missing_with_nuclei = numeric(),
    fraction_nucleus_area_equal_zero = numeric(),
    fraction_nucleus_count_zero = numeric(),
    fraction_nucleus_count_not_one = numeric(),
    fraction_high_percent_control_vs_known_q90 = numeric(),
    cell_area_unique_values = integer(),
    fraction_cell_area_equal_20 = numeric(),
    cell_area_is_informative = logical(),
    nucleus_area_is_informative = logical(),
    cell_to_nucleus_area_ratio_is_informative = logical(),
    transcripts_per_cell_area_is_informative = logical(),
    percent_control_is_informative = logical(), stringsAsFactors = FALSE
  )
  tail_direction <- c(
    nCount_Xenium = "low", nFeature_Xenium = "low",
    cell_area = "low", nucleus_area = "low",
    cell_to_nucleus_area_ratio = "high",
    transcripts_per_cell_area = "low", percent_control = "high"
  )
  known_tail_thresholds <- vapply(qc_metrics, function(metric) {
    probability <- if (tail_direction[[metric]] == "low") 0.10 else 0.90
    unname(stats::quantile(qc_data[[metric]], probability, na.rm = TRUE))
  }, numeric(1))
  metric_informative <- vapply(
    qc_metrics, function(metric) metric_is_informative(qc_data[[metric]]), logical(1)
  )
  threshold_rows <- data.frame(
    sample = sample_name,
    threshold_source = "Known cells in consensus object",
    metric = qc_metrics,
    boundary = ifelse(tail_direction[qc_metrics] == "low", "q10", "q90"),
    value = unname(known_tail_thresholds[qc_metrics]),
    is_informative = unname(metric_informative[qc_metrics]),
    stringsAsFactors = FALSE
  )
  ratio_q10 <- unname(stats::quantile(
    qc_data$cell_to_nucleus_area_ratio, 0.10, na.rm = TRUE
  ))
  threshold_rows <- rbind(
    threshold_rows,
    data.frame(
      sample = sample_name,
      threshold_source = "Known cells in consensus object",
      metric = "cell_to_nucleus_area_ratio",
      boundary = "q10",
      value = ratio_q10,
      is_informative = metric_informative[["cell_to_nucleus_area_ratio"]],
      stringsAsFactors = FALSE
    )
  )
  no_unknown_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No Unknown consensus populations") +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = paste(sample_name, "Unknown-QC audit")) +
    theme_void()
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(output_paths[["plot"]], width = 9, height = 7)
  print(no_unknown_plot)
  grDevices::dev.off()
  write.csv(group_summary, output_paths[["group_summary"]], row.names = FALSE)
  write.csv(effect_summary, output_paths[["statistical_summary"]], row.names = FALSE)
  write.csv(cluster_summary, output_paths[["cluster_summary"]], row.names = FALSE)
  write.csv(threshold_rows, output_paths[["thresholds"]], row.names = FALSE)
  input_info <- file.info(input_path)
  xenium_cells_info <- file.info(xenium_cells_path)
  provenance <- data.frame(
    sample_id = sample_name,
    input_path = normalizePath(input_path, mustWork = TRUE),
    input_size = as.numeric(input_info$size),
    input_mtime = as.numeric(input_info$mtime),
    xenium_cells_path = normalizePath(xenium_cells_path, mustWork = TRUE),
    xenium_cells_size = as.numeric(xenium_cells_info$size),
    xenium_cells_mtime = as.numeric(xenium_cells_info$mtime),
    input_cells = ncol(obj), unknown_cells = 0L, unknown_populations = 0L,
    morphology_cells_matched = length(morphology_index),
    cell_area_matches_raw = all(cell_area_matches_raw),
    cell_area_validation_tolerance = cell_area_rounding_tolerance,
    max_cell_area_absolute_difference = max(cell_area_absolute_difference),
    nucleus_area_available_fraction = mean(is.finite(metadata$nucleus_area)),
    nucleus_count_zero_fraction = mean(metadata$nucleus_count == 0),
    cell_area_is_informative = metric_informative[["cell_area"]],
    nucleus_area_is_informative = metric_informative[["nucleus_area"]],
    percent_control_is_informative = metric_informative[["percent_control"]],
    stringsAsFactors = FALSE
  )
  write.csv(provenance, output_paths[["provenance"]], row.names = FALSE)
  message("No Unknown consensus populations in ", sample_name, "; wrote an empty audit.")
  quit(save = "no", status = 0L)
}

group_rows <- list()
row_index <- 1L
for (group_name in c("Known", "Unknown")) {
  group_data <- qc_data[qc_data$qc_group == group_name, , drop = FALSE]
  for (metric in qc_metrics) {
    group_rows[[row_index]] <- data.frame(
      sample = sample_name,
      qc_group = group_name,
      metric = metric,
      as.list(summarize_values(group_data[[metric]])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    row_index <- row_index + 1L
  }
}
group_summary <- do.call(rbind, group_rows)

known_data <- qc_data[qc_data$qc_group == "Known", , drop = FALSE]
unknown_data <- qc_data[qc_data$qc_group == "Unknown", , drop = FALSE]
cell_area_unique_values <- length(unique(qc_data$cell_area))
fraction_cell_area_equal_20 <- mean(qc_data$cell_area == 20)
cell_area_is_informative <- cell_area_unique_values > 1L &&
  fraction_cell_area_equal_20 < 0.95
tail_direction <- c(
  nCount_Xenium = "low", nFeature_Xenium = "low",
  cell_area = "low", nucleus_area = "low",
  cell_to_nucleus_area_ratio = "high",
  transcripts_per_cell_area = "low", percent_control = "high"
)
known_tail_thresholds <- vapply(qc_metrics, function(metric) {
  probability <- if (tail_direction[[metric]] == "low") 0.10 else 0.90
  unname(stats::quantile(known_data[[metric]], probability, na.rm = TRUE))
}, numeric(1))
ratio_known_q10 <- unname(stats::quantile(
  known_data$cell_to_nucleus_area_ratio, 0.10, na.rm = TRUE
))
metric_informative <- vapply(
  qc_metrics, function(metric) metric_is_informative(qc_data[[metric]]), logical(1)
)

effect_rows <- lapply(qc_metrics, function(metric) {
  unknown_values <- unknown_data[[metric]]
  known_values <- known_data[[metric]]
  unknown_finite <- unknown_values[is.finite(unknown_values)]
  known_finite <- known_values[is.finite(known_values)]
  if (length(unknown_finite) && length(known_finite)) {
    test <- suppressWarnings(stats::wilcox.test(
      unknown_finite, known_finite, exact = FALSE
    ))
    u_statistic <- unname(test$statistic)
    rank_biserial <- 2 * u_statistic /
      (length(unknown_finite) * length(known_finite)) - 1
    wilcoxon_p <- test$p.value
  } else {
    rank_biserial <- NA_real_
    wilcoxon_p <- NA_real_
  }
  threshold <- known_tail_thresholds[[metric]]
  if (!metric_informative[[metric]]) {
    unknown_tail <- rep(NA, length(unknown_values))
    known_tail <- rep(NA, length(known_values))
  } else if (tail_direction[[metric]] == "low") {
    unknown_tail <- unknown_values <= threshold
    known_tail <- known_values <= threshold
  } else {
    unknown_tail <- unknown_values > threshold
    known_tail <- known_values > threshold
  }
  data.frame(
    sample = sample_name,
    metric = metric,
    median_known = finite_median(known_values),
    median_unknown = finite_median(unknown_values),
    median_difference_unknown_minus_known =
      finite_median(unknown_values) - finite_median(known_values),
    rank_biserial_unknown_vs_known = rank_biserial,
    wilcoxon_p = wilcoxon_p,
    known_tail_direction = tail_direction[[metric]],
    known_tail_threshold = threshold,
    known_tail_fraction = tail_fraction(known_tail),
    unknown_tail_fraction = tail_fraction(unknown_tail),
    metric_is_informative = metric_informative[[metric]],
    stringsAsFactors = FALSE
  )
})
effect_summary <- do.call(rbind, effect_rows)
effect_summary$wilcoxon_p_bh <- stats::p.adjust(
  effect_summary$wilcoxon_p, method = "BH"
)

for (metric in qc_metrics) {
  threshold <- known_tail_thresholds[[metric]]
  flag_name <- paste0(metric, "_known_tail")
  if (!metric_informative[[metric]]) {
    qc_data[[flag_name]] <- rep(NA, nrow(qc_data))
  } else if (tail_direction[[metric]] == "low") {
    qc_data[[flag_name]] <- qc_data[[metric]] <= threshold
  } else {
    qc_data[[flag_name]] <- qc_data[[metric]] > threshold
  }
}
if (!metric_informative[["cell_to_nucleus_area_ratio"]]) {
  qc_data$cell_to_nucleus_area_ratio_known_low_tail <- rep(NA, nrow(qc_data))
} else {
  qc_data$cell_to_nucleus_area_ratio_known_low_tail <-
    qc_data$cell_to_nucleus_area_ratio <= ratio_known_q10
}

cluster_rows <- lapply(sort(unique(unknown_data$seurat_clusters)), function(cluster_id) {
  cluster_data <- qc_data[
    qc_data$qc_group == "Unknown" & qc_data$seurat_clusters == cluster_id,
    , drop = FALSE
  ]
  cluster_labels <- unique(cluster_data$consensus_label)
  if (length(cluster_labels) != 1L || !is_unknown_consensus_label(cluster_labels)) {
    stop("Unexpected cluster-to-consensus cardinality for cluster ", cluster_id, ".")
  }
  data.frame(
    sample_id = sample_name,
    seurat_clusters = cluster_id,
    consensus_label = cluster_labels,
    n_cells = nrow(cluster_data),
    fraction_of_sample = nrow(cluster_data) / nrow(qc_data),
    median_nCount_Xenium = median(cluster_data$nCount_Xenium),
    median_nFeature_Xenium = median(cluster_data$nFeature_Xenium),
    median_cell_area = median(cluster_data$cell_area),
    median_nucleus_area = median(cluster_data$nucleus_area),
    median_cell_to_nucleus_area_ratio = finite_median(
      cluster_data$cell_to_nucleus_area_ratio
    ),
    median_transcripts_per_cell_area = median(
      cluster_data$transcripts_per_cell_area
    ),
    median_nucleus_count = median(cluster_data$nucleus_count),
    median_percent_control = median(cluster_data$percent_control),
    fraction_low_nCount_vs_known_q10 = tail_fraction(
      cluster_data$nCount_Xenium_known_tail
    ),
    fraction_low_nFeature_vs_known_q10 = tail_fraction(
      cluster_data$nFeature_Xenium_known_tail
    ),
    fraction_low_cell_area_vs_known_q10 = tail_fraction(
      cluster_data$cell_area_known_tail
    ),
    fraction_low_nucleus_area_vs_known_q10 = tail_fraction(
      cluster_data$nucleus_area_known_tail
    ),
    fraction_low_transcript_density_vs_known_q10 = tail_fraction(
      cluster_data$transcripts_per_cell_area_known_tail
    ),
    fraction_low_cell_to_nucleus_ratio_vs_known_q10 = tail_fraction(
      cluster_data$cell_to_nucleus_area_ratio_known_low_tail
    ),
    fraction_high_cell_to_nucleus_ratio_vs_known_q90 = tail_fraction(
      cluster_data$cell_to_nucleus_area_ratio_known_tail
    ),
    fraction_nucleus_area_missing = mean(!is.finite(cluster_data$nucleus_area)),
    fraction_nucleus_area_missing_with_nuclei = mean(
      !is.finite(cluster_data$nucleus_area) & cluster_data$nucleus_count > 0
    ),
    fraction_nucleus_area_equal_zero = mean(
      is.finite(cluster_data$nucleus_area) & cluster_data$nucleus_area == 0
    ),
    fraction_nucleus_count_zero = mean(cluster_data$nucleus_count == 0),
    fraction_nucleus_count_not_one = mean(cluster_data$nucleus_count != 1),
    fraction_high_percent_control_vs_known_q90 = tail_fraction(
      cluster_data$percent_control_known_tail
    ),
    cell_area_unique_values = cell_area_unique_values,
    fraction_cell_area_equal_20 = fraction_cell_area_equal_20,
    cell_area_is_informative = cell_area_is_informative,
    nucleus_area_is_informative = metric_informative[["nucleus_area"]],
    cell_to_nucleus_area_ratio_is_informative =
      metric_informative[["cell_to_nucleus_area_ratio"]],
    transcripts_per_cell_area_is_informative =
      metric_informative[["transcripts_per_cell_area"]],
    percent_control_is_informative = metric_informative[["percent_control"]],
    stringsAsFactors = FALSE
  )
})
cluster_summary <- do.call(rbind, cluster_rows)

threshold_rows <- data.frame(
  sample = sample_name,
  threshold_source = "Known cells in consensus object",
  metric = qc_metrics,
  boundary = ifelse(tail_direction[qc_metrics] == "low", "q10", "q90"),
  value = unname(known_tail_thresholds[qc_metrics]),
  is_informative = unname(metric_informative[qc_metrics]),
  stringsAsFactors = FALSE
)
threshold_rows <- rbind(
  threshold_rows,
  data.frame(
    sample = sample_name,
    threshold_source = "Known cells in consensus object",
    metric = "cell_to_nucleus_area_ratio",
    boundary = "q10",
    value = ratio_known_q10,
    is_informative = metric_informative[["cell_to_nucleus_area_ratio"]],
    stringsAsFactors = FALSE
  )
)
stored_thresholds <- obj@misc$QC_thresholds
if (!is.null(stored_thresholds) && length(stored_thresholds)) {
  stored_values <- unlist(stored_thresholds, recursive = TRUE, use.names = TRUE)
  stored_values <- suppressWarnings(as.numeric(stored_values))
  stored_names <- names(unlist(stored_thresholds, recursive = TRUE, use.names = TRUE))
  keep <- is.finite(stored_values)
  if (any(keep)) {
    threshold_rows <- rbind(
      threshold_rows,
      data.frame(
        sample = sample_name,
        threshold_source = "object@misc$QC_thresholds",
        metric = stored_names[keep],
        boundary = "stored",
        value = stored_values[keep],
        is_informative = NA,
        stringsAsFactors = FALSE
      )
    )
  }
}

long_data <- do.call(rbind, lapply(qc_metrics, function(metric) {
  data.frame(
    qc_group = factor(qc_data$qc_group, levels = c("Known", "Unknown")),
    seurat_clusters = qc_data$seurat_clusters,
    metric = metric,
    value = qc_data[[metric]],
    stringsAsFactors = FALSE
  )
}))

p_group <- ggplot(long_data, aes(x = qc_group, y = value, fill = qc_group)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Known = "grey70", Unknown = "#D55E00")) +
  labs(
    title = paste(sample_name, "Unknown versus Known retained-cell QC"),
    subtitle = "Cells already passed the original QC filters",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

unknown_clusters <- sort(unique(unknown_data$seurat_clusters))
cluster_plot_data <- long_data[
  long_data$qc_group == "Unknown", , drop = FALSE
]
cluster_plot_data$comparison_group <- paste0(
  "Unknown cluster ", cluster_plot_data$seurat_clusters
)
known_plot_data <- long_data[long_data$qc_group == "Known", , drop = FALSE]
known_plot_data$comparison_group <- "Known cells"
cluster_plot_data <- rbind(known_plot_data, cluster_plot_data)
cluster_levels <- c("Known cells", paste0("Unknown cluster ", unknown_clusters))
cluster_plot_data$comparison_group <- factor(
  cluster_plot_data$comparison_group, levels = cluster_levels
)

p_cluster <- ggplot(
  cluster_plot_data,
  aes(x = comparison_group, y = value, fill = comparison_group)
) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Known cells" = "grey70",
    setNames(rep("#D55E00", length(unknown_clusters)),
             paste0("Unknown cluster ", unknown_clusters))
  )) +
  labs(
    title = paste(sample_name, "QC by Unknown cluster"),
    subtitle = "Known cells provide the within-sample reference distribution",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
grDevices::cairo_pdf(output_paths[["plot"]], width = 9, height = 7)
print(p_group)
print(p_cluster)
grDevices::dev.off()

utils::write.csv(group_summary, output_paths[["group_summary"]], row.names = FALSE)
utils::write.csv(
  effect_summary, output_paths[["statistical_summary"]], row.names = FALSE
)
utils::write.csv(
  cluster_summary, output_paths[["cluster_summary"]], row.names = FALSE
)
utils::write.csv(threshold_rows, output_paths[["thresholds"]], row.names = FALSE)
input_info <- file.info(input_path)
xenium_cells_info <- file.info(xenium_cells_path)
provenance <- data.frame(
  sample_id = sample_name,
  input_path = normalizePath(input_path, mustWork = TRUE),
  input_size = as.numeric(input_info$size),
  input_mtime = as.numeric(input_info$mtime),
  xenium_cells_path = normalizePath(xenium_cells_path, mustWork = TRUE),
  xenium_cells_size = as.numeric(xenium_cells_info$size),
  xenium_cells_mtime = as.numeric(xenium_cells_info$mtime),
  input_cells = ncol(obj),
  unknown_cells = nrow(unknown_data),
  unknown_populations = length(unique(unknown_data$consensus_label)),
  morphology_cells_matched = length(morphology_index),
  cell_area_matches_raw = all(cell_area_matches_raw),
  cell_area_validation_tolerance = cell_area_rounding_tolerance,
  max_cell_area_absolute_difference = max(cell_area_absolute_difference),
  nucleus_area_available_fraction = mean(is.finite(metadata$nucleus_area)),
  nucleus_count_zero_fraction = mean(metadata$nucleus_count == 0),
  cell_area_is_informative = cell_area_is_informative,
  nucleus_area_is_informative = metric_informative[["nucleus_area"]],
  percent_control_is_informative = metric_informative[["percent_control"]],
  stringsAsFactors = FALSE
)
utils::write.csv(provenance, output_paths[["provenance"]], row.names = FALSE)

message("Unknown cells: ", nrow(unknown_data), " / ", nrow(qc_data))
message("Unknown clusters: ", paste(unknown_clusters, collapse = ", "))
message("Wrote:\n- ", paste(output_paths, collapse = "\n- "))
message(
  "Interpretation note: cell-level P values are descriptive because cells within ",
  "a sample/cluster are not independent. Prioritize effect sizes and consistent ",
  "cluster-level QC shifts."
)
