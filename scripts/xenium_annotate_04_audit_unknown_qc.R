#!/usr/bin/env Rscript

# Descriptive QC audit for cells assigned to consensus_label = "Unknown".
# This script does not filter, relabel, or save a modified Seurat object.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
pilot_manifest <- load_resolution2_pilot_manifest(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--pilot-res3", "--pilot-res4", "--pilot-res5", "--weighted-2of3",
  "--list", "--dry-run", "--overwrite"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

pilot_res3 <- "--pilot-res3" %in% args
pilot_res4 <- "--pilot-res4" %in% args
pilot_res5 <- "--pilot-res5" %in% args
pilot_flags <- c(pilot_res3, pilot_res4, pilot_res5)
if (sum(pilot_flags) != 1L) {
  stop("Choose exactly one of --pilot-res3, --pilot-res4, or --pilot-res5.")
}
weighted_2of3 <- "--weighted-2of3" %in% args
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

task_map <- data.frame(
  task_id = seq_len(nrow(pilot_manifest)),
  sample_id = pilot_manifest$sample_id,
  stringsAsFactors = FALSE
)

if ("--list" %in% args) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept TASK_ID.")
  cat("Mode: resolution-", resolution_tag, " pilot\n", sep = "")
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
if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_04_audit_unknown_qc.R ",
    "--pilot-res3|--pilot-res4|--pilot-res5 [--weighted-2of3] ",
    "[--list|--dry-run|--overwrite] [TASK_ID]"
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
annotation_root <- file.path(
  here(config$project$outputs_dir), "xenium", "preprocess", pilot_stage,
  "annotation"
)
input_path <- file.path(
  annotation_root, consensus_stage, "rds",
  paste0(sample_name, "_Consensus_annotated.rds")
)
audit_root <- file.path(annotation_root, consensus_stage, "unknown_qc_audit")
table_dir <- file.path(audit_root, "tables")
plot_dir <- file.path(audit_root, "plots")
output_paths <- c(
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
  plot = file.path(
    plot_dir, paste0(sample_name, "_Unknown_vs_Known_QC.pdf")
  )
)

if (dry_run) {
  cat("Task ID:", task_id, "\n")
  cat("Sample:", sample_name, "\n")
  cat("Input:", input_path, "\n")
  cat("Expected outputs:\n")
  cat(paste0("- ", output_paths), sep = "\n")
  cat("\nDry-run is path-only; no RDS object was loaded.\n")
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

qc_metrics <- c(
  "nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control"
)
required_columns <- c("seurat_clusters", "consensus_label", qc_metrics)
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
if (!any(labels == "Unknown")) {
  stop("No cells have consensus_label = Unknown in ", sample_name, ".")
}
if (!any(labels != "Unknown")) {
  stop("All cells have consensus_label = Unknown in ", sample_name, ".")
}

qc_data <- data.frame(
  cell_id = cell_ids,
  seurat_clusters = as.character(metadata$seurat_clusters),
  consensus_label = labels,
  qc_group = ifelse(labels == "Unknown", "Unknown", "Known"),
  stringsAsFactors = FALSE
)
for (metric in qc_metrics) {
  values <- suppressWarnings(as.numeric(metadata[[metric]]))
  if (any(!is.finite(values))) {
    stop("Non-finite values in ", metric, " for ", sample_name, ".")
  }
  qc_data[[metric]] <- values
}

summarize_values <- function(values) {
  c(
    n = length(values), mean = mean(values), sd = stats::sd(values),
    q10 = unname(stats::quantile(values, 0.10)),
    q25 = unname(stats::quantile(values, 0.25)),
    median = stats::median(values),
    q75 = unname(stats::quantile(values, 0.75)),
    q90 = unname(stats::quantile(values, 0.90))
  )
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
tail_direction <- c(
  nCount_Xenium = "low", nFeature_Xenium = "low",
  cell_area = "low", percent_control = "high"
)
known_tail_thresholds <- vapply(qc_metrics, function(metric) {
  probability <- if (tail_direction[[metric]] == "low") 0.10 else 0.90
  unname(stats::quantile(known_data[[metric]], probability))
}, numeric(1))

effect_rows <- lapply(qc_metrics, function(metric) {
  unknown_values <- unknown_data[[metric]]
  known_values <- known_data[[metric]]
  test <- suppressWarnings(stats::wilcox.test(
    unknown_values, known_values, exact = FALSE
  ))
  u_statistic <- unname(test$statistic)
  rank_biserial <- 2 * u_statistic /
    (length(unknown_values) * length(known_values)) - 1
  threshold <- known_tail_thresholds[[metric]]
  if (tail_direction[[metric]] == "low") {
    unknown_tail <- unknown_values <= threshold
    known_tail <- known_values <= threshold
  } else {
    unknown_tail <- unknown_values >= threshold
    known_tail <- known_values >= threshold
  }
  data.frame(
    sample = sample_name,
    metric = metric,
    median_known = median(known_values),
    median_unknown = median(unknown_values),
    median_difference_unknown_minus_known =
      median(unknown_values) - median(known_values),
    rank_biserial_unknown_vs_known = rank_biserial,
    wilcoxon_p = test$p.value,
    known_tail_direction = tail_direction[[metric]],
    known_tail_threshold = threshold,
    known_tail_fraction = mean(known_tail),
    unknown_tail_fraction = mean(unknown_tail),
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
  if (tail_direction[[metric]] == "low") {
    qc_data[[flag_name]] <- qc_data[[metric]] <= threshold
  } else {
    qc_data[[flag_name]] <- qc_data[[metric]] >= threshold
  }
}

cluster_rows <- lapply(sort(unique(unknown_data$seurat_clusters)), function(cluster_id) {
  cluster_data <- qc_data[
    qc_data$qc_group == "Unknown" & qc_data$seurat_clusters == cluster_id,
    , drop = FALSE
  ]
  cluster_labels <- unique(cluster_data$consensus_label)
  if (length(cluster_labels) != 1L || cluster_labels != "Unknown") {
    stop("Unexpected cluster-to-consensus cardinality for cluster ", cluster_id, ".")
  }
  data.frame(
    sample = sample_name,
    seurat_clusters = cluster_id,
    n_cells = nrow(cluster_data),
    fraction_of_sample = nrow(cluster_data) / nrow(qc_data),
    median_nCount_Xenium = median(cluster_data$nCount_Xenium),
    median_nFeature_Xenium = median(cluster_data$nFeature_Xenium),
    median_cell_area = median(cluster_data$cell_area),
    median_percent_control = median(cluster_data$percent_control),
    fraction_low_nCount_vs_known_q10 = mean(
      cluster_data$nCount_Xenium_known_tail
    ),
    fraction_low_nFeature_vs_known_q10 = mean(
      cluster_data$nFeature_Xenium_known_tail
    ),
    fraction_low_cell_area_vs_known_q10 = mean(
      cluster_data$cell_area_known_tail
    ),
    fraction_high_percent_control_vs_known_q90 = mean(
      cluster_data$percent_control_known_tail
    ),
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
  stringsAsFactors = FALSE
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
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
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
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
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

message("Unknown cells: ", nrow(unknown_data), " / ", nrow(qc_data))
message("Unknown clusters: ", paste(unknown_clusters, collapse = ", "))
message("Wrote:\n- ", paste(output_paths, collapse = "\n- "))
message(
  "Interpretation note: cell-level P values are descriptive because cells within ",
  "a sample/cluster are not independent. Prioritize effect sizes and consistent ",
  "cluster-level QC shifts."
)
