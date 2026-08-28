#!/usr/bin/env Rscript

# Compare weighted-2-of-3 consensus labels and cluster membership across the
# resolution-3, resolution-4, and resolution-5 three-sample pilot objects.
# This is an audit only: it never modifies or replaces consensus labels.

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
sample_manifest <- load_sample_manifest(config)

required_pilot_columns <- c("task_id", "sample_id", "PCW", "age_group")
missing_pilot_columns <- setdiff(required_pilot_columns, names(pilot_manifest))
if (length(missing_pilot_columns)) {
  stop(
    "Pilot manifest is missing required columns: ",
    paste(missing_pilot_columns, collapse = ", ")
  )
}
manifest_counts <- vapply(
  pilot_manifest$sample_id,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(manifest_counts != 1L)) {
  stop("Every pilot sample must map exactly once to config/samples.csv.")
}

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite", "--combine")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
combine_mode <- "--combine" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

task_map <- data.frame(
  task_id = seq_len(nrow(pilot_manifest)),
  sample_id = as.character(pilot_manifest$sample_id),
  PCW = as.character(pilot_manifest$PCW),
  age_group = as.character(pilot_manifest$age_group),
  stringsAsFactors = FALSE
)
if (anyDuplicated(task_map$sample_id) || any(!nzchar(task_map$sample_id))) {
  stop("Pilot manifest contains duplicate or blank sample IDs.")
}

if (list_requested) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args) || combine_mode || dry_run || overwrite) {
    stop("--list must be used alone.")
  }
  cat("Consensus method: weighted_2of3\n")
  cat("Compared resolutions: 3.0, 4.0, 5.0\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

positional_args <- args[!args %in% valid_options]
if (combine_mode && length(positional_args)) {
  stop("--combine does not accept TASK_ID.")
}
if (!combine_mode && length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R ",
    "[--list] [--dry-run|--overwrite] [TASK_ID] or --combine [--dry-run|--overwrite]"
  )
}

output_root <- here(config$project$outputs_dir)
resolution_info <- data.frame(
  resolution = c("3.0", "4.0", "5.0"),
  stage = c(
    "03c_resolution3_pilot", "03d_resolution4_pilot", "03e_resolution5_pilot"
  ),
  stringsAsFactors = FALSE
)
audit_root <- file.path(
  output_root, "xenium", "preprocess", "03f_resolution_consistency",
  "weighted_2of3"
)
table_dir <- file.path(audit_root, "tables")
plot_dir <- file.path(audit_root, "plots")

input_paths_for_sample <- function(sample_id) {
  setNames(
    file.path(
      output_root, "xenium", "preprocess", resolution_info$stage,
      "annotation", "03_consensus_labels_weighted_2of3", "rds",
      paste0(sample_id, "_Consensus_annotated.rds")
    ),
    resolution_info$resolution
  )
}

output_paths_for_sample <- function(sample_id) {
  stem <- paste0(sample_id, "_Res3_Res4_Res5_weighted2of3")
  c(
    cell_labels = file.path(table_dir, paste0(stem, "_cell_labels.csv")),
    agreement = file.path(table_dir, paste0(stem, "_agreement_summary.csv")),
    pairwise_metrics = file.path(table_dir, paste0(stem, "_pairwise_metrics.csv")),
    label_transitions = file.path(table_dir, paste0(stem, "_label_transitions.csv")),
    label_stability = file.path(table_dir, paste0(stem, "_label_stability.csv")),
    cluster_transitions = file.path(table_dir, paste0(stem, "_cluster_transitions.csv")),
    provenance = file.path(table_dir, paste0(stem, "_provenance.csv")),
    report = file.path(plot_dir, paste0(stem, "_consistency_report.pdf")),
    umap = file.path(plot_dir, paste0(stem, "_agreement_UMAP.tif")),
    spatial = file.path(plot_dir, paste0(stem, "_agreement_Spatial.tif"))
  )
}

combined_paths <- c(
  agreement = file.path(table_dir, "All_samples_Res3_Res4_Res5_agreement_summary.csv"),
  pairwise_metrics = file.path(
    table_dir, "All_samples_Res3_Res4_Res5_pairwise_metrics.csv"
  ),
  label_stability = file.path(
    table_dir, "All_samples_Res3_Res4_Res5_label_stability.csv"
  ),
  report = file.path(plot_dir, "All_samples_Res3_Res4_Res5_consistency_report.pdf")
)

refuse_overwrite <- function(paths, description) {
  existing <- paths[file.exists(paths)]
  if (length(existing) && !overwrite) {
    stop(
      "Refusing to overwrite existing ", description, ":\n- ",
      paste(existing, collapse = "\n- "),
      "\nRerun with --overwrite only after reviewing these audit files."
    )
  }
}

read_csv_checked <- function(path, required_columns) {
  if (!file.exists(path)) stop("Missing required audit table: ", path)
  value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_columns <- setdiff(required_columns, names(value))
  if (length(missing_columns)) {
    stop(
      "Table is missing required columns: ", path, "; ",
      paste(missing_columns, collapse = ", ")
    )
  }
  value
}

if (combine_mode) {
  per_sample_paths <- lapply(task_map$sample_id, output_paths_for_sample)
  required_inputs <- c(
    vapply(per_sample_paths, `[[`, character(1), "agreement"),
    vapply(per_sample_paths, `[[`, character(1), "pairwise_metrics"),
    vapply(per_sample_paths, `[[`, character(1), "label_stability")
  )
  if (dry_run) {
    cat("Mode: combine completed per-sample resolution-consistency audits\n")
    cat("Dry-run scope: path inspection only; no CSV or RDS file is loaded.\n")
    write.table(
      data.frame(input = required_inputs, exists = file.exists(required_inputs)),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    write.table(
      data.frame(output = combined_paths, exists = file.exists(combined_paths)),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    quit(
      save = "no",
      status = if (all(file.exists(required_inputs))) 0L else 1L
    )
  }
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Cannot combine; missing per-sample tables:\n- ", paste(missing_inputs, collapse = "\n- "))
  }
  refuse_overwrite(combined_paths, "combined resolution-consistency outputs")

  agreement_all <- do.call(rbind, lapply(per_sample_paths, function(paths) {
    read_csv_checked(paths[["agreement"]], c("sample", "agreement_status", "n_cells"))
  }))
  metrics_all <- do.call(rbind, lapply(per_sample_paths, function(paths) {
    read_csv_checked(
      paths[["pairwise_metrics"]],
      c("sample", "resolution_from", "resolution_to", "label_exact_agreement", "cluster_ARI")
    )
  }))
  label_stability_all <- do.call(rbind, lapply(per_sample_paths, function(paths) {
    read_csv_checked(
      paths[["label_stability"]],
      c("sample", "resolution_from", "resolution_to", "label", "jaccard")
    )
  }))

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(agreement_all, combined_paths[["agreement"]], row.names = FALSE)
  utils::write.csv(metrics_all, combined_paths[["pairwise_metrics"]], row.names = FALSE)
  utils::write.csv(
    label_stability_all, combined_paths[["label_stability"]], row.names = FALSE
  )

  p_agreement <- ggplot(
    agreement_all,
    aes(x = sample, y = fraction_cells, fill = agreement_status)
  ) +
    geom_col() +
    scale_fill_manual(values = c(
      "All 3 agree" = "#009E73", "2 of 3 agree" = "#E69F00",
      "All different" = "#D55E00"
    )) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(
      title = "Consensus-label agreement across resolutions 3, 4, and 5",
      x = NULL, y = "Fraction of cells", fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  metric_long <- rbind(
    data.frame(
      sample = metrics_all$sample,
      comparison = paste(metrics_all$resolution_from, metrics_all$resolution_to, sep = " vs "),
      metric = "Exact label agreement",
      value = metrics_all$label_exact_agreement
    ),
    data.frame(
      sample = metrics_all$sample,
      comparison = paste(metrics_all$resolution_from, metrics_all$resolution_to, sep = " vs "),
      metric = "Cluster ARI",
      value = metrics_all$cluster_ARI
    )
  )
  p_metrics <- ggplot(
    metric_long,
    aes(x = comparison, y = value, color = sample, group = sample)
  ) +
    geom_point(size = 2) +
    geom_line() +
    facet_wrap(~metric) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "Pairwise agreement by sample",
      x = "Resolution comparison", y = "Agreement", color = "Sample"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p_jaccard <- ggplot(
    label_stability_all,
    aes(
      x = paste(resolution_from, resolution_to, sep = " vs "),
      y = jaccard, color = sample
    )
  ) +
    geom_boxplot(aes(group = interaction(sample, resolution_from, resolution_to)),
                 outlier.shape = NA, alpha = 0.15) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "Per-label Jaccard stability",
      x = "Resolution comparison", y = "Jaccard", color = "Sample"
    ) +
    theme_bw(base_size = 10)

  grDevices::cairo_pdf(combined_paths[["report"]], width = 10, height = 7)
  print(p_agreement)
  print(p_metrics)
  print(p_jaccard)
  grDevices::dev.off()
  message("Wrote combined resolution-consistency audit:\n- ", paste(combined_paths, collapse = "\n- "))
  quit(save = "no", status = 0L)
}

task_value <- if (length(positional_args) == 1L) {
  positional_args[[1]]
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

task <- task_map[task_id, , drop = FALSE]
sample_id <- task$sample_id[[1]]
input_paths <- input_paths_for_sample(sample_id)
output_paths <- output_paths_for_sample(sample_id)

if (dry_run) {
  cat("Mode: weighted-2-of-3 consensus consistency across resolutions 3, 4, and 5\n")
  cat("Task:", task_id, "of", nrow(task_map), "\n")
  cat("Sample:", sample_id, "\n")
  cat("Dry-run scope: path inspection only; no RDS object is loaded.\n")
  write.table(
    data.frame(
      resolution = names(input_paths), input = unname(input_paths),
      exists = file.exists(input_paths)
    ),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  write.table(
    data.frame(output = output_paths, exists = file.exists(output_paths)),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  quit(save = "no", status = if (all(file.exists(input_paths))) 0L else 1L)
}

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing weighted-2-of-3 consensus objects:\n- ", paste(missing_inputs, collapse = "\n- "))
}
refuse_overwrite(output_paths, paste0("resolution-consistency outputs for ", sample_id))

extract_metadata <- function(obj, resolution) {
  cell_ids <- Cells(obj)
  metadata <- obj[[]]
  if (anyDuplicated(cell_ids)) {
    stop("Duplicate cell IDs at resolution ", resolution, " for ", sample_id, ".")
  }
  if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Missing or duplicate metadata row names at resolution ", resolution, ".")
  }
  if (!setequal(cell_ids, rownames(metadata))) {
    stop("Cell IDs and metadata row names differ at resolution ", resolution, ".")
  }
  metadata <- metadata[cell_ids, , drop = FALSE]
  required_columns <- c("seurat_clusters", "consensus_label", "consensus_method")
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns)) {
    stop(
      "Resolution-", resolution, " object is missing metadata: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  methods <- unique(as.character(metadata$consensus_method))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  if (!identical(methods, "weighted_2of3")) {
    stop(
      "Resolution-", resolution,
      " object does not contain only consensus_method = weighted_2of3."
    )
  }
  labels <- trimws(as.character(metadata$consensus_label))
  clusters <- trimws(as.character(metadata$seurat_clusters))
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Blank consensus labels at resolution ", resolution, ".")
  }
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    stop("Blank cluster IDs at resolution ", resolution, ".")
  }
  data.frame(
    cell_id = cell_ids, cluster = clusters, label = labels,
    stringsAsFactors = FALSE
  )
}

message("Loading three consensus objects sequentially for ", sample_id, ".")
metadata_by_resolution <- vector("list", length(input_paths))
names(metadata_by_resolution) <- names(input_paths)
reference_umap <- NULL
reference_coords <- NULL
for (resolution in names(input_paths)) {
  current_obj <- readRDS(input_paths[[resolution]])
  metadata_by_resolution[[resolution]] <- extract_metadata(current_obj, resolution)
  if (resolution == "3.0") {
    reference_umap <- Embeddings(current_obj, reduction = "umap")
    reference_coords <- Seurat::GetTissueCoordinates(
      current_obj, image = "fov", full = FALSE
    )
  }
  rm(current_obj)
  invisible(gc())
}
if (is.null(reference_umap) || is.null(reference_coords)) {
  stop("Failed to extract resolution-3 UMAP or spatial coordinates.")
}
reference_cells <- metadata_by_resolution[[1]]$cell_id
for (resolution in names(metadata_by_resolution)[-1]) {
  comparison_cells <- metadata_by_resolution[[resolution]]$cell_id
  if (!setequal(reference_cells, comparison_cells)) {
    stop(
      "Cell sets differ between resolution 3.0 and resolution ", resolution,
      " for ", sample_id, "."
    )
  }
  metadata_by_resolution[[resolution]] <- metadata_by_resolution[[resolution]][
    match(reference_cells, comparison_cells), , drop = FALSE
  ]
  if (!identical(reference_cells, metadata_by_resolution[[resolution]]$cell_id)) {
    stop("Failed to align cell IDs for resolution ", resolution, ".")
  }
}

cell_table <- data.frame(cell_id = reference_cells, stringsAsFactors = FALSE)
for (resolution in names(metadata_by_resolution)) {
  suffix <- sub("\\.0$", "", resolution)
  cell_table[[paste0("cluster_res", suffix)]] <- metadata_by_resolution[[resolution]]$cluster
  cell_table[[paste0("label_res", suffix)]] <- metadata_by_resolution[[resolution]]$label
}

label_matrix <- as.matrix(cell_table[c("label_res3", "label_res4", "label_res5")])
modal_results <- apply(label_matrix, 1L, function(values) {
  counts <- sort(table(values), decreasing = TRUE)
  max_support <- as.integer(counts[[1]])
  modal_label <- if (max_support >= 2L) names(counts)[[1]] else "Resolution_unstable"
  status <- if (max_support == 3L) {
    "All 3 agree"
  } else if (max_support == 2L) {
    "2 of 3 agree"
  } else {
    "All different"
  }
  c(modal_label = modal_label, support_n = max_support, status = status)
})
cell_table$cross_resolution_modal_label <- modal_results["modal_label", ]
cell_table$cross_resolution_support_n <- as.integer(modal_results["support_n", ])
cell_table$agreement_status <- modal_results["status", ]

agreement_levels <- c("All 3 agree", "2 of 3 agree", "All different")
agreement_counts <- table(factor(cell_table$agreement_status, levels = agreement_levels))
agreement_summary <- data.frame(
  sample = sample_id,
  agreement_status = agreement_levels,
  n_cells = as.integer(agreement_counts),
  fraction_cells = as.integer(agreement_counts) / nrow(cell_table),
  stringsAsFactors = FALSE
)

choose2 <- function(x) x * (x - 1) / 2
adjusted_rand_index <- function(x, y) {
  contingency <- table(x, y)
  n <- sum(contingency)
  if (n < 2L) return(NA_real_)
  index <- sum(choose2(contingency))
  row_index <- sum(choose2(rowSums(contingency)))
  column_index <- sum(choose2(colSums(contingency)))
  expected <- row_index * column_index / choose2(n)
  maximum <- (row_index + column_index) / 2
  denominator <- maximum - expected
  if (denominator == 0) return(if (index == maximum) 1 else 0)
  (index - expected) / denominator
}

normalized_mutual_information <- function(x, y) {
  contingency <- table(x, y)
  probabilities <- contingency / sum(contingency)
  row_probabilities <- rowSums(probabilities)
  column_probabilities <- colSums(probabilities)
  nonzero <- probabilities > 0
  expected <- outer(row_probabilities, column_probabilities)
  mutual_information <- sum(
    probabilities[nonzero] * log(probabilities[nonzero] / expected[nonzero])
  )
  entropy_x <- -sum(row_probabilities[row_probabilities > 0] *
                      log(row_probabilities[row_probabilities > 0]))
  entropy_y <- -sum(column_probabilities[column_probabilities > 0] *
                      log(column_probabilities[column_probabilities > 0]))
  denominator <- sqrt(entropy_x * entropy_y)
  if (denominator == 0) return(if (identical(as.character(x), as.character(y))) 1 else 0)
  mutual_information / denominator
}

cohen_kappa <- function(x, y) {
  levels_union <- union(unique(x), unique(y))
  x_factor <- factor(x, levels = levels_union)
  y_factor <- factor(y, levels = levels_union)
  observed <- mean(x_factor == y_factor)
  expected <- sum(
    prop.table(table(x_factor)) * prop.table(table(y_factor))
  )
  if (expected == 1) return(if (observed == 1) 1 else NA_real_)
  (observed - expected) / (1 - expected)
}

resolution_pairs <- list(c("3.0", "4.0"), c("3.0", "5.0"), c("4.0", "5.0"))
pairwise_metrics_rows <- list()
label_transition_rows <- list()
label_stability_rows <- list()
cluster_transition_rows <- list()

for (pair_index in seq_along(resolution_pairs)) {
  resolution_from <- resolution_pairs[[pair_index]][[1]]
  resolution_to <- resolution_pairs[[pair_index]][[2]]
  from_data <- metadata_by_resolution[[resolution_from]]
  to_data <- metadata_by_resolution[[resolution_to]]

  pairwise_metrics_rows[[pair_index]] <- data.frame(
    sample = sample_id,
    resolution_from = resolution_from,
    resolution_to = resolution_to,
    n_cells = nrow(from_data),
    label_exact_agreement = mean(from_data$label == to_data$label),
    label_cohen_kappa = cohen_kappa(from_data$label, to_data$label),
    cluster_ARI = adjusted_rand_index(from_data$cluster, to_data$cluster),
    cluster_NMI = normalized_mutual_information(from_data$cluster, to_data$cluster),
    stringsAsFactors = FALSE
  )

  label_transition <- as.data.frame(
    table(from_label = from_data$label, to_label = to_data$label),
    stringsAsFactors = FALSE
  )
  label_transition <- label_transition[label_transition$Freq > 0L, , drop = FALSE]
  label_transition$n_cells <- label_transition$Freq
  label_transition$Freq <- NULL
  label_transition$fraction_from <- label_transition$n_cells /
    ave(label_transition$n_cells, label_transition$from_label, FUN = sum)
  label_transition$fraction_to <- label_transition$n_cells /
    ave(label_transition$n_cells, label_transition$to_label, FUN = sum)
  label_transition$fraction_sample <- label_transition$n_cells / nrow(from_data)
  label_transition$sample <- sample_id
  label_transition$resolution_from <- resolution_from
  label_transition$resolution_to <- resolution_to
  label_transition_rows[[pair_index]] <- label_transition[c(
    "sample", "resolution_from", "resolution_to", "from_label", "to_label",
    "n_cells", "fraction_from", "fraction_to", "fraction_sample"
  )]

  labels_union <- sort(union(unique(from_data$label), unique(to_data$label)))
  label_stability_rows[[pair_index]] <- do.call(rbind, lapply(labels_union, function(label) {
    in_from <- from_data$label == label
    in_to <- to_data$label == label
    intersection_n <- sum(in_from & in_to)
    union_n <- sum(in_from | in_to)
    data.frame(
      sample = sample_id,
      resolution_from = resolution_from,
      resolution_to = resolution_to,
      label = label,
      n_from = sum(in_from),
      n_to = sum(in_to),
      n_intersection = intersection_n,
      jaccard = if (union_n > 0) intersection_n / union_n else NA_real_,
      retention_from = if (sum(in_from) > 0) intersection_n / sum(in_from) else NA_real_,
      retention_to = if (sum(in_to) > 0) intersection_n / sum(in_to) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  cluster_transition <- as.data.frame(
    table(from_cluster = from_data$cluster, to_cluster = to_data$cluster),
    stringsAsFactors = FALSE
  )
  cluster_transition <- cluster_transition[cluster_transition$Freq > 0L, , drop = FALSE]
  cluster_transition$n_cells <- cluster_transition$Freq
  cluster_transition$Freq <- NULL
  cluster_transition$fraction_from <- cluster_transition$n_cells /
    ave(cluster_transition$n_cells, cluster_transition$from_cluster, FUN = sum)
  cluster_transition$fraction_to <- cluster_transition$n_cells /
    ave(cluster_transition$n_cells, cluster_transition$to_cluster, FUN = sum)
  cluster_transition$sample <- sample_id
  cluster_transition$resolution_from <- resolution_from
  cluster_transition$resolution_to <- resolution_to
  cluster_transition_rows[[pair_index]] <- cluster_transition[c(
    "sample", "resolution_from", "resolution_to", "from_cluster", "to_cluster",
    "n_cells", "fraction_from", "fraction_to"
  )]
}

pairwise_metrics <- do.call(rbind, pairwise_metrics_rows)
label_transitions <- do.call(rbind, label_transition_rows)
label_stability <- do.call(rbind, label_stability_rows)
cluster_transitions <- do.call(rbind, cluster_transition_rows)

provenance <- data.frame(
  sample = sample_id,
  resolution = names(input_paths),
  input_path = unname(input_paths),
  consensus_method = "weighted_2of3",
  cell_count = nrow(cell_table),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = "Audit only; cross_resolution_modal_label does not replace consensus_label",
  stringsAsFactors = FALSE
)

agreement_colors <- c(
  "All 3 agree" = "#009E73", "2 of 3 agree" = "#E69F00",
  "All different" = "#D55E00"
)
cell_table$agreement_status <- factor(
  cell_table$agreement_status, levels = agreement_levels
)

p_agreement <- ggplot(
  agreement_summary,
  aes(x = agreement_status, y = fraction_cells, fill = agreement_status)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = agreement_colors) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = paste(sample_id, "consensus-label agreement"),
    subtitle = "Weighted 2-of-3 labels across resolutions 3, 4, and 5",
    x = NULL, y = "Fraction of cells"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

metric_long <- rbind(
  data.frame(
    comparison = paste(pairwise_metrics$resolution_from,
                       pairwise_metrics$resolution_to, sep = " vs "),
    metric = "Exact label agreement", value = pairwise_metrics$label_exact_agreement
  ),
  data.frame(
    comparison = paste(pairwise_metrics$resolution_from,
                       pairwise_metrics$resolution_to, sep = " vs "),
    metric = "Label Cohen kappa", value = pairwise_metrics$label_cohen_kappa
  ),
  data.frame(
    comparison = paste(pairwise_metrics$resolution_from,
                       pairwise_metrics$resolution_to, sep = " vs "),
    metric = "Cluster ARI", value = pairwise_metrics$cluster_ARI
  ),
  data.frame(
    comparison = paste(pairwise_metrics$resolution_from,
                       pairwise_metrics$resolution_to, sep = " vs "),
    metric = "Cluster NMI", value = pairwise_metrics$cluster_NMI
  )
)
p_metrics <- ggplot(metric_long, aes(x = comparison, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste(sample_id, "pairwise resolution agreement"),
    x = "Resolution comparison", y = "Agreement", fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_label_stability <- ggplot(
  label_stability,
  aes(
    x = paste(resolution_from, resolution_to, sep = " vs "),
    y = jaccard, color = label
  )
) +
  geom_point(size = 2, alpha = 0.8) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste(sample_id, "per-label Jaccard stability"),
    x = "Resolution comparison", y = "Jaccard", color = "Consensus label"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_label_transition <- ggplot(
  label_transitions,
  aes(x = to_label, y = from_label, fill = fraction_from)
) +
  geom_tile() +
  facet_wrap(~paste(resolution_from, resolution_to, sep = " to "), scales = "free") +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    title = paste(sample_id, "consensus-label transitions"),
    x = "Label at higher resolution", y = "Label at lower resolution",
    fill = "Fraction of\nsource label"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))

umap <- reference_umap
if (is.null(rownames(umap)) || anyDuplicated(rownames(umap))) {
  stop("Resolution-3 UMAP lacks unique cell IDs for ", sample_id, ".")
}
umap_indices <- match(reference_cells, rownames(umap))
if (anyNA(umap_indices)) stop("Resolution-3 UMAP is missing object cells.")
umap_data <- data.frame(
  UMAP_1 = umap[umap_indices, 1], UMAP_2 = umap[umap_indices, 2],
  agreement_status = cell_table$agreement_status
)
p_umap <- ggplot(
  umap_data, aes(x = UMAP_1, y = UMAP_2, color = agreement_status)
) +
  geom_point(size = 0.08, alpha = 0.8) +
  scale_color_manual(values = agreement_colors, drop = FALSE) +
  coord_equal() +
  labs(
    title = paste(sample_id, "cross-resolution label agreement"),
    color = NULL
  ) +
  theme_void() +
  theme(legend.position = "bottom")

coords <- reference_coords
required_coordinate_columns <- c("x", "y", "cell")
missing_coordinate_columns <- setdiff(required_coordinate_columns, names(coords))
if (length(missing_coordinate_columns)) {
  stop(
    "Spatial coordinates lack required columns: ",
    paste(missing_coordinate_columns, collapse = ", ")
  )
}
coordinate_cell_ids <- as.character(coords$cell)
if (anyDuplicated(coordinate_cell_ids)) {
  stop("Spatial coordinates contain duplicate cell IDs.")
}
coordinate_indices <- match(coordinate_cell_ids, reference_cells)
if (anyNA(coordinate_indices)) {
  stop("Spatial coordinates contain cell IDs absent from the resolution-3 object.")
}
spatial_data <- data.frame(
  plot_x = coords$y, plot_y = coords$x,
  agreement_status = cell_table$agreement_status[coordinate_indices]
)
p_spatial <- ggplot(
  spatial_data, aes(x = plot_x, y = plot_y, color = agreement_status)
) +
  geom_point(size = 0.03) +
  scale_color_manual(values = agreement_colors, drop = FALSE) +
  coord_fixed() +
  labs(
    title = paste(sample_id, "spatial cross-resolution agreement"),
    color = NULL
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    plot.title = element_text(color = "white"),
    legend.position = "bottom",
    legend.text = element_text(color = "white")
  )

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(cell_table, output_paths[["cell_labels"]], row.names = FALSE)
utils::write.csv(agreement_summary, output_paths[["agreement"]], row.names = FALSE)
utils::write.csv(pairwise_metrics, output_paths[["pairwise_metrics"]], row.names = FALSE)
utils::write.csv(label_transitions, output_paths[["label_transitions"]], row.names = FALSE)
utils::write.csv(label_stability, output_paths[["label_stability"]], row.names = FALSE)
utils::write.csv(cluster_transitions, output_paths[["cluster_transitions"]], row.names = FALSE)
utils::write.csv(provenance, output_paths[["provenance"]], row.names = FALSE)

grDevices::cairo_pdf(output_paths[["report"]], width = 11, height = 8.5)
print(p_agreement)
print(p_metrics)
print(p_label_stability)
print(p_label_transition)
grDevices::dev.off()
Cairo::CairoTIFF(output_paths[["umap"]], width = 8, height = 6, units = "in", res = 600)
print(p_umap)
grDevices::dev.off()
Cairo::CairoTIFF(
  output_paths[["spatial"]], width = 8, height = 8, units = "in", res = 600,
  bg = "black"
)
print(p_spatial)
grDevices::dev.off()

message("Wrote resolution-consistency audit for ", sample_id, ":\n- ",
        paste(output_paths, collapse = "\n- "))
message(
  "The cross-resolution modal label is descriptive only and was not written ",
  "back to any Seurat object."
)
