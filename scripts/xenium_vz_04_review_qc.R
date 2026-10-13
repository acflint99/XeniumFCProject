#!/usr/bin/env Rscript

# Stage 4a reviews merged VZ clusters. Stage 4b records the decision and only
# materializes filtered whole-tissue copies when cells are actually removed.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_vz_rl_unknown_qc.R"))
source(here("scripts", "R", "xenium_vz_rl_qc_routing.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)
cluster_column <- "Xenium_snn_res.0.8"
args <- commandArgs(trailingOnly = TRUE)

remove_args <- grep("^--remove-clusters=", args, value = TRUE)
if (length(remove_args) > 1L) stop("Supply --remove-clusters only once.")
remove_supplied <- length(remove_args) == 1L
remove_value <- if (remove_supplied) sub("^--remove-clusters=", "", remove_args) else ""
exact_options <- c("--qc-only", "--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% exact_options &
                          !grepl("^--remove-clusters=", args)]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!startsWith(args, "--"))) {
  stop("All arguments must be named options. See the usage message in this script.")
}
qc_only <- "--qc-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (qc_only && remove_supplied) stop("Choose --qc-only or --remove-clusters, not both.")

clusters_to_remove <- character()
if (remove_supplied && !identical(tolower(trimws(remove_value)), "none")) {
  clusters_to_remove <- unique(trimws(strsplit(remove_value, ",", fixed = TRUE)[[1]]))
  if (any(!nzchar(clusters_to_remove))) stop("--remove-clusters contains a blank cluster ID.")
}

output_root <- here(config$project$outputs_dir)
merged_path <- file.path(output_root, "xenium", "vz", "03_integrated", "rds", "Xenium_VZ_integrated.rds")
input_dir <- production_paths$consensus_rds_dir
output_dir <- file.path(output_root, "xenium", "vz", "04_qc", "rds")
plot_dir <- file.path(output_root, "xenium", "vz", "04_qc", "plots")
table_dir <- file.path(output_root, "xenium", "vz", "04_qc", "tables")

input_names <- paste0(sample_ids, "_Consensus_annotated.rds")
input_paths <- file.path(input_dir, input_names)
output_names <- paste0(sample_ids, "_VZ_QC.rds")
output_paths <- file.path(output_dir, output_names)
summary_path <- file.path(table_dir, "Xenium_VZ_RawSubcluster_QC_Summary.csv")
cluster_sample_summary_path <- file.path(
  table_dir, "Xenium_VZ_RawSubcluster_BySample_QC_Summary.csv"
)
marker_path <- file.path(table_dir, "Xenium_VZ_RawSubcluster_top10_Markers_Res0.8.csv")
review_path <- file.path(table_dir, "Xenium_VZ_QC_review_manifest.csv")
pdf_path <- file.path(plot_dir, "Xenium_VZ_RawSubcluster_QC_Violins.pdf")
cluster_sample_qc_heatmap_path <- file.path(
  plot_dir, "Xenium_VZ_RawSubcluster_BySample_QC_Heatmaps.pdf"
)
cluster_sample_composition_heatmap_path <- file.path(
  plot_dir, "Xenium_VZ_RawSubcluster_SampleComposition_Heatmap.pdf"
)
removal_path <- file.path(table_dir, "Xenium_VZ_QC_removal_manifest.csv")
unknown_qc_paths <- xenium_vz_rl_unknown_qc_output_paths("VZ", plot_dir, table_dir)
qc_outputs <- c(
  summary_path, cluster_sample_summary_path, marker_path,
  unname(unknown_qc_paths), review_path, pdf_path,
  cluster_sample_qc_heatmap_path, cluster_sample_composition_heatmap_path
)
decision_outputs <- if (remove_supplied && !length(clusters_to_remove)) {
  removal_path
} else {
  c(output_paths, removal_path)
}

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  dry_run_inputs <- if (remove_supplied) {
    c(merged_path, input_paths, qc_outputs)
  } else {
    c(merged_path, input_paths)
  }
  dry_run_outputs <- if (qc_only) qc_outputs else decision_outputs
  compact_dry_run(
    paste0(
      "VZ QC [4",
      if (qc_only) "a review" else if (remove_supplied) "b decision" else " inspection",
      "]"
    ),
    inputs = dry_run_inputs,
    outputs = dry_run_outputs,
    checks = c(no_unexpected_inputs = !length(unexpected_names))
  )
  if (length(unexpected_names)) cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  quit(save = "no", status = 0L)
}

if (!qc_only && !remove_supplied) {
  stop(
    "Choose a QC stage explicitly:\n",
    "- --qc-only to create evidence for review\n",
    "- --remove-clusters=9,25,26 for the reviewed VZ QC decision\n",
    "- --remove-clusters=none after review when no cells should be removed"
  )
}
if (!file.exists(merged_path)) stop("Merged VZ object not found: ", merged_path)

if (qc_only) {
  existing_qc <- qc_outputs[file.exists(qc_outputs)]
  if (length(existing_qc) && !overwrite) {
    stop("Refusing to overwrite existing VZ QC review outputs:\n- ",
         paste(existing_qc, collapse = "\n- "), "\nUse --overwrite only after review.")
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(merged_path)
  if (!cluster_column %in% colnames(obj[[]])) stop("Merged VZ object lacks ", cluster_column, ".")
  if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
    stop("Merged VZ object lacks valid sample-qualified Unknown labels. Rerun stages 2 and 3.")
  }
  if (ncol(obj) == 0L) stop("Merged VZ object contains zero cells.")
  obj <- JoinLayers(obj)
  Idents(obj) <- cluster_column
  set.seed(config$runtime$random_seed)

  cluster_labels <- as.character(obj[[cluster_column, drop = TRUE]])
  sample_labels <- as.character(obj$orig.ident)
  consensus_labels <- as.character(obj$consensus_label)
  if (anyNA(cluster_labels) || any(!nzchar(cluster_labels))) {
    stop("Merged VZ object contains blank ", cluster_column, " values.")
  }
  if (anyNA(sample_labels) || any(!sample_labels %in% sample_ids)) {
    stop("Merged VZ object contains sample IDs absent from config/samples.csv.")
  }
  missing_samples <- setdiff(sample_ids, unique(sample_labels))
  if (length(missing_samples)) {
    stop(
      "Merged VZ object is missing configured samples: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  unknown_cells <- is_unknown_consensus_label(consensus_labels)
  configured_labels <- unlist(
    config$regional_subsets$vz_broad_labels, use.names = FALSE
  )
  unexpected_configured_labels <- setdiff(
    unique(consensus_labels[!unknown_cells]), configured_labels
  )
  if (length(unexpected_configured_labels)) {
    stop(
      "Merged VZ object contains unexpected configured consensus labels: ",
      paste(unexpected_configured_labels, collapse = ", ")
    )
  }

  cluster_levels <- unique(cluster_labels)
  numeric_cluster_levels <- suppressWarnings(as.numeric(cluster_levels))
  if (all(is.finite(numeric_cluster_levels))) {
    cluster_levels <- cluster_levels[order(numeric_cluster_levels)]
  } else {
    cluster_levels <- sort(cluster_levels)
  }

  cluster_sample_metadata <- obj@meta.data %>%
    transmute(
      sample_id = as.character(.data$orig.ident),
      expression_cluster = as.character(.data[[cluster_column]]),
      consensus_label = as.character(.data$consensus_label),
      is_unknown = is_unknown_consensus_label(.data$consensus_label),
      nCount_Xenium = .data$nCount_Xenium,
      nFeature_Xenium = .data$nFeature_Xenium
    )
  if (any(!is.finite(cluster_sample_metadata$nCount_Xenium)) ||
      any(!is.finite(cluster_sample_metadata$nFeature_Xenium))) {
    stop("Merged VZ object contains non-finite QC values.")
  }

  sample_qc_reference <- cluster_sample_metadata %>%
    group_by(.data$sample_id) %>%
    summarise(
      sample_cells = n(),
      sample_median_counts = median(.data$nCount_Xenium),
      sample_median_features = median(.data$nFeature_Xenium),
      .groups = "drop"
    )
  if (any(sample_qc_reference$sample_median_counts <= 0) ||
      any(sample_qc_reference$sample_median_features <= 0)) {
    stop("Sample-level VZ QC medians must be positive.")
  }

  cluster_sample_qc <- cluster_sample_metadata %>%
    group_by(.data$sample_id, .data$expression_cluster) %>%
    summarise(
      cell_count = n(),
      q25_counts = as.numeric(stats::quantile(.data$nCount_Xenium, 0.25)),
      median_counts = median(.data$nCount_Xenium),
      mean_counts = mean(.data$nCount_Xenium),
      q75_counts = as.numeric(stats::quantile(.data$nCount_Xenium, 0.75)),
      q25_features = as.numeric(stats::quantile(.data$nFeature_Xenium, 0.25)),
      median_features = median(.data$nFeature_Xenium),
      mean_features = mean(.data$nFeature_Xenium),
      q75_features = as.numeric(stats::quantile(.data$nFeature_Xenium, 0.75)),
      glia_cells = sum(.data$consensus_label == "Glia"),
      gaba_cells = sum(.data$consensus_label == "GABA"),
      purkinje_cells = sum(.data$consensus_label == "Purkinje"),
      opc_cells = sum(.data$consensus_label == "OPC"),
      unknown_cells = sum(.data$is_unknown),
      .groups = "drop"
    ) %>%
    left_join(sample_qc_reference, by = "sample_id") %>%
    mutate(
      configured_cells = .data$cell_count - .data$unknown_cells,
      fraction_of_sample = .data$cell_count / .data$sample_cells,
      glia_fraction = .data$glia_cells / .data$cell_count,
      gaba_fraction = .data$gaba_cells / .data$cell_count,
      purkinje_fraction = .data$purkinje_cells / .data$cell_count,
      opc_fraction = .data$opc_cells / .data$cell_count,
      unknown_fraction = .data$unknown_cells / .data$cell_count,
      median_counts_ratio_to_sample =
        .data$median_counts / .data$sample_median_counts,
      median_features_ratio_to_sample =
        .data$median_features / .data$sample_median_features,
      log2_median_counts_ratio_to_sample = log2(
        .data$median_counts_ratio_to_sample
      ),
      log2_median_features_ratio_to_sample = log2(
        .data$median_features_ratio_to_sample
      )
    ) %>%
    group_by(.data$expression_cluster) %>%
    mutate(
      cluster_cells = sum(.data$cell_count),
      fraction_of_cluster = .data$cell_count / .data$cluster_cells
    ) %>%
    ungroup() %>%
    arrange(
      factor(.data$sample_id, levels = sample_ids),
      factor(.data$expression_cluster, levels = cluster_levels)
    )
  write.csv(cluster_sample_qc, cluster_sample_summary_path, row.names = FALSE)

  heatmap_data <- cluster_sample_qc %>%
    mutate(
      sample_id = factor(.data$sample_id, levels = rev(sample_ids)),
      expression_cluster = factor(
        .data$expression_cluster, levels = cluster_levels
      )
    )
  heatmap_theme <- theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  counts_ratio_plot <- ggplot(
    heatmap_data,
    aes(
      x = .data$expression_cluster,
      y = .data$sample_id,
      fill = .data$log2_median_counts_ratio_to_sample
    )
  ) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "log2 ratio", na.value = "grey90"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "VZ raw subclusters - within-sample median transcript ratio",
      subtitle = paste(
        "Cluster median nCount_Xenium divided by the median for all",
        "VZ-branch cells from the same sample"
      ),
      x = "Resolution-0.8 expression cluster", y = "Sample"
    ) +
    heatmap_theme
  features_ratio_plot <- ggplot(
    heatmap_data,
    aes(
      x = .data$expression_cluster,
      y = .data$sample_id,
      fill = .data$log2_median_features_ratio_to_sample
    )
  ) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "log2 ratio", na.value = "grey90"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "VZ raw subclusters - within-sample median feature ratio",
      subtitle = paste(
        "Cluster median nFeature_Xenium divided by the median for all",
        "VZ-branch cells from the same sample"
      ),
      x = "Resolution-0.8 expression cluster", y = "Sample"
    ) +
    heatmap_theme
  grDevices::cairo_pdf(cluster_sample_qc_heatmap_path, width = 16, height = 14)
  print(counts_ratio_plot)
  print(features_ratio_plot)
  grDevices::dev.off()

  sample_composition_plot <- ggplot(
    heatmap_data,
    aes(
      x = .data$expression_cluster,
      y = .data$sample_id,
      fill = .data$fraction_of_cluster
    )
  ) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_viridis_c(
      option = "C", trans = "sqrt", name = "Fraction\nof cluster",
      na.value = "grey90"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = "VZ raw subclusters - sample composition",
      subtitle = "Within each expression cluster, fractions sum to 1 across samples",
      x = "Resolution-0.8 expression cluster", y = "Sample"
    ) +
    heatmap_theme
  grDevices::cairo_pdf(
    cluster_sample_composition_heatmap_path, width = 16, height = 14
  )
  print(sample_composition_plot)
  grDevices::dev.off()

  qc_stats <- obj@meta.data %>%
    group_by(.data[[cluster_column]]) %>%
    summarise(cell_count = n(), median_counts = median(nCount_Xenium),
              mean_counts = mean(nCount_Xenium),
              median_features = median(nFeature_Xenium), .groups = "drop")
  write.csv(qc_stats, summary_path, row.names = FALSE)

  grDevices::cairo_pdf(pdf_path, width = 12, height = 8)
  for (feature in c("nCount_Xenium", "nFeature_Xenium")) {
    print(
      VlnPlot(obj, features = feature, group.by = cluster_column, pt.size = 0) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        ylab(feature) + ggtitle(paste("Xenium VZ Subcluster QC -", feature))
    )
  }
  grDevices::dev.off()

  all_markers <- FindAllMarkers(
    obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "wilcox", max.cells.per.ident = 1000, random.seed = config$runtime$random_seed
  )
  top_markers <- all_markers %>% group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE)
  write.csv(top_markers, marker_path, row.names = FALSE)

  message("Generating sample-qualified Unknown composition, QC, and marker evidence...")
  unknown_qc_result <- write_xenium_vz_rl_unknown_qc(
    object = obj,
    branch = "VZ",
    cluster_column = cluster_column,
    sample_ids = sample_ids,
    output_paths = unknown_qc_paths,
    random_seed = config$runtime$random_seed
  )

  merged_info <- file.info(merged_path)
  review_manifest <- data.frame(
    merged_path = merged_path, merged_size = as.numeric(merged_info$size),
    merged_mtime = as.numeric(merged_info$mtime), cells = ncol(obj),
    cluster_column = cluster_column,
    clusters = paste(levels(Idents(obj)), collapse = "|"),
    random_seed = config$runtime$random_seed,
    sample_qualified_unknown_cells = unknown_qc_result$unknown_cells,
    sample_qualified_unknown_populations = unknown_qc_result$unknown_populations,
    unknown_marker_comparisons = unknown_qc_result$marker_comparisons,
    unknown_marker_comparisons_skipped = unknown_qc_result$skipped_marker_comparisons,
    cluster_sample_qc_rows = nrow(cluster_sample_qc),
    cluster_sample_qc_samples = dplyr::n_distinct(cluster_sample_qc$sample_id),
    stringsAsFactors = FALSE
  )
  write.csv(review_manifest, review_path, row.names = FALSE)
  message(
    "VZ QC review complete. Inspect the cluster and sample-qualified Unknown ",
    "composition, QC, and marker evidence before removal."
  )
  quit(save = "no", status = 0L)
}

if (length(missing_names)) stop("Cannot apply VZ QC; missing files:\n- ", paste(missing_names, collapse = "\n- "))
if (length(unexpected_names)) {
  stop("Cannot apply VZ QC; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
missing_review <- qc_outputs[!file.exists(qc_outputs)]
if (length(missing_review)) {
  stop("Run --qc-only and review its outputs first. Missing:\n- ", paste(missing_review, collapse = "\n- "))
}
existing_removal <- decision_outputs[file.exists(decision_outputs)]
if (length(existing_removal) && !overwrite) {
  stop("Refusing to overwrite existing VZ QC decision outputs:\n- ",
       paste(existing_removal, collapse = "\n- "), "\nUse --overwrite only after review.")
}

review <- read.csv(review_path, stringsAsFactors = FALSE)
merged_info <- file.info(merged_path)
review_file_is_current <- nrow(review) == 1L &&
  isTRUE(all.equal(as.numeric(review$merged_size), as.numeric(merged_info$size))) &&
  isTRUE(all.equal(as.numeric(review$merged_mtime), as.numeric(merged_info$mtime)))
if (!review_file_is_current) {
  stop("VZ merged object changed after QC review. Rerun --qc-only before the decision.")
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
if (!length(clusters_to_remove)) {
  removal_manifest <- build_zero_removal_routing_manifest(sample_ids, input_paths)
  inspection <- inspect_regional_qc_manifest(removal_manifest, sample_ids)
  if (!inspection$valid || !inspection$zero_removal) {
    stop(
      "Internal error while building the VZ zero-removal routing manifest:\n- ",
      paste(inspection$errors, collapse = "\n- ")
    )
  }
  write.csv(removal_manifest, removal_path, row.names = FALSE, na = "")
  message(
    "Recorded the VZ no-removal decision for all ", length(sample_ids),
    " samples without copying consensus RDS files: ", removal_path
  )
  quit(save = "no", status = 0L)
}

obj <- readRDS(merged_path)
if (!cluster_column %in% colnames(obj[[]])) stop("Merged VZ object lacks ", cluster_column, ".")
if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
  stop("Merged VZ object lacks valid sample-qualified Unknown labels. Rerun stages 2 and 3.")
}
if (!identical(as.integer(review$cells), as.integer(ncol(obj)))) {
  stop("VZ merged object cell count changed after QC review. Rerun --qc-only before removal.")
}
available_clusters <- unique(as.character(obj[[cluster_column, drop = TRUE]]))
unknown_clusters <- setdiff(clusters_to_remove, available_clusters)
if (length(unknown_clusters)) stop("Requested VZ cluster(s) do not exist: ", paste(unknown_clusters, collapse = ", "))
cells_to_remove <- if (length(clusters_to_remove)) {
  colnames(obj)[as.character(obj[[cluster_column, drop = TRUE]]) %in% clusters_to_remove]
} else character()
rm(obj)
invisible(gc())

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Materialize one sample at a time in the main R process. Parallel PSOCK
# workers previously exhausted R's pointer-protection stack while subsetting
# spatial Seurat objects, even though memory use remained below the request.
rows <- lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  message("Materializing VZ QC sample ", i, "/", length(sample_ids), ": ", sample_id)
  temp_obj <- readRDS(input_paths[[i]])
  source_info <- file.info(input_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(temp_obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  integrated_barcodes <- paste0(sample_id, "_", colnames(temp_obj))
  remove_cells <- integrated_barcodes %in% cells_to_remove
  removed <- sum(remove_cells)
  temp_obj <- subset(temp_obj, cells = colnames(temp_obj)[!remove_cells])
  saveRDS(temp_obj, output_paths[[i]], compress = FALSE)
  output_info <- file.info(output_paths[[i]])
  data.frame(
    sample_id = sample_id,
    input_cells = length(remove_cells),
    removed_cells = removed,
    remaining_cells = ncol(temp_obj),
    clusters_removed = paste(clusters_to_remove, collapse = "|"),
    source_path = normalizePath(input_paths[[i]], mustWork = TRUE),
    mapping_input_path = normalizePath(output_paths[[i]], mustWork = TRUE),
    materialization_mode = "filtered_copy",
    source_size = as.numeric(source_info$size),
    source_mtime = as.numeric(source_info$mtime),
    mapping_input_size = as.numeric(output_info$size),
    mapping_input_mtime = as.numeric(output_info$mtime),
    stringsAsFactors = FALSE
  )
})
removal_manifest <- do.call(rbind, rows)
inspection <- inspect_regional_qc_manifest(removal_manifest, sample_ids)
if (!inspection$valid) {
  stop("VZ removal manifest is invalid:\n- ", paste(inspection$errors, collapse = "\n- "))
}
if (sum(removal_manifest$removed_cells) != length(cells_to_remove)) {
  stop("VZ removal count across individual objects does not match the merged object.")
}
write.csv(removal_manifest, removal_path, row.names = FALSE, na = "")
message(
  "Applied reviewed VZ removal decision to all ", length(sample_ids),
  " samples and wrote manifest: ", removal_path
)
