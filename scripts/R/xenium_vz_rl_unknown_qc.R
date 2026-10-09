# Shared VZ/RL stage-4 evidence for sample-qualified Unknown populations.

xenium_vz_rl_unknown_qc_output_paths <- function(branch, plot_dir, table_dir) {
  branch <- toupper(trimws(as.character(branch)))
  if (length(branch) != 1L || !branch %in% c("VZ", "RL")) {
    stop("branch must be VZ or RL.")
  }

  c(
    composition = file.path(
      table_dir,
      paste0("Xenium_", branch, "_SampleQualifiedUnknown_ClusterComposition.csv")
    ),
    qc_summary = file.path(
      table_dir,
      paste0("Xenium_", branch, "_SampleQualifiedUnknown_QC_Summary.csv")
    ),
    top_markers = file.path(
      table_dir,
      paste0(branch, "_SampleQualifiedUnknown_top10_Markers.csv")
    ),
    marker_status = file.path(
      table_dir,
      paste0(branch, "_SampleQualifiedUnknown_Marker_Status.csv")
    ),
    qc_violins = file.path(
      plot_dir,
      paste0("Xenium_", branch, "_SampleQualifiedUnknown_QC_Violins.pdf")
    )
  )
}

write_xenium_vz_rl_unknown_qc <- function(object,
                                          branch,
                                          cluster_column,
                                          sample_ids,
                                          output_paths,
                                          random_seed) {
  branch <- toupper(trimws(as.character(branch)))
  if (length(branch) != 1L || !branch %in% c("VZ", "RL")) {
    stop("branch must be VZ or RL.")
  }
  required_paths <- c(
    "composition", "qc_summary", "top_markers", "marker_status", "qc_violins"
  )
  if (!all(required_paths %in% names(output_paths))) {
    stop("Sample-qualified Unknown QC output paths are incomplete.")
  }

  required_metadata <- c(
    "orig.ident", "consensus_label", "regional_consensus_label",
    "regional_unknown_label", cluster_column,
    "nCount_Xenium", "nFeature_Xenium"
  )
  missing_metadata <- setdiff(required_metadata, colnames(object[[]]))
  if (length(missing_metadata)) {
    stop(
      branch, " stage-4 object lacks metadata required for Unknown QC: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  metadata <- object[[]]
  metadata$sample_id <- as.character(metadata$orig.ident)
  metadata$expression_cluster <- as.character(metadata[[cluster_column]])
  metadata$regional_consensus_label <- as.character(metadata$regional_consensus_label)
  metadata$regional_unknown_label <- as.character(metadata$regional_unknown_label)
  metadata$is_unknown <- is_unknown_consensus_label(metadata$consensus_label)
  metadata$qc_population <- ifelse(
    metadata$is_unknown,
    metadata$regional_unknown_label,
    "Configured regional cells"
  )
  metadata$population_type <- ifelse(
    metadata$is_unknown, "sample-qualified Unknown", "configured regional"
  )

  if (anyNA(metadata$sample_id) || any(!metadata$sample_id %in% sample_ids)) {
    stop(branch, " stage-4 object contains invalid sample IDs.")
  }
  if (anyNA(metadata$expression_cluster) || any(!nzchar(metadata$expression_cluster))) {
    stop(branch, " stage-4 object contains blank resolution-0.8 clusters.")
  }
  if (!any(metadata$is_unknown)) {
    stop(branch, " stage-4 object contains no Unknown cells.")
  }
  if (any(metadata$is_unknown & metadata$regional_unknown_label == "Non-Unknown") ||
      any(!metadata$is_unknown & metadata$regional_unknown_label != "Non-Unknown")) {
    stop(branch, " regional_unknown_label does not match consensus Unknown status.")
  }
  if (any(!is.finite(metadata$nCount_Xenium)) ||
      any(!is.finite(metadata$nFeature_Xenium))) {
    stop(branch, " stage-4 object contains non-finite QC values.")
  }

  composition <- metadata %>%
    count(
      .data$expression_cluster, .data$sample_id,
      .data$regional_unknown_label, .data$population_type,
      name = "cells"
    ) %>%
    group_by(.data$expression_cluster) %>%
    mutate(
      cluster_cells = sum(.data$cells),
      fraction_of_expression_cluster = .data$cells / .data$cluster_cells
    ) %>%
    ungroup() %>%
    group_by(.data$sample_id, .data$regional_unknown_label) %>%
    mutate(
      population_cells = sum(.data$cells),
      fraction_of_population = .data$cells / .data$population_cells
    ) %>%
    ungroup() %>%
    arrange(.data$expression_cluster, .data$sample_id, .data$regional_unknown_label)
  write.csv(composition, output_paths[["composition"]], row.names = FALSE)

  configured_reference <- metadata %>%
    filter(!.data$is_unknown) %>%
    group_by(.data$sample_id) %>%
    summarise(
      sample_configured_cells = n(),
      sample_configured_median_counts = median(.data$nCount_Xenium),
      sample_configured_median_features = median(.data$nFeature_Xenium),
      .groups = "drop"
    )

  qc_metadata <- metadata %>%
    left_join(configured_reference, by = "sample_id")
  qc_summary <- qc_metadata %>%
    group_by(.data$sample_id, .data$qc_population, .data$population_type) %>%
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
      sample_configured_cells = dplyr::first(.data$sample_configured_cells),
      sample_configured_median_counts = dplyr::first(
        .data$sample_configured_median_counts
      ),
      sample_configured_median_features = dplyr::first(
        .data$sample_configured_median_features
      ),
      fraction_below_sample_configured_median_counts = if (
        is.finite(dplyr::first(.data$sample_configured_median_counts))
      ) {
        mean(.data$nCount_Xenium < dplyr::first(
          .data$sample_configured_median_counts
        ))
      } else NA_real_,
      fraction_below_sample_configured_median_features = if (
        is.finite(dplyr::first(.data$sample_configured_median_features))
      ) {
        mean(.data$nFeature_Xenium < dplyr::first(
          .data$sample_configured_median_features
        ))
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      median_counts_ratio_to_sample_configured = ifelse(
        is.finite(.data$sample_configured_median_counts) &
          .data$sample_configured_median_counts > 0,
        .data$median_counts / .data$sample_configured_median_counts,
        NA_real_
      ),
      median_features_ratio_to_sample_configured = ifelse(
        is.finite(.data$sample_configured_median_features) &
          .data$sample_configured_median_features > 0,
        .data$median_features / .data$sample_configured_median_features,
        NA_real_
      )
    ) %>%
    arrange(
      factor(.data$sample_id, levels = sample_ids),
      .data$population_type, .data$qc_population
    )
  write.csv(qc_summary, output_paths[["qc_summary"]], row.names = FALSE)

  unknown_lookup <- metadata %>%
    filter(.data$is_unknown) %>%
    distinct(.data$sample_id, .data$regional_consensus_label) %>%
    arrange(
      factor(.data$sample_id, levels = sample_ids),
      .data$regional_consensus_label
    )
  unknown_labels <- as.character(unknown_lookup$regional_consensus_label)
  unknown_colors <- stats::setNames(
    grDevices::hcl.colors(length(unknown_labels), "Dynamic"), unknown_labels
  )

  pdf_open <- FALSE
  on.exit({
    if (pdf_open) grDevices::dev.off()
  }, add = TRUE)
  grDevices::cairo_pdf(
    output_paths[["qc_violins"]], width = 18, height = 9
  )
  pdf_open <- TRUE
  samples_with_unknown <- sample_ids[sample_ids %in% unknown_lookup$sample_id]
  for (sample_id in samples_with_unknown) {
    sample_metadata <- qc_metadata[qc_metadata$sample_id == sample_id, , drop = FALSE]
    sample_unknown_levels <- as.character(
      unknown_lookup$regional_consensus_label[unknown_lookup$sample_id == sample_id]
    )
    group_levels <- c("Configured regional cells", sample_unknown_levels)
    sample_metadata$qc_population <- factor(
      sample_metadata$qc_population, levels = group_levels
    )
    sample_colors <- c(
      "Configured regional cells" = "#BDBDBD",
      unknown_colors[sample_unknown_levels]
    )

    for (feature in c("nCount_Xenium", "nFeature_Xenium")) {
      plot <- ggplot(
        sample_metadata,
        aes(x = .data$qc_population, y = .data[[feature]], fill = .data$qc_population)
      ) +
        geom_violin(scale = "width", trim = FALSE, na.rm = TRUE) +
        geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.65) +
        scale_fill_manual(values = sample_colors, drop = FALSE) +
        labs(
          title = paste(branch, sample_id, "sample-qualified Unknown QC"),
          subtitle = paste(
            feature,
            "| Unknown populations compared with configured regional cells from the same sample"
          ),
          x = NULL,
          y = feature
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 55, hjust = 1),
          legend.position = "none",
          plot.title = element_text(face = "bold")
        )
      print(plot)
    }
  }
  grDevices::dev.off()
  pdf_open <- FALSE

  marker_rows <- list()
  status_rows <- list()
  marker_index <- 0L
  status_index <- 0L
  for (i in seq_len(nrow(unknown_lookup))) {
    sample_id <- as.character(unknown_lookup$sample_id[[i]])
    unknown_label <- as.character(unknown_lookup$regional_consensus_label[[i]])
    sample_cells <- rownames(metadata)[metadata$sample_id == sample_id]
    sample_object <- subset(object, cells = sample_cells)
    Idents(sample_object) <- "regional_consensus_label"
    target_cells <- sum(as.character(Idents(sample_object)) == unknown_label)
    reference_cells <- ncol(sample_object) - target_cells

    status_index <- status_index + 1L
    if (target_cells < 3L || reference_cells < 3L) {
      status_rows[[status_index]] <- data.frame(
        sample_id = sample_id,
        regional_unknown_label = unknown_label,
        target_cells = target_cells,
        reference_cells = reference_cells,
        status = "skipped_too_few_cells",
        markers_passing_filters = 0L,
        top_markers_written = 0L,
        stringsAsFactors = FALSE
      )
      rm(sample_object)
      gc(verbose = FALSE)
      next
    }

    set.seed(random_seed)
    markers <- tryCatch(
      FindMarkers(
        sample_object,
        ident.1 = unknown_label,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        test.use = "wilcox",
        max.cells.per.ident = 1000,
        random.seed = random_seed
      ),
      error = function(error) error
    )
    if (inherits(markers, "error")) {
      stop(
        branch, " marker testing failed for ", unknown_label, ": ",
        conditionMessage(markers)
      )
    }

    marker_count <- nrow(markers)
    top_count <- min(10L, marker_count)
    status_rows[[status_index]] <- data.frame(
      sample_id = sample_id,
      regional_unknown_label = unknown_label,
      target_cells = target_cells,
      reference_cells = reference_cells,
      status = if (marker_count) "complete" else "no_markers_passing_filters",
      markers_passing_filters = marker_count,
      top_markers_written = top_count,
      stringsAsFactors = FALSE
    )
    if (marker_count) {
      marker_table <- data.frame(
        gene = rownames(markers), markers,
        row.names = NULL, check.names = FALSE
      )
      effect_column <- if ("avg_log2FC" %in% names(marker_table)) {
        "avg_log2FC"
      } else if ("avg_logFC" %in% names(marker_table)) {
        "avg_logFC"
      } else {
        stop(branch, " marker table lacks an average log-fold-change column.")
      }
      marker_table <- marker_table[
        order(marker_table[[effect_column]], decreasing = TRUE), , drop = FALSE
      ]
      marker_table <- head(marker_table, 10L)
      marker_index <- marker_index + 1L
      marker_rows[[marker_index]] <- cbind(
        data.frame(
          sample_id = sample_id,
          regional_unknown_label = unknown_label,
          comparison = "within_sample_vs_all_other_regional_cells",
          target_cells = target_cells,
          reference_cells = reference_cells,
          stringsAsFactors = FALSE
        ),
        marker_table
      )
    }
    rm(sample_object, markers)
    gc(verbose = FALSE)
  }

  marker_status <- dplyr::bind_rows(status_rows)
  top_markers <- if (length(marker_rows)) {
    dplyr::bind_rows(marker_rows)
  } else {
    data.frame(
      sample_id = character(), regional_unknown_label = character(),
      comparison = character(), target_cells = integer(),
      reference_cells = integer(), gene = character(),
      stringsAsFactors = FALSE
    )
  }
  write.csv(top_markers, output_paths[["top_markers"]], row.names = FALSE)
  write.csv(marker_status, output_paths[["marker_status"]], row.names = FALSE)

  list(
    unknown_cells = sum(metadata$is_unknown),
    unknown_populations = nrow(unknown_lookup),
    marker_comparisons = nrow(marker_status),
    skipped_marker_comparisons = sum(marker_status$status == "skipped_too_few_cells")
  )
}
