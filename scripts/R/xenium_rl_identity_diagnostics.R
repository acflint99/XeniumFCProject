# Descriptive diagnostics only: no filtering, relabeling, or score thresholds.
rl_diagnostic_qc <- function(metadata, group_columns) {
  metadata |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::summarise(
      cells = dplyr::n(),
      q25_counts = unname(stats::quantile(nCount_Xenium, 0.25)),
      median_counts = stats::median(nCount_Xenium),
      q75_counts = unname(stats::quantile(nCount_Xenium, 0.75)),
      median_features = stats::median(nFeature_Xenium),
      q25_features = unname(stats::quantile(nFeature_Xenium, 0.25)),
      q75_features = unname(stats::quantile(nFeature_Xenium, 0.75)),
      .groups = "drop"
    )
}

rl_diagnostic_pairs <- function(counts, metadata, pairs, group_columns) {
  available <- pairs[pairs$gene1 %in% rownames(counts) &
                       pairs$gene2 %in% rownames(counts), , drop = FALSE]
  if (!nrow(available)) stop("No diagnostic marker pairs are available in the panel.")
  rows <- lapply(seq_len(nrow(available)), function(i) {
    a <- as.numeric(counts[available$gene1[i], ]) > 0
    b <- as.numeric(counts[available$gene2[i], ]) > 0
    frame <- metadata[, group_columns, drop = FALSE]
    frame$gene1_detected <- a
    frame$gene2_detected <- b
    result <- frame |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
      dplyr::summarise(
        cells = dplyr::n(),
        gene1_cells = sum(gene1_detected), gene2_cells = sum(gene2_detected),
        both_cells = sum(gene1_detected & gene2_detected),
        neither_cells = sum(!gene1_detected & !gene2_detected),
        .groups = "drop"
      )
    result$gene1 <- available$gene1[i]
    result$gene2 <- available$gene2[i]
    result$gene1_fraction <- result$gene1_cells / result$cells
    result$gene2_fraction <- result$gene2_cells / result$cells
    result$both_fraction <- result$both_cells / result$cells
    result$both_fraction_among_gene1 <- ifelse(
      result$gene1_cells > 0, result$both_cells / result$gene1_cells, NA_real_
    )
    result
  })
  dplyr::bind_rows(rows)
}

rl_diagnostic_reference <- function(frame, group_columns) {
  denominator <- frame |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::summarise(cluster_cells = dplyr::n(), .groups = "drop")
  result <- frame |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_columns, "predicted_label")))) |>
    dplyr::summarise(
      cells = dplyr::n(), median_score = stats::median(prediction_score),
      q25_score = unname(stats::quantile(prediction_score, 0.25)),
      q75_score = unname(stats::quantile(prediction_score, 0.75)),
      .groups = "drop"
    ) |>
    dplyr::left_join(denominator, by = group_columns)
  if (anyNA(result$cluster_cells)) stop("Reference denominator join failed.")
  result$fraction_of_cluster <- result$cells / result$cluster_cells
  result
}

rl_diagnostic_match_cells <- function(merged_cell_ids, reference_cell_ids, sample_id) {
  if (anyDuplicated(merged_cell_ids) || anyDuplicated(reference_cell_ids)) {
    stop("Duplicate cell IDs prevent one-to-one reference matching for ", sample_id)
  }
  # Stage 2 uses RenameCells(add.cell.id = sample_id). Reproduce that exact
  # prefix; do not strip suffixes, guess IDs, or match by row position.
  prefixed_ids <- paste0(sample_id, "_", reference_cell_ids)
  index <- match(merged_cell_ids, prefixed_ids)
  if (anyNA(index) || anyDuplicated(index)) {
    stop("Reference cell matching is incomplete or not one-to-one for ", sample_id)
  }
  index
}
