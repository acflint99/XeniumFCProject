# Sample-aware positive and anti-marker evidence for regional subclusters.

get_seurat_assay_layer <- function(object, assay, layer) {
  available <- SeuratObject::Layers(object[[assay]])
  if (!layer %in% available) {
    stop(
      "Assay '", assay, "' lacks the exact '", layer,
      "' layer after joining. Available layers: ",
      paste(available, collapse = ", "), "."
    )
  }
  SeuratObject::LayerData(object[[assay]], layer = layer)
}

two_sided_sign_p <- function(positive, negative) {
  n <- positive + negative
  result <- rep(NA_real_, length(n))
  eligible <- n > 0L
  result[eligible] <- pmin(
    1,
    2 * stats::pbinom(
      pmin(positive[eligible], negative[eligible]),
      size = n[eligible], prob = 0.5
    )
  )
  result
}

normalize_replicate_marker_region <- function(region) {
  region <- trimws(as.character(region))
  if (length(region) != 1L || is.na(region) || !nzchar(region)) {
    stop("region must be one nonblank value: VZ, RL, or Residual.")
  }
  canonical <- c(VZ = "VZ", RL = "RL", RESIDUAL = "Residual")
  key <- toupper(region)
  if (!key %in% names(canonical)) {
    stop("region must be VZ, RL, or Residual.")
  }
  unname(canonical[[key]])
}

build_replicate_marker_evidence <- function(region,
                                            input_path,
                                            cluster_column,
                                            output_dir,
                                            overwrite = FALSE,
                                            min_cells_per_group = 10L,
                                            top_n = 50L) {
  region <- normalize_replicate_marker_region(region)
  min_cells_per_group <- as.integer(min_cells_per_group)
  top_n <- as.integer(top_n)
  if (is.na(min_cells_per_group) || min_cells_per_group < 2L ||
      is.na(top_n) || top_n < 1L) {
    stop("Marker-evidence cell and ranking limits are invalid.")
  }

  output_paths <- c(
    all = file.path(output_dir, paste0("Xenium_", region, "_ReplicateAware_Marker_Evidence.csv")),
    positive = file.path(output_dir, paste0("Xenium_", region, "_Positive_Marker_Candidates.csv")),
    anti = file.path(output_dir, paste0("Xenium_", region, "_AntiMarker_Candidates.csv")),
    coverage = file.path(output_dir, paste0("Xenium_", region, "_Marker_Sample_Coverage.csv")),
    provenance = file.path(output_dir, paste0("Xenium_", region, "_ReplicateAware_Marker_Provenance.csv"))
  )
  existing <- output_paths[file.exists(output_paths)]
  if (length(existing) && !overwrite) {
    stop(
      "Refusing to overwrite replicate-aware marker output(s):\n- ",
      paste(existing, collapse = "\n- "),
      "\nRerun with --overwrite only after reviewing the existing files."
    )
  }
  if (!file.exists(input_path)) stop("Marker-evidence input not found: ", input_path)

  object <- readRDS(input_path)
  required_metadata <- c("orig.ident", cluster_column)
  missing <- setdiff(required_metadata, names(object[[]]))
  if (length(missing)) stop("Marker-evidence input lacks: ", paste(missing, collapse = ", "))
  if (anyDuplicated(Seurat::Cells(object))) stop("Marker-evidence input has duplicate cell IDs.")
  samples <- trimws(as.character(object$orig.ident))
  clusters <- trimws(as.character(object[[cluster_column, drop = TRUE]]))
  if (anyNA(samples) || any(!nzchar(samples)) ||
      anyNA(clusters) || any(!nzchar(clusters))) {
    stop("Marker-evidence sample or cluster metadata contains blank values.")
  }
  assay <- if ("Xenium" %in% Seurat::Assays(object)) "Xenium" else Seurat::DefaultAssay(object)
  preserved_cell_ids <- Seurat::Cells(object)
  object[[assay]] <- SeuratObject::JoinLayers(object[[assay]])
  if (!identical(Seurat::Cells(object), preserved_cell_ids)) {
    stop("Joining marker-evidence assay layers changed cell IDs or their order.")
  }
  normalized <- get_seurat_assay_layer(object, assay, "data")
  counts <- get_seurat_assay_layer(object, assay, "counts")
  if (!identical(dim(normalized), dim(counts)) ||
      !identical(rownames(normalized), rownames(counts)) ||
      !identical(colnames(normalized), Seurat::Cells(object)) ||
      !identical(colnames(counts), Seurat::Cells(object))) {
    stop("Marker-evidence assay layers do not align one-to-one with cells and genes.")
  }

  group_keys <- paste(samples, clusters, sep = "\u001f")
  group_levels <- unique(group_keys)
  group_index <- match(group_keys, group_levels)
  design <- Matrix::sparseMatrix(
    i = seq_along(group_index), j = group_index, x = 1,
    dims = c(length(group_index), length(group_levels)),
    dimnames = list(colnames(normalized), group_levels)
  )
  group_cells <- as.integer(Matrix::colSums(design))
  group_parts <- strsplit(group_levels, "\u001f", fixed = TRUE)
  group_table <- data.frame(
    group_index = seq_along(group_levels),
    sample_id = vapply(group_parts, `[[`, character(1), 1L),
    cluster = vapply(group_parts, `[[`, character(1), 2L),
    n_cells = group_cells,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(paste(group_table$sample_id, group_table$cluster))) {
    stop("Sample/cluster aggregation produced duplicate groups.")
  }

  expression_sum <- normalized %*% design
  detection_sum <- (counts > 0) %*% design
  genes <- rownames(normalized)
  cluster_levels <- unique(clusters)
  cluster_numeric <- suppressWarnings(as.numeric(cluster_levels))
  cluster_levels <- cluster_levels[order(
    is.na(cluster_numeric), cluster_numeric, cluster_levels
  )]
  evidence <- vector("list", length(cluster_levels))
  coverage <- vector("list", length(cluster_levels))

  for (cluster_position in seq_along(cluster_levels)) {
    cluster_id <- cluster_levels[[cluster_position]]
    target_groups <- which(group_table$cluster == cluster_id)
    candidate_samples <- group_table$sample_id[target_groups]
    eligible_samples <- candidate_samples[vapply(candidate_samples, function(sample_id) {
      target <- which(group_table$sample_id == sample_id & group_table$cluster == cluster_id)
      background <- which(group_table$sample_id == sample_id & group_table$cluster != cluster_id)
      length(target) == 1L && group_cells[target] >= min_cells_per_group &&
        length(background) > 0L && sum(group_cells[background]) >= min_cells_per_group
    }, logical(1))]
    eligible_samples <- unique(eligible_samples)
    coverage[[cluster_position]] <- data.frame(
      cluster = cluster_id,
      samples_with_cluster = length(unique(candidate_samples)),
      eligible_samples = length(eligible_samples),
      min_cells_per_target_and_background = min_cells_per_group,
      eligible_sample_ids = paste(sort(eligible_samples), collapse = "|"),
      stringsAsFactors = FALSE
    )
    if (!length(eligible_samples)) next

    expression_difference <- matrix(
      NA_real_, nrow = length(genes), ncol = length(eligible_samples),
      dimnames = list(genes, eligible_samples)
    )
    detection_difference <- expression_difference
    target_detection <- expression_difference
    background_detection <- expression_difference

    for (sample_position in seq_along(eligible_samples)) {
      sample_id <- eligible_samples[[sample_position]]
      target <- which(group_table$sample_id == sample_id & group_table$cluster == cluster_id)
      background <- which(group_table$sample_id == sample_id & group_table$cluster != cluster_id)
      target_n <- group_cells[target]
      background_n <- sum(group_cells[background])
      target_expression <- as.numeric(expression_sum[, target, drop = FALSE]) / target_n
      background_expression <- as.numeric(Matrix::rowSums(
        expression_sum[, background, drop = FALSE]
      )) / background_n
      target_pct <- as.numeric(detection_sum[, target, drop = FALSE]) / target_n
      background_pct <- as.numeric(Matrix::rowSums(
        detection_sum[, background, drop = FALSE]
      )) / background_n
      expression_difference[, sample_position] <- target_expression - background_expression
      detection_difference[, sample_position] <- target_pct - background_pct
      target_detection[, sample_position] <- target_pct
      background_detection[, sample_position] <- background_pct
    }

    positive_n <- Matrix::rowSums(expression_difference > 0, na.rm = TRUE)
    negative_n <- Matrix::rowSums(expression_difference < 0, na.rm = TRUE)
    nonzero_n <- positive_n + negative_n
    sign_p <- two_sided_sign_p(positive_n, negative_n)
    median_effect <- apply(expression_difference, 1L, stats::median, na.rm = TRUE)
    median_detection_difference <- apply(
      detection_difference, 1L, stats::median, na.rm = TRUE
    )
    median_target_detection <- apply(
      target_detection, 1L, stats::median, na.rm = TRUE
    )
    median_background_detection <- apply(
      background_detection, 1L, stats::median, na.rm = TRUE
    )
    direction <- ifelse(median_effect >= 0, "positive_marker", "anti_marker")
    direction_support_n <- ifelse(direction == "positive_marker", positive_n, negative_n)
    direction_support_fraction <- direction_support_n / length(eligible_samples)
    cluster_evidence <- data.frame(
      region = region,
      cluster = cluster_id,
      gene = genes,
      candidate_type = direction,
      n_eligible_samples = length(eligible_samples),
      n_nonzero_sample_effects = nonzero_n,
      direction_support_n = direction_support_n,
      direction_support_fraction = direction_support_fraction,
      median_sample_mean_log_expression_difference = median_effect,
      median_target_detection_fraction = median_target_detection,
      median_background_detection_fraction = median_background_detection,
      median_detection_fraction_difference = median_detection_difference,
      sample_level_sign_p = sign_p,
      sample_level_sign_p_adj = stats::p.adjust(sign_p, method = "BH"),
      stringsAsFactors = FALSE
    )
    cluster_evidence$evidence_rank_score <- (
      cluster_evidence$direction_support_fraction *
        abs(cluster_evidence$median_sample_mean_log_expression_difference) *
        (1 + abs(cluster_evidence$median_detection_fraction_difference))
    )
    evidence[[cluster_position]] <- cluster_evidence
  }

  all_evidence <- do.call(rbind, evidence[!vapply(evidence, is.null, logical(1))])
  if (is.null(all_evidence) || !nrow(all_evidence)) {
    stop("No cluster had enough within-sample target and background cells for marker evidence.")
  }
  all_evidence <- all_evidence[order(
    match(all_evidence$cluster, cluster_levels),
    all_evidence$candidate_type,
    -all_evidence$evidence_rank_score,
    all_evidence$gene
  ), , drop = FALSE]
  rownames(all_evidence) <- NULL

  top_candidates <- function(type) {
    subset <- all_evidence[all_evidence$candidate_type == type, , drop = FALSE]
    if (!nrow(subset)) {
      subset$rank_within_cluster <- integer()
      return(subset)
    }
    groups <- split(seq_len(nrow(subset)), subset$cluster)
    indices <- unlist(lapply(groups, function(index) {
      index[order(-subset$evidence_rank_score[index], subset$gene[index])][
        seq_len(min(top_n, length(index)))
      ]
    }), use.names = FALSE)
    result <- subset[indices, , drop = FALSE]
    result$rank_within_cluster <- ave(
      result$evidence_rank_score, result$cluster,
      FUN = function(x) rank(-x, ties.method = "first")
    )
    result
  }
  positive <- top_candidates("positive_marker")
  anti <- top_candidates("anti_marker")
  coverage_table <- do.call(rbind, coverage)
  input_info <- file.info(input_path)
  provenance <- data.frame(
    region = region,
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_rds = normalizePath(input_path, mustWork = TRUE),
    input_size = as.numeric(input_info$size),
    input_mtime = as.numeric(input_info$mtime),
    cells = ncol(object),
    genes = length(genes),
    assay = assay,
    assay_layers_joined_for_marker_evidence = TRUE,
    sample_column = "orig.ident",
    cluster_column = cluster_column,
    min_cells_per_target_and_background = min_cells_per_group,
    top_candidates_per_cluster_and_direction = top_n,
    inference_unit = "biological sample",
    evidence_test = "two-sided sign test across within-sample target-minus-background effects",
    cell_level_p_values_used_for_inference = FALSE,
    stringsAsFactors = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_evidence, output_paths[["all"]], row.names = FALSE)
  write.csv(positive, output_paths[["positive"]], row.names = FALSE)
  write.csv(anti, output_paths[["anti"]], row.names = FALSE)
  write.csv(coverage_table, output_paths[["coverage"]], row.names = FALSE)
  write.csv(provenance, output_paths[["provenance"]], row.names = FALSE)
  invisible(output_paths)
}

# A review can generate absent evidence, reuse matching evidence, or stop before
# touching partial/stale files. This inspection never loads an analysis object.
plan_replicate_marker_review <- function(paths, provenance_matches, overwrite = FALSE) {
  info <- file.info(paths)
  complete <- all(file.exists(paths) & !is.na(info$size) & info$size > 0)
  if (complete && isTRUE(provenance_matches)) return("reuse")
  if (!any(file.exists(paths))) return("build")
  if (isTRUE(overwrite)) return("rebuild")
  "blocked"
}
