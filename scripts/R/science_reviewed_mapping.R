# Pure metadata validation/mapping; no filtering or Seurat object writes.
validate_science_review <- function(review) {
  required <- c("cluster", "broad_label", "expected_cells", "action", "top20_genes")
  if (!nrow(review) || !all(required %in% names(review))) stop("Incomplete Science review table.")
  for (column in required) {
    if (anyNA(review[[column]]) || any(!nzchar(trimws(as.character(review[[column]]))))) {
      stop("Missing reviewed Science values: ", column)
    }
  }
  if (anyDuplicated(review$cluster) || any(!is.finite(review$expected_cells)) ||
      any(review$expected_cells <= 0 | review$expected_cells != floor(review$expected_cells))) {
    stop("Duplicate reviewed clusters or invalid expected counts.")
  }
  removed_labels <- c("Brainstem", "RBC", "Other")
  expected_action <- ifelse(review$broad_label %in% removed_labels, "remove", "retain")
  if (!identical(as.character(review$action), expected_action)) stop("Science removal actions differ from approved labels.")
  if (any(review$broad_label %in% c("Cycling", "Pericytes", "Endothelial"))) stop("Unapproved legacy broad labels in review.")
  invisible(TRUE)
}

map_science_reviewed_metadata <- function(metadata, review) {
  validate_science_review(review)
  required <- c("science_cluster_res2", "sample_id")
  if (!nrow(metadata) || !all(required %in% names(metadata)) ||
      is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) stop("Invalid Science cell metadata.")
  for (column in required) {
    if (anyNA(metadata[[column]]) || any(!nzchar(trimws(as.character(metadata[[column]]))))) {
      stop("Missing Science grouping values: ", column)
    }
  }
  clusters <- as.character(metadata$science_cluster_res2)
  if (!setequal(unique(clusters), as.character(review$cluster))) stop("Reviewed/source Science cluster sets differ.")
  index <- match(clusters, as.character(review$cluster))
  counts <- as.integer(table(factor(clusters, levels = review$cluster)))
  if (!identical(counts, as.integer(review$expected_cells))) stop("Science cluster cell counts differ from reviewed run.")
  # The unique lookup guarantees each cell maps to exactly one reviewed label.
  data.frame(cell_id = rownames(metadata), sample_id = as.character(metadata$sample_id),
    source_cluster = clusters, broad_label = as.character(review$broad_label[index]),
    action = as.character(review$action[index]), stringsAsFactors = FALSE)
}

check_science_reviewed_markers <- function(markers, review) {
  if (!all(c("cluster", "gene", "marker_rank") %in% names(markers)) ||
      anyNA(markers[, c("cluster", "gene", "marker_rank")]) ||
      anyDuplicated(markers[, c("cluster", "marker_rank")]) ||
      !setequal(as.character(markers$cluster), as.character(review$cluster))) {
    stop("Invalid Science source marker table.")
  }
  for (i in seq_len(nrow(review))) {
    rows <- markers[as.character(markers$cluster) == as.character(review$cluster[i]), ]
    rows <- rows[order(rows$marker_rank), ]
    if (!identical(as.integer(rows$marker_rank), 1:20) ||
        paste(rows$gene, collapse = ";") != review$top20_genes[i]) {
      stop("Source marker evidence differs from approved review for cluster ", review$cluster[i])
    }
  }
  invisible(TRUE)
}


# Harmonize labels when ingesting an older published/reviewed source. Preserve
# original metadata separately; only the standardized reference identity changes.
harmonize_vascular_reference_labels <- function(labels) {
  labels <- as.character(labels)
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Reference labels contain missing/blank values.")
  labels[labels %in% c("Endothelial", "Pericytes")] <- "Vascular"
  labels
}
