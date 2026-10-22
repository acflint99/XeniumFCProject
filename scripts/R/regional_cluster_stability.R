# Small, dependency-free calculations for the regional resolution review.
# Jaccard matching follows Hennig (2007); comparisons use retained cells only.
regional_stability_tree <- function(assignments, prefix = "Xenium_snn_res.") {
  # A string guide ("edge_colourbar") can fail ggplot2's guide lookup when
  # ggraph is loaded as a namespace but not attached. Supply the guide object
  # directly; preserve clustree's mappings, colours, edges and node positions.
  clustree::clustree(assignments, prefix = prefix) +
    ggplot2::guides(edge_colour = ggraph::guide_edge_colourbar())
}

regional_stability_check_plot <- function() {
  assignments <- data.frame(
    Xenium_snn_res.0.3 = rep(0:1, each = 6),
    Xenium_snn_res.0.5 = rep(0:3, each = 3),
    Xenium_snn_res.0.8 = rep(0:5, each = 2), check.names = FALSE
  )
  path <- tempfile(pattern = "regional_clustree_check_", fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  ggplot2::ggsave(path, regional_stability_tree(assignments), width = 8, height = 6)
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Synthetic clustree rendering did not produce a nonempty PDF.")
  }
  message("Synthetic clustree PDF rendering passed (12 artificial cells; no project RDS loaded).")
  invisible(TRUE)
}

regional_stability_sample <- function(sample_ids, fraction, seed) {
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      length(fraction) != 1L || !is.finite(fraction) || fraction <= 0 || fraction >= 1) {
    stop("Sampling requires nonblank sample IDs and a fraction strictly between 0 and 1.")
  }
  set.seed(seed)
  indices <- split(seq_along(sample_ids), sample_ids)
  sort(unlist(lapply(indices, function(i) {
    # sample.int avoids sample(x) treating a length-one numeric x as 1:x.
    i[sample.int(length(i), max(1L, floor(length(i) * fraction)))]
  }), use.names = FALSE))
}

regional_stability_jaccard <- function(reference, recovered, all_clusters) {
  if (!length(reference) || length(reference) != length(recovered) ||
      anyNA(reference) || anyNA(recovered) || any(!nzchar(reference)) ||
      any(!nzchar(recovered)) || anyDuplicated(all_clusters) ||
      !all(reference %in% all_clusters)) stop("Invalid aligned cluster labels.")
  overlap <- table(factor(reference, levels = all_clusters), as.character(recovered))
  reference_n <- rowSums(overlap)
  recovered_n <- colSums(overlap)
  unions <- outer(reference_n, recovered_n, "+") - overlap
  score <- overlap / unions
  best <- max.col(score, ties.method = "first")
  result <- data.frame(
    cluster = all_clusters, retained_cells = as.integer(reference_n),
    matched_cluster = colnames(overlap)[best],
    intersection = as.integer(overlap[cbind(seq_along(best), best)]),
    matched_cells = as.integer(recovered_n[best]),
    jaccard = as.numeric(score[cbind(seq_along(best), best)]),
    stringsAsFactors = FALSE
  )
  absent <- reference_n == 0L
  result$matched_cluster[absent] <- NA_character_
  result$intersection[absent] <- NA_integer_
  result$matched_cells[absent] <- NA_integer_
  result$jaccard[absent] <- NA_real_
  result
}
