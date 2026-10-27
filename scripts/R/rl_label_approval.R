# Explicit label decisions are separate from automatically generated evidence.
# Source file metadata, review contents and ordered cell-to-cluster assignments
# bind approval to one fit without hashing/reloading the full production RDS.
rl_partition_md5 <- function(cell_ids, cluster_ids) {
  cell_ids <- as.character(cell_ids)
  cluster_ids <- as.character(cluster_ids)
  if (!length(cell_ids) || length(cell_ids) != length(cluster_ids) ||
      anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids) ||
      anyNA(cluster_ids) || any(!nzchar(cluster_ids))) {
    stop("Invalid RL cell IDs or cluster assignments for label approval.")
  }
  path <- tempfile("rl_partition_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(cell_ids = cell_ids, cluster_ids = cluster_ids), path,
          compress = FALSE, version = 2)
  unname(tools::md5sum(path))
}

rl_label_approval_template <- function(input, review_path, cell_ids, cluster_ids,
                                       cluster_column) {
  info <- file.info(input)
  if (!file.exists(input) || is.na(info$size) || info$size <= 0 ||
      !file.exists(review_path)) stop("Missing RL source/review for label approval.")
  partition_md5 <- rl_partition_md5(cell_ids, cluster_ids)
  clusters <- sort(unique(as.character(cluster_ids)))
  data.frame(
    cluster_id = clusters, subcluster_label = "", approved = FALSE,
    review_notes = "", input_rds = normalizePath(input, mustWork = TRUE),
    input_size = as.numeric(info$size), input_mtime = as.numeric(info$mtime),
    cells = length(cell_ids), cluster_column = cluster_column,
    source_partition_md5 = partition_md5,
    review_manifest_md5 = unname(tools::md5sum(review_path)),
    stringsAsFactors = FALSE
  )
}

validate_rl_label_approval <- function(approval_path, input, review_path,
                                      cluster_column, labels = NULL,
                                      cell_ids = NULL, cluster_ids = NULL) {
  if (!all(file.exists(c(approval_path, input, review_path)))) {
    stop("Missing RL label approval or source/review. Complete the fresh label ",
         "approval template after scientific review: ", approval_path)
  }
  approval <- read.csv(approval_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("cluster_id", "subcluster_label", "approved", "review_notes",
                "input_rds", "input_size", "input_mtime", "cells", "cluster_column",
                "source_partition_md5", "review_manifest_md5")
  if (!nrow(approval) || !all(required %in% names(approval))) {
    stop("RL label approval has an incomplete schema.")
  }
  approval$cluster_id <- trimws(as.character(approval$cluster_id))
  approval$subcluster_label <- trimws(as.character(approval$subcluster_label))
  approved <- suppressWarnings(as.logical(approval$approved))
  text_fields <- setdiff(required, c("input_size", "input_mtime", "cells", "approved"))
  if (anyNA(approved) || any(!approved) || anyDuplicated(approval$cluster_id) ||
      any(vapply(approval[text_fields], function(x)
        anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1)))) {
    stop("Every RL cluster requires a unique, explicitly approved label and review notes.")
  }
  review <- read.csv(review_path, stringsAsFactors = FALSE)
  if (nrow(review) != 1L ||
      !all(c("cells", "clusters", "cluster_column") %in% names(review)) ||
      anyNA(review[c("cells", "clusters", "cluster_column")]) ||
      !identical(as.character(review$cluster_column), cluster_column)) {
    stop("RL review manifest is invalid for label approval.")
  }
  expected_clusters <- strsplit(as.character(review$clusters), "|", fixed = TRUE)[[1]]
  info <- file.info(input)
  numeric_matches <- function(values, expected, tolerance = 0) {
    values <- suppressWarnings(as.numeric(values))
    length(expected) == 1L && is.finite(expected) && all(is.finite(values)) &&
      all(abs(values - as.numeric(expected)) <= tolerance)
  }
  if (!setequal(approval$cluster_id, expected_clusters) ||
      any(approval$input_rds != normalizePath(input, mustWork = TRUE)) ||
      any(approval$cluster_column != cluster_column) ||
      !numeric_matches(approval$input_size, as.numeric(info$size)) ||
      !numeric_matches(approval$input_mtime, as.numeric(info$mtime), tolerance = 1e-5) ||
      !numeric_matches(approval$cells, as.numeric(review$cells)) ||
      any(approval$review_manifest_md5 != unname(tools::md5sum(review_path))) ||
      length(unique(approval$source_partition_md5)) != 1L) {
    stop("RL label approval is stale or does not match the source/review; ",
         "approve a fresh template instead of reusing integer cluster IDs.")
  }
  if (xor(is.null(cell_ids), is.null(cluster_ids))) stop("Supply both RL cell and cluster IDs.")
  if (!is.null(cell_ids) &&
      (length(cell_ids) != as.numeric(review$cells) ||
       !setequal(as.character(cluster_ids), expected_clusters) ||
       any(approval$source_partition_md5 != rl_partition_md5(cell_ids, cluster_ids)))) {
    stop("RL approved cell-to-cluster assignments differ from the current fit.")
  }
  if (!is.null(labels)) {
    if (is.null(names(labels)) || anyDuplicated(names(labels)) || anyNA(labels) ||
        any(!nzchar(as.character(labels))) || !setequal(names(labels), approval$cluster_id) ||
        !identical(unname(as.character(labels[approval$cluster_id])), approval$subcluster_label)) {
      stop("The configured RL label map differs from the explicitly approved labels.")
    }
  }
  invisible(approval)
}
