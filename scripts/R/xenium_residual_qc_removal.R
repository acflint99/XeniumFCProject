# Reviewed exclusions are bound to the pre-QC object, not a reusable cluster number.
residual_qc_removal_paths <- function(config) {
  root <- here::here(config$project$outputs_dir, "xenium", "residual")
  qc <- file.path(root, "04_qc_filtered")
  list(
    decision = resolve_config_path("config/reviewed_residual_qc_removals.csv", config),
    source = file.path(root, "03_integrated", "rds", "Xenium_Residual_Integrated.rds"),
    qc_cells = file.path(root, "04_qc_review", "tables", "Xenium_Residual_QC_PerCell.csv"),
    qc_provenance = file.path(root, "04_qc_review", "tables", "Xenium_Residual_QC_Provenance.csv"),
    filtered = file.path(qc, "rds", "Xenium_Residual_QC_Filtered.rds"),
    removed = file.path(qc, "tables", "Xenium_Residual_QC_Removed_Cells.csv"),
    manifest = file.path(qc, "tables", "Xenium_Residual_QC_Removal_Manifest.csv"),
    provenance = file.path(qc, "tables", "Xenium_Residual_QC_Removal_Provenance.csv"),
    markers = file.path(root, "04_post_qc_subcluster_review", "tables"),
    original_markers = file.path(root, "04_subcluster_review", "tables")
  )
}

residual_qc_file_matches <- function(path, size, mtime) {
  info <- file.info(path)
  isTRUE(as.numeric(info$size) == as.numeric(size)) &&
    isTRUE(abs(as.numeric(info$mtime) - as.numeric(mtime)) <= 1e-4)
}

read_residual_qc_decision <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  required <- c("cluster", "cluster_column", "expected_input_cells", "expected_removed_cells",
                "source_size", "source_mtime", "review_notes", "approved")
  if (!all(required %in% names(d)) || !nrow(d)) stop("Incomplete residual QC decision.")
  for (name in c("cluster", "cluster_column", "review_notes")) {
    if (anyNA(d[[name]]) || any(!nzchar(trimws(as.character(d[[name]]))))) stop("Blank decision ", name)
  }
  if (anyDuplicated(d$cluster) || anyNA(as.logical(d$approved)) || any(!as.logical(d$approved))) {
    stop("Residual QC decisions must be unique and explicitly approved.")
  }
  for (name in c("expected_input_cells", "expected_removed_cells", "source_size", "source_mtime")) {
    x <- suppressWarnings(as.numeric(d[[name]]))
    if (any(!is.finite(x)) || any(x <= 0)) stop("Invalid decision ", name)
    d[[name]] <- x
  }
  for (name in c("expected_input_cells", "expected_removed_cells")) {
    if (any(d[[name]] != round(d[[name]]))) stop("Noninteger decision cell count.")
  }
  if (length(unique(d$cluster_column)) != 1L ||
      length(unique(d$expected_input_cells)) != 1L ||
      length(unique(d$source_size)) != 1L || length(unique(d$source_mtime)) != 1L ||
      sum(d$expected_removed_cells) >= d$expected_input_cells[[1]]) {
    stop("Inconsistent source or cell counts in residual QC decision.")
  }
  d$cluster <- as.character(d$cluster)
  d
}

load_residual_qc_state <- function(config) {
  p <- residual_qc_removal_paths(config)
  empty <- data.frame(cell_id = character(), sample_id = character(), cluster = character(),
                       consensus_label = character())
  if (!file.exists(p$decision)) {
    if (any(file.exists(unlist(p[c("filtered", "removed", "manifest", "provenance")])))) {
      stop("Residual QC outputs exist without their approved decision.")
    }
    return(list(active = FALSE, input = p$source, markers = p$original_markers,
                removed = empty, paths = p, decision_md5 = ""))
  }
  d <- read_residual_qc_decision(p$decision)
  required <- unlist(p[c("source", "filtered", "removed", "manifest", "provenance")])
  if (any(!file.exists(required))) stop("Approved residual QC removal must be applied first. Missing: ",
                                      paste(required[!file.exists(required)], collapse = ", "))
  provenance <- read.csv(p$provenance, stringsAsFactors = FALSE)
  if (nrow(provenance) != 1L ||
      !all(c("decision_md5", "filtered_size", "filtered_mtime", "source_size", "source_mtime") %in% names(provenance))) {
    stop("Invalid residual QC removal provenance.")
  }
  digest <- unname(tools::md5sum(p$decision))
  if (!identical(provenance$decision_md5, digest) ||
      !residual_qc_file_matches(p$source, d$source_size[[1]], d$source_mtime[[1]]) ||
      !residual_qc_file_matches(p$source, provenance$source_size, provenance$source_mtime) ||
      !residual_qc_file_matches(p$filtered, provenance$filtered_size, provenance$filtered_mtime)) {
    stop("Residual QC decision or source/filtered object changed; review and regenerate QC outputs.")
  }
  removed <- read.csv(p$removed, stringsAsFactors = FALSE, na.strings = character())
  if (!all(names(empty) %in% names(removed)) || anyDuplicated(removed$cell_id) ||
      anyNA(removed$cell_id) || any(!nzchar(removed$cell_id)) ||
      nrow(removed) != sum(d$expected_removed_cells)) stop("Invalid removed-cell manifest.")
  for (i in seq_len(nrow(d))) {
    if (sum(as.character(removed$cluster) == d$cluster[[i]]) != d$expected_removed_cells[[i]]) {
      stop("Removed-cell cluster counts differ from approved decision.")
    }
  }
  ids <- as.character(load_sample_manifest(config)$sample_id)
  if (any(!removed$sample_id %in% ids) || any(!vapply(seq_len(nrow(removed)), function(i) {
    startsWith(removed$cell_id[[i]], paste0(removed$sample_id[[i]], "_"))
  }, logical(1)))) stop("Removed-cell sample IDs/prefixes are invalid.")
  m <- read.csv(p$manifest, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "input_cells", "removed_cells", "remaining_cells") %in% names(m)) ||
      nrow(m) != length(ids) || anyDuplicated(m$sample_id) || !setequal(m$sample_id, ids)) {
    stop("Residual QC removal manifest must map one-to-one to the active sample manifest.")
  }
  m <- m[match(ids, m$sample_id), , drop = FALSE]
  for (name in c("input_cells", "removed_cells", "remaining_cells")) {
    x <- as.numeric(m[[name]])
    if (any(!is.finite(x)) || any(x < 0 | x != round(x))) stop("Invalid QC manifest cell counts.")
  }
  observed <- as.integer(table(factor(removed$sample_id, levels = ids)))
  if (any(m$input_cells - m$removed_cells != m$remaining_cells) ||
      any(m$removed_cells != observed) || sum(m$input_cells) != d$expected_input_cells[[1]]) {
    stop("QC removal sample totals do not reconcile.")
  }
  list(active = TRUE, input = p$filtered, markers = p$markers, removed = removed,
       manifest = m, paths = p, decision_md5 = digest)
}

validate_residual_qc_object <- function(object, state) {
  if (!state$active) return(invisible(TRUE))
  stored <- object@misc$residual_qc_removal
  if (is.null(stored) || !identical(stored$decision_md5, state$decision_md5) ||
      any(colnames(object) %in% state$removed$cell_id) ||
      ncol(object) != sum(state$manifest$remaining_cells)) {
    stop("Residual object does not carry the current approved QC exclusion.")
  }
  invisible(TRUE)
}

residual_qc_exclusions_for_sample <- function(cell_ids, consensus, sample_id, state) {
  if (length(cell_ids) != length(consensus) || anyNA(cell_ids) ||
      any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) stop("Invalid base cell IDs for ", sample_id)
  rows <- state$removed[state$removed$sample_id == sample_id, , drop = FALSE]
  qualified <- paste0(sample_id, "_", cell_ids)
  index <- match(rows$cell_id, qualified)
  if (anyNA(index) ||
      !identical(as.character(consensus[index]), as.character(rows$consensus_label)) ||
      any(!is_unknown_consensus_label(consensus[index]))) {
    stop(sample_id, ": approved residual QC exclusions are absent or have changed consensus labels.")
  }
  qualified %in% rows$cell_id
}

# A separate Stage-5 source prevents old filtered partitions from being reused
# as though they had been refitted after exclusions.
residual_post_qc_paths <- function(config) {
  root <- here::here(config$project$outputs_dir, "xenium", "residual", "05_post_qc")
  list(input = file.path(root, "rds", "Xenium_Residual_postQC.rds"),
       provenance = file.path(root, "tables", "Xenium_Residual_postQC_manifest.csv"),
       coverage = file.path(root, "tables", "Xenium_Residual_postQC_sample_counts.csv"),
       markers = file.path(root, "marker_review", "tables"))
}

load_residual_analysis_state <- function(config) {
  state <- load_residual_qc_state(config)
  p <- residual_post_qc_paths(config)
  required <- unlist(p[c("input", "provenance", "coverage")])
  if (any(!file.exists(required))) stop("Run residual post-QC processing before review/labels/mapping. Missing: ",
                                        paste(required[!file.exists(required)], collapse = ", "))
  m <- read.csv(p$provenance, stringsAsFactors = FALSE, na.strings = character(),
                colClasses = c(decision_md5 = "character", sample_manifest_md5 = "character", coverage_md5 = "character"))
  fields <- c("input_rds", "input_size", "input_mtime", "output_rds", "output_size", "output_mtime",
              "decision_md5", "cells_after_qc", "resolutions", "random_seed", "sample_manifest_md5", "coverage_md5")
  if (nrow(m) != 1L || !all(fields %in% names(m)) ||
      !identical(as.character(m$decision_md5), state$decision_md5) ||
      !identical(normalizePath(m$input_rds, mustWork = TRUE), normalizePath(state$input, mustWork = TRUE)) ||
      !identical(normalizePath(m$output_rds, mustWork = TRUE), normalizePath(p$input, mustWork = TRUE)) ||
      !residual_qc_file_matches(state$input, m$input_size, m$input_mtime) ||
      !residual_qc_file_matches(p$input, m$output_size, m$output_mtime) ||
      !identical(as.integer(m$random_seed), as.integer(config$runtime$random_seed)) ||
      !identical(as.character(m$sample_manifest_md5), unname(tools::md5sum(resolve_config_path(config$manifests$samples, config)))) ||
      !identical(as.character(m$coverage_md5), unname(tools::md5sum(p$coverage)))) {
    stop("Residual post-QC provenance is stale; regenerate Stage 5 before review.")
  }
  ids <- as.character(load_sample_manifest(config)$sample_id)
  coverage <- read.csv(p$coverage, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "remaining_cells") %in% names(coverage)) ||
      !identical(as.character(coverage$sample_id), ids) || anyDuplicated(coverage$sample_id)) stop("Invalid post-QC sample map.")
  counts <- suppressWarnings(as.numeric(coverage$remaining_cells))
  if (any(!is.finite(counts)) || any(counts < 0 | counts != floor(counts)) ||
      !isTRUE(sum(counts) == as.numeric(m$cells_after_qc)) ||
      (state$active && !identical(as.integer(counts), as.integer(state$manifest$remaining_cells)))) {
    stop("Post-QC counts differ from the reviewed manifest.")
  }
  state$qc_input <- state$input
  state$input <- p$input
  state$markers <- p$markers
  state$post_qc_paths <- p
  state$post_qc_manifest <- m
  state$expected_counts <- setNames(as.integer(counts), ids)
  state
}
