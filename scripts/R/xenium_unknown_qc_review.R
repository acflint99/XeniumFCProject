# Shared validation and routing helpers for the post-consensus Unknown-QC gate.

unknown_qc_review_key_columns <- function() {
  c("sample_id", "consensus_label")
}

unknown_qc_review_key <- function(data) {
  required <- unknown_qc_review_key_columns()
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Unknown-QC data lack key column(s): ", paste(missing, collapse = ", "))
  }
  values <- lapply(data[required], function(value) trimws(as.character(value)))
  if (any(vapply(values, function(value) anyNA(value) || any(!nzchar(value)), logical(1)))) {
    stop("Unknown-QC keys contain missing or blank values.")
  }
  do.call(paste, c(values, sep = "\u001f"))
}

validate_unknown_qc_review <- function(review, sample_ids, expected_keys = NULL,
                                       require_complete = TRUE) {
  required <- c(unknown_qc_review_key_columns(), "decision", "review_notes")
  missing <- setdiff(required, names(review))
  if (length(missing)) {
    stop("Reviewed Unknown-QC manifest lacks column(s): ", paste(missing, collapse = ", "))
  }
  if (!nrow(review)) {
    if (require_complete) {
      stop(
        "Reviewed Unknown-QC manifest has no decisions. Copy the generated review ",
        "template to the configured manifest and fill decision plus review_notes."
      )
    }
    return(invisible(review))
  }

  keys <- unknown_qc_review_key(review)
  if (anyDuplicated(keys)) stop("Reviewed Unknown-QC manifest contains duplicate keys.")
  if (any(!review$sample_id %in% sample_ids)) {
    stop("Reviewed Unknown-QC manifest contains sample IDs absent from config/samples.csv.")
  }
  if (any(!is_unknown_consensus_label(review$consensus_label))) {
    stop("Reviewed Unknown-QC manifest contains a non-Unknown consensus label.")
  }

  decisions <- tolower(trimws(as.character(review$decision)))
  notes <- trimws(as.character(review$review_notes))
  if (require_complete && (anyNA(decisions) || any(!decisions %in% c("keep", "remove")))) {
    stop("Every reviewed Unknown-QC decision must be 'keep' or 'remove'.")
  }
  if (require_complete && (anyNA(notes) || any(!nzchar(notes)))) {
    stop("Every reviewed Unknown-QC row must have nonblank review_notes.")
  }
  if (!is.null(expected_keys) && !setequal(keys, expected_keys)) {
    stop(
      "Reviewed Unknown-QC decisions do not exactly match the current audit template. ",
      "Missing decisions: ", length(setdiff(expected_keys, keys)),
      "; unused decisions: ", length(setdiff(keys, expected_keys)), "."
    )
  }
  invisible(review)
}

unknown_qc_routing_columns <- function() {
  c(
    "sample_id", "source_path", "mapping_input_path", "materialization_mode",
    "reviewed_qc_path", "reviewed_qc_size", "reviewed_qc_mtime",
    "source_size", "source_mtime", "mapping_input_size", "mapping_input_mtime",
    "input_cells", "removed_unknown_cells", "remaining_cells",
    "removed_unknown_labels"
  )
}

inspect_unknown_qc_routing_manifest <- function(manifest, sample_ids) {
  errors <- character()
  required <- unknown_qc_routing_columns()
  if (is.null(manifest)) return(list(valid = FALSE, errors = "manifest is missing"))
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    return(list(valid = FALSE, errors = paste0(
      "missing routing column(s): ", paste(missing, collapse = ", ")
    )))
  }
  if (nrow(manifest) != length(sample_ids) ||
      !identical(as.character(manifest$sample_id), sample_ids)) {
    errors <- c(errors, "sample rows must exactly match config/samples.csv order")
  }
  if (anyDuplicated(manifest$source_path) || anyDuplicated(manifest$mapping_input_path)) {
    errors <- c(errors, "routing paths must be one-to-one with sample IDs")
  }
  modes <- as.character(manifest$materialization_mode)
  if (anyNA(modes) || any(!modes %in% c("reuse_consensus", "filtered_copy"))) {
    errors <- c(errors, "materialization_mode must be reuse_consensus or filtered_copy")
  }
  numeric_columns <- c(
    "reviewed_qc_size", "reviewed_qc_mtime",
    "source_size", "source_mtime", "mapping_input_size", "mapping_input_mtime",
    "input_cells", "removed_unknown_cells", "remaining_cells"
  )
  numeric_values <- lapply(
    manifest[numeric_columns], function(value) suppressWarnings(as.numeric(value))
  )
  if (any(vapply(numeric_values, function(value) any(!is.finite(value)), logical(1)))) {
    errors <- c(errors, "routing numeric provenance and cell counts must be finite")
  } else {
    if (any(numeric_values$reviewed_qc_size <= 0 |
            numeric_values$source_size <= 0 |
            numeric_values$mapping_input_size <= 0)) {
      errors <- c(errors, "routing file sizes must be positive")
    }
    if (any(numeric_values$input_cells - numeric_values$removed_unknown_cells !=
            numeric_values$remaining_cells)) {
      errors <- c(errors, "routing cell counts are inconsistent")
    }
  }
  removed_labels <- as.character(manifest$removed_unknown_labels)
  if (anyNA(removed_labels) || any(!nzchar(removed_labels))) {
    errors <- c(errors, "removed_unknown_labels contains missing or blank values")
  }
  reuse <- !is.na(modes) & modes == "reuse_consensus"
  if (any(reuse) && any(
    suppressWarnings(as.numeric(manifest$removed_unknown_cells[reuse])) != 0 |
      removed_labels[reuse] != "none"
  )) {
    errors <- c(errors, "reuse_consensus rows must record zero removals")
  }
  review_paths <- as.character(manifest$reviewed_qc_path)
  if (anyNA(review_paths) || any(!nzchar(review_paths)) ||
      length(unique(review_paths)) != 1L) {
    errors <- c(errors, "all rows must record the same reviewed_qc_path")
  }
  list(valid = !length(errors), errors = unique(errors))
}

load_unknown_qc_routing_manifest <- function(config = load_pipeline_config()) {
  output_root <- here::here(config$project$outputs_dir)
  path <- file.path(
    output_root, "xenium", "annotation", "04_unknown_qc_reviewed", "tables",
    "Xenium_Unknown_QC_routing_manifest.csv"
  )
  if (!file.exists(path)) stop("Unknown-QC routing manifest not found: ", path)
  list(
    path = path,
    manifest = read.csv(path, stringsAsFactors = FALSE, na.strings = character())
  )
}

resolve_unknown_qc_input_paths <- function(config = load_pipeline_config(),
                                           sample_ids = load_sample_manifest(config)$sample_id) {
  loaded <- load_unknown_qc_routing_manifest(config)
  inspection <- inspect_unknown_qc_routing_manifest(loaded$manifest, sample_ids)
  if (!inspection$valid) {
    stop("Unknown-QC routing manifest is invalid:\n- ", paste(inspection$errors, collapse = "\n- "))
  }
  paths <- as.character(loaded$manifest$mapping_input_path)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Unknown-QC routed input(s) are missing:\n- ", paste(missing, collapse = "\n- "))
  info <- file.info(paths)
  expected_size <- suppressWarnings(as.numeric(loaded$manifest$mapping_input_size))
  expected_mtime <- suppressWarnings(as.numeric(loaded$manifest$mapping_input_mtime))
  unchanged <- vapply(seq_along(paths), function(index) {
    isTRUE(all.equal(expected_size[[index]], as.numeric(info$size[[index]]))) &&
      isTRUE(all.equal(expected_mtime[[index]], as.numeric(info$mtime[[index]])))
  }, logical(1))
  if (any(!unchanged)) {
    stop(
      "Unknown-QC routed input provenance changed for sample(s): ",
      paste(sample_ids[!unchanged], collapse = ", "), ". Repeat the review/apply gate."
    )
  }
  stats::setNames(paths, sample_ids)
}
