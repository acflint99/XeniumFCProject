# Helpers for reviewing and applying VZ/RL annotations to whole-tissue
# sample-qualified Unknown populations.

unknown_relabel_key_columns <- function() {
  c("sample_id", "consensus_label", "VZ_subcluster", "RL_subcluster")
}

unknown_relabel_key <- function(data) {
  key_columns <- unknown_relabel_key_columns()
  missing_columns <- setdiff(key_columns, names(data))
  if (length(missing_columns)) {
    stop(
      "Unknown-relabel data lack key column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  values <- lapply(data[key_columns], function(x) trimws(as.character(x)))
  if (any(vapply(values, function(x) anyNA(x) || any(!nzchar(x)), logical(1)))) {
    stop("Unknown-relabel keys contain missing or blank values.")
  }
  do.call(paste, c(values, sep = "\u001f"))
}

validate_unknown_relabel_review <- function(review,
                                            sample_ids,
                                            allowed_final_labels = NULL,
                                            require_complete_labels = TRUE) {
  required_columns <- c(
    unknown_relabel_key_columns(), "final_label", "review_notes"
  )
  missing_columns <- setdiff(required_columns, names(review))
  if (length(missing_columns)) {
    stop(
      "Reviewed Unknown-relabel manifest lacks column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (!nrow(review)) {
    if (require_complete_labels) {
      stop(
        "Reviewed Unknown-relabel manifest has no decisions. Copy the generated ",
        "review template to the configured manifest and fill final_label."
      )
    }
    return(invisible(review))
  }

  key <- unknown_relabel_key(review)
  if (anyDuplicated(key)) {
    stop("Reviewed Unknown-relabel manifest contains duplicate annotation-pair keys.")
  }
  if (any(!review$sample_id %in% sample_ids)) {
    stop(
      "Reviewed Unknown-relabel manifest contains sample IDs absent from ",
      "config/samples.csv: ",
      paste(unique(review$sample_id[!review$sample_id %in% sample_ids]), collapse = ", ")
    )
  }
  if (any(!is_unknown_consensus_label(review$consensus_label))) {
    stop("Reviewed Unknown-relabel manifest contains a non-Unknown consensus label.")
  }

  final_labels <- trimws(as.character(review$final_label))
  if (require_complete_labels && (anyNA(final_labels) || any(!nzchar(final_labels)))) {
    stop("Every reviewed Unknown annotation pair must have a nonblank final_label.")
  }
  review_notes <- trimws(as.character(review$review_notes))
  if (require_complete_labels && (anyNA(review_notes) || any(!nzchar(review_notes)))) {
    stop("Every reviewed Unknown annotation pair must have nonblank review_notes.")
  }
  if (!is.null(allowed_final_labels)) {
    nonblank <- !is.na(final_labels) & nzchar(final_labels)
    invalid <- unique(final_labels[nonblank & !final_labels %in% allowed_final_labels])
    if (length(invalid)) {
      stop(
        "Reviewed Unknown final_label values are outside the approved broad/subcluster ",
        "taxonomy: ", paste(invalid, collapse = ", ")
      )
    }
  }
  invisible(review)
}

validate_unknown_comparison_provenance <- function(manifest,
                                                   sample_ids,
                                                   expected_inputs) {
  required_columns <- c(
    "sample_id", "input_type", "input_path", "input_size", "input_mtime"
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns)) {
    stop(
      "Unknown-comparison manifest lacks provenance column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  required_types <- c("consensus", "VZ", "RL")
  expected_rows <- expand.grid(
    sample_id = sample_ids,
    input_type = required_types,
    stringsAsFactors = FALSE
  )
  expected_rows <- expected_rows[order(
    match(expected_rows$sample_id, sample_ids),
    match(expected_rows$input_type, required_types)
  ), , drop = FALSE]
  observed_rows <- manifest[, c("sample_id", "input_type"), drop = FALSE]
  observed_rows <- observed_rows[order(
    match(observed_rows$sample_id, sample_ids),
    match(observed_rows$input_type, required_types)
  ), , drop = FALSE]
  rownames(expected_rows) <- rownames(observed_rows) <- NULL
  if (!identical(observed_rows, expected_rows)) {
    stop(
      "Unknown-comparison manifest must contain exactly one consensus, VZ, and RL ",
      "provenance row for every configured sample."
    )
  }

  expected_lookup <- expected_inputs
  expected_lookup$key <- paste(expected_lookup$sample_id, expected_lookup$input_type, sep = "\u001f")
  manifest_key <- paste(manifest$sample_id, manifest$input_type, sep = "\u001f")
  expected_index <- match(manifest_key, expected_lookup$key)
  if (anyNA(expected_index)) stop("Unknown-comparison manifest contains unexpected inputs.")
  expected_paths <- normalizePath(
    expected_lookup$input_path[expected_index], mustWork = TRUE
  )
  recorded_paths <- normalizePath(as.character(manifest$input_path), mustWork = FALSE)
  if (!identical(recorded_paths, expected_paths)) {
    stop("Unknown-comparison input paths no longer match the active workflow inputs.")
  }

  current_info <- file.info(expected_paths)
  current_size <- as.numeric(current_info$size)
  current_mtime <- as.numeric(current_info$mtime)
  recorded_size <- suppressWarnings(as.numeric(manifest$input_size))
  recorded_mtime <- suppressWarnings(as.numeric(manifest$input_mtime))
  unchanged <- vapply(seq_along(expected_paths), function(i) {
    isTRUE(all.equal(recorded_size[[i]], current_size[[i]])) &&
      isTRUE(all.equal(recorded_mtime[[i]], current_mtime[[i]]))
  }, logical(1))
  if (any(!unchanged)) {
    changed <- paste(
      manifest$sample_id[!unchanged], manifest$input_type[!unchanged], sep = "/"
    )
    stop(
      "VZ/RL comparison is stale for changed input(s): ",
      paste(changed, collapse = ", "),
      ". Rerun the comparison before applying reviewed labels."
    )
  }
  invisible(TRUE)
}
