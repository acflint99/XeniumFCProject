# Shared Stage-4 decision-manifest and Stage-7 input-routing helpers for the
# VZ and RL branches.

regional_qc_base_manifest_columns <- function() {
  c(
    "sample_id", "input_cells", "removed_cells", "remaining_cells",
    "clusters_removed"
  )
}

regional_qc_routing_manifest_columns <- function() {
  c(
    "source_path", "mapping_input_path", "materialization_mode",
    "source_size", "source_mtime", "mapping_input_size",
    "mapping_input_mtime"
  )
}

inspect_regional_qc_manifest <- function(manifest, sample_ids) {
  errors <- character()
  base_columns <- regional_qc_base_manifest_columns()
  routing_columns <- regional_qc_routing_manifest_columns()

  if (is.null(manifest)) {
    return(list(valid = FALSE, errors = "manifest is missing"))
  }
  missing_base <- setdiff(base_columns, names(manifest))
  if (length(missing_base)) {
    return(list(
      valid = FALSE,
      errors = paste0("missing base column(s): ", paste(missing_base, collapse = ", "))
    ))
  }
  if (nrow(manifest) != length(sample_ids)) {
    errors <- c(errors, "row count does not match config/samples.csv")
  } else if (!identical(as.character(manifest$sample_id), sample_ids)) {
    errors <- c(errors, "sample rows are not in config/samples.csv order")
  }

  removed_cells <- suppressWarnings(as.numeric(manifest$removed_cells))
  if (any(!is.finite(removed_cells)) || any(removed_cells < 0)) {
    errors <- c(errors, "removed_cells must contain finite nonnegative counts")
  }
  cluster_decisions <- as.character(manifest$clusters_removed)
  if (anyNA(cluster_decisions) || any(!nzchar(cluster_decisions))) {
    errors <- c(errors, "clusters_removed contains missing or blank values")
  }
  if (length(unique(cluster_decisions)) != 1L) {
    errors <- c(errors, "clusters_removed is inconsistent across samples")
  }

  routing_present <- routing_columns %in% names(manifest)
  if (any(routing_present) && !all(routing_present)) {
    errors <- c(errors, "routing columns are only partially present")
  }
  has_routing <- all(routing_present)

  input_cells <- suppressWarnings(as.numeric(manifest$input_cells))
  remaining_cells <- suppressWarnings(as.numeric(manifest$remaining_cells))
  if (has_routing) {
    modes <- as.character(manifest$materialization_mode)
    allowed_modes <- c("reuse_original", "filtered_copy")
    if (anyNA(modes) || any(!modes %in% allowed_modes)) {
      errors <- c(errors, "materialization_mode must be reuse_original or filtered_copy")
    }
    source_paths <- as.character(manifest$source_path)
    mapping_paths <- as.character(manifest$mapping_input_path)
    if (anyNA(source_paths) || any(!nzchar(source_paths)) ||
        anyNA(mapping_paths) || any(!nzchar(mapping_paths))) {
      errors <- c(errors, "routing paths contain missing or blank values")
    }
    if (anyDuplicated(source_paths) || anyDuplicated(mapping_paths)) {
      errors <- c(errors, "routing paths must be one-to-one with sample IDs")
    }
    provenance_columns <- c(
      "source_size", "source_mtime", "mapping_input_size", "mapping_input_mtime"
    )
    provenance_values <- lapply(
      provenance_columns,
      function(column) suppressWarnings(as.numeric(manifest[[column]]))
    )
    if (any(vapply(provenance_values, function(x) any(!is.finite(x)), logical(1)))) {
      errors <- c(errors, "routing file size/mtime provenance must be finite")
    }
    size_values <- provenance_values[c(1L, 3L)]
    if (any(vapply(size_values, function(x) any(x <= 0), logical(1)))) {
      errors <- c(errors, "routing file sizes must be positive")
    }

    reuse_rows <- !is.na(modes) & modes == "reuse_original"
    filtered_rows <- !is.na(modes) & modes == "filtered_copy"
    if (any(reuse_rows)) {
      if (any(removed_cells[reuse_rows] != 0, na.rm = TRUE) ||
          any(cluster_decisions[reuse_rows] != "none", na.rm = TRUE)) {
        errors <- c(errors, "reuse_original rows must record no removals")
      }
      if (any(source_paths[reuse_rows] != mapping_paths[reuse_rows], na.rm = TRUE)) {
        errors <- c(errors, "reuse_original rows must route to their source paths")
      }
      reuse_counts_present <- is.finite(input_cells[reuse_rows]) |
        is.finite(remaining_cells[reuse_rows])
      if (any(reuse_counts_present) &&
          any(
            !is.finite(input_cells[reuse_rows]) |
              !is.finite(remaining_cells[reuse_rows]) |
              input_cells[reuse_rows] != remaining_cells[reuse_rows]
          )) {
        errors <- c(errors, "reuse_original counts must both be missing or unchanged")
      }
    }
    if (any(filtered_rows)) {
      filtered_counts_valid <-
        is.finite(input_cells[filtered_rows]) &
        is.finite(remaining_cells[filtered_rows]) &
        input_cells[filtered_rows] >= 0 &
        remaining_cells[filtered_rows] >= 0 &
        input_cells[filtered_rows] - removed_cells[filtered_rows] ==
          remaining_cells[filtered_rows]
      if (any(!filtered_counts_valid)) {
        errors <- c(errors, "filtered_copy rows contain inconsistent cell counts")
      }
    }
  } else {
    legacy_counts_valid <-
      is.finite(input_cells) & is.finite(remaining_cells) &
      input_cells >= 0 & remaining_cells >= 0 &
      input_cells - removed_cells == remaining_cells
    if (any(!legacy_counts_valid)) {
      errors <- c(errors, "legacy manifest contains inconsistent cell counts")
    }
    modes <- rep("legacy_filtered_copy", nrow(manifest))
  }

  zero_removal <- !length(errors) &&
    identical(unique(cluster_decisions), "none") &&
    all(removed_cells == 0) &&
    if (has_routing) {
      all(modes == "reuse_original")
    } else {
      all(input_cells == remaining_cells)
    }

  list(
    valid = !length(errors),
    errors = unique(errors),
    has_routing = has_routing,
    legacy = !has_routing,
    zero_removal = zero_removal,
    removed_cells = removed_cells,
    materialization_mode = modes
  )
}

build_zero_removal_routing_manifest <- function(sample_ids, source_paths) {
  if (length(source_paths) != length(sample_ids)) {
    stop("Source-path count does not match the configured sample count.")
  }
  missing_paths <- source_paths[!file.exists(source_paths)]
  if (length(missing_paths)) {
    stop("Cannot route missing consensus input(s):\n- ", paste(missing_paths, collapse = "\n- "))
  }
  normalized_sources <- normalizePath(source_paths, mustWork = TRUE)
  source_info <- file.info(normalized_sources)
  if (anyNA(source_info$size) || anyNA(source_info$mtime)) {
    stop("Could not read file provenance for every consensus input.")
  }
  data.frame(
    sample_id = sample_ids,
    input_cells = NA_integer_,
    removed_cells = 0L,
    remaining_cells = NA_integer_,
    clusters_removed = "none",
    source_path = normalized_sources,
    mapping_input_path = normalized_sources,
    materialization_mode = "reuse_original",
    source_size = as.numeric(source_info$size),
    source_mtime = as.numeric(source_info$mtime),
    mapping_input_size = as.numeric(source_info$size),
    mapping_input_mtime = as.numeric(source_info$mtime),
    stringsAsFactors = FALSE
  )
}

resolve_regional_qc_mapping_inputs <- function(manifest,
                                               sample_ids,
                                               legacy_paths) {
  inspection <- inspect_regional_qc_manifest(manifest, sample_ids)
  if (!isTRUE(inspection$valid)) {
    stop(
      "QC decision manifest is invalid:\n- ",
      paste(inspection$errors, collapse = "\n- ")
    )
  }
  if (inspection$has_routing) {
    input_paths <- as.character(manifest$mapping_input_path)
    modes <- as.character(manifest$materialization_mode)
    expected_size <- suppressWarnings(as.numeric(manifest$mapping_input_size))
    expected_mtime <- suppressWarnings(as.numeric(manifest$mapping_input_mtime))
  } else {
    if (length(legacy_paths) != length(sample_ids)) {
      stop("Legacy input-path count does not match the configured sample count.")
    }
    input_paths <- legacy_paths
    modes <- rep("legacy_filtered_copy", length(sample_ids))
    expected_size <- expected_mtime <- rep(NA_real_, length(sample_ids))
  }

  missing_paths <- input_paths[!file.exists(input_paths)]
  if (length(missing_paths)) {
    stop("QC-routed mapping input(s) are missing:\n- ", paste(missing_paths, collapse = "\n- "))
  }
  if (inspection$has_routing) {
    current_info <- file.info(input_paths)
    current_size <- as.numeric(current_info$size)
    current_mtime <- as.numeric(current_info$mtime)
    changed <- !vapply(seq_along(input_paths), function(i) {
      isTRUE(all.equal(expected_size[[i]], current_size[[i]])) &&
        isTRUE(all.equal(expected_mtime[[i]], current_mtime[[i]]))
    }, logical(1))
    if (any(changed)) {
      stop(
        "QC-routed mapping input(s) changed after the Stage-4 decision:\n- ",
        paste(input_paths[changed], collapse = "\n- ")
      )
    }
  }

  list(
    input_paths = input_paths,
    materialization_mode = modes,
    inspection = inspection
  )
}
