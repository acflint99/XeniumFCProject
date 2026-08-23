# Read-only configuration utilities for the fetal cerebellum Xenium pipeline.

find_pipeline_config <- function() {
  candidates <- c(
    here::here("config", "config.yml"),
    here::here("scripts", "config", "config.yml")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Could not find config/config.yml or scripts/config/config.yml.")
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

load_pipeline_config <- function(config_path = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required to read config.yml.")
  }
  if (is.null(config_path)) config_path <- find_pipeline_config()
  config_path <- normalizePath(config_path, mustWork = TRUE)
  config <- yaml::read_yaml(config_path)
  attr(config, "config_path") <- config_path
  attr(config, "config_dir") <- dirname(config_path)
  config
}

resolve_config_path <- function(path, config) {
  candidates <- c(
    here::here(path),
    file.path(attr(config, "config_dir"), basename(path))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(normalizePath(existing[[1]], mustWork = TRUE))
  }
  candidates[[1]]
}

load_sample_manifest <- function(config = load_pipeline_config()) {
  path <- resolve_config_path(config$manifests$samples, config)
  if (!file.exists(path)) stop("Sample manifest not found: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
}

load_slide_manifest <- function(config = load_pipeline_config()) {
  path <- resolve_config_path(config$manifests$slides, config)
  if (!file.exists(path)) stop("Slide manifest not found: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
}

get_sample_record <- function(sample_id,
                              samples = NULL,
                              config = load_pipeline_config()) {
  if (is.null(samples)) samples <- load_sample_manifest(config)
  matches <- which(samples$sample_id == sample_id)
  if (length(matches) != 1L) {
    stop("Expected exactly one manifest row for sample '", sample_id,
         "'; found ", length(matches), ".")
  }
  samples[matches, , drop = FALSE]
}

resolve_sample_paths <- function(sample_record, config = load_pipeline_config()) {
  if (nrow(sample_record) != 1L) stop("sample_record must contain exactly one row.")
  stats_file <- sample_record$cell_stats_file[[1]]
  if (!nzchar(stats_file)) stats_file <- config$inputs$cerebellum_cell_stats_file

  input_dir <- here::here(
    config$project$data_dir,
    sample_record$input_directory[[1]]
  )

  list(
    sample_id = sample_record$sample_id[[1]],
    input_dir = input_dir,
    cell_stats_file = stats_file,
    cell_stats_path = file.path(input_dir, stats_file)
  )
}

validate_pipeline_config <- function(config = load_pipeline_config(), check_files = FALSE) {
  errors <- character()
  add_error <- function(...) errors <<- c(errors, paste0(...))

  required_sections <- c(
    "project", "manifests", "inputs", "qc", "initial_clustering",
    "label_transfer", "regional_subsets", "plotting"
  )
  missing_sections <- setdiff(required_sections, names(config))
  if (length(missing_sections)) {
    add_error("Missing config sections: ", paste(missing_sections, collapse = ", "))
  }

  samples <- load_sample_manifest(config)
  slides <- load_slide_manifest(config)

  required_sample_columns <- c(
    "sample_id", "input_directory", "input_layout", "cell_stats_file"
  )
  missing_sample_columns <- setdiff(required_sample_columns, names(samples))
  if (length(missing_sample_columns)) {
    add_error("Missing sample columns: ", paste(missing_sample_columns, collapse = ", "))
  }
  if (anyDuplicated(samples$sample_id)) add_error("Duplicate sample IDs found.")
  if (any(!nzchar(samples$sample_id))) add_error("Blank sample IDs found.")
  if (any(!samples$input_layout %in% c("single", "multi_sample_slide"))) {
    add_error("input_layout must be 'single' or 'multi_sample_slide'.")
  }

  multi <- samples[samples$input_layout == "multi_sample_slide", , drop = FALSE]
  if (any(!nzchar(multi$cell_stats_file))) {
    add_error("Every multi-sample-slide component must have a cell_stats_file.")
  }

  if (anyDuplicated(slides$input_directory)) add_error("Duplicate slide directories found.")
  slide_members <- strsplit(slides$biological_samples, "|", fixed = TRUE)
  declared_counts <- suppressWarnings(as.integer(slides$n_biological_samples))
  observed_counts <- lengths(slide_members)
  if (!identical(declared_counts, as.integer(observed_counts))) {
    add_error("One or more slide sample counts do not match biological_samples.")
  }
  mapped_samples <- unlist(slide_members, use.names = FALSE)
  if (!setequal(mapped_samples, multi$sample_id)) {
    add_error("Multi-sample manifest rows and slide mappings do not contain the same samples.")
  }
  slide_lookup <- setNames(slides$input_directory, slides$input_directory)
  if (any(!multi$input_directory %in% names(slide_lookup))) {
    add_error("A multi-sample component refers to an unknown slide directory.")
  }

  numeric_checks <- c(
    mad_multiplier = config$qc$mad_multiplier,
    min_transcripts = config$qc$min_transcripts,
    min_features = config$qc$min_features,
    min_cell_area = config$qc$min_cell_area,
    max_control_percent = config$qc$max_control_percent,
    clustering_dimensions = config$initial_clustering$dimensions,
    clustering_resolution = config$initial_clustering$resolution,
    transfer_dimensions = config$label_transfer$dimensions,
    regional_integration_dimensions = config$regional_subsets$integration_dimensions_default,
    regional_post_qc_dimensions = config$regional_subsets$post_qc_dimensions_observed,
    regional_neighbor_k = config$regional_subsets$neighbor_k_post_qc
  )
  if (any(!is.finite(numeric_checks)) || any(numeric_checks <= 0)) {
    add_error("QC, clustering, and transfer numeric settings must be finite and positive.")
  }

  vz_labels <- unlist(config$regional_subsets$vz_broad_labels, use.names = FALSE)
  rl_labels <- unlist(config$regional_subsets$rl_broad_labels, use.names = FALSE)
  if (!length(vz_labels) || any(!nzchar(vz_labels)) || anyDuplicated(vz_labels)) {
    add_error("regional_subsets.vz_broad_labels must contain unique, nonblank labels.")
  }
  if (!length(rl_labels) || any(!nzchar(rl_labels)) || anyDuplicated(rl_labels)) {
    add_error("regional_subsets.rl_broad_labels must contain unique, nonblank labels.")
  }

  if (isTRUE(check_files)) {
    for (i in seq_len(nrow(samples))) {
      paths <- resolve_sample_paths(samples[i, , drop = FALSE], config)
      if (!dir.exists(paths$input_dir)) {
        add_error("Missing input directory for ", paths$sample_id, ": ", paths$input_dir)
        next
      }
      if (!file.exists(paths$cell_stats_path)) {
        add_error("Missing cell-stat CSV for ", paths$sample_id, ": ", paths$cell_stats_path)
        next
      }
      header <- tryCatch(
        names(read.csv(
          paths$cell_stats_path,
          nrows = 0,
          comment.char = "#",
          check.names = FALSE,
          stringsAsFactors = FALSE
        )),
        error = function(e) character()
      )
      if (!"Cell ID" %in% header) {
        add_error("Cell-stat CSV lacks 'Cell ID' for ", paths$sample_id, ".")
      }
    }
  }

  if (length(errors)) {
    stop("Configuration validation failed:\n- ", paste(errors, collapse = "\n- "))
  }

  message(
    "Configuration valid: ", nrow(samples), " biological samples; ",
    nrow(slides), " multi-sample slide directories; file checks ",
    if (isTRUE(check_files)) "enabled." else "disabled."
  )
  invisible(list(config = config, samples = samples, slides = slides))
}
