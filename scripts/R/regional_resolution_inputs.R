# Resolve current post-QC analysis sources without loading RDS files.
# All branches refit after removals or reuse when no cells are removed.
regional_resolution_paths <- function(config, region) {
  root <- here::here(config$project$outputs_dir, "xenium", tolower(region), "05_post_qc")
  p <- list(input = file.path(root, "rds", paste0("Xenium_", region, "_postQC.rds")),
            provenance = file.path(root, "tables", paste0("Xenium_", region, "_postQC_manifest.csv")),
            parent = file.path(root, "diagnostics", "resolution_stability"))
  p$inputs <- c(p$input, p$provenance)
  if (region == "Residual") {
    q <- residual_qc_removal_paths(config)
    post <- residual_post_qc_paths(config)
    p$coverage <- post$coverage
    p$inputs <- c(p$input, p$provenance, p$coverage, q$source,
      if (file.exists(q$decision)) unlist(q[c("decision", "filtered", "removed", "manifest", "provenance")]))
  }
  p
}

regional_resolution_manifest <- function(paths, config, region) {
  manifest <- read.csv(paths$provenance, stringsAsFactors = FALSE)
  if (region != "Residual") return(list(manifest = manifest, state = NULL, expected_counts = NULL))
  state <- load_residual_analysis_state(config)
  if (!identical(normalizePath(paths$input, mustWork = TRUE), normalizePath(state$input, mustWork = TRUE))) {
    stop("Residual post-QC input changed after resolving the diagnostic source.")
  }
  list(manifest = state$post_qc_manifest, state = state, expected_counts = state$expected_counts)
}

regional_resolution_check_coverage <- function(sample_id, samples, expected_counts = NULL) {
  if (anyNA(sample_id) || any(!nzchar(sample_id)) || any(!sample_id %in% samples)) {
    stop("Missing or unexpected sample IDs in the resolution-review input.")
  }
  expected_samples <- samples
  if (!is.null(expected_counts)) {
    if (!identical(names(expected_counts), samples) || anyNA(expected_counts) ||
        !identical(as.integer(table(factor(sample_id, levels = samples))), unname(expected_counts))) {
      stop("Residual per-sample counts differ from the current coverage/QC manifest.")
    }
    expected_samples <- samples[expected_counts > 0L]
  }
  if (!setequal(unique(sample_id), expected_samples)) stop("Object coverage differs from the expected nonempty samples.")
  invisible(TRUE)
}
