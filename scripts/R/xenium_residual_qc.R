# Descriptive residual QC. These functions never filter cells or assign labels.

residual_qc_numeric <- function(x, name, allow_missing = FALSE) {
  value <- suppressWarnings(as.numeric(as.character(x)))
  invalid <- is.infinite(value) | (!allow_missing & !is.finite(value)) |
    (!is.na(x) & is.na(value)) | (!is.na(value) & value < 0)
  if (any(invalid)) stop("Invalid or non-finite ", name, " values.")
  value
}

residual_qc_read_morphology <- function(path, raw_ids, retained_area) {
  required <- c("cell_id", "cell_area", "nucleus_area", "nucleus_count")
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  if (!all(required %in% header)) stop("Missing morphology columns in ", path)
  classes <- rep("NULL", length(header))
  classes[match(required, header)] <- c("character", "numeric", "numeric", "numeric")
  raw <- utils::read.csv(path, colClasses = classes, check.names = FALSE,
                         stringsAsFactors = FALSE)
  if (anyNA(raw$cell_id) || any(!nzchar(raw$cell_id)) || anyDuplicated(raw$cell_id)) {
    stop("Missing or duplicate raw cell IDs in ", path)
  }
  if (anyDuplicated(raw_ids)) stop("Duplicate residual raw cell IDs for ", path)
  index <- match(raw_ids, raw$cell_id)
  if (anyNA(index)) stop("Residual cell IDs absent from raw morphology: ", path)
  raw <- raw[index, , drop = FALSE]
  if (!identical(raw_ids, raw$cell_id)) stop("Morphology join changed cell order.")
  raw$cell_area <- residual_qc_numeric(raw$cell_area, "cell_area")
  raw$nucleus_area <- residual_qc_numeric(raw$nucleus_area, "nucleus_area", TRUE)
  raw$nucleus_count <- residual_qc_numeric(raw$nucleus_count, "nucleus_count")
  if (any(raw$cell_area <= 0) ||
      any(abs(raw$nucleus_count - round(raw$nucleus_count)) > 1e-8)) {
    stop("Morphology requires positive cell area and integer nucleus counts.")
  }
  # Same rounding tolerance as the production Unknown-QC audit.
  difference <- abs(retained_area - raw$cell_area)
  if (any(!is.finite(difference)) || any(difference > 0.00501)) {
    stop("Retained/raw cell areas differ beyond 0.00501 square micrometers: ", path)
  }
  attr(raw, "max_cell_area_difference") <- max(difference)
  raw
}

residual_qc_summarize <- function(cells, metrics, min_cells_per_group = 10L) {
  required <- c("cell_id", "sample_id", "cluster", "nucleus_area",
                "nucleus_count", "cell_area", metrics)
  if (!all(required %in% names(cells)) || !nrow(cells)) {
    stop("QC cells are empty or lack required columns.")
  }
  if (anyNA(cells$cell_id) || any(!nzchar(cells$cell_id)) ||
      anyDuplicated(cells$cell_id)) stop("Missing or duplicate QC cell IDs.")
  for (column in c("sample_id", "cluster")) {
    if (anyNA(cells[[column]]) || any(!nzchar(cells[[column]]))) {
      stop("Blank QC ", column, ".")
    }
  }
  for (metric in metrics) {
    if (!is.numeric(cells[[metric]]) || any(is.infinite(cells[[metric]])) ||
        any(cells[[metric]] < 0, na.rm = TRUE)) stop("Invalid QC metric: ", metric)
  }
  safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) stats::median(x) else NA_real_
  }
  describe <- function(x) {
    finite <- x[is.finite(x)]
    q <- if (length(finite)) stats::quantile(finite, c(.1, .25, .5, .75, .9),
                                           names = FALSE) else rep(NA_real_, 5L)
    data.frame(n_cells = length(x), n_measured = length(finite),
               n_missing = length(x) - length(finite),
               p10 = q[1], p25 = q[2], median = q[3], p75 = q[4], p90 = q[5])
  }
  clusters <- unique(cells$cluster)
  number <- suppressWarnings(as.numeric(clusters))
  clusters <- clusters[order(is.na(number), number, clusters)]
  pooled <- list(); paired <- list()
  for (cluster in clusters) {
    target <- cells$cluster == cluster
    for (metric in metrics) {
      pooled[[length(pooled) + 1L]] <- cbind(
        data.frame(cluster = cluster, metric = metric), describe(cells[[metric]][target])
      )
      for (sample_id in unique(cells$sample_id[target])) {
        in_sample <- cells$sample_id == sample_id
        x <- cells[[metric]][in_sample & target]
        y <- cells[[metric]][in_sample & !target]
        a <- describe(x); b <- describe(y)
        sample_median <- safe_median(cells[[metric]][in_sample])
        eligible <- a$n_measured >= min_cells_per_group &&
          b$n_measured >= min_cells_per_group
        paired[[length(paired) + 1L]] <- cbind(
          data.frame(sample_id = sample_id, cluster = cluster, metric = metric), a,
          data.frame(fraction_of_sample = length(x) / sum(in_sample),
                     fraction_of_cluster = length(x) / sum(target),
                     sample_median = sample_median,
                     median_ratio_to_sample = if (is.finite(sample_median) && sample_median > 0)
                       a$median / sample_median else NA_real_,
                     background_cells = length(y), background_measured = b$n_measured,
                     background_median = b$median,
                     median_difference = a$median - b$median,
                     median_ratio = if (is.finite(b$median) && b$median > 0)
                       a$median / b$median else NA_real_,
                     eligible = eligible)
        )
      }
    }
  }
  paired <- do.call(rbind, paired)
  balanced <- list()
  for (cluster in clusters) for (metric in metrics) {
    rows <- paired[paired$cluster == cluster & paired$metric == metric, , drop = FALSE]
    used <- rows[rows$eligible, , drop = FALSE]
    balanced[[length(balanced) + 1L]] <- data.frame(
      cluster = cluster, metric = metric,
      samples_with_cluster = nrow(rows), eligible_samples = nrow(used),
      min_measured_target_and_background = min_cells_per_group,
      eligible_sample_ids = paste(used$sample_id, collapse = "|"),
      median_sample_target_median = safe_median(used$median),
      median_sample_background_median = safe_median(used$background_median),
      median_sample_difference = safe_median(used$median_difference),
      median_sample_ratio = safe_median(used$median_ratio),
      samples_with_lower_target_median = sum(used$median_difference < 0),
      samples_with_equal_target_median = sum(used$median_difference == 0),
      samples_with_higher_target_median = sum(used$median_difference > 0)
    )
  }
  segmentation <- do.call(rbind, lapply(clusters, function(cluster) {
    do.call(rbind, lapply(unique(cells$sample_id[cells$cluster == cluster]), function(id) {
      d <- cells[cells$cluster == cluster & cells$sample_id == id, , drop = FALSE]
      data.frame(sample_id = id, cluster = cluster, n_cells = nrow(d),
        fraction_nucleus_count_zero = mean(d$nucleus_count == 0),
        fraction_nucleus_count_not_one = mean(d$nucleus_count != 1),
        fraction_nucleus_area_missing = mean(!is.finite(d$nucleus_area)),
        fraction_nucleus_area_zero = mean(d$nucleus_area == 0, na.rm = TRUE),
        fraction_missing_area_with_nucleus = mean(
          !is.finite(d$nucleus_area) & d$nucleus_count > 0),
        cell_area_unique_values = length(unique(d$cell_area)))
    }))
  }))
  # An all-missing area has no defined zero fraction, rather than NaN.
  segmentation$fraction_nucleus_area_zero[
    is.nan(segmentation$fraction_nucleus_area_zero)] <- NA_real_
  list(pooled = do.call(rbind, pooled), paired = paired,
       balanced = do.call(rbind, balanced), segmentation = segmentation)
}
