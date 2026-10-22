# Atlas reporting helpers. No filtering, annotation, or analysis-object writes.
atlas_read_csv <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  if (anyDuplicated(names(x))) stop("Duplicate CSV column names: ", path)
  x
}

atlas_keys <- function(x, key, description) {
  if (!all(key %in% names(x))) stop(description, " lacks: ", paste(setdiff(key, names(x)), collapse = ", "))
  if (!nrow(x)) stop(description, " is empty.")
  for (k in key) if (anyNA(x[[k]]) || any(!nzchar(trimws(as.character(x[[k]]))))) {
    stop(description, " has missing/blank ", k)
  }
  if (anyDuplicated(x[key])) stop(description, " has duplicate keys: ", paste(key, collapse = ", "))
  invisible(TRUE)
}

atlas_align <- function(x, ids, description) {
  atlas_keys(x, "sample_id", description)
  if (!setequal(x$sample_id, ids)) stop(description, " does not exactly cover the active sample manifest.")
  x[match(ids, x$sample_id), , drop = FALSE]
}

atlas_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wx")
  on.exit(close(con))
  write.csv(x, con, row.names = FALSE, na = "")
}

atlas_provenance <- function(paths) {
  info <- file.info(paths)
  small <- file.exists(paths) & !is.na(info$size) & info$size <= 10 * 1024^2
  hashes <- rep(NA_character_, length(paths))
  hashes[small] <- unname(tools::md5sum(paths[small]))
  data.frame(path = paths, exists = file.exists(paths), bytes = info$size,
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    md5 = hashes, hash_policy = "MD5 only for files <=10 MiB; large objects use path/size/mtime", stringsAsFactors = FALSE)
}

atlas_check_metadata <- function(x, ids) {
  x <- atlas_align(x, ids, "Atlas metadata")
  required <- c("donor_id", "donor_verified", "PCW", "section_id", "section_plane",
    "slide_id", "run_id", "panel_id", "xoa_version", "segmentation_method")
  if (!all(required %in% names(x))) stop("Atlas metadata lacks: ", paste(setdiff(required, names(x)), collapse = ", "))
  if (any(!is.na(x$donor_verified) & !x$donor_verified %in% c("yes", "no", "provisional"))) stop("donor_verified must be yes/no/provisional or blank.")
  verified <- !is.na(x$donor_verified) & x$donor_verified == "yes"
  if (any(verified & (is.na(x$donor_id) | !nzchar(trimws(x$donor_id))))) stop("Verified donors require donor_id.")
  pcw <- suppressWarnings(as.numeric(x$PCW))
  if (any(!is.na(x$PCW) & (!is.finite(pcw) | pcw <= 0))) stop("PCW must be positive numeric or blank.")
  x$PCW <- pcw
  if (any(verified & !is.na(x$section_id))) {
    k <- x[verified & !is.na(x$section_id), c("donor_id", "section_id"), drop = FALSE]
    if (anyDuplicated(k)) stop("Duplicate donor/section pairs; resolve before reporting independent sections.")
  }
  grouped <- !is.na(x$donor_id)
  for (donor in unique(x$donor_id[grouped])) {
    ages <- unique(x$PCW[which(grouped & x$donor_id == donor & !is.na(x$PCW))])
    if (length(ages) > 1L) stop("Conflicting PCW values within verified donor: ", donor)
  }
  x
}

atlas_qc_flags <- function(metadata, thresholds) {
  required <- c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control")
  if (!all(required %in% names(metadata))) stop("Missing QC metadata.")
  if (any(!vapply(metadata[required], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(metadata[required])))) stop("Non-numeric/non-finite QC metadata.")
  limits <- c("nCount_lower", "nCount_upper", "nFeature_lower", "nFeature_upper",
    "area_lower", "area_upper", "max_control_percent")
  if (!all(limits %in% names(thresholds)) || any(!is.finite(unlist(thresholds[limits])))) stop("Missing/invalid stored thresholds.")
  if (thresholds$nCount_lower > thresholds$nCount_upper ||
      thresholds$nFeature_lower > thresholds$nFeature_upper || thresholds$area_lower > thresholds$area_upper) stop("Reversed stored QC limits.")
  with(metadata, data.frame(low_counts = nCount_Xenium < thresholds$nCount_lower,
    high_counts = nCount_Xenium > thresholds$nCount_upper,
    low_features = nFeature_Xenium < thresholds$nFeature_lower,
    high_features = nFeature_Xenium > thresholds$nFeature_upper,
    low_area = cell_area < thresholds$area_lower,
    high_area = cell_area > thresholds$area_upper,
    high_controls = percent_control > thresholds$max_control_percent))
}

atlas_quantiles <- function(metadata, sample_id, stage) {
  metrics <- intersect(c("nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control"), names(metadata))
  do.call(rbind, lapply(metrics, function(metric) {
    value <- metadata[[metric]]
    if (!is.numeric(value) || any(!is.finite(value))) stop("Invalid metric: ", metric)
    q <- stats::quantile(value, c(.05, .25, .5, .75, .95), names = FALSE)
    data.frame(sample_id = sample_id, stage = stage, metric = metric, n_cells = length(value),
      q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5])
  }))
}

atlas_metadata <- function(path, derive_qc = FALSE) {
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) stop("Expected Seurat RDS: ", path)
  cells <- colnames(object)
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) stop("Invalid cell IDs: ", path)
  metadata <- object[[]]
  if (!identical(rownames(metadata), cells)) stop("Cell/metadata alignment failed: ", path)
  if (derive_qc) {
    # Cropped RDS precedes qc_xenium(): percent_control is not yet stored there.
    # Reproduce the original calculation without changing stored production data.
    if (!"Xenium" %in% names(object@assays)) stop("Cropped object lacks Xenium assay.")
    pull_counts <- function(assay) {
      layers <- grep("^counts($|\\.)", SeuratObject::Layers(object[[assay]]), value = TRUE)
      if (length(layers) != 1L || layers != "counts") stop("Expected one unsplit counts layer for ", assay)
      # GetAssayData preserves dimensions for legacy one-feature control assays.
      counts <- SeuratObject::GetAssayData(object[[assay]], layer = "counts")
      if (!identical(colnames(counts), cells)) stop("Counts/cell alignment mismatch: ", assay)
      values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
      if (any(!is.finite(values) | values < 0 | values != floor(values))) stop("Invalid raw counts: ", assay)
      counts
    }
    counts <- pull_counts("Xenium")
    metadata$nCount_Xenium <- Matrix::colSums(counts)
    metadata$nFeature_Xenium <- Matrix::colSums(counts > 0)
    controls <- rep(0, length(cells))
    for (assay in intersect(c("BlankCodeword", "ControlCodeword", "ControlProbe", "GenomicControl"), names(object@assays))) {
      controls <- controls + Matrix::colSums(pull_counts(assay))
    }
    metadata$percent_control <- 100 * controls / (metadata$nCount_Xenium + controls)
    metadata$percent_control[is.na(metadata$percent_control)] <- 0
  }
  list(metadata = metadata, thresholds = object@misc$QC_thresholds,
    assays = names(object@assays))
}

atlas_validate_extracted <- function(retention, populations, failures, ids) {
  atlas_keys(retention, c("sample_id", "stage"), "Extracted retention")
  atlas_keys(populations, c("sample_id", "label"), "Extracted populations")
  for (x in list(retention, populations, failures)) {
    if (any(!x$sample_id %in% ids) || any(!is.finite(x$n_cells) | x$n_cells < 0 | x$n_cells != floor(x$n_cells))) stop("Invalid extracted cell counts/sample IDs.")
  }
  for (id in ids) {
    r <- retention[retention$sample_id == id, ]
    if (!setequal(r$stage, c("cropped", "initial_qc", "final"))) stop("Partial retention stages for ", id)
    n <- setNames(r$n_cells, r$stage)
    if (n["cropped"] <= 0 || n["initial_qc"] > n["cropped"] || n["final"] > n["initial_qc"]) stop("Nonmonotonic cell retention for ", id)
    if (sum(populations$n_cells[populations$sample_id == id]) != n["final"]) stop("Final population total mismatch: ", id)
    f <- failures[failures$sample_id == id & failures$kind == "exclusive_combination", ]
    if (!nrow(f) || anyDuplicated(f$reason) || sum(f$n_cells) != n["cropped"] ||
        sum(f$n_cells[f$reason == "pass"]) != n["initial_qc"]) stop("Failure-pattern reconciliation failed: ", id)
  }
  invisible(TRUE)
}
