# Descriptive diagnostics only: no identities, filtering, or input objects change.
rl_bg_genes <- c("EOMES", "SLC1A3", "BCAN", "SOX2", "NES", "PTPRZ1", "FABP7",
                 "SOX9", "TNC", "HOPX", "GDF10", "ATP1A2", "OTX2", "RSPO3",
                 "TRPC3", "CACNA2D1", "ATOH1", "MKI67", "TOP2A")
rl_bg_partners <- c("SLC1A3", "BCAN", "SOX2", "NES", "PTPRZ1", "FABP7")
rl_bg_target <- "BG Progenitors"
rl_bg_suffixes <- c("Cells.csv.gz", "Marker_Availability.csv", "Marker_Detection.csv",
                    "Coexpression.csv", "Spatial_Distances.csv", "Summary.csv",
                    "Spatial_Markers.pdf", "Provenance.csv")

rl_bg_extract_counts <- function(object, genes, cells) {
  assay <- object[["Xenium"]]
  layers <- SeuratObject::Layers(assay, search = "^counts($|\\.)")
  if (!length(layers)) stop("Xenium assay has no raw counts layers.")
  # Missing genes stay unavailable (NA), never interpreted as zero detection.
  result <- matrix(NA_real_, length(genes), length(cells), dimnames = list(genes, cells))
  seen <- integer(length(cells))
  for (layer in layers) {
    layer_cells <- intersect(cells, SeuratObject::Cells(assay, layer = layer))
    if (!length(layer_cells)) next
    available <- intersect(genes, SeuratObject::Features(assay, layer = layer))
    index <- match(layer_cells, cells)
    seen[index] <- seen[index] + 1L
    if (!length(available)) next
    block <- SeuratObject::LayerData(assay, layer = layer,
                                    features = available, cells = layer_cells)
    block <- as.matrix(block)
    if (!setequal(colnames(block), layer_cells) ||
        !setequal(rownames(block), available) || any(!is.finite(block)) ||
        any(block < 0) || any(block != floor(block))) {
      stop("Invalid raw counts or count-layer cell/gene alignment: ", layer)
    }
    result[available, colnames(block)] <- block[available, , drop = FALSE]
  }
  if (any(seen != 1L)) stop("Every RL cell must occur in exactly one counts layer.")
  result
}

rl_bg_pair_summary <- function(frame, partners = rl_bg_partners) {
  rows <- list()
  for (label in unique(frame$RL_subcluster)) {
    z <- frame[frame$RL_subcluster == label, , drop = FALSE]
    for (partner in partners) for (threshold in c(1L, 2L)) {
      x <- z[["EOMES"]]; y <- z[[partner]]
      measured <- !is.na(x) & !is.na(y)
      a <- x[measured] >= threshold; b <- y[measured] >= threshold
      n <- sum(measured)
      rows[[length(rows) + 1L]] <- data.frame(
        RL_subcluster = label, gene1 = "EOMES", gene2 = partner,
        minimum_transcripts = threshold, cells = nrow(z), measured_cells = n,
        status = if (n == nrow(z)) "available" else if (n == 0L) "unavailable" else "partial_availability",
        both_positive = if (n) sum(a & b) else NA_integer_,
        EOMES_only = if (n) sum(a & !b) else NA_integer_,
        partner_only = if (n) sum(!a & b) else NA_integer_,
        neither = if (n) sum(!a & !b) else NA_integer_,
        both_fraction = if (n) mean(a & b) else NA_real_,
        partner_fraction_among_EOMES_positive = if (sum(a)) mean(b[a]) else NA_real_,
        EOMES_fraction_among_partner_positive = if (sum(b)) mean(a[b]) else NA_real_
      )
    }
  }
  do.call(rbind, rows)
}

rl_bg_spatial_distances <- function(frame) {
  target <- frame[frame$RL_subcluster == rl_bg_target, , drop = FALSE]
  groups <- list(RL_SVZ = "RL SVZ", UBCs = c("Immature UBCs", "Mature UBCs"),
                 Granule = c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs",
                             "Migrating GCs", "Mature GCs"))
  # Native centroid coordinates, not anatomical boundaries or inferred lineage.
  out <- lapply(names(groups), function(group) {
    comparison <- frame[frame$RL_subcluster %in% groups[[group]], , drop = FALSE]
    if (!nrow(target)) return(NULL)
    distances <- if (nrow(comparison)) {
      as.numeric(RANN::nn2(as.matrix(comparison[, c("x", "y")]),
                          as.matrix(target[, c("x", "y")]), k = 1)$nn.dists)
    } else rep(NA_real_, nrow(target))
    data.frame(cell_id = target$cell_id, comparison = group,
               comparison_cells = nrow(comparison), nearest_distance = distances,
               status = if (nrow(comparison)) "available" else "comparison_absent")
  })
  out <- do.call(rbind, out)
  if (is.null(out)) out <- data.frame(cell_id = character(), comparison = character(),
    comparison_cells = integer(), nearest_distance = numeric(), status = character())
  out
}

rl_bg_analyze <- function(input, sample_id, task_id, output_paths, config_path, manifest_path) {
  for (package in c("Seurat", "Matrix", "ggplot2", "RANN", "scattermore")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing dependency: ", package, ". Check the project renv environment; no packages installed.")
    }
  }
  before <- file.info(input)
  message(Sys.time(), " Loading ", sample_id)
  obj <- readRDS(input)
  cells <- colnames(obj); meta <- obj[[]]
  required <- c("RL_subcluster", "PCW", "nCount_Xenium", "nFeature_Xenium")
  if (anyDuplicated(cells) || !identical(rownames(meta), cells) ||
      !all(required %in% names(meta)) || !"Xenium" %in% SeuratObject::Assays(obj)) {
    stop("Invalid cell IDs, required metadata, or Xenium assay: ", sample_id)
  }
  if ("orig.ident" %in% names(meta) &&
      (anyNA(meta$orig.ident) || !all(as.character(meta$orig.ident) == sample_id))) {
    stop("Object sample metadata does not match manifest task: ", sample_id)
  }
  labels <- as.character(meta$RL_subcluster)
  if (any(!is.na(labels) & !nzchar(trimws(labels)))) stop("Blank RL labels found.")
  if (anyNA(meta[, required[-1], drop = FALSE]) ||
      any(!is.finite(meta$nCount_Xenium)) || any(meta$nCount_Xenium <= 0) ||
      any(!is.finite(meta$nFeature_Xenium)) || any(meta$nFeature_Xenium <= 0) ||
      length(unique(as.character(meta$PCW))) != 1L || any(!nzchar(as.character(meta$PCW)))) {
    stop("Missing/invalid sample age or QC metadata.")
  }
  coords <- Seurat::GetTissueCoordinates(obj, type = "Xenium")
  if (!all(c("cell", "x", "y") %in% names(coords)) ||
      anyDuplicated(coords$cell) || !setequal(as.character(coords$cell), cells)) {
    stop("Spatial centroids must cover every object cell exactly once.")
  }
  coords <- coords[match(cells, coords$cell), , drop = FALSE]
  if (!is.numeric(coords$x) || !is.numeric(coords$y) ||
      any(!is.finite(coords$x)) || any(!is.finite(coords$y))) stop("Invalid spatial coordinates.")
  rl <- which(!is.na(labels))
  if (!length(rl)) stop("No mapped RL cells: ", sample_id)
  background <- data.frame(x = coords$y, y = coords$x)
  counts <- rl_bg_extract_counts(obj, rl_bg_genes, cells[rl])
  frame <- data.frame(cell_id = cells[rl], sample_id = sample_id,
    PCW = as.character(meta$PCW[rl]), RL_subcluster = labels[rl],
    x = coords$y[rl], y = coords$x[rl], nCount_Xenium = meta$nCount_Xenium[rl],
    nFeature_Xenium = meta$nFeature_Xenium[rl], stringsAsFactors = FALSE)
  frame <- cbind(frame, as.data.frame(t(counts)))
  availability <- data.frame(gene = rl_bg_genes,
    measured_cells = rowSums(!is.na(counts)), total_rl_cells = length(rl))
  availability$status <- ifelse(availability$measured_cells == 0L, "unavailable",
    ifelse(availability$measured_cells == length(rl), "available", "partial_availability"))
  rm(obj, meta, coords, counts); invisible(gc())
  coexpression <- rl_bg_pair_summary(frame)
  detection <- do.call(rbind, lapply(unique(frame$RL_subcluster), function(label) {
    z <- frame[frame$RL_subcluster == label, , drop = FALSE]
    do.call(rbind, lapply(rl_bg_genes, function(gene) {
      v <- z[[gene]]; n <- sum(!is.na(v))
      data.frame(RL_subcluster = label, gene = gene, cells = nrow(z), measured_cells = n,
        detected_ge1 = if (n) sum(v >= 1, na.rm = TRUE) else NA_integer_,
        fraction_ge1 = if (n) mean(v >= 1, na.rm = TRUE) else NA_real_,
        fraction_ge2 = if (n) mean(v >= 2, na.rm = TRUE) else NA_real_,
        mean_raw_counts = if (n) mean(v, na.rm = TRUE) else NA_real_)
    }))
  }))
  distances <- rl_bg_spatial_distances(frame)
  target <- frame[frame$RL_subcluster == rl_bg_target, , drop = FALSE]
  summary <- data.frame(sample_id = sample_id, task_id = task_id, PCW = frame$PCW[1],
    object_cells = nrow(background), mapped_rl_cells = nrow(frame), BG_cells = nrow(target),
    status = if (nrow(target)) "BG_present" else "BG_absent",
    coordinate_units = "native Xenium centroid units; not calibrated here")
  for (group in c("RL_SVZ", "UBCs", "Granule")) {
    v <- distances$nearest_distance[distances$comparison == group]
    summary[[paste0("median_distance_to_", group)]] <- if (any(is.finite(v))) median(v, na.rm = TRUE) else NA_real_
  }
  add_sample <- function(z) {
    z$sample_id <- rep(sample_id, nrow(z)); z$task_id <- rep(task_id, nrow(z)); z
  }
  dir.create(dirname(output_paths[1]), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(output_paths[["Cells.csv.gz"]], "wt")
  tryCatch(write.csv(frame, connection, row.names = FALSE), finally = close(connection))
  for (name in c("Marker_Availability", "Marker_Detection", "Coexpression", "Spatial_Distances", "Summary")) {
    z <- switch(name, Marker_Availability = availability, Marker_Detection = detection,
      Coexpression = coexpression, Spatial_Distances = distances, Summary = summary)
    write.csv(add_sample(z), output_paths[[paste0(name, ".csv")]], row.names = FALSE)
  }
  # Plot all cells without random sampling; same orientation as production reports.
  base <- function(title) {
    ggplot2::ggplot() + scattermore::geom_scattermore(data = background,
      ggplot2::aes(x, y), colour = "grey85", pointsize = 1, pixels = c(1600, 1600)) +
      ggplot2::coord_fixed() + ggplot2::theme_void() +
      ggplot2::labs(title = paste(sample_id, title), subtitle = "Centroids; display x = source y, display y = source x")
  }
  points <- function(z, mapping) scattermore::geom_scattermore(data = z,
    mapping = mapping, pointsize = 4, pixels = c(1600, 1600))
  grDevices::pdf(output_paths[["Spatial_Markers.pdf"]], width = 11, height = 8.5)
  tryCatch({
    print(base("BG Progenitors only") +
      points(target, ggplot2::aes(x, y)) + ggplot2::labs(caption = paste(nrow(target), "BG cells; black foreground")))
    context <- frame
    context$group <- ifelse(context$RL_subcluster == rl_bg_target, rl_bg_target,
      ifelse(context$RL_subcluster == "RL SVZ", "RL SVZ",
        ifelse(context$RL_subcluster %in% c("Immature UBCs", "Mature UBCs"), "UBCs", "Other RL")))
    context <- context[order(context$group == rl_bg_target), , drop = FALSE]
    print(base("RL context") + points(context, ggplot2::aes(x, y, colour = group)) +
      ggplot2::scale_colour_manual(values = c("BG Progenitors" = "#8A2BE2", "RL SVZ" = "#D60050",
                                            "UBCs" = "#FF8C00", "Other RL" = "#70B59A")))
    for (partner in rl_bg_partners) {
      if (!nrow(target)) {
        print(base(paste("BG: EOMES +", partner)) +
          ggplot2::labs(caption = "No BG Progenitors cells in this sample; no co-expression estimate."))
        next
      }
      z <- target
      a <- z$EOMES >= 1; b <- z[[partner]] >= 1
      z$detection <- ifelse(is.na(a) | is.na(b), "Unavailable",
        ifelse(a & b, "Both", ifelse(a, "EOMES only", ifelse(b, "Partner only", "Neither"))))
      print(base(paste("BG: EOMES +", partner)) + points(z, ggplot2::aes(x, y, colour = detection)) +
        ggplot2::scale_colour_manual(values = c("Both" = "#D60050", "EOMES only" = "#0072B2",
          "Partner only" = "#E69F00", "Neither" = "grey45", "Unavailable" = "black")) +
        ggplot2::labs(caption = "Positive = at least 1 raw transcript; descriptive detection, not an identity threshold"))
    }
    print(ggplot2::ggplot(detection, ggplot2::aes(gene, RL_subcluster)) +
      ggplot2::geom_point(ggplot2::aes(size = fraction_ge1, colour = log1p(mean_raw_counts)), na.rm = TRUE) +
      ggplot2::scale_colour_viridis_c() + ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
      ggplot2::labs(title = paste(sample_id, "raw marker detection across mapped RL identities"),
        size = "Fraction >=1", colour = "log1p mean\nraw counts", x = NULL, y = NULL,
        caption = "Missing genes are unavailable; mean counts are descriptive and not library-size normalized."))
  }, finally = grDevices::dev.off())
  after <- file.info(input)
  if (!identical(before$size, after$size) || !identical(before$mtime, after$mtime)) {
    stop("Input changed during diagnostics; outputs are incomplete/unverified.")
  }
  provenance <- data.frame(sample_id = sample_id, task_id = task_id,
    input_path = normalizePath(input), input_size = before$size,
    input_mtime = as.numeric(before$mtime), target_label = rl_bg_target,
    genes = paste(rl_bg_genes, collapse = "|"), detection_thresholds = "raw counts >=1 and >=2",
    coordinate_transform = "display x=source y; display y=source x; Euclidean distance unchanged",
    distance_comparisons = "BG to RL SVZ; BG to Immature+Mature UBCs; BG to granule identities",
    interpretation = "descriptive; samples are not necessarily independent donors; distances depend on cell density and section geometry",
    config_path = normalizePath(config_path), config_md5 = unname(tools::md5sum(config_path)),
    manifest_path = normalizePath(manifest_path), manifest_md5 = unname(tools::md5sum(manifest_path)),
    helper_md5 = unname(tools::md5sum(here::here("scripts", "R", "rl_bg_spatial_diagnostics.R"))),
    R_version = as.character(getRversion()),
    package_versions = paste(vapply(c("Seurat", "SeuratObject", "RANN", "scattermore"),
      function(p) paste(p, as.character(utils::packageVersion(p))), character(1)), collapse = "; "),
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE)
  if (any(!file.exists(output_paths[names(output_paths) != "Provenance.csv"]))) stop("Incomplete diagnostic outputs.")
  # Written last: a completion record, not a modified analysis manifest.
  write.csv(provenance, output_paths[["Provenance.csv"]], row.names = FALSE)
  message("Completed ", sample_id, ": ", nrow(target), " BG cells; input unchanged.")
  invisible(summary)
}
