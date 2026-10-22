# Residual Stage-7 paths, metadata checks, plotting, and lightweight audit.

residual_report_paths <- function(root, sample_id) {
  prefix <- file.path(root, sample_id, paste0(sample_id, "_Residual"))
  setNames(paste0(prefix, c("_Global_Spatial_pt0.6.tif", "_Global_Spatial_pt0.85.tif",
                           "_Faceted_Clusters_pt0.6.tif", "_Spatial_Report.pdf",
                           "_Markers_DotPlot.tif", "_Markers_DotPlot.pdf", "_Marker_Availability.csv",
                           "_cluster_counts.csv", "_Report_Status.csv")),
           c("global06", "global085", "facets", "pdf", "dot_tiff", "dot_pdf", "markers", "counts", "status"))
}

residual_shared_marker_panels <- function() {
  # Exact production panel contracts, not raw exploratory review panels.
  panels <- list(
    "eCN" = rl_other_markers[["eCN"]],
    "BG" = vz_compact_lineage_markers[["Glia"]],
    "Immature UBCs" = rl_production_marker_panel("Immature UBCs")[["UBC"]],
    # Preserve the residual UBC display genes after the separate VZ cluster
    # formerly bearing this label was reassigned to the granule lineage.
    "Mature UBCs" = c("CBLN1", "RELN"),
    "GABA Progenitors" = rl_other_markers[["GABA Progenitors"]],
    "Microglia" = rl_other_markers[["Microglia"]],
    "OPCs" = vz_compact_lineage_markers[["OPCs"]]
  )
  attr(panels, "marker_sources") <- c(
    "eCN" = "RL Other: eCN", "BG" = "VZ production: Glia",
    "Immature UBCs" = "RL production: UBC", "Mature UBCs" = "VZ production: UBCs",
    "GABA Progenitors" = "RL Other: GABA Progenitors", "Microglia" = "RL Other: Microglia",
    "OPCs" = "VZ production: OPCs"
  )
  panels
}

residual_report_marker_sets <- function(review, shared_panels, allow_repeated_genes = TRUE) {
  if (!all(c("final_label", "positive_marker_support", "approved") %in% names(review)) ||
      !nrow(review) || anyNA(as.logical(review$approved)) || any(!as.logical(review$approved))) {
    stop("Residual report markers require the approved label-review table.")
  }
  sets <- lapply(split(as.character(review$positive_marker_support), review$final_label), function(values) {
    genes <- unique(trimws(unlist(strsplit(values, ";", fixed = TRUE))))
    genes <- genes[nzchar(genes)]
    if (!length(genes) || any(!grepl("^[A-Za-z0-9_.-]+$", genes))) stop("Invalid reviewed residual marker list.")
    head(genes, 3L)
  })
  sources <- setNames(rep("Approved residual review: positive_marker_support", length(sets)), names(sets))
  for (label in names(sets)) {
    if (label %in% names(shared_panels)) {
      sets[[label]] <- shared_panels[[label]]
      sources[[label]] <- attr(shared_panels, "marker_sources")[[label]]
    } else if (label %in% union(vz_subcluster_order, rl_subcluster_order)) {
      stop("Shared residual identity requires an explicit VZ/RL production marker mapping: ", label)
    }
    if (!length(sets[[label]]) || anyNA(sets[[label]]) || any(!nzchar(sets[[label]])) ||
        anyDuplicated(sets[[label]])) stop("Invalid residual marker panel for ", label)
  }
  if (!allow_repeated_genes && anyDuplicated(unlist(sets, use.names = FALSE))) {
    stop("Residual marker groups contain repeated genes; define nonduplicated shared groups before plotting.")
  }
  attr(sets, "marker_sources") <- sources
  sets
}

residual_marker_panel_signature <- function(marker_sets) {
  sources <- attr(marker_sets, "marker_sources")
  lines <- vapply(names(marker_sets), function(label) {
    paste(label, sources[[label]], paste(marker_sets[[label]], collapse = ";"), sep = "\t")
  }, character(1))
  path <- tempfile("residual_marker_panel_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(lines, path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

check_residual_plot_metadata <- function(metadata, expected, approved_labels) {
  required <- c("residual_subcluster", "residual_qc_status", "PCW")
  if (!all(required %in% names(metadata)) || !nrow(metadata) ||
      is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Residual plotting requires nonempty, unique cell metadata and: ", paste(required, collapse = ", "))
  }
  labels <- as.character(metadata$residual_subcluster)
  status <- as.character(metadata$residual_qc_status)
  if (any(!is.na(labels) & !nzchar(trimws(labels))) ||
      anyNA(status) || any(!status %in% c("retained", "removed", "not_residual")) ||
      !identical(!is.na(labels), status == "retained")) {
    stop("Residual labels and approved QC statuses are inconsistent.")
  }
  if (any(!labels[!is.na(labels)] %in% approved_labels)) stop("Mapped residual labels differ from the approved taxonomy.")
  if (nrow(metadata) != expected$total_cells ||
      sum(status == "retained") != expected$mapped_residual_cells ||
      sum(status == "removed") != expected$removed_residual_qc_cells) {
    stop("Sample cell counts differ from the residual mapping manifest.")
  }
  pcw <- unique(trimws(as.character(metadata$PCW)))
  if (length(pcw) != 1L || anyNA(pcw) || !nzchar(pcw) ||
      !identical(pcw, trimws(as.character(expected$PCW)))) stop("Sample PCW differs from its mapping manifest.")
  list(labels = labels, retained = status == "retained", removed = status == "removed",
       PCW = pcw, residual_cells = sum(status == "retained"))
}

generate_residual_spatial_report <- function(sample_id, input_path, paths, expected, approved_labels,
                                             marker_sets, review_md5) {
  # Grid unit measurements must not open an unsolicited Rplots.pdf in the project.
  grDevices::pdf(file = NULL)
  measurement_device <- grDevices::dev.cur()
  on.exit({
    if (measurement_device %in% grDevices::dev.list()) grDevices::dev.off(measurement_device)
  }, add = TRUE)
  input_info <- file.info(input_path)
  object <- readRDS(input_path)
  checked <- check_residual_plot_metadata(object[[]], expected, approved_labels)
  present <- intersect(master_subcluster_order, unique(checked$labels[checked$retained]))
  if (length(setdiff(checked$labels[checked$retained], present))) stop("Residual plot order is incomplete.")
  palette <- subcluster_palette[present]
  if (anyNA(palette)) stop("Residual identities lack palette colors.")
  spatial <- prepare_spatial_plot_data(object, "residual_subcluster", present, sample_id)
  if (!setequal(spatial$cell, colnames(object)) ||
      any(!is.finite(spatial$plot_x)) || any(!is.finite(spatial$plot_y))) {
    stop(sample_id, " spatial coordinates do not cover every mapped object cell exactly once with finite values.")
  }
  # Only the already-approved residual QC exclusions are hidden from the plot.
  spatial <- spatial[!checked$removed[match(spatial$cell, rownames(object[[]]))], , drop = FALSE]
  if (!nrow(spatial)) stop(sample_id, " has no retained spatial context cells.")
  global <- function(size, rasterize = FALSE) {
    if (checked$residual_cells > 0L) {
      return(build_layered_spatial_plot(spatial, palette,
        paste(sample_id, "Residual subclusters"), targets = present,
        foreground_size = size, rasterize = rasterize))
    }
    ggplot2::ggplot() + spatial_point_layer(spatial,
      ggplot2::aes(x = plot_x, y = plot_y), point_size = 0.5, color = "gray50",
      rasterize = rasterize) + ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::ggtitle(paste(sample_id, "No retained residual cells")) + spatial_black_theme()
  }
  groups <- if (length(present)) unname(split(present, ceiling(seq_along(present) / 3L))) else list()
  facet_plot <- function(rasterize = FALSE) {
    if (length(present)) return(build_faceted_spatial_rows(spatial, palette, groups, rasterize))
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No retained residual cells in this sample", color = "white", size = 6) + spatial_black_theme()
  }
  global06 <- global(0.6)
  global085 <- global(0.85)
  facets <- facet_plot()
  requested <- marker_sets[present]
  availability <- data.frame(residual_subcluster = as.character(rep(names(requested), lengths(requested))),
                              gene = as.character(unlist(requested, use.names = FALSE)), stringsAsFactors = FALSE)
  availability$marker_source <- as.character(attr(marker_sets, "marker_sources")[availability$residual_subcluster])
  if (checked$residual_cells > 0L) {
    if (!"Xenium" %in% names(object@assays)) stop("Residual DotPlot requires the Xenium assay.")
    availability$present <- availability$gene %in% rownames(object[["Xenium"]])
    features <- lapply(requested, function(genes) intersect(genes, rownames(object[["Xenium"]])))
    if (any(!lengths(features))) stop("No reviewed DotPlot marker is present for: ",
                                     paste(names(features)[!lengths(features)], collapse = ", "))
    dot_object <- subset(object, cells = rownames(object[[]])[checked$retained])
    Seurat::DefaultAssay(dot_object) <- "Xenium"
    layers <- SeuratObject::Layers(dot_object[["Xenium"]], search = "^data($|\\.)")
    if (!length(layers)) stop("Mapped Xenium assay lacks normalized data; no normalization is performed by the plotting stage.")
    if (length(layers) > 1L) dot_object <- SeuratObject::JoinLayers(dot_object, assay = "Xenium", layers = "^data($|[.])", new = "data")
    dot_object$residual_subcluster <- factor(dot_object$residual_subcluster, levels = rev(present))
    Seurat::Idents(dot_object) <- "residual_subcluster"
    dot <- workflow_dotplot(dot_object, features = features, assay = "Xenium",
                            idents = rev(present), cols = c("lightgrey", "red"),
                            dot.scale = 6, cluster.idents = FALSE, scale = length(present) > 1L) +
      Seurat::RotatedAxis() + ggplot2::ggtitle(paste(sample_id, "Residual markers"))
    dot_dimensions <- workflow_dotplot_dimensions(dot)
    dot_display <- patchwork::wrap_elements(full = workflow_dotplot_grob(dot, strip_colors = palette))
  } else {
    availability$present <- logical()
    dot_display <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No retained residual cells: marker DotPlot not applicable", size = 5) + ggplot2::theme_void()
    dot_dimensions <- c(width = 12, height = 4)
  }
  report <- patchwork::wrap_plots(
    A = global(0.6, TRUE) + ggplot2::labs(subtitle = "Residual point size 0.6") +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(color = "white")),
    B = global(0.85, TRUE) + ggplot2::labs(subtitle = "Residual point size 0.85") +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(color = "white")),
    C = facet_plot(TRUE), D = dot_display, design = "AB\nCC\nDD", heights = c(1, 1.4, 0.8)
  ) + patchwork::plot_annotation(
    title = paste0(sample_id, " | PCW ", gsub("[^0-9.]", "", checked$PCW),
                   " | Residual cells: ", checked$residual_cells),
    subtitle = "Grey: other retained cells. Approved residual QC exclusions omitted. Facets: residual cells only.",
    theme = ggplot2::theme(plot.background = ggplot2::element_rect(fill = "black", color = NA),
      plot.title = ggplot2::element_text(color = "white", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "white", hjust = 0.5))
  )
  dir.create(dirname(paths[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
  save_cairo_tiff(global06, paths[["global06"]], 12, 10)
  save_cairo_tiff(global085, paths[["global085"]], 12, 10)
  save_cairo_tiff(facets, paths[["facets"]], 15, max(5, 4 * length(groups)))
  save_cairo_tiff(dot_display, paths[["dot_tiff"]], dot_dimensions[["width"]],
                   dot_dimensions[["height"]], background = "white")
  ggplot2::ggsave(paths[["dot_pdf"]], dot_display, device = grDevices::cairo_pdf,
                 width = dot_dimensions[["width"]], height = dot_dimensions[["height"]], bg = "white")
  write.csv(availability, paths[["markers"]], row.names = FALSE)
  ggplot2::ggsave(paths[["pdf"]], report, device = grDevices::cairo_pdf,
                 width = 24, height = 18, units = "in", bg = "black", limitsize = FALSE)
  if (pdf_page_count(paths[["pdf"]]) != 1L) stop("Residual sample report must have exactly one page.")
  counts <- data.frame(residual_subcluster = present,
                        cell_count = as.integer(table(factor(checked$labels, levels = present))))
  counts$percentage <- if (checked$residual_cells) 100 * counts$cell_count / checked$residual_cells else numeric()
  if (sum(counts$cell_count) != checked$residual_cells) stop("Residual plot counts do not cover all mapped residual cells.")
  write.csv(counts, paths[["counts"]], row.names = FALSE)
  source_info <- file.info(input_path)
  if (!identical(as.numeric(input_info$size), as.numeric(source_info$size)) ||
      !identical(as.numeric(input_info$mtime), as.numeric(source_info$mtime))) stop("Mapped RDS changed during reporting.")
  status <- data.frame(sample_id = sample_id, input_path = input_path,
    input_size = as.numeric(source_info$size), input_mtime = as.numeric(source_info$mtime),
    total_cells = ncol(object), mapped_residual_cells = checked$residual_cells,
    removed_residual_qc_cells = sum(checked$removed), plotted_cells = nrow(spatial),
    PCW = checked$PCW, review_md5 = review_md5,
    marker_panel_md5 = residual_marker_panel_signature(marker_sets),
    status = if (checked$residual_cells) "plotted" else "no_retained_residual_cells",
    stringsAsFactors = FALSE)
  for (key in names(paths)[names(paths) != "status"]) {
    status[[paste0(key, "_path")]] <- paths[[key]]
    status[[paste0(key, "_bytes")]] <- as.numeric(file.info(paths[[key]])$size)
  }
  write.csv(status, paths[["status"]], row.names = FALSE)
  invisible(status)
}

audit_residual_spatial_reports <- function(root, mapping, input_paths, review_md5, marker_panel_md5) {
  rows <- lapply(seq_len(nrow(mapping)), function(i) {
    paths <- residual_report_paths(root, as.character(mapping$sample_id[[i]]))
    if (any(!file.exists(paths))) stop("Incomplete residual report outputs for ", mapping$sample_id[[i]])
    status <- read.csv(paths[["status"]], stringsAsFactors = FALSE)
    info <- file.info(input_paths[[i]])
    required <- c("sample_id", "input_path", "input_size", "input_mtime", "total_cells",
                  "mapped_residual_cells", "removed_residual_qc_cells", "review_md5", "marker_panel_md5", "status")
    if (nrow(status) != 1L || !all(required %in% names(status)) ||
        !identical(as.character(status$sample_id), as.character(mapping$sample_id[[i]])) ||
        !identical(as.character(status$input_path), input_paths[[i]]) ||
        !identical(as.character(status$review_md5), review_md5) ||
        !identical(as.character(status$marker_panel_md5), marker_panel_md5) ||
        !isTRUE(status$input_size == as.numeric(info$size)) ||
        !isTRUE(abs(status$input_mtime - as.numeric(info$mtime)) <= 1e-4) ||
        !isTRUE(status$total_cells == mapping$total_cells[[i]]) ||
        !isTRUE(status$mapped_residual_cells == mapping$mapped_residual_cells[[i]]) ||
        !isTRUE(status$removed_residual_qc_cells == mapping$removed_residual_qc_cells[[i]])) {
      stop("Residual report provenance is incomplete or stale for ", mapping$sample_id[[i]])
    }
    for (key in names(paths)[names(paths) != "status"]) {
      if (!isTRUE(status[[paste0(key, "_path")]] == paths[[key]]) ||
          !isTRUE(status[[paste0(key, "_bytes")]] == as.numeric(file.info(paths[[key]])$size))) {
        stop("Residual report output changed for ", mapping$sample_id[[i]], ": ", key)
      }
    }
    counts <- read.csv(paths[["counts"]], stringsAsFactors = FALSE)
    if (!all(c("residual_subcluster", "cell_count", "percentage") %in% names(counts)) ||
        anyDuplicated(counts$residual_subcluster) || anyNA(counts$cell_count) ||
        any(counts$cell_count <= 0) || sum(counts$cell_count) != status$mapped_residual_cells) {
      stop("Invalid residual report counts for ", mapping$sample_id[[i]])
    }
    if (pdf_page_count(paths[["pdf"]]) != 1L) stop("Residual sample PDF is not one page.")
    status
  })
  do.call(rbind, rows)
}
