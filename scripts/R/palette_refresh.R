# Plot-only snapshots: stored metadata, reductions and assay values are read-only.
refresh_corrected_broad_paths <- function(output_root, sample_ids) {
  stages <- file.path(output_root, "xenium", "vz_rl",
                      c("02_corrected_broad", "01_corrected_broad"), "rds")
  candidates <- lapply(stages, function(stage) {
    file.path(stage, paste0(sample_ids, "_VZ_RL_QC_CorrectedBroad.rds"))
  })
  found <- vapply(candidates, function(paths) sum(file.exists(paths)), integer(1))
  populated <- which(found > 0L)
  if (length(populated) > 1L) {
    stop("Corrected-broad RDS files occur in both documented layouts; resolve duplicate/partial sets before plotting:\n- ",
         paste0(stages[populated], " (", found[populated], "/", length(sample_ids), ")", collapse = "\n- "))
  }
  chosen <- if (length(populated)) populated[[1]] else 1L
  message("Corrected-broad layout: ", stages[[chosen]], " (", found[[chosen]], "/", length(sample_ids), ")")
  if (!length(populated)) message("Also checked: ", stages[[2]])
  # Keep an incomplete set incomplete; never fill gaps from another stage.
  stats::setNames(candidates[[chosen]], sample_ids)
}

refresh_display_labels <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & x == "Microglia"] <- "Immune"
  x
}

refresh_normalized_layers <- function(assay) {
  # Legacy Assay treats `search` as exact names; Assay5 supports patterns.
  # Enumerate first so both classes select the same stored normalized layers.
  grep("^data($|[.])", SeuratObject::Layers(assay, search = NA), value = TRUE)
}

refresh_marker_panels <- function(kind, identities, review = NULL) {
  raw_identities <- identities
  raw_identities[raw_identities == "Immune"] <- "Microglia"
  rl_identities <- intersect(rl_subcluster_order, raw_identities)
  panels <- switch(kind,
    broad = markers,
    vz = active_vz_dotplot_markers(raw_identities),
    rl = rl_production_marker_panel(rl_identities),
    residual = residual_report_marker_sets(review, residual_shared_marker_panels(), allow_repeated_genes = TRUE),
    combined = c(active_vz_dotplot_markers(raw_identities),
                 if (length(rl_identities)) rl_production_marker_panel(rl_identities) else list()),
    stop("Unknown marker-panel kind: ", kind)
  )
  if (kind == "residual") panels <- panels[refresh_display_labels(names(panels)) %in% identities]
  # Combined panels share genes. Retain their first group and exact within-group order.
  if (kind == "combined") {
    unique_panels <- list()
    seen <- character()
    for (i in seq_along(panels)) {
      features <- setdiff(panels[[i]], seen)
      if (length(features)) {
        group <- names(panels)[[i]]
        unique_panels[[group]] <- unique(c(unique_panels[[group]], features))
        seen <- c(seen, features)
      }
    }
    panels <- unique_panels
  }
  names(panels) <- refresh_display_labels(names(panels))
  # One displayed copy per gene is required by workflow_dotplot(). Keep the
  # first configured group, with an explicit audit rather than a silent drop.
  seen <- character()
  omitted <- character()
  for (i in seq_along(panels)) {
    duplicate <- panels[[i]][panels[[i]] %in% seen]
    if (length(duplicate)) omitted <- c(omitted, paste(names(panels)[[i]], duplicate, sep = ":"))
    panels[[i]] <- setdiff(panels[[i]], seen)
    seen <- c(seen, panels[[i]])
  }
  panels <- panels[lengths(panels) > 0L]
  attr(panels, "omitted_duplicate_markers") <- omitted
  panels
}

refresh_plot_object <- function(object, input_path, output_dir, label_column,
                                kind, title, spatial = FALSE, review = NULL) {
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Snapshot directory is not empty: ", output_dir, ". Choose a new --tag; existing snapshots are protected.")
  }
  metadata <- object[[]]
  cells <- colnames(object)
  if (!length(cells) || anyDuplicated(cells) || !identical(rownames(metadata), cells)) {
    stop("Empty object, duplicate cell IDs, or metadata/cell order mismatch: ", input_path)
  }
  if (!label_column %in% names(metadata)) stop("Missing metadata: ", label_column, " in ", input_path)
  raw_labels <- as.character(metadata[[label_column]])
  # Mapped regional objects deliberately retain nonbranch cells as NA background.
  allow_background <- spatial && kind %in% c("vz", "rl", "residual")
  missing_labels <- is.na(raw_labels) | !nzchar(trimws(raw_labels))
  if (any(missing_labels) && !allow_background) stop("Unexpected missing plotting labels: ", input_path)
  display_labels <- refresh_display_labels(raw_labels)
  observed <- unique(display_labels[!missing_labels])
  preferred <- switch(kind, broad = celltype_order, vz = vz_subcluster_order,
                      rl = rl_subcluster_order, residual = master_subcluster_order,
                      combined = master_subcluster_order)
  preferred <- unique(refresh_display_labels(preferred))
  if (kind == "broad") {
    identity_order <- order_consensus_levels(observed, celltype_order)
    colors <- resolve_cluster_colors(identity_order)
  } else {
    unsupported <- setdiff(observed, c(preferred, observed[is_unknown_consensus_label(observed)]))
    if (length(unsupported)) stop(
      "Unregistered plotting identities: ", paste(unsupported, collapse = ", "),
      ". Input: ", input_path,
      if (kind == "vz" && "Mature UBCs" %in% unsupported) {
        ". This VZ object contains a label retired from the current VZ taxonomy. Verify that reviewed-label application and sample mapping have completed; a plot refresh will not relabel these cells."
      } else ""
    )
    identity_order <- c(intersect(preferred, observed), setdiff(observed, preferred))
    colors <- subcluster_palette[identity_order]
    unknown <- is_unknown_consensus_label(identity_order)
    colors[unknown] <- resolve_cluster_colors(identity_order[unknown])
    if (anyNA(colors)) stop("Missing subcluster colors in ", input_path)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- character()
  report_plots <- list()
  save_plot <- function(plot, stem, width, height, background = "white", report = TRUE) {
    pdf <- file.path(output_dir, paste0(stem, ".pdf"))
    tif <- file.path(output_dir, paste0(stem, ".tif"))
    ggplot2::ggsave(pdf, plot, device = grDevices::cairo_pdf, width = width, height = height, bg = background)
    grDevices::tiff(tif, width = width, height = height, units = "in", res = 600,
                   type = "cairo", compression = "lzw", bg = background)
    tryCatch(print(plot), finally = grDevices::dev.off())
    outputs <<- c(outputs, pdf, tif)
    if (report) report_plots[[stem]] <<- plot
  }
  if (!length(observed)) {
    # Explicitly valid for a mapped regional sample with no cells in that branch.
    placeholder <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = paste(title, "- no mapped cells"))
    save_plot(placeholder, "No_mapped_cells", 9, 5)
  } else {
    assay <- if (kind == "broad" && label_column == "clusters_refined") Seurat::DefaultAssay(object) else "Xenium"
    if (!assay %in% names(object@assays)) stop("Missing assay ", assay, " in ", input_path)
    Seurat::DefaultAssay(object) <- assay
    layers <- refresh_normalized_layers(object[[assay]])
    if (!length(layers)) stop("Normalized data missing; refresh never normalizes assays: ", input_path)
    if (length(layers) > 1L) {
      object <- SeuratObject::JoinLayers(object, assay = assay, layers = "^data($|[.])", new = "data")
    }
    if (".palette_refresh_label" %in% names(metadata)) stop("Temporary plotting column already exists.")
    labels_for_plot <- display_labels
    labels_for_plot[missing_labels] <- ".background"
    object$.palette_refresh_label <- factor(labels_for_plot, levels = c(rev(identity_order), ".background"))
    Seurat::Idents(object) <- ".palette_refresh_label"
    panels <- refresh_marker_panels(kind, observed, review)
    if (!length(panels) || any(!lengths(panels))) stop("No valid marker panels for ", title)
    availability <- data.frame(group = rep(names(panels), lengths(panels)),
                               gene = unlist(panels, use.names = FALSE), stringsAsFactors = FALSE)
    availability$present <- availability$gene %in% rownames(object[[assay]])
    duplicate_path <- file.path(output_dir, "Duplicate_marker_display_notes.txt")
    notes <- attr(panels, "omitted_duplicate_markers")
    writeLines(c("Duplicate genes retain their first configured group; source marker lists are unchanged.",
                 if (length(notes)) notes else "None"), duplicate_path)
    outputs <- c(outputs, duplicate_path)
    # Production VZ panels have an exact feature contract. References/broad
    # panels retain their established present-gene intersection behavior.
    if (kind == "vz" && any(!availability$present)) {
      stop("Required VZ display markers missing: ", paste(availability$gene[!availability$present], collapse = ", "))
    }
    features <- lapply(panels, intersect, y = rownames(object[[assay]]))
    if (kind == "broad" && any(!lengths(features))) {
      message("Broad marker groups with no assay genes omitted (recorded in Marker_availability.csv): ",
              paste(names(features)[!lengths(features)], collapse = ", "))
      features <- features[lengths(features) > 0L]
    }
    if (!length(features)) stop("No configured markers occur in the selected assay: ", input_path)
    if (any(!lengths(features))) stop("No available markers for: ", paste(names(features)[!lengths(features)], collapse = ", "))
    marker_path <- file.path(output_dir, "Marker_availability.csv")
    utils::write.csv(availability, marker_path, row.names = FALSE)
    outputs <- c(outputs, marker_path)
    dot <- workflow_dotplot(object, features = features, assay = assay,
                            idents = rev(identity_order), cluster.idents = FALSE,
                            scale = length(observed) > 1L) + ggplot2::ggtitle(paste(title, "markers"))
    dimensions <- workflow_dotplot_dimensions(dot)
    strips <- if (kind == "vz") vz_dotplot_strip_colors else NULL
    display <- patchwork::wrap_elements(full = workflow_dotplot_grob(dot, strip_colors = strips))
    save_plot(display, "Marker_DotPlot", dimensions[["width"]], dimensions[["height"]])
    reductions <- grep("umap", Seurat::Reductions(object), value = TRUE, ignore.case = TRUE)
    if (kind %in% c("broad", "vz", "rl", "combined") && !spatial && !length(reductions)) {
      stop("Expected a saved UMAP reduction in ", input_path)
    }
    for (reduction in reductions) {
      umap <- Seurat::DimPlot(object, reduction = reduction, group.by = ".palette_refresh_label",
                              cols = c(colors, ".background" = "#D9D9D9"), raster = TRUE,
                              raster.dpi = c(1200L, 1200L)) + ggplot2::ggtitle(paste(title, reduction))
      save_plot(umap, paste0("UMAP_", reduction), 11, 8)
    }
    counts <- data.frame(identity = identity_order,
                         cells = as.integer(table(factor(display_labels, levels = identity_order))))
    counts$percent_of_labeled <- 100 * counts$cells / sum(counts$cells)
    counts_path <- file.path(output_dir, "Cluster_counts.csv")
    utils::write.csv(counts, counts_path, row.names = FALSE)
    outputs <- c(outputs, counts_path)
    counts$identity <- factor(counts$identity, levels = rev(identity_order))
    # Retain both established output stems for consumers of existing snapshots;
    # both figures now display abundance, while the CSV retains raw counts.
    for (metric in c("cells", "percent_of_labeled")) {
      bar <- ggplot2::ggplot(counts, ggplot2::aes(x = identity, y = percent_of_labeled, fill = identity)) +
        ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(values = colors, guide = "none") +
        ggplot2::scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
        ggplot2::labs(title = title, x = NULL, y = "Abundance (% of labeled cells)")
      save_plot(bar, paste0("Counts_", metric), 10, max(5, 2 + .28 * length(observed)))
    }
    if (spatial) {
      spatial_data <- prepare_spatial_plot_data(object, ".palette_refresh_label", identity_order, title)
      if (nrow(spatial_data) != length(cells) || !setequal(spatial_data$cell, cells)) {
        stop("Spatial coordinates must cover every input cell exactly once: ", input_path)
      }
      if (any(!is.finite(spatial_data$plot_x) | !is.finite(spatial_data$plot_y))) stop("Non-finite spatial coordinates.")
      global <- build_layered_spatial_plot(spatial_data, colors, title, targets = identity_order)
      save_plot(global, "Global_Spatial", 12, 9, "black", report = FALSE)
      report_plots$Global_Spatial <- build_layered_spatial_plot(spatial_data, colors, title,
                                                              targets = identity_order, rasterize = TRUE)
      rows <- split(identity_order, ceiling(seq_along(identity_order) / 4))
      facets <- build_faceted_spatial_rows(spatial_data, colors, rows)
      save_plot(facets, "Faceted_Spatial", 16, max(5, 4 * length(rows)), "black", report = FALSE)
      report_plots$Faceted_Spatial <- build_faceted_spatial_rows(spatial_data, colors, rows, rasterize = TRUE)
    }
  }
  report_path <- file.path(output_dir, "Plot_Report.pdf")
  grDevices::cairo_pdf(report_path, width = 18, height = 12, onefile = TRUE)
  tryCatch(for (plot in report_plots) print(plot), finally = grDevices::dev.off())
  outputs <- c(outputs, report_path)
  info <- file.info(outputs)
  if (anyNA(info$size) || any(info$size <= 0)) stop("Partial refresh output set for ", input_path)
  if (!identical(colnames(object), cells)) stop("Cell IDs changed during plot refresh.")
  # Every source metadata column must remain unchanged; the temporary plot
  # column and assay-layer view are never saved to an RDS.
  if (!identical(object[[]][names(metadata)], metadata)) stop("Source metadata changed during plot refresh.")
  data.frame(input_rds = normalizePath(input_path), input_size = file.info(input_path)$size,
             input_mtime = as.numeric(file.info(input_path)$mtime), label_column = label_column,
             cells = length(cells), plotted_cells = sum(!missing_labels),
             output_path = normalizePath(outputs), output_bytes = info$size, stringsAsFactors = FALSE)
}
