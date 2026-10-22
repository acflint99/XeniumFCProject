# Plot-only rendering for the saved combined object. No embeddings or labels are
# recomputed, and display subsets are never saved back to the processed RDS.
nonspatial_plot_spec <- function() {
  stems <- c(
    batch = "Xenium_Batch_Comp_UMAP",
    consensus = "Xenium_Merged_ConsensusLabel_UMAP",
    corrected = "Xenium_Merged_CorrectedBroad_UMAP",
    subclusters = "Xenium_Merged_Subclusters_UMAP",
    dot = "Xenium_Merged_CorrectedBroad_DotPlot_markers",
    purkinje = "Xenium_Purkinje_Specific_DotPlot",
    heatmap = "Xenium_Merged_Cluster_Top5_Heatmap",
    signaling = "Xenium_Signaling_Genes_dual_VlnPlot"
  )
  do.call(rbind, lapply(names(stems), function(family) {
    do.call(rbind, lapply(c(FALSE, TRUE), function(no_unknown) {
      stem <- stems[[family]]
      if (no_unknown) {
        stem <- if (grepl("_UMAP$", stem)) sub("_UMAP$", "_NoUnknown_UMAP", stem) else
          paste0(stem, "_NoUnknown")
      }
      extensions <- if (family %in% c("batch", "consensus", "corrected", "subclusters")) "tif" else c("tif", "pdf")
      data.frame(family = family, no_unknown = no_unknown,
                 filename = paste0(stem, ".", extensions), stringsAsFactors = FALSE)
    }))
  }))
}

# Unresolved is a heterogeneous retained population, not a claimed cell type.
# It may have no positive markers under the unchanged criteria. Keep its cells
# as heatmap columns and document the absent gene block. Other missing groups
# still indicate an unexpected/incomplete marker result and remain errors.
nonspatial_heatmap_marker_summary <- function(top, broad_levels) {
  if (!nrow(top)) stop("No positive corrected-broad heatmap markers found.")
  missing <- setdiff(broad_levels, as.character(top$cluster))
  unexpected <- setdiff(missing, "Unresolved")
  if (length(unexpected)) {
    stop("No top heatmap markers for: ", paste(unexpected, collapse = ", "))
  }
  if ("Unresolved" %in% missing) {
    message("HEATMAP NOTE: Unresolved has no qualifying positive markers; ",
            "retaining its cells and displaying markers selected for other groups. ",
            "Marker thresholds are unchanged.")
  }
  genes <- lapply(broad_levels, function(label) unique(as.character(
    top$gene[as.character(top$cluster) == label]
  )))
  data.frame(
    identity = broad_levels, markers_selected = lengths(genes),
    genes = vapply(genes, paste, collapse = ";", FUN.VALUE = character(1)),
    status = ifelse(lengths(genes) > 0L, "markers_selected", "no_qualifying_positive_markers_cells_retained"),
    stringsAsFactors = FALSE
  )
}

# Strict display order: fail instead of silently dropping an unregistered label.
nonspatial_plot_levels <- function(labels, known_order) {
  ordered <- order_consensus_levels(labels, known_order)
  missing <- setdiff(ordered[!is_unknown_consensus_label(ordered)], known_order)
  if (length(missing)) stop("Unregistered plot identity: ", paste(missing, collapse = ", "))
  ordered
}

nonspatial_label_umap <- function(object, column, label_order, colors, title, seed) {
  labels <- as.character(object[[column, drop = TRUE]])
  levels_present <- label_order[label_order %in% labels]
  # Explicit scale breaks govern the legend even when metadata was stored as text.
  p <- Seurat::DimPlot(
    object, reduction = "umap.harmony", group.by = column, label = FALSE,
    raster = TRUE, pt.size = 0.5, alpha = 0.8, cols = colors,
    shuffle = TRUE, seed = seed
  ) + ggplot2::scale_color_manual(values = colors, breaks = levels_present, name = NULL) +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(12, 24, 12, 24)) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 4, byrow = TRUE,
                                                override.aes = list(size = 3, alpha = 1))) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.16)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.16))
  xy <- Seurat::Embeddings(object, "umap.harmony")[colnames(object), 1:2, drop = FALSE]
  centers <- do.call(rbind, lapply(levels_present, function(label) {
    keep <- labels == label
    data.frame(x = stats::median(xy[keep, 1]), y = stats::median(xy[keep, 2]), label = label)
  }))
  p + ggrepel::geom_label_repel(
    data = centers, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
    seed = seed, size = if (column == "comb_subcluster") 3.5 else 4,
    fill = "white", color = "black", box.padding = 0.65, point.padding = 0.3,
    force = 4, max.overlaps = Inf, max.iter = 30000, max.time = 10,
    min.segment.length = 0
  )
}

render_nonspatial_merge_plots <- function(object, plot_dir, sample_ids, seed, processed_path) {
  required <- c("orig.ident", "consensus_label", "corrected_broad_label", "comb_subcluster", "VZ_subcluster")
  missing <- setdiff(required, colnames(object[[]]))
  if (length(missing)) stop("Processed object lacks plot metadata: ", paste(missing, collapse = ", "))
  if (anyDuplicated(colnames(object)) || !identical(rownames(object[[]]), colnames(object))) {
    stop("Processed cell IDs must be unique and aligned to metadata.")
  }
  for (column in required[1:4]) {
    labels <- as.character(object[[column, drop = TRUE]])
    if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Missing/blank plot labels in ", column)
  }
  if (!setequal(as.character(object$orig.ident), sample_ids)) stop("Plot samples differ from config/samples.csv.")
  for (reduction in c("umap.unintegrated", "umap.harmony")) {
    xy <- Seurat::Embeddings(object, reduction)
    if (!identical(rownames(xy), colnames(object)) || ncol(xy) < 2L || any(!is.finite(xy[, 1:2]))) {
      stop("Invalid or misaligned UMAP coordinates in ", reduction)
    }
  }
  if (!requireNamespace("Cairo", quietly = TRUE)) stop("Cairo is required in the project renv environment.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required in the project renv environment.")
  Seurat::DefaultAssay(object) <- "Xenium"
  broad_levels <- nonspatial_plot_levels(object$corrected_broad_label, celltype_order)
  consensus_levels <- nonspatial_plot_levels(object$consensus_label, celltype_order)
  # Place subclusters within the same broad-family sequence, preserving the
  # established within-family order. Unknowns always appear last.
  fine_observed <- nonspatial_plot_levels(object$comb_subcluster, master_subcluster_order)
  fine_broad <- map_corrected_broad_labels(fine_observed)
  fine_levels <- fine_observed[order(match(fine_broad, celltype_order), seq_along(fine_observed))]
  fine_colors <- stats::setNames(unname(subcluster_palette[fine_levels]), fine_levels)
  unknown_fine <- fine_levels[is_unknown_consensus_label(fine_levels)]
  fine_colors[unknown_fine] <- resolve_cluster_colors(unknown_fine)
  if (anyNA(fine_colors)) stop("Missing subcluster colors: ", paste(names(fine_colors)[is.na(fine_colors)], collapse = ", "))
  object$corrected_broad_label <- factor(as.character(object$corrected_broad_label), levels = broad_levels)
  object$consensus_label <- factor(as.character(object$consensus_label), levels = consensus_levels)
  object$comb_subcluster <- factor(as.character(object$comb_subcluster), levels = fine_levels)
  object$orig.ident <- factor(as.character(object$orig.ident), levels = sample_ids)
  broad_colors <- resolve_cluster_colors(broad_levels)
  sample_colors <- stats::setNames(scales::hue_pal()(length(sample_ids)), sample_ids)
  specs <- nonspatial_plot_spec()
  records <- list()
  save_plot <- function(plot, family, no_unknown, width, height, cells, column, note = "") {
    rows <- specs[specs$family == family & specs$no_unknown == no_unknown, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      path <- file.path(plot_dir, rows$filename[[i]])
      if (endsWith(path, ".tif")) {
        Cairo::CairoTIFF(path, width = width, height = height, units = "in", res = 600)
        tryCatch({
          if (inherits(plot, "grob")) grid::grid.draw(plot) else print(plot)
        }, finally = grDevices::dev.off())
      } else {
        ggplot2::ggsave(path, plot, device = grDevices::cairo_pdf, width = width, height = height, limitsize = FALSE)
      }
      records[[length(records) + 1L]] <<- data.frame(
        filename = rows$filename[[i]], family = family, no_unknown = no_unknown,
        grouping = column, cells = cells, total_cells = ncol(object),
        processed_path = processed_path, seed = seed, note = note, stringsAsFactors = FALSE
      )
    }
  }
  display_object <- function(column, no_unknown) {
    if (!no_unknown) return(object)
    cells <- colnames(object)[!is_unknown_consensus_label(object[[column, drop = TRUE]])]
    if (!length(cells)) stop("No non-Unknown cells for ", column)
    if (length(cells) == ncol(object)) return(object)
    subset(object, cells = cells)
  }
  suffix <- function(no_unknown) if (no_unknown) " (Unknown excluded)" else " (with Unknowns)"

  # Shared batch key below two equal-aspect panels; no per-panel side legends.
  for (no_unknown in c(FALSE, TRUE)) {
    display <- display_object("corrected_broad_label", no_unknown)
    batch_panel <- function(reduction, title) {
      Seurat::DimPlot(display, reduction = reduction, group.by = "orig.ident", raster = TRUE,
                      cols = sample_colors, shuffle = TRUE, seed = seed) +
        ggplot2::scale_color_manual(values = sample_colors, breaks = sample_ids, drop = FALSE, name = "Sample") +
        ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
        ggplot2::theme_classic(base_size = 14) + ggplot2::coord_fixed() +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                    override.aes = list(size = 3, alpha = 1)))
    }
    p <- (batch_panel("umap.unintegrated", "Pre-Harmony (PCA UMAP)") +
            batch_panel("umap.harmony", "Post-Harmony (Harmony UMAP)")) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(title = paste0("Batch comparison", suffix(no_unknown)))
    p <- p & ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 11))
    save_plot(p, "batch", no_unknown, 16, 12, ncol(display), "orig.ident",
              "NoUnknown uses corrected_broad_label; sample order follows manifest.")
    rm(display, p)
  }

  umap_columns <- c(consensus = "consensus_label", corrected = "corrected_broad_label", subclusters = "comb_subcluster")
  umap_titles <- c(consensus = "Original Broad Labels (Consensus)", corrected = "Corrected Broad Labels", subclusters = "Subcluster Labels")
  for (family in names(umap_columns)) {
    column <- umap_columns[[family]]
    label_order <- switch(family, consensus = consensus_levels, corrected = broad_levels, subclusters = fine_levels)
    colors <- if (family == "subclusters") fine_colors else resolve_cluster_colors(label_order)
    for (no_unknown in c(FALSE, TRUE)) {
      display <- display_object(column, no_unknown)
      p <- nonspatial_label_umap(display, column, label_order, colors,
                                paste0("Xenium Merged: ", umap_titles[[family]], suffix(no_unknown)), seed)
      save_plot(p, family, no_unknown, if (family == "subclusters") 16 else 12,
                if (family == "subclusters") 14 else 11, ncol(display), column,
                "NoUnknown uses the displayed label column; saved embeddings retained.")
      rm(display, p)
    }
  }

  # Marker panels and rows follow the same broad order. Shared genes retain
  # their first configured group, as in the existing workflow helper.
  marker_groups <- markers[intersect(celltype_order, names(markers))]
  marker_groups <- lapply(marker_groups, intersect, y = rownames(object))
  marker_groups <- marker_groups[lengths(marker_groups) > 0L]
  if (!length(marker_groups)) stop("No broad marker genes are present.")
  for (no_unknown in c(FALSE, TRUE)) {
    display <- display_object("corrected_broad_label", no_unknown)
    present <- broad_levels[broad_levels %in% as.character(display$corrected_broad_label)]
    Seurat::Idents(display) <- factor(display$corrected_broad_label, levels = rev(present))
    p <- workflow_dotplot(display, features = marker_groups, assay = "Xenium") +
      ggplot2::labs(title = paste0("Corrected Broad Marker Expression", suffix(no_unknown)), x = NULL, y = NULL)
    dimensions <- workflow_dotplot_dimensions(p)
    colored <- workflow_dotplot_grob(p, broad_colors = cluster_colors)
    strip_colors <- attr(colored, "workflow_strip_colors")
    if (!length(strip_colors) || anyNA(strip_colors) ||
        !identical(unname(strip_colors), unname(cluster_colors[names(strip_colors)]))) {
      stop("Broad marker-section colors do not match their registered cluster labels.")
    }
    save_plot(colored, "dot", no_unknown, dimensions[["width"]], dimensions[["height"]],
              ncol(display), "corrected_broad_label")
    rm(display, p, colored)
  }

  pc_clusters <- c("Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs")
  pc_genes <- intersect(c("FOXP1", "ITPR1", "COL5A1", "VSTM2L", "NDNF", "EBF2", "CALB1",
                          "TRPC3", "NEFL", "ETV1", "PCDH10", "EBF1", "BCL11A", "RORB", "EN1"), rownames(object))
  pc_cells <- colnames(object)[as.character(object$VZ_subcluster) %in% pc_clusters]
  if (!length(pc_cells) || !length(pc_genes)) stop("Purkinje cells or marker genes are absent.")
  display <- subset(object, cells = pc_cells)
  Seurat::Idents(display) <- factor(display$VZ_subcluster, levels = rev(pc_clusters))
  p <- workflow_dotplot(display, features = pc_genes, assay = "Xenium") +
    ggplot2::labs(title = "Purkinje Cell Lineage Marker Expression", x = NULL, y = NULL)
  dimensions <- workflow_dotplot_dimensions(p)
  for (no_unknown in c(FALSE, TRUE)) {
    save_plot(p, "purkinje", no_unknown, dimensions[["width"]], dimensions[["height"]],
              ncol(display), "VZ_subcluster", "Identical pair: this panel contains only the four named Purkinje identities.")
  }
  rm(display, p)
  invisible(gc())

  # Reuse the existing marker-test settings with the requested corrected broad
  # grouping. Compute once; paired heatmaps share known genes and sampled cells.
  Seurat::Idents(object) <- factor(object$corrected_broad_label, levels = broad_levels)
  set.seed(seed)
  heatmap_markers <- Seurat::FindAllMarkers(
    object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 500, random.seed = seed
  )
  if (!nrow(heatmap_markers)) stop("No positive corrected-broad heatmap markers found.")
  top <- heatmap_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(n = 5, order_by = avg_log2FC, with_ties = FALSE) |>
    dplyr::ungroup()
  marker_summary <- nonspatial_heatmap_marker_summary(top, broad_levels)
  write.csv(marker_summary, file.path(plot_dir, "Xenium_Merged_Heatmap_marker_summary.csv"), row.names = FALSE)
  marker_note <- if (any(marker_summary$markers_selected == 0L)) {
    "Unresolved has no qualifying positive markers; its cells are retained."
  } else ""
  top <- top[order(match(as.character(top$cluster), broad_levels), -top$avg_log2FC), , drop = FALSE]
  top_genes <- unique(as.character(top$gene))
  temp <- Seurat::ScaleData(object, features = top_genes, verbose = FALSE)
  set.seed(seed)
  temp <- subset(temp, downsample = 100)
  for (no_unknown in c(FALSE, TRUE)) {
    display <- if (no_unknown) subset(temp, cells = colnames(temp)[!is_unknown_consensus_label(temp$corrected_broad_label)]) else temp
    present <- broad_levels[broad_levels %in% as.character(display$corrected_broad_label)]
    display$corrected_broad_label <- factor(display$corrected_broad_label, levels = present)
    Seurat::Idents(display) <- "corrected_broad_label"
    genes <- unique(as.character(top$gene[as.character(top$cluster) %in% present]))
    p <- Seurat::DoHeatmap(
      display, features = rev(genes), group.by = "corrected_broad_label",
      group.colors = unname(broad_colors[present]), size = 3.5, angle = 45,
      draw.lines = TRUE, raster = FALSE, combine = FALSE
    )[[1]] + ggplot2::scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
      ggplot2::guides(color = "none") +
      ggplot2::labs(
        title = paste0("Top 5 Markers per Corrected Broad Cluster", suffix(no_unknown)),
        subtitle = if (nzchar(marker_note)) marker_note else NULL
      )
    save_plot(p, "heatmap", no_unknown, 16, 12, ncol(display), "corrected_broad_label",
              paste("At most 100 cells per corrected broad identity; shared full-cohort marker tests and scaling.", marker_note))
    rm(display, p)
  }
  rm(temp)
  invisible(gc())

  kit_genes <- intersect(unique(c("APP", "TNFRSF21", "CADM3", "CADM4", "NECTIN3", "CNTN2", "L1CAM",
                                 "CXCL12", "CXCR4", "EFNB2", "EPHA4", "GJA1", "NCAM1", "NRXN2",
                                 "CLSTN1", "ADGRL1", "DAG1", "LRRTM1", "RELN", "VLDLR", "CD99", "NTF3", "NTRK2")), rownames(object))
  kit_idents <- intersect(celltype_order, c("Granule", "Purkinje", "Glia", "GABA", "OPC", "UBC"))
  kit_idents <- intersect(kit_idents, broad_levels)
  if (!length(kit_genes) || !length(kit_idents)) stop("Signaling-gene violin inputs are absent.")
  kit_cells <- colnames(object)[as.character(object$corrected_broad_label) %in% kit_idents]
  display <- subset(object, cells = kit_cells)
  Seurat::Idents(display) <- factor(display$corrected_broad_label, levels = kit_idents)
  p <- Seurat::VlnPlot(display, features = kit_genes, pt.size = 0, ncol = 3,
                       cols = broad_colors[kit_idents]) &
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10))
  for (no_unknown in c(FALSE, TRUE)) {
    save_plot(p, "signaling", no_unknown, 20, 20, ncol(display), "corrected_broad_label",
              "Identical pair: the selected six broad identities contain no Unknowns.")
  }
  manifest <- do.call(rbind, records)
  if (!setequal(manifest$filename, specs$filename) || anyDuplicated(manifest$filename) ||
      any(!file.exists(file.path(plot_dir, specs$filename))) ||
      any(file.info(file.path(plot_dir, specs$filename))$size <= 0)) {
    stop("Incomplete combined plot output set.")
  }
  write.csv(manifest, file.path(plot_dir, "Xenium_Merged_plot_manifest.csv"), row.names = FALSE)
  message("Saved ", nrow(manifest), " combined plot files plus their provenance manifest.")
  invisible(manifest)
}
