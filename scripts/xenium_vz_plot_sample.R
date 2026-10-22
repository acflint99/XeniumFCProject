# Per-sample VZ spatial and marker reports.
library(Seurat)
library(ggplot2)
library(here)
library(dplyr)

source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))

generate_spatial_reports <- function(sample_id, input_rds_path, base_plot_dir) {
  message("--- Processing Sample: ", sample_id, " ---")
  on.exit(while (!is.null(dev.list())) dev.off(), add = TRUE)

  sample_plot_dir <- file.path(base_plot_dir, sample_id)
  if (!dir.exists(sample_plot_dir)) dir.create(sample_plot_dir, recursive = TRUE)
  output_stem <- paste0(sample_id, "_VZ")

  obj <- readRDS(input_rds_path)
  required_metadata <- c("VZ_subcluster", "consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop(sample_id, " lacks required VZ plotting metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!"Xenium" %in% names(obj@assays)) stop(sample_id, " lacks the Xenium assay.")
  DefaultAssay(obj) <- "Xenium"

  observed_vz_labels <- unique(as.character(obj$VZ_subcluster))
  observed_vz_labels <- observed_vz_labels[!is.na(observed_vz_labels) & nzchar(observed_vz_labels)]
  unsupported_vz_labels <- setdiff(observed_vz_labels, vz_subcluster_order)
  missing_vz_colors <- setdiff(observed_vz_labels, names(vz_palette))
  if (length(unsupported_vz_labels) || length(missing_vz_colors)) {
    stop(
      sample_id, " has VZ labels not fully registered in the shared taxonomy. Missing from order: ",
      if (length(unsupported_vz_labels)) paste(unsupported_vz_labels, collapse = ", ") else "none",
      "; missing from palette: ",
      if (length(missing_vz_colors)) paste(missing_vz_colors, collapse = ", ") else "none"
    )
  }
  sample_vz_order <- vz_subcluster_order[vz_subcluster_order %in% observed_vz_labels]
  if (!length(sample_vz_order)) stop(sample_id, " contains no mapped VZ subclusters.")
  obj$VZ_subcluster <- factor(as.character(obj$VZ_subcluster), levels = sample_vz_order)

  stats <- obj@meta.data %>%
    group_by(VZ_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    mutate(percentage = (cell_count / sum(cell_count)) * 100)

  spatial_data <- prepare_spatial_plot_data(
    object = obj,
    cluster_column = "VZ_subcluster",
    cluster_order = sample_vz_order,
    sample_id = sample_id
  )

  p_global <- build_layered_spatial_plot(
    spatial_data,
    vz_palette,
    paste(sample_id, "VZ Subclusters"),
    targets = sample_vz_order
  )
  vz_lineage_clusters <- intersect(vz_lineage_subclusters, sample_vz_order)
  p_global_no_other <- build_layered_spatial_plot(
    spatial_data,
    vz_palette,
    paste(sample_id, "VZ Lineages"),
    targets = vz_lineage_clusters
  )
  highlight_definitions <- list(
    Purkinje_Lineage = c(
      "VZPs", "Maturing PCs", "Early-born PCs",
      "Late-born PCs", "Patterning PCs"
    ),
    Glia_Lineage = c(
      "VZPs", "RG Progenitors", "Cycling Glia",
      "BG", "Astrocytes/Ependyma", "Pre-OPCs", "OPCs"
    ),
    GABA_Lineage = c(
      "VZPs", "GABA Progenitors", "GABA Interneurons", "Golgi Cells", "MLIs", "iCN"
    )
  )
  highlight_plots <- lapply(names(highlight_definitions), function(plot_name) {
    build_layered_spatial_plot(
      spatial_data,
      vz_palette,
      paste(sample_id, plot_name),
      targets = highlight_definitions[[plot_name]]
    )
  })
  names(highlight_plots) <- names(highlight_definitions)

  lineage_rows <- list(
    c(
      "VZPs", "Maturing PCs", "Early-born PCs",
      "Late-born PCs", "Patterning PCs"
    ),
    c("GABA Progenitors", "GABA Interneurons", "Golgi Cells", "MLIs", "iCN"),
    c("RG Progenitors", "Cycling Glia", "BG", "Astrocytes/Ependyma", "Pre-OPCs", "OPCs")
  )
  other_rows <- list(
    c("RL VZ", "Differentiating GCs", "eCN"),
    "Unresolved"
  )
  p_facet <- build_faceted_spatial_rows(spatial_data, vz_palette, lineage_rows)
  p_other <- build_faceted_spatial_rows(spatial_data, vz_palette, other_rows)

  # Pre-rasterize only the spatial panels embedded in the PDF report. The
  # standalone TIFFs use direct points so their nominal sizes remain exact.
  p_global_report <- build_layered_spatial_plot(
    spatial_data,
    vz_palette,
    paste(sample_id, "VZ Subclusters"),
    targets = sample_vz_order,
    rasterize = TRUE
  )
  p_global_no_other_report <- build_layered_spatial_plot(
    spatial_data,
    vz_palette,
    paste(sample_id, "VZ Lineages"),
    targets = vz_lineage_clusters,
    rasterize = TRUE
  )
  p_facet_report <- build_faceted_spatial_rows(
    spatial_data,
    vz_palette,
    lineage_rows,
    rasterize = TRUE
  )

  plotted_cells <- Cells(obj)[!is.na(obj$VZ_subcluster)]
  if (!length(plotted_cells)) stop(sample_id, " contains no mapped VZ cells.")
  dotplot_identity <- as.character(obj$VZ_subcluster)
  dotplot_identity[is.na(dotplot_identity)] <- ".not_plotted"
  obj$subcluster_rev <- factor(
    dotplot_identity,
    levels = c(rev(sample_vz_order), ".not_plotted")
  )
  Idents(obj) <- "subcluster_rev"
  xenium_features <- rownames(obj[["Xenium"]])
  requested_compact_vz_markers <- unlist(
    active_vz_dotplot_markers(sample_vz_order), use.names = FALSE
  )
  missing_compact_vz_markers <- setdiff(
    requested_compact_vz_markers, xenium_features
  )
  if (length(missing_compact_vz_markers)) {
    stop(
      sample_id,
      " is missing compact VZ DotPlot marker genes from the Xenium assay: ",
      paste(missing_compact_vz_markers, collapse = ", ")
    )
  }
  p_dot <- workflow_dotplot(
    obj,
    features = active_vz_dotplot_markers(sample_vz_order),
    assay = "Xenium",
    idents = rev(sample_vz_order),
    cols = c("lightgrey", "red"),
    dot.scale = 6,
    cluster.idents = FALSE,
    scale = length(sample_vz_order) > 1L
  ) +
    RotatedAxis() +
    theme(
      axis.text.x = element_text(size = 12, face = "italic"),
      axis.text.y = element_text(size = 12, face = "plain")
    ) +
    ggtitle(paste(sample_id, "VZ Markers"))

  dotplot_dimensions <- workflow_dotplot_dimensions(p_dot)
  p_dot_display <- patchwork::wrap_elements(full = workflow_dotplot_grob(
    p_dot, strip_colors = vz_dotplot_strip_colors
  ))

  p_report <- build_spatial_report_page(
    p_global_report,
    p_global_no_other_report,
    p_facet_report,
    p_dot_display,
    sample_id
  )

  # Build every plot before writing files so a DotPlot error cannot leave a partial report.
  write.csv(
    stats,
    file.path(sample_plot_dir, paste0(output_stem, "_cluster_counts.csv")),
    row.names = FALSE
  )
  save_cairo_tiff(
    p_global,
    file.path(sample_plot_dir, paste0(output_stem, "_Global_Spatial_pt0.75.tif")),
    width = 10,
    height = 8
  )
  save_cairo_tiff(
    p_global_no_other,
    file.path(sample_plot_dir, paste0(output_stem, "_Global_Spatial_NoOther_pt0.75.tif")),
    width = 10,
    height = 8
  )
  highlight_files <- c(
    Purkinje_Lineage = "_Purkinje_Spatial_pt0.75.tif",
    Glia_Lineage = "_Glia_Spatial_pt0.75.tif",
    GABA_Lineage = "_GABA_Spatial_pt0.75.tif"
  )
  for (plot_name in names(highlight_files)) {
    save_cairo_tiff(
      highlight_plots[[plot_name]],
      file.path(sample_plot_dir, paste0(output_stem, highlight_files[[plot_name]])),
      width = 10,
      height = 8
    )
  }
  save_cairo_tiff(
    p_facet,
    file.path(sample_plot_dir, paste0(output_stem, "_Faceted_Clusters_pt0.6.tif")),
    width = 18,
    height = 12
  )
  save_cairo_tiff(
    p_other,
    file.path(sample_plot_dir, paste0(output_stem, "_Faceted_Other_pt0.6.tif")),
    width = 14,
    height = 8
  )
  save_cairo_tiff(
    p_dot_display,
    file.path(sample_plot_dir, paste0(output_stem, "_Markers_DotPlot.tif")),
    width = dotplot_dimensions[["width"]],
    height = dotplot_dimensions[["height"]],
    background = "white"
  )
  ggplot2::ggsave(
    file.path(sample_plot_dir, paste0(output_stem, "_Markers_DotPlot.pdf")),
    plot = p_dot_display,
    device = grDevices::cairo_pdf,
    width = dotplot_dimensions[["width"]],
    height = dotplot_dimensions[["height"]],
    bg = "white"
  )
  save_spatial_report_pdf(
    p_report,
    file.path(sample_plot_dir, paste0(output_stem, "_Spatial_Report.pdf"))
  )

  message("Finished sample: ", sample_id)
  rm(
    obj, p_global, p_global_no_other, p_facet, p_other, p_dot, p_dot_display, p_report,
    p_global_report, p_global_no_other_report, p_facet_report, spatial_data, stats
  )
  gc()
}

options(bitmapType = "cairo")
