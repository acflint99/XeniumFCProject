# Per-sample RL spatial and marker reports.
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
  output_stem <- paste0(sample_id, "_RL")

  obj <- readRDS(input_rds_path)
  required_metadata <- c("RL_subcluster", "consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop(sample_id, " lacks required RL plotting metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!"Xenium" %in% names(obj@assays)) stop(sample_id, " lacks the Xenium assay.")
  DefaultAssay(obj) <- "Xenium"

  observed_rl_labels <- unique(as.character(obj$RL_subcluster))
  observed_rl_labels <- observed_rl_labels[!is.na(observed_rl_labels) & nzchar(observed_rl_labels)]
  unsupported_rl_labels <- setdiff(observed_rl_labels, rl_subcluster_order)
  missing_rl_colors <- setdiff(observed_rl_labels, names(rl_palette))
  if (length(unsupported_rl_labels) || length(missing_rl_colors)) {
    stop(
      sample_id, " has RL labels not fully registered in the shared taxonomy. Missing from order: ",
      if (length(unsupported_rl_labels)) paste(unsupported_rl_labels, collapse = ", ") else "none",
      "; missing from palette: ",
      if (length(missing_rl_colors)) paste(missing_rl_colors, collapse = ", ") else "none"
    )
  }
  sample_rl_order <- rl_subcluster_order[rl_subcluster_order %in% observed_rl_labels]
  if (!length(sample_rl_order)) stop(sample_id, " contains no mapped RL subclusters.")
  obj$RL_subcluster <- factor(as.character(obj$RL_subcluster), levels = sample_rl_order)

  stats <- obj@meta.data %>%
    group_by(RL_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    mutate(percentage = (cell_count / sum(cell_count)) * 100)

  spatial_data <- prepare_spatial_plot_data(
    object = obj,
    cluster_column = "RL_subcluster",
    cluster_order = sample_rl_order,
    sample_id = sample_id
  )

  p_global <- build_layered_spatial_plot(
    spatial_data,
    rl_palette,
    paste(sample_id, "RL Subclusters"),
    targets = sample_rl_order
  )
  rl_lineage_clusters <- c(
    "RL SVZ", "Immature UBCs", "Mature UBCs", "Prolif GCPs", "Maturing GCPs",
    "Differentiating GCs", "Migrating GCs", "Mature GCs"
  )
  p_global_no_other <- build_layered_spatial_plot(
    spatial_data,
    rl_palette,
    paste(sample_id, "RL Lineages"),
    targets = rl_lineage_clusters
  )
  highlight_definitions <- list(
    RL_Lineage = "RL SVZ",
    Granule_Lineage = c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
    `Granule&RL_Lineage` = c("RL SVZ", "Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
    UBC_Lineage = c("RL SVZ", "Immature UBCs", "Mature UBCs")
  )
  highlight_plots <- lapply(names(highlight_definitions), function(plot_name) {
    build_layered_spatial_plot(
      spatial_data,
      rl_palette,
      paste(sample_id, plot_name),
      targets = highlight_definitions[[plot_name]]
    )
  })
  names(highlight_plots) <- names(highlight_definitions)

  lineage_rows <- list(
    c("RL SVZ", "Immature UBCs", "Mature UBCs"),
    c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs")
  )
  other_rows <- list(c("BG", "BG Progenitors", "Maturing PCs", "eCN", "GABA Progenitors", "Microglia", "Endothelial"))
  p_facet <- build_faceted_spatial_rows(spatial_data, rl_palette, lineage_rows)
  p_other <- build_faceted_spatial_rows(spatial_data, rl_palette, other_rows)

  # Pre-rasterize only the spatial panels embedded in the PDF report. The
  # standalone TIFFs use direct points so their nominal sizes remain exact.
  p_global_report <- build_layered_spatial_plot(
    spatial_data,
    rl_palette,
    paste(sample_id, "RL Subclusters"),
    targets = sample_rl_order,
    rasterize = TRUE
  )
  p_global_no_other_report <- build_layered_spatial_plot(
    spatial_data,
    rl_palette,
    paste(sample_id, "RL Lineages"),
    targets = rl_lineage_clusters,
    rasterize = TRUE
  )
  p_facet_report <- build_faceted_spatial_rows(
    spatial_data,
    rl_palette,
    lineage_rows,
    rasterize = TRUE
  )

  plotted_cells <- Cells(obj)[!is.na(obj$RL_subcluster)]
  if (!length(plotted_cells)) stop(sample_id, " contains no mapped RL cells.")
  dotplot_identity <- as.character(obj$RL_subcluster)
  dotplot_identity[is.na(dotplot_identity)] <- ".not_plotted"
  obj$subcluster_rev <- factor(
    dotplot_identity,
    levels = c(rev(sample_rl_order), ".not_plotted")
  )
  Idents(obj) <- "subcluster_rev"
  xenium_features <- rownames(obj[["Xenium"]])
  existing_rl_markers <- rl_production_marker_panel(sample_rl_order)
  missing_rl_markers <- setdiff(unlist(existing_rl_markers, use.names = FALSE), xenium_features)
  if (length(missing_rl_markers)) {
    stop(
      sample_id, " Xenium assay lacks requested RL marker gene(s): ",
      paste(missing_rl_markers, collapse = ", ")
    )
  }

  p_dot <- workflow_dotplot(
    obj,
    features = existing_rl_markers,
    assay = "Xenium",
    idents = rev(sample_rl_order),
    cols = c("lightgrey", "red"),
    dot.scale = 6,
    cluster.idents = FALSE,
    scale = length(sample_rl_order) > 1L
  ) +
    RotatedAxis() +
    theme(
      axis.text.x = element_text(size = 12, face = "italic"),
      axis.text.y = element_text(size = 12, face = "plain")
    ) +
    ggtitle(paste(sample_id, "RL Markers"))

  dotplot_dimensions <- workflow_dotplot_dimensions(p_dot)
  # A wrapped grob keeps the broad-colored strips in TIFF, PDF, and the
  # patchwork report while remaining printable by the existing save helper.
  p_dot <- patchwork::wrap_elements(full = workflow_dotplot_grob(p_dot))

  p_report <- build_spatial_report_page(
    p_global_report,
    p_global_no_other_report,
    p_facet_report,
    p_dot,
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
    RL_Lineage = "_Lineage_Spatial_pt0.75.tif",
    Granule_Lineage = "_Granule_Spatial_pt0.75.tif",
    `Granule&RL_Lineage` = "_Granule&RL_Spatial_pt0.75.tif",
    UBC_Lineage = "_UBC_Spatial_pt0.75.tif"
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
    height = 10
  )
  save_cairo_tiff(
    p_other,
    file.path(sample_plot_dir, paste0(output_stem, "_Faceted_Other_pt0.6.tif")),
    width = 18,
    height = 6
  )
  save_cairo_tiff(
    p_dot,
    file.path(sample_plot_dir, paste0(output_stem, "_Markers_DotPlot.tif")),
    width = dotplot_dimensions[["width"]],
    height = dotplot_dimensions[["height"]],
    background = "white"
  )
  ggplot2::ggsave(
    file.path(sample_plot_dir, paste0(output_stem, "_Markers_DotPlot.pdf")),
    plot = p_dot,
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
    obj, p_global, p_global_no_other, p_facet, p_other, p_dot, p_report,
    p_global_report, p_global_no_other_report, p_facet_report, spatial_data, stats
  )
  gc()
}

options(bitmapType = "cairo")
