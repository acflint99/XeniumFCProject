# Display scales only: normalization, identities, marker sets and expression
# summaries remain controlled by the calling analysis. Scaled expression is
# gene-wise within each plotted dataset, not an absolute cross-sample measure.
broad_dotplot_col_min <- -2.5
broad_dotplot_col_max <- 2.5
broad_dotplot_dot_min <- 0
broad_dotplot_dot_max <- 100
broad_dotplot_dot_scale <- 6

standardize_workflow_dotplot <- function(plot, scaled = TRUE) {
  plot +
    ggplot2::scale_color_gradient(
      low = "lightgrey", high = "red",
      limits = c(broad_dotplot_col_min, broad_dotplot_col_max),
      breaks = c(broad_dotplot_col_min, 0, broad_dotplot_col_max),
      oob = scales::squish,
      name = if (scaled) "Average Expression (scaled)" else "Average Expression (log-normalized)"
    ) +
    ggplot2::scale_radius(
      range = c(0, broad_dotplot_dot_scale),
      limits = c(broad_dotplot_dot_min, broad_dotplot_dot_max),
      breaks = c(0, 25, 50, 75, 100),
      name = "Percent Expressed"
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, direction = "horizontal", nrow = 1,
                                  title.position = "top"),
      color = ggplot2::guide_colorbar(order = 2, direction = "horizontal",
                                     title.position = "top")
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 12, face = "italic", angle = 45,
                                         hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 12, face = "plain"),
      axis.title = ggplot2::element_text(size = 13),
      strip.text = ggplot2::element_text(size = 12, face = "plain"),
      strip.text.x = ggplot2::element_text(size = 12, face = "plain", angle = 90),
      panel.spacing.x = grid::unit(0.08, "in"),
      legend.position = "top", legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 12, face = "plain"),
      legend.text = ggplot2::element_text(size = 11),
      legend.spacing.x = grid::unit(0.12, "in")
    )
}

# Fix Seurat's percent limits before building the plot: otherwise it can cap
# percent values to each sample's observed maximum before ggplot sees them.
workflow_dotplot <- function(...) {
  args <- list(...)
  args$col.min <- broad_dotplot_col_min
  args$col.max <- broad_dotplot_col_max
  args$dot.min <- broad_dotplot_dot_min / 100
  args$dot.scale <- broad_dotplot_dot_scale
  args$scale.min <- broad_dotplot_dot_min
  args$scale.max <- broad_dotplot_dot_max
  # Add named marker groups with the supported ggplot2 facet API. This
  # preserves Seurat's marker ordering without its legacy facets= argument.
  marker_groups <- NULL
  if (is.list(args$features) || !is.null(names(args$features))) {
    groups <- args$features
    if (is.null(names(groups)) || anyNA(names(groups)) || any(!nzchar(names(groups)))) {
      stop("DotPlot marker groups must all have nonblank names.")
    }
    features <- unname(unlist(groups, use.names = FALSE))
    if (anyDuplicated(features)) stop("DotPlot marker genes must be unique across groups.")
    marker_groups <- stats::setNames(rep(names(groups), lengths(groups)), features)
    args$features <- features
  }
  plot <- do.call(Seurat::DotPlot, args)
  if (!is.null(marker_groups)) {
    plot$data$feature.groups <- factor(
      unname(marker_groups[as.character(plot$data$features.plot)]),
      levels = unique(unname(marker_groups))
    )
    plot <- plot + ggplot2::facet_grid(
      cols = ggplot2::vars(feature.groups), scales = "free_x", space = "free_x"
    )
  }
  standardize_workflow_dotplot(
    plot, scaled = !identical(args$scale, FALSE) && length(unique(plot$data$id)) > 1L
  )
}

# Inches: compact gene/identity spacing plus room for labels and top legends.
# Derive both TIFF and PDF dimensions from the same rendered plot data.
workflow_dotplot_dimensions <- function(plot) {
  feature_col <- if ("features.plot" %in% names(plot$data)) "features.plot" else "feature"
  identity_col <- if ("id" %in% names(plot$data)) "id" else "cell_type"
  c(width = max(9, 3 + 0.28 * length(unique(plot$data[[feature_col]]))),
    height = max(4, 2.8 + 0.28 * length(unique(plot$data[[identity_col]]))))
}
