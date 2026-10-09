# Shared plotting helpers for the RL and VZ per-sample spatial reports.

spatial_point_layer <- function(
    data,
    mapping,
    point_size,
    color = NULL,
    rasterize = FALSE,
    raster_dpi = 600,
    raster_pixels = c(1800, 1800)) {
  layer_args <- list(
    data = data,
    mapping = mapping,
    inherit.aes = FALSE
  )
  if (!is.null(color)) layer_args$color <- color

  # The standalone spatial TIFFs are already raster devices. Drawing their
  # points directly preserves ggplot2's requested size at the final 600-dpi
  # resolution and avoids enlarging a lower-resolution raster layer. Only the
  # combined PDF report needs pre-rasterized layers to keep file size bounded.
  if (!rasterize) {
    return(do.call(
      ggplot2::geom_point,
      c(layer_args, list(size = point_size, shape = 16, stroke = 0))
    ))
  }

  # scattermore is an imported Seurat dependency and is therefore the most
  # reliable raster backend on the production R environment.
  if (requireNamespace("scattermore", quietly = TRUE)) {
    return(do.call(
      scattermore::geom_scattermore,
      c(
        layer_args,
        list(
          pointsize = max(1L, round(point_size * 5)),
          pixels = raster_pixels
        )
      )
    ))
  }
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    return(do.call(
      ggrastr::geom_point_rast,
      c(layer_args, list(size = point_size, raster.dpi = raster_dpi))
    ))
  }
  stop(
    "A raster spatial backend is required. Install the project-locked ",
    "Seurat dependency 'scattermore' or make 'ggrastr' available through renv."
  )
}

prepare_spatial_plot_data <- function(object, cluster_column, cluster_order, sample_id) {
  coords <- Seurat::GetTissueCoordinates(object, type = "Xenium")
  required_columns <- c("cell", "x", "y")
  missing_columns <- setdiff(required_columns, colnames(coords))
  if (length(missing_columns)) {
    stop(
      sample_id, " spatial coordinates lack required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  coordinate_cells <- match(as.character(coords$cell), colnames(object))
  if (anyNA(coordinate_cells)) {
    stop(sample_id, " spatial coordinates do not align one-to-one with object cells.")
  }
  if (anyDuplicated(coords$cell)) {
    stop(sample_id, " spatial coordinates contain duplicate cell identifiers.")
  }

  data.frame(
    cell = as.character(coords$cell),
    plot_x = coords$y,
    plot_y = coords$x,
    cluster = factor(
      as.character(object[[cluster_column]][coordinate_cells, 1]),
      levels = cluster_order
    ),
    stringsAsFactors = FALSE
  )
}

spatial_black_theme <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = "white"),
      legend.background = ggplot2::element_rect(fill = "black", color = NA),
      legend.key = ggplot2::element_rect(fill = "black", color = NA),
      panel.background = ggplot2::element_rect(fill = "black", color = "black"),
      plot.background = ggplot2::element_rect(fill = "black", color = "black"),
      plot.title = ggplot2::element_text(color = "white", hjust = 0.5),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
    )
}

build_layered_spatial_plot <- function(
    spatial_data,
    palette,
    title,
    targets = names(palette),
    background_color = "gray50",
    background_size = 0.5,
    foreground_size = 0.75,
    rasterize = FALSE) {
  observed_targets <- targets[targets %in% unique(as.character(spatial_data$cluster))]
  foreground <- spatial_data[
    !is.na(spatial_data$cluster) & spatial_data$cluster %in% observed_targets,
    , drop = FALSE
  ]
  background <- spatial_data[
    is.na(spatial_data$cluster) | !spatial_data$cluster %in% observed_targets,
    , drop = FALSE
  ]
  if (!nrow(foreground)) stop("No cells were available for spatial plot: ", title)

  ggplot2::ggplot() +
    spatial_point_layer(
      data = background,
      mapping = ggplot2::aes(x = plot_x, y = plot_y),
      color = background_color,
      point_size = background_size,
      rasterize = rasterize
    ) +
    spatial_point_layer(
      data = foreground,
      mapping = ggplot2::aes(x = plot_x, y = plot_y, color = cluster),
      point_size = foreground_size,
      rasterize = rasterize
    ) +
    ggplot2::scale_color_manual(
      values = palette[observed_targets],
      limits = observed_targets,
      drop = FALSE
    ) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::ggtitle(title) +
    spatial_black_theme() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 2.3, alpha = 1))
    )
}

build_spatial_cluster_panel <- function(
    spatial_data,
    cluster_name,
    palette,
    x_limits,
    y_limits,
    foreground_size = 0.6,
    rasterize = FALSE) {
  foreground <- spatial_data[
    !is.na(spatial_data$cluster) & spatial_data$cluster == cluster_name,
    , drop = FALSE
  ]

  plot <- ggplot2::ggplot()
  if (nrow(foreground)) {
    plot <- plot + spatial_point_layer(
      data = foreground,
      mapping = ggplot2::aes(x = plot_x, y = plot_y),
      color = unname(palette[[cluster_name]]),
      point_size = foreground_size,
      rasterize = rasterize
    )
  }

  plot +
    ggplot2::coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    ggplot2::ggtitle(cluster_name) +
    spatial_black_theme() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(color = "white", hjust = 0.5, size = 11)
    )
}

build_faceted_spatial_rows <- function(
    spatial_data,
    palette,
    row_groups,
    rasterize = FALSE) {
  requested_clusters <- unlist(row_groups, use.names = FALSE)
  missing_colors <- setdiff(requested_clusters, names(palette))
  if (length(missing_colors)) {
    stop("Faceted spatial labels are missing palette colors: ", paste(missing_colors, collapse = ", "))
  }
  if (!length(requested_clusters) || anyDuplicated(requested_clusters)) {
    stop("Faceted spatial row groups must contain unique cluster labels.")
  }

  x_limits <- range(spatial_data$plot_x, finite = TRUE)
  y_limits <- range(spatial_data$plot_y, finite = TRUE)
  row_plots <- lapply(row_groups, function(cluster_row) {
    panels <- lapply(cluster_row, function(cluster_name) {
      build_spatial_cluster_panel(
        spatial_data = spatial_data,
        cluster_name = cluster_name,
        palette = palette,
        x_limits = x_limits,
        y_limits = y_limits,
        rasterize = rasterize
      )
    })
    patchwork::wrap_plots(panels, ncol = length(cluster_row))
  })
  patchwork::wrap_plots(row_plots, ncol = 1) &
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "black", color = "black"))
}

save_cairo_tiff <- function(plot, filename, width, height, background = "black") {
  Cairo::CairoTIFF(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = 600,
    bg = background
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(filename)
}

build_spatial_report_page <- function(
    global_plot,
    global_no_other_plot,
    faceted_plot,
    marker_plot,
    sample_id) {
  report_design <- "
AABB
CCCC
DDDD
"
  patchwork::wrap_plots(
    A = global_plot,
    B = global_no_other_plot,
    C = faceted_plot,
    D = marker_plot,
    design = report_design,
    heights = c(1, 1.25, 0.75)
  ) +
    patchwork::plot_annotation(
      title = paste(sample_id, "spatial and marker report"),
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.title = ggplot2::element_text(color = "black", hjust = 0.5)
      )
    )
}

save_spatial_report_pdf <- function(plot, filename) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = grDevices::cairo_pdf,
    width = 24,
    height = 18,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

pdf_page_count <- function(path, pdfinfo_command = Sys.which("pdfinfo")) {
  if (!nzchar(pdfinfo_command)) stop("The Poppler command 'pdfinfo' is unavailable on PATH.")
  info_output <- suppressWarnings(system2(
    pdfinfo_command,
    args = shQuote(path),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(info_output, "status")
  if (!is.null(status) && status != 0L) {
    stop("pdfinfo failed for ", path, ":\n", paste(info_output, collapse = "\n"))
  }
  page_line <- grep("^Pages:[[:space:]]+[0-9]+", info_output, value = TRUE)
  if (length(page_line) != 1L) stop("Could not read a unique page count from: ", path)
  as.integer(sub("^Pages:[[:space:]]+", "", page_line))
}

build_stage8_pcw_page_map <- function(config, samples) {
  required_manifest_columns <- "sample_id"
  if (!all(required_manifest_columns %in% names(samples))) {
    stop("Sample manifest must contain: ", paste(required_manifest_columns, collapse = ", "))
  }
  if (nrow(samples) != 34L || anyDuplicated(samples$sample_id)) {
    stop(
      "The current Stage-8 combined report requires exactly 34 unique samples ",
      "from config/samples.csv; observed ", nrow(samples), "."
    )
  }

  metadata_path <- resolve_config_path(config$manifests$sample_metadata, config)
  if (!file.exists(metadata_path)) stop("Sample metadata not found: ", metadata_path)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "The 'readxl' package is required to order Stage-8 report pages by PCW. ",
      "Restore it through the project renv environment before running this report."
    )
  }
  sample_metadata <- readxl::read_excel(metadata_path)
  required_metadata_columns <- c("sample", "PCW")
  if (!all(required_metadata_columns %in% names(sample_metadata))) {
    stop("Sample metadata must contain: ", paste(required_metadata_columns, collapse = ", "))
  }

  manifest_order <- seq_len(nrow(samples))
  base_sample_ids <- sub("_\\d+$", "", samples$sample_id)
  metadata_matches <- lapply(base_sample_ids, function(sample_id) {
    which(as.character(sample_metadata$sample) == sample_id)
  })
  bad_match_counts <- which(lengths(metadata_matches) != 1L)
  if (length(bad_match_counts)) {
    stop(
      "Expected exactly one PCW metadata match for each biological sample; failures:\n- ",
      paste0(
        samples$sample_id[bad_match_counts], " (base ID ",
        base_sample_ids[bad_match_counts], "): ",
        lengths(metadata_matches)[bad_match_counts], " matches",
        collapse = "\n- "
      )
    )
  }

  pcw_values <- vapply(
    metadata_matches,
    function(i) trimws(as.character(sample_metadata$PCW[[i]])),
    character(1)
  )
  if (anyNA(pcw_values) || any(!nzchar(pcw_values))) {
    stop(
      "Blank or missing PCW metadata for: ",
      paste(samples$sample_id[is.na(pcw_values) | !nzchar(pcw_values)], collapse = ", ")
    )
  }
  pcw_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pcw_values)))
  if (any(!is.finite(pcw_numeric))) {
    stop(
      "Could not parse numeric PCW values for: ",
      paste(samples$sample_id[!is.finite(pcw_numeric)], collapse = ", ")
    )
  }

  page_order <- order(pcw_numeric, manifest_order)
  data.frame(
    page = seq_len(nrow(samples)),
    sample_id = samples$sample_id[page_order],
    PCW = pcw_values[page_order],
    PCW_numeric = pcw_numeric[page_order],
    manifest_order = manifest_order[page_order],
    stringsAsFactors = FALSE
  )
}

combine_stage8_sample_reports <- function(
    page_map,
    report_root,
    branch,
    dry_run = FALSE,
    overwrite = FALSE) {
  branch <- toupper(branch)
  if (!branch %in% c("RL", "VZ")) stop("Combined-report branch must be RL or VZ.")
  required_page_columns <- c("page", "sample_id", "PCW", "PCW_numeric", "manifest_order")
  missing_page_columns <- setdiff(required_page_columns, names(page_map))
  if (length(missing_page_columns)) {
    stop("Stage-8 page map is missing: ", paste(missing_page_columns, collapse = ", "))
  }
  if (nrow(page_map) != 34L || anyDuplicated(page_map$sample_id)) {
    stop(
      "The current Stage-8 combined report requires exactly 34 unique samples ",
      "from config/samples.csv; observed ", nrow(page_map), "."
    )
  }
  if (!identical(page_map$page, seq_len(nrow(page_map))) ||
      any(!is.finite(page_map$PCW_numeric)) ||
      is.unsorted(page_map$PCW_numeric)) {
    stop("Stage-8 page map is not in valid ascending numeric PCW order.")
  }
  sample_ids <- as.character(page_map$sample_id)
  if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

  input_paths <- file.path(
    report_root,
    sample_ids,
    paste0(sample_ids, "_", branch, "_Spatial_Report.pdf")
  )
  output_pdf <- file.path(
    report_root,
    paste0("Xenium_", branch, "_AllSamples_Spatial_Report.pdf")
  )
  output_manifest <- file.path(
    report_root,
    paste0("Xenium_", branch, "_AllSamples_Spatial_Report_manifest.csv")
  )

  pdf_tools <- Sys.which(c("pdfinfo", "pdfunite"))
  missing_tools <- names(pdf_tools)[!nzchar(pdf_tools)]
  if (length(missing_tools)) {
    stop(
      "Required Poppler command(s) are unavailable: ",
      paste(missing_tools, collapse = ", "),
      ". Load the Cheaha Poppler module or make these commands available on PATH; ",
      "do not install packages inside the production job."
    )
  }

  missing_inputs <- input_paths[!file.exists(input_paths)]
  if (length(missing_inputs)) {
    stop(
      branch, " combined report is missing ", length(missing_inputs),
      " manifest-defined input PDF(s):\n- ",
      paste(missing_inputs, collapse = "\n- ")
    )
  }
  input_info <- file.info(input_paths)
  invalid_inputs <- input_paths[is.na(input_info$size) | input_info$size < 10000L]
  if (length(invalid_inputs)) {
    stop(
      branch, " combined report has empty or implausibly small input PDF(s):\n- ",
      paste(invalid_inputs, collapse = "\n- ")
    )
  }

  input_pages <- vapply(
    input_paths,
    pdf_page_count,
    integer(1),
    pdfinfo_command = pdf_tools[["pdfinfo"]]
  )
  non_single_page <- input_paths[input_pages != 1L]
  if (length(non_single_page)) {
    stop(
      branch, " Stage-8 inputs must each contain exactly one page:\n- ",
      paste(
        paste0(non_single_page, " (", input_pages[input_pages != 1L], " pages)"),
        collapse = "\n- "
      )
    )
  }

  manifest <- data.frame(
    page_map[, required_page_columns, drop = FALSE],
    input_pdf = normalizePath(input_paths, mustWork = TRUE),
    input_bytes = input_info$size,
    input_modified = format(input_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )

  if (dry_run) {
    cat(
      "DRY-RUN PASS | ", branch, " combined Stage-8 report | inputs=34/34 | ",
      "pages=34 | order=numeric PCW, then config/samples.csv | output=", output_pdf, "\n",
      sep = ""
    )
    write.table(
      page_map[, required_page_columns, drop = FALSE],
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    return(invisible(list(output_pdf = output_pdf, manifest = manifest)))
  }

  existing_outputs <- c(output_pdf, output_manifest)[file.exists(c(output_pdf, output_manifest))]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite existing combined report output(s):\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing them."
    )
  }

  dir.create(report_root, recursive = TRUE, showWarnings = FALSE)
  temporary_pdf <- tempfile(
    pattern = paste0(".", tolower(branch), "_combined_"),
    tmpdir = report_root,
    fileext = ".pdf"
  )
  temporary_manifest <- tempfile(
    pattern = paste0(".", tolower(branch), "_combined_manifest_"),
    tmpdir = report_root,
    fileext = ".csv"
  )
  on.exit(unlink(c(temporary_pdf, temporary_manifest), force = TRUE), add = TRUE)

  merge_output <- suppressWarnings(system2(
    pdf_tools[["pdfunite"]],
    args = c(shQuote(input_paths), shQuote(temporary_pdf)),
    stdout = TRUE,
    stderr = TRUE
  ))
  merge_status <- attr(merge_output, "status")
  if (!is.null(merge_status) && merge_status != 0L) {
    stop("pdfunite failed:\n", paste(merge_output, collapse = "\n"))
  }
  if (!file.exists(temporary_pdf) || file.info(temporary_pdf)$size < 10000L) {
    stop("pdfunite did not create a valid temporary combined PDF.")
  }

  combined_pages <- pdf_page_count(temporary_pdf, pdf_tools[["pdfinfo"]])
  if (combined_pages != length(sample_ids)) {
    stop(
      "Combined ", branch, " report has ", combined_pages,
      " pages; expected ", length(sample_ids), "."
    )
  }
  manifest$output_pdf <- normalizePath(output_pdf, mustWork = FALSE)
  manifest$combined_page_count <- combined_pages
  write.csv(manifest, temporary_manifest, row.names = FALSE)

  if (!file.rename(temporary_pdf, output_pdf)) {
    stop("Could not atomically place the combined PDF at: ", output_pdf)
  }
  if (!file.rename(temporary_manifest, output_manifest)) {
    stop("Combined PDF was created, but its manifest could not be placed at: ", output_manifest)
  }

  final_pages <- pdf_page_count(output_pdf, pdf_tools[["pdfinfo"]])
  if (final_pages != 34L) stop("Final combined PDF failed its 34-page validation.")
  message("Created 34-page ", branch, " combined Stage-8 report: ", output_pdf)
  message("Created page manifest: ", output_manifest)
  invisible(list(output_pdf = output_pdf, output_manifest = output_manifest))
}
