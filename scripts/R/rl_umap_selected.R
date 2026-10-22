# Render every cell from the user-selected, previously fitted comparison.
# Load the approved object for its original DimPlot format; never write it.
render_selected_rl_umap <- function(config, input, comparison_dir,
                                    dry_run = FALSE, overwrite = FALSE, object = NULL) {
  selected <- "nn30_md0.05"
  coordinates <- file.path(comparison_dir, paste0("Xenium_RL_", selected, "_coordinates.rds"))
  source_manifest <- file.path(comparison_dir, "Xenium_RL_UMAP_Comparison_manifest.csv")
  base <- here::here(config$project$outputs_dir, "xenium", "rl", "06_subclusters")
  plot_dir <- file.path(base, "plots")
  table_dir <- file.path(base, "tables")
  full_tiff <- file.path(plot_dir, "Xenium_RL_Subcluster_UMAP.tif")
  filtered_tiff <- file.path(plot_dir, "Xenium_RL_Subcluster_rmBG,Ep,PC_UMAP.tif")
  pdf_path <- file.path(plot_dir, "Xenium_RL_Subcluster_UMAP_AllCells.pdf")
  receipt <- file.path(table_dir, "Xenium_RL_Selected_UMAP_manifest.csv")
  inputs <- c(input, coordinates, source_manifest)
  outputs <- c(full_tiff, filtered_tiff, pdf_path, receipt)
  if (dry_run) {
    ok <- compact_dry_run("Selected RL UMAP all-cell rendering (paths only)", inputs, outputs)
    if (!ok) stop("Selected UMAP preflight failed; no objects were loaded.")
    return(invisible(outputs))
  }
  if (!all(file.exists(inputs))) stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = "; "))
  if (any(file.exists(outputs)) && !overwrite) {
    stop("Existing selected UMAP plots/manifests are protected. Use --overwrite after review.")
  }
  m <- read.csv(source_manifest, stringsAsFactors = FALSE)
  required <- c("comparison", "n_neighbors", "min_dist", "spread", "metric", "seed",
                "input", "input_bytes", "input_mtime", "n_cells", "coordinates", "dimensions")
  if (!all(required %in% names(m))) stop("Incomplete comparison provenance.")
  m <- m[m$comparison == selected & !is.na(m$comparison), , drop = FALSE]
  info <- file.info(input)
  if (nrow(m) != 1L || anyNA(m[, required]) ||
      m$n_neighbors != 30L || m$min_dist != 0.05 || m$spread != 1 || m$metric != "cosine" ||
      m$seed != config$runtime$random_seed ||
      m$dimensions != paste(seq_len(config$regional_subsets$post_qc_dimensions_observed), collapse = "|") ||
      !identical(normalizePath(m$input, mustWork = TRUE), normalizePath(input, mustWork = TRUE)) ||
      !identical(normalizePath(m$coordinates, mustWork = TRUE), normalizePath(coordinates, mustWork = TRUE)) ||
      m$input_bytes != info$size || m$input_mtime != as.character(info$mtime)) {
    stop("Selected comparison parameters or source-object provenance no longer match. Rebuild/review comparisons.")
  }
  d <- readRDS(coordinates)
  if (!is.data.frame(d) || !all(c("cell_id", "RL_subcluster", "UMAP_1", "UMAP_2") %in% names(d))) {
    stop("Malformed coordinate table.")
  }
  if (!nrow(d) || nrow(d) != m$n_cells || anyNA(d$cell_id) ||
      any(!nzchar(d$cell_id)) || anyDuplicated(d$cell_id) ||
      anyNA(d$RL_subcluster) || any(!nzchar(d$RL_subcluster)) ||
      !is.numeric(d$UMAP_1) || !is.numeric(d$UMAP_2) ||
      any(!is.finite(d$UMAP_1)) || any(!is.finite(d$UMAP_2))) {
    stop("Invalid/missing/duplicate cells, labels, or coordinates.")
  }
  if (!requireNamespace("Cairo", quietly = TRUE)) stop("Missing project plotting dependency: Cairo.")
  suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
  source(here::here("scripts", "color_palette.R"), local = TRUE)
  if (any(!d$RL_subcluster %in% rl_subcluster_order)) stop("Unregistered RL labels.")
  obj <- if (is.null(object)) readRDS(input) else object
  if (!inherits(obj, "Seurat") || !"RL_subcluster" %in% colnames(obj[[]]) ||
      !"umap_clean" %in% Reductions(obj) ||
      !identical(colnames(obj), d$cell_id) ||
      !identical(as.character(obj$RL_subcluster), as.character(d$RL_subcluster))) {
    stop("Selected coordinates do not exactly match the approved object's cell IDs/order and labels.")
  }
  # Replace coordinates in memory only, retaining the original reduction key.
  coords <- as.matrix(d[, c("UMAP_1", "UMAP_2")])
  rownames(coords) <- d$cell_id
  colnames(coords) <- colnames(Embeddings(obj, "umap_clean"))
  obj[["umap_clean"]] <- CreateDimReducObject(
    embeddings = coords, key = Key(obj[["umap_clean"]]),
    assay = DefaultAssay(obj[["umap_clean"]])
  )
  active_order <- rl_subcluster_order[rl_subcluster_order %in% d$RL_subcluster]
  obj$RL_subcluster <- factor(as.character(obj$RL_subcluster), levels = active_order)
  Idents(obj) <- "RL_subcluster"
  plot_cells <- function(x, title) {
    xy <- Embeddings(x, "umap_clean")
    centers <- stats::aggregate(data.frame(label_x = xy[, 1], label_y = xy[, 2]),
                               list(RL_subcluster = as.character(x$RL_subcluster)), stats::median)
    # Match Stage 6 exactly: default point size/order/rasterization, boxed
    # median labels, original axis names, title, theme and legend text size.
    DimPlot(x, reduction = "umap_clean", label = FALSE, cols = rl_palette) +
      geom_label(data = centers,
        aes(x = label_x, y = label_y, label = RL_subcluster),
        inherit.aes = FALSE, size = 4, fill = "white", color = "black",
        label.padding = grid::unit(0.15, "lines"),
        label.r = grid::unit(0.12, "lines"), linewidth = 0.25) +
      ggtitle(title) + theme(legend.text = element_text(size = 8))
  }
  full <- plot_cells(obj, "Xenium Merged RL Subclusters")
  excluded <- c("BG", "BG Progenitors", "Maturing PCs", "eCN", "GABA Progenitors", "Microglia", "Endothelial")
  retained_ids <- d$cell_id[!d$RL_subcluster %in% excluded]
  if (!length(retained_ids)) stop("No cells remain for the established lineage-only display.")
  lineage <- subset(obj, cells = retained_ids)
  filtered <- plot_cells(lineage, "Xenium Merged RL Subclusters (Filtered)")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  save_tiff <- function(p, path) {
    Cairo::CairoTIFF(path, width = 12, height = 9, units = "in", res = 600)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  save_tiff(full, full_tiff)
  save_tiff(filtered, filtered_tiff)
  save_pdf <- function() {
    grDevices::cairo_pdf(pdf_path, width = 12, height = 9)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(full)
    print(filtered)
  }
  save_pdf()
  write.csv(data.frame(m, selected_on = as.character(Sys.time()),
    all_cells_displayed = nrow(d), lineage_cells_displayed = ncol(lineage),
    full_tiff = full_tiff, filtered_tiff = filtered_tiff, pdf = pdf_path), receipt, row.names = FALSE)
  message("Rendered all ", nrow(d), " cells with saved nn30/min.dist=0.05 coordinates: ", full_tiff)
  invisible(outputs)
}
