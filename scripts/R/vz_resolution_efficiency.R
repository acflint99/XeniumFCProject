# Runtime-only helpers. Never write an analysis object or change cell membership.
vz_comparison_slim <- function(object) {
  ids <- colnames(object)
  metadata <- object[[]]
  coords <- SeuratObject::Embeddings(object[["umap_clean"]])
  SeuratObject::DefaultAssay(object) <- "Xenium"
  # JoinLayers is performed by the driver before this step; retain both layers
  # unchanged for the existing marker and DotPlot calculations.
  if (!all(c("counts", "data") %in% SeuratObject::Layers(object[["Xenium"]]))) {
    stop("Slimming requires joined Xenium counts and normalized data layers.")
  }
  object <- Seurat::DietSeurat(object, assays = "Xenium", layers = c("counts", "data"),
    dimreducs = "umap_clean", graphs = NULL, misc = FALSE)
  if (!identical(colnames(object), ids) || !identical(object[[]], metadata) ||
      !identical(SeuratObject::Embeddings(object[["umap_clean"]]), coords) ||
      length(object@graphs) != 0L) stop("Slimming changed cells, metadata, UMAP or retained graphs.")
  object
}

vz_comparison_validate_checkpoint <- function(checkpoint, key, columns, anchors) {
  if (!is.list(checkpoint) || !identical(checkpoint$key, key) ||
      !is.data.frame(checkpoint$assignments) ||
      !all(c("cell_id", "sample_id", columns) %in% names(checkpoint$assignments))) {
    stop("Partition checkpoint provenance/schema differs from this input or implementation.")
  }
  assignments <- checkpoint$assignments
  if (!identical(assignments$cell_id, key$cell_ids) ||
      !identical(assignments$sample_id, key$sample_ids) || anyDuplicated(assignments$cell_id) ||
      !identical(assignments[[columns[1]]], as.character(anchors[[1]])) ||
      !identical(assignments[[columns[4]]], as.character(anchors[[2]]))) {
    stop("Partition checkpoint changed IDs/order, sample mapping, or anchor assignments.")
  }
  for (column in columns) {
    values <- suppressWarnings(as.numeric(assignments[[column]]))
    if (anyNA(values) || any(!is.finite(values)) || any(values != floor(values)) ||
        length(unique(values)) < 2L) stop("Invalid checkpoint partition: ", column)
  }
  invisible(TRUE)
}
