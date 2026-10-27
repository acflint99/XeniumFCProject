# Refit/reuse residual representations, preserving all original metadata values.
# Post-QC assignments receive their own prefix; numeric cluster IDs are local
# to a fitted partition and must receive fresh biological review after refitting.
residual_postqc_representations <- function(object, config, no_removals) {
  ids <- colnames(object)
  metadata <- object[[]]
  resolutions <- c(0.5, 0.8, 1)
  columns <- paste0("Xenium_snn_res.", resolutions)
  post_columns <- paste0("residual_postQC_", columns)
  if (!length(ids) || anyDuplicated(ids) || !identical(rownames(metadata), ids) ||
      !all(c("orig.ident", columns) %in% names(metadata)) || any(post_columns %in% names(metadata))) {
    stop("Invalid original residual metadata or already processed input.")
  }
  dimensions <- as.integer(config$regional_subsets$integration_dimensions_default)
  neighbors <- as.integer(config$regional_subsets$neighbor_k_post_qc)
  seed <- as.integer(config$runtime$random_seed)
  set.seed(seed)
  Seurat::DefaultAssay(object) <- "Xenium"
  if (!no_removals) {
    object[["Xenium"]] <- SeuratObject::JoinLayers(object[["Xenium"]])
    object <- Seurat::NormalizeData(object, verbose = FALSE)
    object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    if (length(Seurat::VariableFeatures(object)) <= dimensions) stop("Too few variable genes for PCA.")
    object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
    object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = dimensions,
                           seed.use = seed, verbose = FALSE)
    object <- harmony::RunHarmony(object, group.by.vars = "orig.ident", dims.use = seq_len(dimensions),
                                  reduction.save = "harmony", verbose = TRUE)
    object <- Seurat::RunUMAP(object, reduction = "harmony", dims = seq_len(dimensions),
      reduction.name = "umap_harmony", n.neighbors = neighbors, seed.use = seed, verbose = FALSE)
    object <- Seurat::FindNeighbors(object, reduction = "harmony", dims = seq_len(dimensions),
                                    k.param = neighbors, verbose = FALSE)
    # Preserve the existing residual FindClusters defaults (including its seed).
    object <- Seurat::FindClusters(object, resolution = resolutions, verbose = FALSE)
  }
  if (!all(c("pca", "harmony", "umap_harmony") %in% Seurat::Reductions(object)) ||
      !"Xenium_snn" %in% SeuratObject::Graphs(object)) stop("Incomplete post-QC representations/graph.")
  for (i in seq_along(columns)) {
    object[[post_columns[i]]] <- object[[columns[i], drop = TRUE]]
    object[[columns[i]]] <- metadata[[columns[i]]]
  }
  if ("seurat_clusters" %in% names(metadata)) object$seurat_clusters <- metadata$seurat_clusters
  Seurat::Idents(object) <- post_columns[2]
  if (!identical(colnames(object), ids) || !identical(object[[]][, names(metadata), drop = FALSE], metadata)) {
    stop("Post-QC processing changed cell IDs/order or original metadata.")
  }
  for (reduction in c("pca", "harmony", "umap_harmony")) {
    coordinates <- SeuratObject::Embeddings(object[[reduction]])
    if (!identical(rownames(coordinates), ids) || any(!is.finite(coordinates))) stop("Invalid post-QC coordinates.")
  }
  graph <- object[["Xenium_snn"]]
  if (!identical(rownames(graph), ids) || !identical(colnames(graph), ids) || any(!is.finite(graph@x))) {
    stop("Post-QC graph is not aligned with retained cells.")
  }
  object
}
