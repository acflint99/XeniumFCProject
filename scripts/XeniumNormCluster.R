library(Seurat)
library(ggplot2)
library(patchwork)
library(randomcoloR)
library(future)
library(igraph)
library(dplyr)
library(here)

process_xenium_clusters <- function(xenium_obj, sample_name) {
  
  ## =========================================================
  ## 0. Setup
  ## =========================================================
  set.seed(42)
  plan("sequential")  # avoid future export issues
  
  DefaultAssay(xenium_obj) <- "Xenium"
  assay_name <- DefaultAssay(xenium_obj)
  message("Using assay: ", assay_name)
  
  ## =========================================================
  ## 1. Normalization & Preprocessing
  ## =========================================================
  scale_factor <- median(xenium_obj[[paste0("nCount_", assay_name)]][,1])
  
  xenium_obj <- NormalizeData(xenium_obj, scale.factor = scale_factor)
  xenium_obj <- FindVariableFeatures(xenium_obj, selection.method = "vst", nfeatures = 2000)
  xenium_obj <- ScaleData(xenium_obj)
  
  xenium_obj <- RunPCA(xenium_obj, npcs = 50, verbose = FALSE, seed.use = 42)
  
  # Remove scale.data to save memory
  LayerData(xenium_obj[[assay_name]], layer = "scale.data") <- NULL
  message("Removed scale.data layer from assay: ", assay_name)
  
  ## =========================================================
  ## 2. UMAP + Neighbor Graph
  ## =========================================================
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:50, seed.use = 42, verbose = FALSE)
  xenium_obj <- FindNeighbors(xenium_obj, dims = 1:50, nn.method = "annoy")
  
  graph_name <- paste0(assay_name, "_snn")
  if (!graph_name %in% names(xenium_obj@graphs)) {
    stop("SNN graph not found. Available graphs: ", paste(names(xenium_obj@graphs), collapse = ", "))
  }
  message("Using graph: ", graph_name)
  
  umap_plot <- DimPlot(xenium_obj, reduction = "umap") +
    ggtitle(paste0(sample_name, " - UMAP Projection (dims = 1:50)"))
  
  ## =========================================================
  ## 3. Optimal Resolution via Modularity
  ## =========================================================
  resolutions <- seq(0.2, 1.2, by = 0.2)
  mod_scores <- numeric(length(resolutions))
  
  snn_sparse <- xenium_obj@graphs[[graph_name]]
  subset_size <- min(5000, ncol(xenium_obj))
  set.seed(42)
  subset_cells <- sample(seq_len(ncol(xenium_obj)), subset_size)
  snn_sub <- snn_sparse[subset_cells, subset_cells]
  
  for (i in seq_along(resolutions)) {
    tmp_obj <- FindClusters(
      xenium_obj,
      resolution = resolutions[i],
      graph.name = graph_name,
      algorithm = 1,
      n.iter = 50,
      random.seed = 42,
      verbose = FALSE
    )
    
    g <- graph_from_adjacency_matrix(snn_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
    clusters_sub <- as.numeric(Idents(tmp_obj))[subset_cells]
    mod_scores[i] <- modularity(g, clusters_sub)
  }
  
  best_res <- resolutions[which.max(mod_scores)]
  message("Optimal clustering resolution (highest modularity): ", best_res)
  
  mod_plot <- ggplot(data.frame(resolution = resolutions, modularity = mod_scores),
                     aes(x = resolution, y = modularity)) +
    geom_line() + geom_point(size = 2) +
    geom_vline(xintercept = best_res, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = resolutions) +
    labs(title = paste(sample_name, "- Modularity vs Resolution"),
         x = "Resolution", y = "Modularity") +
    theme_minimal()
  
  ## =========================================================
  ## 4. Final Clustering
  ## =========================================================
  xenium_obj <- FindClusters(
    xenium_obj,
    resolution = best_res,
    graph.name = graph_name,
    algorithm = 1,
    n.iter = 100,
    random.seed = 42,
    verbose = FALSE
  )
  
  umap_cluster_plot <- DimPlot(xenium_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste(sample_name, "- UMAP Clusters (Res =", best_res, ")"))
  
  ## =========================================================
  ## 5. Export UMAP + Modularity PDF
  ## =========================================================
  pdf_file <- file.path(here("outputs"), paste0(sample_name, "_ClusterOptim_UMAP_Plots.pdf"))
  pdf(pdf_file, width = 8, height = 6)
  print(umap_plot)
  print(mod_plot)
  print(umap_cluster_plot)
  dev.off()
  message("Saved PDF plots: ", pdf_file)
  
  ## =========================================================
  ## 6. Spatial Cluster Plots
  ## =========================================================
  xenium_obj$seurat_clusters <- factor(Idents(xenium_obj))
  clusters <- levels(xenium_obj$seurat_clusters)
  n_clusters <- length(clusters)
  
  # Report number of clusters
  message("Number of clusters generated: ", n_clusters)
  
  cluster_colors <- distinctColorPalette(n_clusters)
  names(cluster_colors) <- clusters
  
  global_cluster_plot <- ImageDimPlot(
    object = xenium_obj,
    fov = "fov",
    group.by = "seurat_clusters",
    cols = cluster_colors,
    size = 0.75
  ) + scale_y_reverse() +
    ggtitle(paste(sample_name, "- Raw Clusters"))
  
  coords <- GetTissueCoordinates(xenium_obj)
  metadata <- xenium_obj@meta.data
  plot_data <- cbind(coords, cluster = as.character(metadata$seurat_clusters))
  
  facet_cluster_plot <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.2) +
    facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "black", color = NA),
          plot.background = element_rect(fill = "black", color = NA),
          strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
          plot.title = element_text(color = "white", hjust = 0.5, size = 14)) +
    ggtitle(paste(sample_name, "- Raw Clusters"))
  
  pdf_file_spatial <- file.path(here("outputs"), paste0(sample_name, "_RawClustersSpatialPlots.pdf"))
  pdf(pdf_file_spatial, width = 12, height = 8)
  print(global_cluster_plot)
  print(facet_cluster_plot)
  dev.off()
  message("Saved spatial plots PDF: ", pdf_file_spatial)
  
  ## =========================================================
  ## 7. Save Final Object
  ## =========================================================
  rds_file <- paste0(sample_name, "_CB_QC_cluster.rds")
  saveRDS(xenium_obj, rds_file)
  message("Saved processed Seurat object: ", rds_file)
  
  return(xenium_obj)
}