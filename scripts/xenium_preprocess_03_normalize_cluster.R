# =========================================================
# Required Libraries
# =========================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(randomcoloR)
library(future)
library(igraph)
library(dplyr)
library(here)
library(Cairo) # Required for reliable headless cluster TIFF writing

# =========================================================
# Xenium Clustering Function
# =========================================================
process_xenium_clusters <- function(xenium_obj, sample_name) {
  
  ## =========================================================
  ## 0. Setup
  ## =========================================================
  set.seed(42)
  plan("sequential")  # avoid future export issues within the parallel workers
  
  plot_output <- here("outputs", "xenium", "preprocess", "03_clustered", "plots")
  if (!dir.exists(plot_output)) dir.create(plot_output)
  rds_output <- here("outputs", "xenium", "preprocess", "03_clustered", "rds")
  if (!dir.exists(rds_output)) dir.create(rds_output)
  
  DefaultAssay(xenium_obj) <- "Xenium"
  assay_name <- DefaultAssay(xenium_obj)
  cat(paste0("Using assay: ", assay_name, ":\n"))
  
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
  cat(paste0("Using graph: ", graph_name, ":\n"))
  
  umap_plot <- DimPlot(xenium_obj, reduction = "umap") +
    ggtitle(paste0(sample_name, " - UMAP Projection (dims = 1:50)"))
  
  ## =========================================================
  ## 4. Final Clustering
  ## =========================================================
  xenium_obj <- FindClusters(
    xenium_obj,
    resolution = 1.5,
    graph.name = graph_name,
    algorithm = 1,
    n.iter = 100,
    random.seed = 42,
    verbose = FALSE
  )
  
  umap_cluster_plot <- DimPlot(xenium_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste(sample_name, "- UMAP Clusters (Res = 1.5)"))
  
  ## =========================================================
  ## 5. Export UMAP + Modularity TIFF (Using CairoTIFF)
  ## =========================================================
  CairoTIFF(file.path(plot_output, paste0(sample_name, "_UMAP.tif")), width = 10, height = 10, units = "in", res = 600)
  print(umap_plot)
  dev.off()
  
  CairoTIFF(file.path(plot_output, paste0(sample_name, "_RawCluster_UMAP.tif")), width = 10, height = 10, units = "in", res = 600)
  print(umap_cluster_plot)
  dev.off()
  
  ## =========================================================
  ## 6. Spatial Cluster Plots
  ## =========================================================
  xenium_obj$seurat_clusters <- factor(Idents(xenium_obj))
  clusters <- levels(xenium_obj$seurat_clusters)
  n_clusters <- length(clusters)
  
  # Report number of clusters
  cat(paste0("Number of clusters generated: ", n_clusters, ":\n"))
  
  cluster_colors <- distinctColorPalette(n_clusters)
  names(cluster_colors) <- clusters
  
  global_cluster_plot <- ImageDimPlot(
    object = xenium_obj,
    fov = "fov",
    group.by = "seurat_clusters",
    cols = cluster_colors,
    size = 0.75
) +
    ggtitle(paste(sample_name, "- Raw Clusters"))
  
  # This line removes the stroke from the actual points in the plot
  global_cluster_plot$layers[[1]]$aes_params$stroke <- 0
  
  # Export Global Plot
  tif_file_global <- file.path(plot_output, paste0(sample_name, "_GlobalRawClustersSpatialPlot.tif"))
  CairoTIFF(tif_file_global, width = 10, height = 10, units = "in", res = 600)
  print(global_cluster_plot)
  dev.off()
  
  coords <- GetTissueCoordinates(xenium_obj)
  metadata <- xenium_obj@meta.data
  plot_data <- cbind(coords, cluster = as.character(metadata$seurat_clusters))
  
  facet_cluster_plot <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.2) +
    facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "black", color = NA),
          plot.background = element_rect(fill = "black", color = NA),
          strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
          plot.title = element_text(color = "white", hjust = 0.5, size = 14)) +
    ggtitle(paste(sample_name, "- Raw Clusters"))
  
  
  # Export Facet Plot
  tif_file_facet <- file.path(plot_output, paste0(sample_name, "_FacetRawClustersSpatialPlot.tif"))
  CairoTIFF(tif_file_facet, width = 12, height = 10, units = "in", res = 600)
  print(facet_cluster_plot)
  dev.off()
  
  message("Saved spatial plots TIFs")
  
  ## =========================================================
  ## 7. Save Final Object
  ## =========================================================
  rds_file <- file.path(rds_output, paste0(sample_name, "_CB_QC_cluster.rds"))
  saveRDS(xenium_obj, rds_file)
  message("Saved processed Seurat object: ", rds_file)
  
  return(xenium_obj)
}
