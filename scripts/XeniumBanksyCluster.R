library(Seurat)
library(SeuratWrappers) # Required for RunBanksy wrapper
library(Banksy)         # Required for the BANKSY algorithm
library(ggplot2)
library(patchwork)
library(randomcoloR)
library(future)
library(igraph)
library(dplyr)
library(here)

process_xenium_clusters <- function(xenium_obj, sample_name, lambda = 0.1, k_geom = 10) {
  
  ## =========================================================
  ## 0. Setup
  ## =========================================================
  set.seed(42)
  plan("sequential")
  
  DefaultAssay(xenium_obj) <- "Xenium"
  assay_name <- DefaultAssay(xenium_obj)
  cat(paste0("Using assay: ", assay_name, ":\n"))
  
  ## =========================================================
  ## 1. Normalization & BANKSY Setup
  ## =========================================================
  message("Running Log Normalization...")
  
  # REPLACE SCTransform with standard Log-Normalization
  xenium_obj <- NormalizeData(xenium_obj, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000, 
                              verbose = FALSE)
  
  # Identify Variable Features (For Xenium, using all genes is often best)
  xenium_obj <- FindVariableFeatures(xenium_obj, 
                                     selection.method = "vst", 
                                     nfeatures = 3000, 
                                     verbose = FALSE)
  
  # Scale the data (Required for PCA)
  # Scale only the 3000 Variable Features for Banksy input
  xenium_obj <- ScaleData(xenium_obj, 
                          assay = "Xenium", 
                          features = VariableFeatures(xenium_obj), 
                          verbose = FALSE)
  
  # Update assay variable (remains "Xenium" for standard normalization)
  # DefaultAssay(xenium_obj) is already "Xenium"
  
  # Run BANKSY
  message("Running BANKSY...")
  xenium_obj <- RunBanksy(
    xenium_obj,
    lambda = lambda,     
    k_geom = k_geom,     
    assay = assay_name,
    slot = "data",       # Use the log-normalized data slot
    features = "variable",
    verbose = FALSE
  )
  
  # BANKSY creates a new assay automatically. Set as default.
  DefaultAssay(xenium_obj) <- "BANKSY"
  message("Scaling BANKSY assay...")
  xenium_obj <- ScaleData(xenium_obj, assay = "BANKSY", verbose = FALSE)
  
  # Ensure VariableFeatures are set for the BANKSY assay specifically
  VariableFeatures(xenium_obj) <- rownames(xenium_obj)
  
  message("BANKSY assay created. Running PCA on augmented features...")
  
  # Increase npcs and use all Banksy features
  xenium_obj <- RunPCA(xenium_obj, 
                       assay = "BANKSY", 
                       features = VariableFeatures(xenium_obj), 
                       npcs = 50, 
                       verbose = FALSE, 
                       seed.use = 42)
  
  # Note: No need to NULL out scale.data here like we did for SCT residuals
  # unless you are extremely tight on memory.
  
  ## =========================================================
  ## 2. UMAP + Neighbor Graph
  ## =========================================================
  # Using 50 dims and higher neighbors to resolve the "hairball"
  message("Running UMAP and Finding Neighbors...")
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:50, n.neighbors = 50, min.dist = 0.1, seed.use = 42, verbose = FALSE)
  xenium_obj <- FindNeighbors(xenium_obj, dims = 1:50, nn.method = "annoy")
  
  # The graph should now dynamically name itself based on the active BANKSY assay
  graph_name <- "BANKSY_snn"
  if (!graph_name %in% names(xenium_obj@graphs)) {
    graph_name <- paste0(DefaultAssay(xenium_obj), "_snn")
    if (!graph_name %in% names(xenium_obj@graphs)) {
      stop("SNN graph not found. Available graphs: ", paste(names(xenium_obj@graphs), collapse = ", "))
    }
  }
  cat(paste0("Using graph: ", graph_name, ":\n"))
  
  umap_plot <- DimPlot(xenium_obj, reduction = "umap") +
    ggtitle(paste0(sample_name, " - BANKSY UMAP Projection"))
  
  ## =========================================================
  ## 3 & 4. Fixed High-Resolution Clustering
  ## =========================================================
  # Skip the modularity loop for consistency across samples
  fixed_res <- 1 
  
  message(paste0("Running final clustering at fixed resolution: ", fixed_res))
  xenium_obj <- FindClusters(
    xenium_obj,
    resolution = fixed_res,
    graph.name = graph_name,
    algorithm = 1,
    n.iter = 100,
    random.seed = 42,
    verbose = FALSE
  )
  
  umap_cluster_plot <- DimPlot(xenium_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste(sample_name, "- BANKSY Clusters (Res =", fixed_res, ")"))
  
  png_file_umap <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_BanksyCluster1_UMAP_Plot.png"))
  png(png_file_umap, width = 10, height = 10, units = "in", res = 300)
  print(umap_cluster_plot)
  dev.off()
  
  ## =========================================================
  ## 6. Spatial Cluster Plots
  ## =========================================================
  xenium_obj$seurat_clusters <- factor(Idents(xenium_obj))
  clusters <- levels(xenium_obj$seurat_clusters)
  n_clusters <- length(clusters)
  
  cat(paste0("Number of clusters generated: ", n_clusters, ":\n"))
  
  cluster_colors <- distinctColorPalette(n_clusters)
  names(cluster_colors) <- clusters
  
  global_cluster_plot <- ImageDimPlot(
    object = xenium_obj,
    fov = "fov",
    group.by = "seurat_clusters",
    cols = cluster_colors,
    size = 0.75,
  ) + 
    ggtitle(paste(sample_name, "- BANKSY Raw Clusters"))
  
  global_cluster_plot$layers[[1]]$aes_params$stroke <- 0
  
  png_file_global <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_GlobalRawBanksyClusters1SpatialPlot.png"))
  png(png_file_global, width = 10, height = 10, units = "in", res = 300)
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
    ggtitle(paste(sample_name, "- BANKSY Raw Clusters"))
  
  
  png_file_facet <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_FacetRawBanksyClusters1SpatialPlot.png"))
  png(png_file_facet, width = 12, height = 10, units = "in", res = 300)
  print(facet_cluster_plot)
  dev.off()
  
  message("Saved spatial plots PNGs")
  
  ## =========================================================
  ## 7. Save Final Object
  ## =========================================================
  rds_file <- file.path(here("outputs", "XeniumBanksyRDS"), paste0(sample_name, "_CB_QC_Bcluster.rds"))
  saveRDS(xenium_obj, rds_file)
  message("Saved processed Seurat object: ", rds_file)
  
  return(xenium_obj)
}