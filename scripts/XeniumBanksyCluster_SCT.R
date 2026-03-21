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

process_xenium_clusters <- function(xenium_obj, sample_name, lambda = 0.2, k_geom = 15) {
  
  ## =========================================================
  ## 0. Setup
  ## =========================================================
  set.seed(42)
  plan("sequential")  # avoid future export issues
  
  DefaultAssay(xenium_obj) <- "Xenium"
  assay_name <- DefaultAssay(xenium_obj)
  cat(paste0("Using assay: ", assay_name, ":\n"))
  
  ## =========================================================
  ## 1. Normalization & BANKSY Setup
  ## =========================================================
  message("Running SCTransform...")
  
  # SCTransform replaces NormalizeData, FindVariableFeatures, and ScaleData
  xenium_obj <- SCTransform(xenium_obj, assay = assay_name, clip.range = c(-10, 10), verbose = FALSE)
  
  # Update our active assay variable to SCT
  assay_name <- "SCT"
  DefaultAssay(xenium_obj) <- assay_name
  
  # Run BANKSY to compute the spatial neighborhood augmented matrix
  message("Running BANKSY...")
  xenium_obj <- RunBanksy(
    xenium_obj,
    lambda = lambda,     
    k_geom = k_geom,     
    assay = assay_name,
    slot = "data",       # Use the normalized data slot from SCT
    features = "variable",
    verbose = FALSE
  )
  
  # BANKSY creates a new assay automatically. We set it as default.
  DefaultAssay(xenium_obj) <- "BANKSY"
  message("BANKSY assay created. Running PCA on augmented features...")
  
  # Run PCA on the BANKSY assay 
  xenium_obj <- RunPCA(xenium_obj, assay = "BANKSY", features = rownames(xenium_obj), npcs = 50, verbose = FALSE, seed.use = 42)
  
  # Free up memory (SCTransform stores residuals in scale.data which can be large)
  LayerData(xenium_obj[[assay_name]], layer = "scale.data") <- NULL
  message("Removed scale.data layer from assay: ", assay_name)
  
  ## =========================================================
  ## 2. UMAP + Neighbor Graph
  ## =========================================================
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:50, seed.use = 42, verbose = FALSE)
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
  cat(paste0("Optimal clustering resolution (highest modularity): ", best_res, ":\n"))
  
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
    ggtitle(paste(sample_name, "- BANKSY Clusters (Res =", best_res, ")"))
  
  ## =========================================================
  ## 5. Export UMAP + Modularity PDF
  ## =========================================================
  pdf_file <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_BanksyClusterOptim_UMAP_Plots.pdf"))
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
  
  png_file_global <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_GlobalRawBanksyClustersSpatialPlot.png"))
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
  
  
  png_file_facet <- file.path(here("outputs", "XeniumBanksyPlots"), paste0(sample_name, "_FacetRawBanksyClustersSpatialPlot.png"))
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