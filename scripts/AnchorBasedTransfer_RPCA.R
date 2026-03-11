library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(future)

source(here("scripts", "color_palette.R"))


# ----------------------------
# Main function
# ----------------------------
run_label_transfer <- function(sample_name,
                               reference_name,  # Keep the name for file paths
                               pred_score_thresh = 0.6,
                               output_root = "outputs") {
  
  # ----------------------------
  # 1. Load Xenium datasets
  # ----------------------------
  reference_file <- here(output_root, "SingleCellRDS", paste0(reference_name, "_newClusters_newUMAPv2_5k.rds"))
  xenium_file <- here(output_root, "XeniumRDS", paste0(sample_name, "_CB_QC_cluster.rds"))
  
  reference <- readRDS(reference_file)
  xenium <- readRDS(xenium_file)
  
  # ----------------------------
  # 2. Restrict to shared genes
  # ----------------------------
  shared_genes <- intersect(rownames(reference), rownames(xenium))
  reference <- subset(reference, features = shared_genes)
  xenium <- subset(xenium, features = shared_genes)
  cat("Number of shared genes:", length(shared_genes), "\n")
  
  # 1. Downsample the dominant reference to balance the clusters
  reference_balanced <- subset(reference, downsample = 1000)
  
  # 2. Re-calculate Variable Features on the balanced reference
  reference_balanced <- FindVariableFeatures(reference_balanced, nfeatures = 2000)
  transfer_features <- intersect(VariableFeatures(reference_balanced), rownames(xenium))
  
  # 3. Scale the data (Crucial for PCA/RPCA)
  # Only scale the features we are using for the transfer to save time
  reference_balanced <- ScaleData(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium <- ScaleData(xenium, features = transfer_features, verbose = FALSE)
  
  # 4. Run PCA on both
  reference_balanced <- RunPCA(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium <- RunPCA(xenium, features = transfer_features, verbose = FALSE)
  
  # ----------------------------
  # 3. Find transfer anchors
  # ----------------------------
  anchors <- FindTransferAnchors(
    reference = reference_balanced,
    query = xenium,
    normalization.method = "LogNormalize", 
    reference.reduction = "pca",
    reduction = "rpca",         # Switching to RPCA to prevent dominance
    features = transfer_features,
    dims = 1:30,                # Drop to 30 to reduce noise
    k.anchor = 5                # Default is 5, keep it local
  )
  
  # ----------------------------
  # 4. Transfer cell type labels
  # ----------------------------
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_balanced$clusters_refined,
    dims = 1:30
  )
  xenium <- AddMetaData(xenium, predictions)
  cat("Finished anchor transfer", "\n")
  # ----------------------------
  # 5. Filter low-confidence calls
  # ----------------------------
  xenium$high_conf <- xenium$prediction.score.max > pred_score_thresh
  
  # ----------------------------
  # 6. Cluster-level majority voting
  # ----------------------------
  majority_labels <- xenium@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1])
  
  xenium$cluster_majority <- majority_labels$cluster_majority[
    match(xenium$seurat_clusters, majority_labels$seurat_clusters)
  ]
  
  # ----------------------------
  # 7. Cluster-level weighted voting
  # ----------------------------
  weighted_labels <- xenium@meta.data %>%
    group_by(seurat_clusters, predicted.id) %>%
    summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    slice_max(score_sum, n = 1) %>%
    select(seurat_clusters, cluster_weighted = predicted.id)
  
  xenium$cluster_weighted <- weighted_labels$cluster_weighted[
    match(xenium$seurat_clusters, weighted_labels$seurat_clusters)
  ]
  
  # ----------------------------
  # 8. Save comparison tables
  # ----------------------------
  tables_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_Tables"))
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison <- xenium@meta.data %>%
    select(seurat_clusters, cluster_majority, cluster_weighted) %>%
    distinct() %>%
    # Convert to numeric if it's a factor, then sort
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(
    comparison,
    file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted_comparison.csv")),
    row.names = FALSE
  )
  
  comparison_table <- table(xenium$seurat_clusters, xenium$predicted.id)
  write.csv(
    comparison_table,
    file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")),
    row.names = FALSE
  )
  cat("Saved cluster comparison tables to", tables_dir, "\n")
  
  # ----------------------------
  # . Lock factor order
  # ----------------------------
  validate_palette(unique(xenium$cluster_weighted))
  
  xenium$cluster_weighted <- factor(
    xenium$cluster_weighted,
    levels = celltype_order
  )
  
  # ----------------------------
  # 9. Spatial visualization
  # ----------------------------
  plots_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_Plots"))
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ImageDimPlot majority
  png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_majority_GlobalSpatialPlot.png")),
      width = 6*300, height = 6*300, res = 300)
  print(
    ImageDimPlot(xenium, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) +
    ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Majority)"))
  )
  dev.off()
  
  # ImageDimPlot weighted
  png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_GlobalSpatialPlot.png")),
      width = 6*300, height = 6*300, res = 300)
  print(
  ImageDimPlot(xenium, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) +
    ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Weighted)"))
  )
  dev.off()
  
  # ----------------------------
  # Spatial ggplot
  # ----------------------------
  coords <- GetTissueCoordinates(xenium, type = "Xenium")
  metadata <- xenium@meta.data
  plot_data <- cbind(coords, cluster = as.character(metadata$cluster_weighted))
  
  # Make cluster a factor with your desired order
  plot_data$cluster <- factor(plot_data$cluster, levels = celltype_order)
  
  png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_FacetSpatialPlot.png")),
      width = 12*300, height = 8*300, res = 300)
  print(
  ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.2) +
    facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
      plot.title = element_text(color = "white", hjust = 0.5, size = 14)
    ) +
    ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Weighted)"))
  )
  dev.off()
  
  # ----------------------------
  # UMAP with majority labels
  # ----------------------------
  png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_UMAP.png")),
      width = 8*300, height = 6*300, res = 300)
  print(DimPlot(xenium, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors))
  dev.off()
  
  # ----------------------------
  # DotPlot for marker expression
  # ----------------------------
  markers <- c(
    "MKI67", "LTBP1", "OTX2","EOMES", "ATOH1", "PAX6", "NEUROD1", "RELN",
    "PRDM13", "DLL1", "ASCL1",
    "FOXP2", "CALB1", "DAB1", "PAX2", "GAD1", "GAD2",
    "SOX9", "ADCY2", "PDGFRA", "OLIG1", "FOXC1", "SLC7A11", "CLDN5", "PECAM1", "P2RY12"
  )
  
  Idents(xenium) <- factor(xenium$cluster_weighted, levels = rev(celltype_order))
  
  png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_DotPlot_markers.png")),
      width = 10*300, height = 6*300, res = 300)
  print(
  DotPlot(xenium, features = markers) +
    RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(paste(sample_name, "-", reference_name, "- Marker Expression by Cluster (Weighted)"))
  )
  dev.off()
  cat("Saved plots to", plots_dir, "\n")
  
  # ----------------------------
  # Save annotated object
  # ----------------------------
    output_file <- here(output_root, paste0("Xenium", reference_name, "ABT_RDS"), paste0(sample_name, "_", reference_name, "_annotated.rds"))
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(xenium, output_file)
    cat("Annotated Xenium object saved to", output_file, "\n")
    rm(xenium, reference, anchors) # Explicitly clear large objects
    gc() # Force garbage collection
  
  # Return True
  return(TRUE)
}