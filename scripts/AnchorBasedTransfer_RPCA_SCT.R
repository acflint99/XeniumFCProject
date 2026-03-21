options(device = function(...) {
  pdf(NULL)
})

library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(future)

source(here("scripts", "color_palette.R"))

# ----------------------------
# Main function
# ----------------------------
run_label_transfer_sct <- function(sample_name,
                                   reference_name,
                                   pred_score_thresh = 0.6,
                                   output_root = "outputs") {
  # Set a global option to use a null device for all ggplot-based exports
  options(ggplot2.as_ggplot = TRUE)
  
  # Force R to use a null device for ghost plots
  options(bitmapType='cairo') 
  pdf(NULL) 
  on.exit(while (!is.null(dev.list())) dev.off()) # Clean up when function ends
  
  plots_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_BanksySCT_Plots"))
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_BanksySCT_Tables"))
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  # 1. Load Datasets
  reference_file <- here(output_root, "SingleCellRDS", paste0(reference_name, "_newClusters_5k_SCT_newUMAPv2.rds"))
  xenium_file <- here(output_root, "XeniumBanksySCT_RDS", paste0(sample_name, "_CB_QC_Bcluster.rds"))
  
  reference <- readRDS(reference_file)
  xenium <- readRDS(xenium_file)
  
  # 2. Preparation for SCT-based transfer
  DefaultAssay(reference) <- "SCT"
  DefaultAssay(xenium) <- "SCT"
  
  shared_genes <- intersect(rownames(reference), rownames(xenium))
  reference_balanced <- subset(reference, downsample = 1000)
  transfer_features <- intersect(VariableFeatures(reference_balanced, assay = "SCT"), shared_genes)
  
  # Ensure the scale.data slot is populated for RPCA
  reference_balanced <- ScaleData(reference_balanced, features = transfer_features, assay = "SCT", verbose = FALSE)
  xenium <- ScaleData(xenium, features = transfer_features, assay = "SCT", verbose = FALSE)
  
  # Now RunPCA will have the data it needs
  reference_balanced <- RunPCA(reference_balanced, features = transfer_features, assay = "SCT", verbose = FALSE)
  xenium <- RunPCA(xenium, features = transfer_features, assay = "SCT", verbose = FALSE)
  
  cat("Number of shared transfer features:", length(transfer_features), "\n")
  if(length(transfer_features) < 100) warning("Very few shared features found for transfer!")
  
  # 3. Find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference_balanced,
    query = xenium,
    normalization.method = "SCT",
    reference.reduction = "pca",
    reduction = "rpca", 
    features = transfer_features,
    dims = 1:30,
    k.anchor = 30,
    k.score = 50,           # Added: Helps filter anchors in dense spatial regions
    max.features = 200,     # Added: Limits to top features for RPCA projection speed/stability
    approx.pca = TRUE
  )
  
  # 4. Transfer cell type labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_balanced$clusters_refined,
    dims = 1:30,
    store.weights = TRUE
  )
  
  xenium <- AddMetaData(xenium, predictions)
  cat("Finished anchor transfer. Median prediction score:", median(xenium$prediction.score.max), "\n")
  
  # NEW: Save Prediction Score Histogram
  p_hist <- ggplot(xenium@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste(sample_name, "Prediction Scores"),
         subtitle = paste("Threshold:", pred_score_thresh),
         x = "Max Prediction Score", y = "Cell Count")
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_prediction_scores_hist.png")),
         plot = p_hist, width = 6, height = 4, dpi = 300)
  
  # 5. Filter and Voting
  # Check if the threshold is too aggressive
  pass_count <- sum(xenium$prediction.score.max > pred_score_thresh)
  cat("Cells passing threshold:", pass_count, "out of", nrow(xenium@meta.data), "\n")
  
  xenium$high_conf <- xenium$prediction.score.max > pred_score_thresh
  
  majority_labels <- xenium@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1], .groups = "drop")
  
  xenium$cluster_majority <- majority_labels$cluster_majority[match(xenium$seurat_clusters, majority_labels$seurat_clusters)]
  
  weighted_labels <- xenium@meta.data %>%
    group_by(seurat_clusters, predicted.id) %>%
    summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    slice_max(score_sum, n = 1) %>%
    select(seurat_clusters, cluster_weighted = predicted.id)
  
  xenium$cluster_weighted <- weighted_labels$cluster_weighted[match(xenium$seurat_clusters, weighted_labels$seurat_clusters)]
  
  # ----------------------------
  # 8. Save Tables
  # ----------------------------
  comparison <- xenium@meta.data %>%
    select(seurat_clusters, cluster_majority, cluster_weighted) %>%
    distinct() %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(comparison, file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted_comparison.csv")), row.names = FALSE)
  
  # ----------------------------
  # 9. Visualization with ggsave
  # ----------------------------
  validate_palette(unique(xenium$cluster_weighted))
  xenium$cluster_weighted <- factor(xenium$cluster_weighted, levels = celltype_order)
  
  # ImageDimPlot majority
  p1 <- ImageDimPlot(xenium, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) +
    ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Majority)"))
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_majority_GlobalSpatialPlot.png")),
         plot = p1, width = 6, height = 6, dpi = 300)
  
  # ImageDimPlot weighted
  p2 <- ImageDimPlot(xenium, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) +
    ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Weighted)"))
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_GlobalSpatialPlot.png")),
         plot = p2, width = 6, height = 6, dpi = 300)
  
  # Spatial Facet Plot
  coords <- GetTissueCoordinates(xenium, type = "Xenium")
  plot_data <- cbind(coords, cluster = factor(as.character(xenium$cluster_weighted), levels = celltype_order))
  
  p3 <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.2) +
    facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
      plot.title = element_text(color = "white", hjust = 0.5, size = 14)
    ) +
    ggtitle(paste(sample_name, "Facet (Weighted)"))
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_FacetSpatialPlot.png")),
         plot = p3, width = 12, height = 8, dpi = 300)
  
  # UMAP
  p4 <- DimPlot(xenium, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_UMAP.png")),
         plot = p4, width = 8, height = 6, dpi = 300)
  
  # DotPlot
  markers <- list(
    "RL" = c("MKI67", "LTBP1", "OTX2"),
    "UBC" = c("EOMES"),
    "Granule" = c("ATOH1", "PAX6", "NEUROD1", "RELN"),
    #"VZ" = c("PRDM13", "DLL1", "ASCL1"),
    "Purkinje" = c("FOXP2", "CALB1", "DAB1"),
    "GABA" = c("PAX2", "GAD1", "GAD2"),
    "Glia" = c("SOX9", "ADCY2"),
    "OPC" = c("PDGFRA", "OLIG1"),
    "Meninges" = c("FOXC1", "SLC7A11"),
    "Endothelial" = c("CLDN5", "PECAM1"),
    "Immune" = c("P2RY12")
  )
  
  
  Idents(xenium) <- factor(xenium$cluster_weighted, levels = rev(celltype_order))
  
  p5 <- DotPlot(xenium, features = markers, assay = "SCT") +
    RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(sample_name, "Marker Expression"))
  
  ggsave(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_DotPlot_markers.png")),
         plot = p5, width = 10, height = 6, dpi = 300)
  
  cat("Saved plots to", plots_dir, "\n")
  
  # 10. Save annotated object
  output_file <- here(output_root, paste0("Xenium", reference_name, "ABT_BanksySCT_RDS"), paste0(sample_name, "_", reference_name, "_annotated.rds"))
  if(!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
  saveRDS(xenium, output_file)
  
  cat("Annotated Xenium object saved to", output_file, "\n")
  rm(xenium, reference, anchors)
  gc()
  
  return(TRUE)
}