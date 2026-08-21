# ----------------------------
# Setup & Libraries
# ----------------------------
library(future)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)
library(tidyr)
library(Cairo)

# Increase the global object size limit
options(future.globals.maxSize = 400 * 1024^3)

# =========================================================
# Sample Mapping Logic (Slurm Array Support)
# =========================================================
sample_list <- c(
  # "FB328_1_X_G",
  #  "FB330_1_X_G",
  # "FB78_X_G",
  # "GZFB5_X_G",
  #  "GZFB_12_X_G_1",
  # "GZFB_12_X_G_2",
  # "GZFB_12_X_G_3",
  #  "GZFB_12_X_G_4",
  # "GZFB_12_X_G_5",
  # "GZFB_1_X_G",
  # "GZFB_9_X_G_1",
  # "GZFB_9_X_G_2",
  # "GZFB_9_X_G_3",
  # "GZFB_20_X_G_1",
  # "GZFB_20_X_G_11",
  # "GZFB_20_X_G_12"
  # "GZFB_20_X_G_2",
  # "GZFB_20_X_G_3",
  # "GZFB_20_X_G_4",
  # "GZFB_22_X_G_6",
  # "GZFB_20_X_G_9",
  # "GZFB_20_X_G_7",
  # "GZFB_20_X_G_5",
  # "GZFB_20_X_G_10",
  # "GZFB_20_X_G_8",
  # "GZFB_20_X_G_6",
  # "GZFB_22_X_G_4",
  # "GZFB_22_X_G_2",
  # "GZFB_22_X_G_5",
  # "GZFB_22_X_G_3",
  # "GZFB_22_X_G_1"
   "GZFB4_X_G",
   "FB124_X_G",
   "FB198_X_G"
  
)

# Get the Task ID from Slurm (e.g., --array=1-16)
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

if (is.na(task_id) || task_id < 1 || task_id > length(sample_list)) {
  stop("Error: Task ID is out of bounds or not provided.")
}

current_sample <- sample_list[task_id]
message("### Processing Annotation for Sample [", task_id, "]: ", current_sample, " ###")

# =========================================================
# Annotation Function
# =========================================================
annotate_xenium_from_ref <- function(xenium_obj, sample_name, reference_name = "Science") {
  
  ## ----------------------------
  ## 0. Setup & Paths
  ## ----------------------------
  set.seed(42)
  plan("sequential") 
  
  # Load custom palette and ordering from your specific script
  source(here("scripts", "color_palette.R")) 
  
  pred_score_thresh <- 0.4
  output_root       <- "outputs"
  
  plots_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Plots"))
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Tables"))
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  ## ----------------------------
  ## 1. Load Reference
  ## ----------------------------
  ref_path <- here(output_root, paste0(reference_name, "RDS"), paste0(reference_name, "_newClusters_newUMAPv2_5k.rds"))
  if (!file.exists(ref_path)) stop("Reference file not found: ", ref_path)
  reference <- readRDS(ref_path)
  
  ## ----------------------------
  ## 2. Preparation & Anchor Finding
  ## ----------------------------
  # Set reference assay safely
  if ("RNA" %in% Assays(reference)) {
    DefaultAssay(reference) <- "RNA"
  } else {
    # Fallback to whatever assay exists (e.g., "originalexp")
    DefaultAssay(reference) <- Assays(reference)[1] 
  }
  
  # Set xenium assay safely
  if ("Xenium" %in% Assays(xenium_obj)) {
    DefaultAssay(xenium_obj) <- "Xenium"
  } else {
    DefaultAssay(xenium_obj) <- Assays(xenium_obj)[1]
  }
  
  reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  shared_genes <- intersect(rownames(reference), rownames(xenium_obj))
  
  # Downsample reference for speed and balance
  # Set the active identity to your cell types BEFORE downsampling
  Idents(reference) <- "clusters_refined"
  reference_balanced <- subset(reference, downsample = 1000)
  transfer_features  <- intersect(VariableFeatures(reference_balanced), shared_genes)
  
  reference_balanced <- ScaleData(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium_obj         <- ScaleData(xenium_obj, features = transfer_features, verbose = FALSE)
  
  reference_balanced <- RunPCA(reference_balanced, features = transfer_features, verbose = FALSE)
  xenium_obj         <- RunPCA(xenium_obj, features = transfer_features, verbose = FALSE)
  
  cat("Shared features for transfer:", length(transfer_features), "\n")
  
  anchors <- FindTransferAnchors(
    reference = reference_balanced,
    query = xenium_obj,
    normalization.method = "LogNormalize", 
    reduction = "rpca", 
    features = transfer_features,
    dims = 1:30,
    k.anchor = 20,
    k.score = 30,
    approx.pca = TRUE
  )
  
  ## ----------------------------
  ## 3. Transfer & Voting Logic
  ## ----------------------------
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_balanced$clusters_refined,
    dims = 1:30,
    store.weights = TRUE
  )
  
  xenium_obj <- AddMetaData(xenium_obj, predictions)
  xenium_obj$high_conf <- xenium_obj$prediction.score.max > pred_score_thresh
  
  xenium_obj <- RunUMAP(xenium_obj, dims = 1:30, reduction = "pca")
  
  # -------------------
  # Thresholded Majority Voting
  # -------------------
  majority_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1],
      mean_max_score = mean(prediction.score.max), # Calculate cluster average confidence
      .groups = "drop"
    ) %>%
    mutate(
      cluster_majority = ifelse(mean_max_score >= pred_score_thresh, cluster_majority, "Unknown")
    )
  
  xenium_obj$cluster_majority <- majority_labels$cluster_majority[match(xenium_obj$seurat_clusters, majority_labels$seurat_clusters)]
  
  # -------------------
  # Thresholded Rigorous Weighted Voting 
  # -------------------
  score_cols <- grep("^prediction\\.score\\.", colnames(xenium_obj@meta.data), value = TRUE)
  score_cols <- score_cols[!score_cols %in% c("prediction.score.max", "prediction.score.id")]
  
  # Get cluster sizes to normalize the total probability sums
  cluster_counts <- xenium_obj@meta.data %>% count(seurat_clusters, name = "n_cells")
  
  weighted_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), sum), .groups = "drop") %>%
    pivot_longer(cols = -seurat_clusters, names_to = "cluster_weighted", values_to = "total_score") %>%
    group_by(seurat_clusters) %>%
    slice_max(total_score, n = 1, with_ties = FALSE) %>%
    left_join(cluster_counts, by = "seurat_clusters") %>%
    mutate(
      cluster_weighted = sub("^prediction\\.score\\.", "", cluster_weighted),
      avg_prob = total_score / n_cells, # Calculate the mean probability per cell for the winning class
      cluster_weighted = ifelse(avg_prob >= pred_score_thresh, cluster_weighted, "Unknown")
    )
  
  xenium_obj$cluster_weighted <- weighted_labels$cluster_weighted[match(xenium_obj$seurat_clusters, weighted_labels$seurat_clusters)]
  
  ## ----------------------------
  ## 4. Export Comparison Table
  ## ----------------------------
  comparison <- xenium_obj@meta.data %>%
    select(seurat_clusters, cluster_majority, cluster_weighted) %>%
    distinct() %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(comparison, file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted.csv")), row.names = FALSE)
  
  ## ----------------------------
  ## 4c. Export Cluster-by-CellType Contingency Table
  ## ----------------------------
  # Note: Retains `predicted.id` so you can see the raw background noise/variance that caused a cluster to become "Unknown"
  count_matrix <- xenium_obj@meta.data %>%
    count(seurat_clusters, predicted.id) %>%
    pivot_wider(names_from = predicted.id, values_from = n, values_fill = list(n = 0)) %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(count_matrix, 
            file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")), 
            row.names = FALSE)
  
  message("Contingency matrix exported for: ", sample_name)
  
  ## ----------------------------
  ## 5. Visualizations
  ## ----------------------------
  
  # Safely add "Unknown" to your palette if it isn't already in color_palette.R
  if (!"Unknown" %in% celltype_order) {
    celltype_order <- c(celltype_order, "Unknown")
    if (!"Unknown" %in% names(cluster_colors)) {
      cluster_colors["Unknown"] <- "grey80" 
    }
  }
  
  # Validate palette (exclude "Unknown" from the strict check so it doesn't crash)
  clusters_to_validate <- setdiff(unique(c(xenium_obj$cluster_weighted, xenium_obj$cluster_majority)), "Unknown")
  validate_palette(clusters_to_validate)
  
  xenium_obj$cluster_weighted <- factor(xenium_obj$cluster_weighted, levels = celltype_order)
  xenium_obj$cluster_majority <- factor(xenium_obj$cluster_majority, levels = celltype_order)
  
  # --- UMAPs ---
  p_umap_wei <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterWeighted_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_wei)
  dev.off()
  
  p_umap_maj <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_majority", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterMajority_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_maj)
  dev.off()
  
  # --- Histogram of scores ---
  p_hist <- ggplot(xenium_obj@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
    theme_minimal() + labs(title = paste(sample_name, "Scores"), x = "Max Prediction Score")
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.tif")), width = 6, height = 4, units = "in", res = 600)
  print(p_hist)
  dev.off()
  ggsave(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.pdf")), p_hist, device = "pdf", width = 6, height = 4)
  
  # --- Spatial Plots ---
  p_wei <- ImageDimPlot(xenium_obj, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Weighted"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterWeighted.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_wei)
  dev.off() # Fixed missing parentheses
  
  p_maj <- ImageDimPlot(xenium_obj, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Majority")) # Fixed title
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterMajority.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_maj)
  dev.off()
  
  # --- Spatial Facet Plots ---
  coords <- GetTissueCoordinates(xenium_obj) 
  
  # Weighted Facet
  plot_data_wei <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_weighted), levels = celltype_order))
  p_facet_wei <- ggplot(plot_data_wei, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = 0.1) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterWeighted.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_wei)
  dev.off()
  
  # Majority Facet
  plot_data_maj <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_majority), levels = celltype_order))
  p_facet_maj <- ggplot(plot_data_maj, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = 0.1) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterMajority.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_maj)
  dev.off()
  
  # --- DotPlots ---
  existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium_obj)))
  
  # Weighted DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_weighted, levels = rev(celltype_order))
  p_dot_wei <- DotPlot(xenium_obj, features = existing_markers, assay = "Xenium") + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "red") +
    ggtitle(paste(sample_name, "Markers (Weighted)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_wei)
  dev.off()
  ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.pdf")), p_dot_wei, device = "pdf", width = 10, height = 6)
  
  # Majority DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_majority, levels = rev(celltype_order))
  p_dot_maj <- DotPlot(xenium_obj, features = existing_markers, assay = "Xenium") + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "red") +
    ggtitle(paste(sample_name, "Markers (Majority)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_maj)
  dev.off()
  ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.pdf")), p_dot_maj, device = "pdf", width = 10, height = 6)
  
  ## ----------------------------
  ## 6. Save & Return
  ## ----------------------------
  output_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_RDS"))
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(xenium_obj, file.path(output_dir, paste0(sample_name, "_", reference_name, "_annotated.rds")))
  message("Successfully annotated and saved: ", sample_name)
  
  return(xenium_obj)
}

# =========================================================
# Execution
# =========================================================
input_file <- here("outputs", "Xenium_Res1.5_RDS", paste0(current_sample, "_CB_QC_cluster.rds"))

if (file.exists(input_file)) {
  seu <- readRDS(input_file)
  annotate_xenium_from_ref(seu, current_sample)
} else {
  stop("Input file not found: ", input_file)
}