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

# Increase the global object size limit
options(future.globals.maxSize = 400 * 1024^3)

# =========================================================
# Sample Mapping Logic (Slurm Array Support)
# =========================================================
sample_list <- c(
  # "GZFB4_X_G",
  # "FB124_X_G",
  # "FB198_X_G",
  # "FB328_1_X_G",
  # "FB330_1_X_G",
  # "FB78_X_G",
  # "GZFB5_X_G",
  "GZFB_12_X_G_1",
  "GZFB_12_X_G_2",
  "GZFB_12_X_G_3",
  "GZFB_12_X_G_4",
  "GZFB_12_X_G_5",
  "GZFB_1_X_G",
  "GZFB_9_X_G_1",
  "GZFB_9_X_G_2",
  "GZFB_9_X_G_3"
  
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
annotate_xenium_from_ref <- function(xenium_obj, sample_name, reference_name = "Aldinger") {
  
  ## ----------------------------
  ## 0. Setup & Paths
  ## ----------------------------
  set.seed(42)
  plan("sequential") # Switch to multisession if local
  #plan("multisession", workers = 8)
  
  # Load custom palette and ordering from your specific script
  # Note: Ensure this file exists relative to the execution path
  source(here("scripts", "color_palette.R")) 
  
  pred_score_thresh <- 0.6
  output_root       <- "outputs"
  
  plots_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Plots"))
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Tables"))
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  ## ----------------------------
  ## 1. Load Reference
  ## ----------------------------
  ref_path <- here(output_root, "SingleCellRDS", paste0(reference_name, "_newClusters_newUMAPv2_5k.rds"))
  if (!file.exists(ref_path)) stop("Reference file not found: ", ref_path)
  reference <- readRDS(ref_path)
  
  ## ----------------------------
  ## 2. Preparation & Anchor Finding
  ## ----------------------------
  DefaultAssay(reference) <- "RNA"
  DefaultAssay(xenium_obj) <- "Xenium"
  
  reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  shared_genes <- intersect(rownames(reference), rownames(xenium_obj))
  
  # Downsample reference for speed and balance
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
  
  # Majority Voting
  majority_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1], .groups = "drop")
  
  xenium_obj$cluster_majority <- majority_labels$cluster_majority[match(xenium_obj$seurat_clusters, majority_labels$seurat_clusters)]
  
  # Weighted Voting (by prediction scores)
  weighted_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters, predicted.id) %>%
    summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    slice_max(score_sum, n = 1) %>%
    select(seurat_clusters, cluster_weighted = predicted.id)
  
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
  # Create the cross-tabulation (Clusters x Predicted ID)
  count_matrix <- xenium_obj@meta.data %>%
    group_by(seurat_clusters, cluster_weighted) %>%
    tally() %>%
    pivot_wider(names_from = cluster_weighted, values_from = n, values_fill = 0) %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  # Export the CSV
  write.csv(count_matrix, 
            file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")), 
            row.names = FALSE)
  
  message("Contingency matrix exported for: ", sample_name)
  
  ## ----------------------------
  ## 5. Visualizations
  ## ----------------------------
  # Validate palette and factor levels (from color_palette.R)
  validate_palette(unique(xenium_obj$cluster_weighted))
  xenium_obj$cluster_weighted <- factor(xenium_obj$cluster_weighted, levels = celltype_order)
  
  p_umap<- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)
  ggsave(here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_UMAP.tif")), p_umap, device = "tiff", width = 8, height = 6, dpi = 600, compression = "lzw")
  
  # Histogram of scores
  p_hist <- ggplot(xenium_obj@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
    theme_minimal() + labs(title = paste(sample_name, "Scores"), x = "Max Prediction Score")
  
  ggsave(here(plots_dir, paste0(sample_name, "_prediction_scores_hist.tif")), p_hist, device = "tiff", width = 6, height = 4, dpi = 600, compression = "lzw")
  
  # Spatial Plots (Majority and Weighted)
  p_maj <- ImageDimPlot(xenium_obj, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Majority"))
  
  p_wei <- ImageDimPlot(xenium_obj, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Weighted"))
  
  ggsave(here(plots_dir, paste0(sample_name, "_Spatial_Majority.tif")), p_maj, device = "tiff", width = 8, height = 8, dpi = 600, compression = "lzw")
  ggsave(here(plots_dir, paste0(sample_name, "_Spatial_Weighted.tif")), p_wei, device = "tiff", width = 8, height = 8, dpi = 600, compression = "lzw")
  
  # Spatial Facet Plot
  coords <- GetTissueCoordinates(xenium_obj) 
  plot_data <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_weighted), levels = celltype_order))
  
  p_facet <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = 0.1) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  ggsave(here(plots_dir, paste0(sample_name, "_FacetSpatial.tif")), p_facet, device = "tiff", width = 12, height = 8, dpi = 600, compression = "lzw")
  
  # DotPlot
  existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium_obj)))
  Idents(xenium_obj) <- factor(xenium_obj$cluster_weighted, levels = rev(celltype_order))
  
  p_dot <- DotPlot(xenium_obj, features = existing_markers, assay = "Xenium") + 
    RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "red") +
    ggtitle(paste(sample_name, "Markers"))
  
  ggsave(here(plots_dir, paste0(sample_name, "_DotPlot.tif")), p_dot, device = "tiff", width = 10, height = 6, dpi = 600, compression = "lzw")
  
  ## ----------------------------
  ## 6. Save & Return
  ## ----------------------------
  output_dir <- here(output_root, paste0("Xenium_Res1.5", reference_name, "ABT_RDS"))
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