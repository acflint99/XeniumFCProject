# ----------------------------
# Setup & Libraries
# ----------------------------
options(device = function(...) {
  pdf(NULL)
})

library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(future)

# Increase the global object size limit to 100GB
options(future.globals.maxSize = 100 * 1024^3)

# Load custom palette and ordering
source(here("scripts", "color_palette.R"))

# ----------------------------
# Manual Variables (Set these for troubleshooting)
# ----------------------------
sample_name       <- "FB78_X_G" # e.g., "Sample1"
reference_name    <- "Aldinger"    # e.g., "FetalCerebellum"
pred_score_thresh <- 0.6
output_root       <- "outputs"

# ----------------------------
# 1. Directory Setup
# ----------------------------
plots_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_Banksy_Plots"))
if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

tables_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_Banksy_Tables"))
if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

# ----------------------------
# 2. Load Datasets
# ----------------------------
reference_file <- here(output_root, "SingleCellRDS", paste0(reference_name, "_newClusters_newUMAPv2_5k.rds"))
xenium_file    <- here(output_root, "XeniumBanksyRDS", paste0(sample_name, "_CB_QC_Bcluster.rds"))

reference <- readRDS(reference_file)
xenium    <- readRDS(xenium_file)

# ----------------------------
# 3. Preparation & Anchor Finding
# ----------------------------
DefaultAssay(reference) <- "RNA"
DefaultAssay(xenium)    <- "Xenium"

reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
shared_genes <- intersect(rownames(reference), rownames(xenium))
reference_balanced <- subset(reference, downsample = 1000)
transfer_features  <- intersect(VariableFeatures(reference_balanced), shared_genes)

reference_balanced <- ScaleData(reference_balanced, features = transfer_features, verbose = FALSE)
xenium             <- ScaleData(xenium, features = transfer_features, verbose = FALSE)

reference_balanced <- RunPCA(reference_balanced, features = transfer_features, verbose = FALSE)
xenium             <- RunPCA(xenium, features = transfer_features, verbose = FALSE)

cat("Shared features:", length(transfer_features), "\n")

anchors <- FindTransferAnchors(
  reference = reference_balanced,
  query = xenium,
  normalization.method = "LogNormalize", 
  reduction = "rpca", 
  features = transfer_features,
  dims = 1:30,
  k.anchor = 20,
  k.score = 30,
  approx.pca = TRUE
)

# ----------------------------
# 4. Transfer & Metadata
# ----------------------------
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference_balanced$clusters_refined,
  dims = 1:30,
  store.weights = TRUE
)

xenium <- AddMetaData(xenium, predictions)
cat("Median prediction score:", median(xenium$prediction.score.max), "\n")

# ----------------------------
# 5. Filtering and Voting
# ----------------------------
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
# 6. Saving Tables
# ----------------------------
comparison <- xenium@meta.data %>%
  select(seurat_clusters, cluster_majority, cluster_weighted) %>%
  distinct() %>%
  arrange(as.numeric(as.character(seurat_clusters)))

write.csv(comparison, file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted_comparison.csv")), row.names = FALSE)

# ----------------------------
# 7. Visualization Section
# ----------------------------
validate_palette(unique(xenium$cluster_weighted))
xenium$cluster_weighted <- factor(xenium$cluster_weighted, levels = celltype_order)

# --- Plot 1: Histogram ---
p_hist <- ggplot(xenium@meta.data, aes(x = prediction.score.max)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = paste(sample_name, "LogNorm Prediction Scores"),
       subtitle = paste("Threshold:", pred_score_thresh),
       x = "Max Prediction Score", y = "Cell Count")

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_prediction_scores_hist.png")),
    width = 6, height = 4, units = "in", res = 300)
print(p_hist)
dev.off()

# --- Plot 2: ImageDimPlot majority ---
p1 <- ImageDimPlot(xenium, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) +
  ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Majority)"))

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_clustermajority_GlobalSpatialPlot.png")),
    width = 6, height = 6, units = "in", res = 300)
print(p1)
dev.off()

# --- Plot 3: ImageDimPlot weighted ---
p2 <- ImageDimPlot(xenium, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) +
  ggtitle(paste(sample_name, "-", reference_name, "- Clusters (Weighted)"))

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_clusterweighted_GlobalSpatialPlot.png")), 
    width = 6, height = 6, units = "in", res = 300)
print(p2)
dev.off()

# --- Plot 4: Spatial Facet Plot ---
# Check available image names if this fails: print(names(xenium@images))
coords <- GetTissueCoordinates(xenium) 
plot_data <- cbind(coords, cluster = factor(as.character(xenium$cluster_weighted), levels = celltype_order))

p3 <- ggplot(plot_data, aes(x = x, y = y, color = cluster)) + # Swapped to standard x/y
  geom_point(size = 0.1) +
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

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_clusterweighted_FacetSpatialPlot.png")), 
    width = 12, height = 8, units = "in", res = 300)
print(p3)
dev.off()

# --- Plot 5: UMAP ---
# Ensure UMAP was run on the xenium object previously, or run it here:
# xenium <- RunUMAP(xenium, dims = 1:30, verbose = FALSE)
p4 <- DimPlot(xenium, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_clusterweighted_UMAP.png")), 
    width = 8, height = 6, units = "in", res = 300)
print(p4)
dev.off()

# --- Plot 6: DotPlot ---
markers <- list(
  "RL" = c("MKI67", "LTBP1", "OTX2"),
  "UBC" = c("EOMES"),
  "Granule" = c("ATOH1", "PAX6", "NEUROD1", "RELN"),
  "Purkinje" = c("FOXP2", "CALB1", "DAB1"),
  "GABA" = c("PAX2", "GAD1", "GAD2"),
  "Glia" = c("SOX9", "ADCY2"),
  "OPC" = c("PDGFRA", "OLIG1"),
  "Meninges" = c("FOXC1", "SLC7A11"),
  "Endothelial" = c("CLDN5", "PECAM1"),
  "Immune" = c("P2RY12")
)

# Filter markers to only those present in the Xenium assay
existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium)))

Idents(xenium) <- factor(xenium$cluster_weighted, levels = rev(celltype_order))

p5 <- DotPlot(xenium, features = existing_markers, assay = "Xenium") + 
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(sample_name, "Marker Expression (LogNorm)"))

png(filename = here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_DotPlot_markers.png")), 
    width = 10, height = 6, units = "in", res = 300)
print(p5)
dev.off()

# ----------------------------
# 8. Save Annotated Object
# ----------------------------
output_dir <- here(output_root, paste0("Xenium", reference_name, "ABT_BanksyLogNorm_RDS"))
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, paste0(sample_name, "_", reference_name, "_annotated.rds"))

saveRDS(xenium, output_file)
cat("Annotated Xenium object saved to:", output_file, "\n")