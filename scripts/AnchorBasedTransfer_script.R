# Clear the environment
rm(list = ls())

# ============================
# Anchor-based label transfer for Xenium (50 PCs) + Cluster-level voting
# ============================

library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(future)
library(randomcoloR)

# allow up to 32 GB per future (worker)
options(future.globals.maxSize = 32 * 1024^3)
plan("multisession", workers = 8)

# ----------------------------
# 1. Load processed RDS objects
# ----------------------------
reference_file <- here("outputs", "SingleCellRDS", "Aldinger_newClusters_newUMAPv2_5k.rds")
reference <- readRDS(reference_file)

xenium_file <- here("outputs", "XeniumRDS", "GZFB5_X_G_CB_QC_cluster.rds")
xenium <- readRDS(xenium_file)

# ----------------------------
# 2. Restrict to shared genes
# ----------------------------
shared_genes <- intersect(rownames(reference), rownames(xenium))
reference <- subset(reference, features = shared_genes)
xenium <- subset(xenium, features = shared_genes)
cat("Number of shared genes:", length(shared_genes), "\n")

# ----------------------------
# 3. Find transfer anchors using 50 PCs
# ----------------------------
anchors <- FindTransferAnchors(
  reference = reference,
  query = xenium,
  normalization.method = "LogNormalize",
  dims = 1:50
)

# ----------------------------
# 4. Transfer cell type labels using 50 PCs
# ----------------------------
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$clusters_refined,  
  dims = 1:50
)
xenium <- AddMetaData(xenium, predictions)

# ----------------------------
# 5. Optional: filter low-confidence calls
# ----------------------------
xenium$high_conf <- xenium$prediction.score.max > 0.6

# ----------------------------
# 6. Quick sanity check
# ----------------------------
comparison_table <- table(xenium$seurat_clusters, xenium$predicted.id)
print(comparison_table)

# ----------------------------
# 7. Cluster-level majority voting
# ----------------------------
majority_labels <- xenium@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1])

xenium$cluster_majority <- majority_labels$cluster_majority[match(xenium$seurat_clusters, majority_labels$seurat_clusters)]

# ----------------------------
# 8. Cluster-level weighted voting
# ----------------------------
weighted_labels <- xenium@meta.data %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(score_sum, n = 1) %>%
  select(seurat_clusters, cluster_weighted = predicted.id)

xenium$cluster_weighted <- weighted_labels$cluster_weighted[match(xenium$seurat_clusters, weighted_labels$seurat_clusters)]

# ----------------------------
# 9. Compare majority vs weighted labels
# ----------------------------
comparison <- xenium@meta.data %>%
  select(seurat_clusters, cluster_majority, cluster_weighted) %>%
  distinct()

print(comparison)
num_diff <- sum(comparison$cluster_majority != comparison$cluster_weighted)
cat("Number of clusters where majority vs weighted labels differ:", num_diff, "\n")

# ----------------------------
# 10. Save tables
# ----------------------------
tables_dir <- here("outputs", "XeniumABTTables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(comparison, file = here(tables_dir, "GZFB5_Aldinger_cluster_majority_vs_weighted_comparison.csv"), row.names = FALSE)
write.csv(comparison_table, file = here(tables_dir, "GZFB5_Aldinger_cluster_prediction_cellcounts.csv"), row.names = FALSE)


# ----------------------------
# 11. Spatial visualization (PNG version)
# ----------------------------
plots_dir <- here("outputs", "XeniumABTPlots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ImageDimPlot majority
png(filename = here(plots_dir, "GZFB5_Aldinger_cluster_majority_GlobalSpatialPlot.png"),
    width = 6*300, height = 6*300, res = 300)
ImageDimPlot(xenium, group.by = "cluster_majority", size = 0.75) + ggtitle("GZFB5 - Aldinger - Clusters (Majority)")
dev.off()

# ImageDimPlot weighted
png(filename = here(plots_dir, "GZFB5_Aldinger_cluster_weighted_GlobalSpatialPlot.png"),
    width = 6*300, height = 6*300, res = 300)
ImageDimPlot(xenium, group.by = "cluster_weighted", size = 0.75) + ggtitle("GZFB5 - Aldinger - Clusters (Weighted)")
dev.off()

# Spatial ggplot
clusters <- unique(xenium$cluster_weighted)
n_clusters <- length(clusters)
cluster_colors <- distinctColorPalette(n_clusters)
names(cluster_colors) <- clusters

coords <- GetTissueCoordinates(xenium)
metadata <- xenium@meta.data
plot_data <- cbind(coords, cluster = as.character(metadata$cluster_weighted))

png(filename = here(plots_dir, "GZFB5_Aldinger_cluster_weighted_FacetSpatialPlot.png"),
    width = 12*300, height = 8*300, res = 300)
ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
  geom_point(size = 0.2) +
  facet_wrap(~cluster) +
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
  ggtitle("GZFB5 - Aldinger - Clusters (Weighted)")
dev.off()

# UMAP with majority labels
png(filename = here(plots_dir, "GZFB5_Aldinger_cluster_weighted_UMAP.png"),
    width = 8*300, height = 6*300, res = 300)
DimPlot(xenium, reduction = "umap", label = TRUE, group.by = "cluster_weighted")
dev.off()

# ----------------------------
# 12. DotPlot for marker expression (PNG version)
# ----------------------------
markers <- c(
  "FOXP2", "CALB1", "DAB1", "PAX2", "GAD1", "GAD2", "MKI67", "LTBP1", "OTX2",
  "EOMES", "ATOH1", "PAX6", "NEUROD1", "RELN", "PRDM13", "DLL1", "ASCL1",
  "SOX9", "ADCY2", "PDGFRA", "OLIG1", "FOXC1", "SLC7A11", "CLDN5", "PECAM1", "P2RY12"
)

celltype_order <- c(
  "Purkinje", "GABA", "RL", "UBC", "Granule", "VZ", 
  "Glia", "OPC", "Meninges", "Endothelial", "Immune"
)

Idents(xenium) <- factor(
  xenium$cluster_weighted,
  levels = rev(celltype_order)
)

png(filename = here(plots_dir, "GZFB5_Aldinger_cluster_weighted_DotPlot_markers.png"),
    width = 10*300, height = 6*300, res = 300)
DotPlot(xenium, features = markers) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("GZFB5 - Aldinger - Marker Expression by Cluster (Weighted)")
dev.off()

# ----------------------------
# 13. Save annotated Xenium object
# ----------------------------
output_file <- here("outputs", "XeniumABTRDS", "GZFB5_Aldinger_annotated.rds")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(xenium, output_file)
cat("Annotated Xenium object saved to", output_file, "\n")