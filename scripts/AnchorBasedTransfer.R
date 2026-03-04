# ============================
# Anchor-based label transfer for Xenium (50 PCs) + Cluster-level voting
# ============================

library(Seurat)
library(ggplot2)
library(here)
library(dplyr)

library(future)

# allow up to 32 GB per future (worker)
options(future.globals.maxSize = 32 * 1024^3)

# use all 8 CPUs
plan("multisession", workers = 8)

# ----------------------------
# 1. Load processed RDS objects
# ----------------------------

reference_file <- here("data", "Aldinger_filtered_5kgenes_newUMAP.rds")
reference <- readRDS(reference_file)

xenium_file <- here("outputs", "GZFB5_X_G_CB_QC_cluster.rds")
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
  refdata = reference$clusters_refined,  #changed to test
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
print(table(xenium$seurat_clusters, xenium$predicted.id))

# ----------------------------
# 7. Cluster-level majority voting
# ----------------------------
majority_labels <- xenium@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1])

# Map back to Xenium metadata
xenium$cluster_majority <- majority_labels$cluster_majority[match(xenium$seurat_clusters, majority_labels$seurat_clusters)]

# ----------------------------
# 8. Cluster-level weighted voting (by prediction.score.max)
# ----------------------------
weighted_labels <- xenium@meta.data %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(score_sum, n = 1) %>%
  select(seurat_clusters, cluster_weighted = predicted.id)

# Map back to Xenium metadata
xenium$cluster_weighted <- weighted_labels$cluster_weighted[match(xenium$seurat_clusters, weighted_labels$seurat_clusters)]

# ----------------------------
# 9. Compare majority vs weighted labels
# ----------------------------
comparison <- xenium@meta.data %>%
  select(seurat_clusters, cluster_majority, cluster_weighted) %>%
  distinct()

print(comparison)

# Optional: how many clusters differ between methods
num_diff <- sum(comparison$cluster_majority != comparison$cluster_weighted)
cat("Number of clusters where majority vs weighted labels differ:", num_diff, "\n")

# ----------------------------
# 10. Spatial visualization
# ----------------------------
ImageDimPlot(xenium, group.by = "cluster_majority", size = 0.75) +
  ggtitle("Cluster Majority Voting Labels")

ImageDimPlot(xenium, group.by = "cluster_weighted", size = 0.75) +
  ggtitle("Cluster Weighted Voting Labels")

clusters <- levels(xenium$cluster_weighted)
n_clusters <- length(clusters)

cluster_colors <- distinctColorPalette(n_clusters)
names(cluster_colors) <- clusters

coords <- GetTissueCoordinates(xenium)
metadata <- xenium@meta.data
plot_data <- cbind(coords, cluster = as.character(metadata$cluster_weighted))

ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
  geom_point(size = 0.2) +
  facet_wrap(~cluster) +
#  scale_color_manual(values = cluster_colors) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA),
        strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
        plot.title = element_text(color = "white", hjust = 0.5, size = 14)) +
  ggtitle(paste("GZFB5 - Clusters"))

DimPlot(xenium, reduction = "umap", label = TRUE, group.by = "cluster_majority")

#check key cell type markers----
markers <- c(
  "EBF2",      # Purkinje
  "ROR2",
  "FOXP2",
  "CALB1",
  "DAB1",
  "PAX2",      # GABA
  "GAD1",
  "GAD2",
  "PRDM13",      # VZ
  "DLL1",
  "TFAP2A",    # VZ
  "ATOH1",     # GCP / RL
  "PAX6",
  "MKI67",
  "WNT2B",     # RL
  "INHBB",     
  "LTBP1",
  "OTX2",
  "EOMES",     # UBC
  "NEUROD1",   # GN
  "RELN",
  "SLC1A2",     # Glia
  "SOX9",
  "ADCY2",
  "SOX2",
  "GATA1",     # RBC
  "PHOX2B",    # Brainstem
  "HOXB4",
  "LMX1B",
  "NEUROD2",
  "HTR2C",     # Choroid
  "FOXJ1",     # Ependymal
  "P2RY12",    # Microglia
  "PDGFRA",    # OPC
  "OLIG1",
  "PDGFRB",    # Pericytes
  "RGS5",
  "CSPG4",
  "CLDN5",     # Endothelial
  "PECAM1",
  "FOXC1",      # Meninges
  "SLC7A11"
)

# markers <- c(
#   "LAMA2", "COL3A1", "COL1A2", "SLC7A11"
# )

DotPlot(
  xenium,
  features = markers,
  group.by = "cluster_majority"  # change if needed
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

Idents(xenium) <- xenium$cluster_majority

markers_Brainstem <- FindMarkers(xenium, ident.1 = "Brainstem", only.pos = TRUE, max.cells.per.ident = 2000, test.use = "LR")

# ----------------------------
# 11. Save annotated Xenium object
# ----------------------------
#saveRDS(xenium, here("outputs", "xenium_annotated_cluster_labels.rds"))
cat("Annotated Xenium object with cluster labels saved to outputs/xenium_annotated_cluster_labels.rds\n")