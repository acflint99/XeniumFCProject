# Clear the environment
rm(list = ls())

library(here)
library(Seurat)

# Source the function
source(here("scripts", "XeniumCropCerebellum.R"))

# Load and crop a sample
xenium_cereb <- XeniumCropCerebellum("GZFB5_X_G")

# Source the function
source(here("scripts", "XeniumQC.R"))

# Filter out low quality cells
xenium_cereb_QC <- qc_xenium(xenium_cereb, sample_name = "GZFB5_X_G")











library(Seurat)
library(pracma)    # for numeric derivative (elbow)
library(cluster)   # for silhouette()
library(ggplot2)
library(patchwork)
library(RColorBrewer)

## ------------------------------
## 1. Normalization and preprocessing
## ------------------------------

xenium_cereb <- NormalizeData(xenium_cereb, scale.factor = median(xenium_cereb$nCount_Xenium))
xenium_cereb <- FindVariableFeatures(xenium_cereb, selection.method = "vst", nfeatures = 2000)
xenium_cereb <- ScaleData(xenium_cereb)
xenium_cereb <- RunPCA(xenium_cereb, verbose = FALSE)

## ------------------------------
## 2. Automatic PCA elbow detection
## ------------------------------

pca_sd <- xenium_cereb[["pca"]]@stdev
pca_var <- pca_sd^2
pca_var_ratio <- pca_var / sum(pca_var)
second_deriv <- diff(diff(pca_var_ratio))
elbow_pc <- which.min(second_deriv) + 1
if(elbow_pc > 50) elbow_pc <- 50
message("Optimal number of PCs selected: ", elbow_pc)

# Elbow plot
elbow_plot <- ElbowPlot(xenium_cereb, ndims = 50) +
  geom_vline(xintercept = elbow_pc, linetype = "dashed", color = "red") +
  ggtitle("PCA Elbow Plot")

## ------------------------------
## 3. UMAP embedding
## ------------------------------

xenium_cereb <- RunUMAP(xenium_cereb, dims = 1:elbow_pc, verbose = FALSE)
umap_plot <- DimPlot(xenium_cereb, reduction = "umap", label = FALSE) +
  ggtitle("UMAP Projection")

## ------------------------------
## 4. Optimal clustering resolution via silhouette
## ------------------------------

resolutions <- seq(0.2, 1.2, by = 0.2)
sil_scores <- numeric(length(resolutions))
embeddings <- Embeddings(xenium_cereb, reduction = "pca")[, 1:elbow_pc]

for (i in seq_along(resolutions)) {
  res <- resolutions[i]
  xenium_cereb <- FindClusters(xenium_cereb, resolution = res, verbose = FALSE)
  clusters <- as.numeric(Idents(xenium_cereb))
  sil <- silhouette(clusters, dist(embeddings))
  sil_scores[i] <- mean(sil[, 3])
}

best_idx <- which.max(sil_scores)
best_res <- resolutions[best_idx]
message("Optimal clustering resolution: ", best_res)

# Plot silhouette scores
sil_plot <- ggplot(data.frame(resolution = resolutions, silhouette = sil_scores),
                   aes(x = resolution, y = silhouette)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = best_res, linetype = "dashed", color = "red") +
  labs(title = "Silhouette Score vs Clustering Resolution",
       x = "Resolution", y = "Average Silhouette Score") +
  theme_minimal()

## ------------------------------
## 5. Run final clustering
## ------------------------------

xenium_cereb <- FindClusters(xenium_cereb, resolution = best_res)
umap_cluster_plot <- DimPlot(xenium_cereb, reduction = "umap", label = TRUE) +
  ggtitle(paste("UMAP with Clusters (Resolution =", best_res, ")"))

## ------------------------------
## 6. Export all plots to a PDF
## ------------------------------

pdf("Xenium_Clustering_plots.pdf", width = 8, height = 6)
print(elbow_plot)
print(umap_plot)
print(sil_plot)
print(umap_cluster_plot)
dev.off()


## ------------------------------
## 1. Ensure cluster identities are factors
## ------------------------------

xenium_cereb$seurat_clusters <- factor(Idents(xenium_cereb))

## ------------------------------
## 2. Global cluster plot using ImageFeaturePlot
## ------------------------------

# Global: all clusters on all FOVs together
global_cluster_plot <- ImageFeaturePlot(
  object = xenium_cereb,
  fov = "fov",  # field-of-view identifier in your object
  features = "seurat_clusters",
  cols = brewer.pal(n = max(as.numeric(xenium_cereb$seurat_clusters)), name = "Set3")
) + 
  scale_y_reverse() +
  ggtitle("Global Spatial Plot of Clusters")

## ------------------------------
## 3. Faceted cluster plot (one panel per cluster)
## ------------------------------

# Faceted plot: split by cluster
facet_cluster_plot <- ImageFeaturePlot(
  object = xenium_cereb,
  fov = "fov",
  features = "seurat_clusters",
  cols = brewer.pal(n = max(as.numeric(xenium_cereb$seurat_clusters)), name = "Set3")
) +
  scale_y_reverse() +
  ggtitle("Spatial Plot Faceted by Cluster") +
  facet_wrap(~ seurat_clusters)

## ------------------------------
## 4. Export both plots to a single PDF
## ------------------------------

pdf("Xenium_Cluster_SpatialPlots.pdf", width = 12, height = 8)
print(global_cluster_plot)
print(facet_cluster_plot)
dev.off()

message("PDF export complete: Xenium_Clusters_SpatialPlots.pdf")

message("All plots saved to Xenium_analysis_plots.pdf")
