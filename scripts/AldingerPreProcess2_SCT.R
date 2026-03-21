# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
#library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(here) # Load here

# --- Load Data ---
# Since this path starts with /data/, it is likely outside your project folder.
# We'll keep it absolute unless it's inside your project directory.
Aldinger <- readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat_updated.rds")

# Efficiently clear scale data if present
Aldinger[["RNA"]]@scale.data <- matrix()

p3 <- DimPlot(Aldinger, reduction = "umap", group.by = "figure_clusters") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))
p3

# --- Proportions and Filtering ---
cluster_counts <- table(Idents(Aldinger))
df_clusters <- as.data.frame(cluster_counts)
colnames(df_clusters) <- c("Cluster", "Count")
df_clusters$Proportion <- df_clusters$Count / sum(df_clusters$Count)

Aldinger_filtered <- subset(Aldinger, idents = c("19-Ast/Ependymal", "21-BS Choroid/Ependymal", "16-Pericytes", "17-Brainstem", "20-Choroid"), invert = TRUE)

# --- Rename clusters ---
Aldinger_filtered@meta.data$clusters_refined <- dplyr::recode(
  Aldinger_filtered@meta.data$figure_clusters,
  "01-PC" = "Purkinje", "02-RL" = "RL", "03-GCP" = "Granule", "04-GN" = "Granule",
  "05-eCN/UBC" = "UBC", "06-iCN" = "GABA", "07-PIP" = "GABA", "08-BG" = "Glia",
  "09-Ast" = "Glia", "10-Glia" = "Glia", "11-OPC" = "OPC", "12-Committed OPC" = "OPC",
  "13-Endothelial" = "Endothelial", "14-Microglia" = "Immune", "15-Meninges" = "Meninges",
  "18-MLI" = "GABA"
)

Idents(Aldinger_filtered) <- "clusters_refined"

# --- Redo PCA & UMAP using SCTransform ---
Aldinger_filtered <- SCTransform(Aldinger_filtered, 
                                 assay = "RNA", 
                                 new.assay.name = "SCT",
                                 verbose = TRUE)

Aldinger_filtered <- RunPCA(Aldinger_filtered, verbose = FALSE)
Aldinger_filtered <- FindNeighbors(Aldinger_filtered, dims = 1:50)
Aldinger_filtered <- RunUMAP(Aldinger_filtered, dims = 1:50)

# --- Save Plots and Objects using here() ---

# Plot 1: Initial Clusters (from earlier in the script logic)
ggsave(here("outputs", "SingleCellPlots", "Aldinger_origclusters_UMAP.pdf"), 
       plot = p3, width = 7, height = 6)

# Plot 2: New SCT Clusters
p2 <- DimPlot(Aldinger_filtered, reduction = "umap", group.by = "clusters_refined")
ggsave(here("outputs", "SingleCellPlots", "AldingerUMAP_newClusters_SCT_newUMAPv1.pdf"), 
       plot = p2, width = 7, height = 6)

# --- Clean up and Save ---
assay_names <- names(Aldinger_filtered@assays)
for (assay in assay_names) {
  Aldinger_filtered[[assay]]@scale.data <- matrix()
}

# Save RDS
saveRDS(Aldinger_filtered, 
        here("outputs", "SingleCellRDS", "Aldinger_newClusters_SCT_newUMAPv1.rds"))