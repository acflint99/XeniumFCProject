# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)

Aldinger <- readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat_updated.rds")

Aldinger[["RNA"]]@scale.data <- matrix()

p3 <- DimPlot(Aldinger, reduction = "umap", group.by = "figure_clusters")+
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))
p3
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/AldingerUMAP_origclusters.pdf", plot = p3, width = 7, height = 6)

#check proportions of clusters----
cluster_counts <- table(Idents(Aldinger))

cluster_props <- prop.table(cluster_counts)

df_clusters <- as.data.frame(cluster_counts)
colnames(df_clusters) <- c("Cluster", "Count")
df_clusters$Proportion <- df_clusters$Count / sum(df_clusters$Count)

ggplot(df_clusters, aes(x = Cluster, y = Proportion)) +
  geom_bar(stat = "identity") +
  ylab("Proportion of cells") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#remove clusters----
Aldinger_filtered <- subset(Aldinger, idents = c("19-Ast/Ependymal", "21-BS Choroid/Ependymal", "16-Pericytes", "17-Brainstem", "20-Choroid"), invert = TRUE)


#rename clusters with simple/harmonized names----
Aldinger_filtered@meta.data$clusters_refined <- Aldinger_filtered@meta.data$figure_clusters

Aldinger_filtered@meta.data$clusters_refined <- dplyr::recode(
  Aldinger_filtered@meta.data$clusters_refined,
  "01-PC" = "Purkinje",
  "02-RL" = "RL",
  "03-GCP" = "Granule",
  "04-GN" = "Granule",
  "05-eCN/UBC" = "UBC",
  "06-iCN" = "GABA",
  "07-PIP" = "GABA",
  "08-BG" = "Glia",
  "09-Ast" = "Glia",
  "10-Glia" = "Glia",
  "11-OPC" = "OPC",
  "12-Committed OPC" = "OPC",
  "13-Endothelial" = "Endothelial",
  "14-Microglia" = "Immune",
  "15-Meninges" = "Meninges",
  "18-MLI" = "GABA"
)


# plot UMAP again with new cluster labels----
p <- DimPlot(Aldinger_filtered, reduction = "umap", group.by = "clusters_refined")
p
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots//AldingerUMAP_newClusters.pdf", plot = p, width = 7, height = 6)

saveRDS(Aldinger_filtered, "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters.rds")

#redo PCA & UMAP----
# Switch to raw RNA assay
DefaultAssay(Aldinger_filtered) <- "RNA"
# 1️⃣ Normalize data
Aldinger_filtered <- NormalizeData(Aldinger_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
Aldinger_filtered <- FindVariableFeatures(Aldinger_filtered, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
Aldinger_filtered <- ScaleData(Aldinger_filtered, features = rownames(Aldinger_filtered))

# 4️⃣ Run PCA
Aldinger_filtered <- RunPCA(Aldinger_filtered, features = VariableFeatures(Aldinger_filtered))

# 5️⃣ Find neighbors
Aldinger_filtered <- FindNeighbors(Aldinger_filtered, dims = 1:50)

# retain previous cluster identities
Idents(Aldinger_filtered) <- "clusters_refined"

# 7️⃣ Run UMAP
Aldinger_filtered <- RunUMAP(Aldinger_filtered, dims = 1:50)

Aldinger_filtered[["RNA"]]@scale.data <- matrix()

# 8️⃣ Plot UMAP
p2 <- DimPlot(Aldinger_filtered, reduction = "umap", group.by = "clusters_refined")
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/AldingerUMAP_newClusters_newUMAPv1.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(Aldinger_filtered@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  Aldinger_filtered[[assay]]@scale.data <- matrix()
}


saveRDS(Aldinger_filtered, "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters_newUMAPv1.rds")

