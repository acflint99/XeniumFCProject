# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(here)
library(Cairo)

# Source the color palette ----
# Ensure color_palette.R is saved in your project root, or update this path accordingly.
source(here("scripts", "color_palette.R"))

# Define output directories
plot_path <- here("outputs", "AldingerPlots")
RDS_path <- here("outputs", "AldingerRDS")

# Ensure directories exist before saving any files
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
if (!dir.exists(RDS_path)) dir.create(RDS_path, recursive = TRUE)

Aldinger <- readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat_updated.rds")

Aldinger[["RNA"]]@scale.data <- matrix()

p3 <- DimPlot(Aldinger, reduction = "umap", group.by = "figure_clusters") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))

CairoTIFF(filename = file.path(plot_path, "AldingerUMAP_origclusters.tiff"), 
          width = 8, height = 6, units = "in", res = 600)
print(p3)
dev.off()

# check proportions of clusters----
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

# remove clusters----
Aldinger_filtered <- subset(Aldinger, idents = c("19-Ast/Ependymal", "21-BS Choroid/Ependymal", "16-Pericytes", "17-Brainstem", "20-Choroid"), invert = TRUE)

# rename clusters with simple/harmonized names----
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

# Apply celltype order from color_palette.R as factor levels
Aldinger_filtered$clusters_refined <- factor(
  Aldinger_filtered$clusters_refined, 
  levels = intersect(celltype_order, unique(Aldinger_filtered$clusters_refined))
)
Idents(Aldinger_filtered) <- "clusters_refined"

# plot UMAP again with new cluster labels and color_palette.R colors ----
p <- DimPlot(Aldinger_filtered, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "AldingerUMAP_newClusters.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p)
dev.off()

saveRDS(Aldinger_filtered, file.path(RDS_path, "Aldinger_newClusters.rds"))

# redo PCA & UMAP----
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

# retain previous cluster identities & order
Aldinger_filtered$clusters_refined <- factor(
  Aldinger_filtered$clusters_refined,
  levels = intersect(celltype_order, unique(Aldinger_filtered$clusters_refined))
)
Idents(Aldinger_filtered) <- "clusters_refined"

# 7️⃣ Run UMAP
Aldinger_filtered <- RunUMAP(Aldinger_filtered, dims = 1:50)

Aldinger_filtered[["RNA"]]@scale.data <- matrix()

# 8️⃣ Plot UMAP with color_palette.R colors
p2 <- DimPlot(Aldinger_filtered, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "AldingerUMAP_newClusters_newUMAPv1.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

# remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(Aldinger_filtered@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  Aldinger_filtered[[assay]]@scale.data <- matrix()
}

saveRDS(Aldinger_filtered, file.path(RDS_path, "Aldinger_newClusters_newUMAPv1.rds"))

# Verify all expected outputs exist
expected_files <- c(
  file.path(plot_path, "AldingerUMAP_origclusters.tiff"),
  file.path(plot_path, "AldingerUMAP_newClusters.tiff"),
  file.path(RDS_path, "Aldinger_newClusters.rds"),
  file.path(plot_path, "AldingerUMAP_newClusters_newUMAPv1.tiff"),
  file.path(RDS_path, "Aldinger_newClusters_newUMAPv1.rds")
)

missing_files <- expected_files[!file.exists(expected_files)]

if (length(missing_files) > 0) {
  stop("The following expected files were not created:\n", paste(missing_files, collapse = "\n"))
} else {
  message("All output files were successfully created and verified!")
}