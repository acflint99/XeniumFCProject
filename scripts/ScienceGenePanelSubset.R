# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(Cairo)

# Source the color palette ----
source(here("scripts", "color_palette.R"))

# Define output directories
plot_path <- here("outputs", "SciencePlots")
RDS_path <- here("outputs", "ScienceRDS")

# Ensure directories exist before saving any files
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
if (!dir.exists(RDS_path)) dir.create(RDS_path, recursive = TRUE)

# Load the FC dataset ----
Science = readRDS(file.path(RDS_path, "Science_newClusters_newUMAPv1.rds"))

xenium_genes <- readRDS(here("inputs", "xenium_5k_genes.rds"))

genes_present <- intersect(xenium_genes, rownames(Science))

ScienceSubset <- subset(
  Science,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

# Apply celltype order from color_palette.R as factor levels
ScienceSubset$clusters_refined <- factor(
  ScienceSubset$clusters_refined,
  levels = intersect(celltype_order, unique(ScienceSubset$clusters_refined))
)

# check UMAP
p <- DimPlot(ScienceSubset, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_FC_newClusters_5k.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p)
dev.off()

# export new Seurat object as .rds ----
saveRDS(ScienceSubset, file = file.path(RDS_path, "Science_newClusters_UMAPv1_5k.rds"))


## redo PCA & UMAP ----
# 1️⃣ Normalize data
ScienceSubset_newUMAP <- NormalizeData(ScienceSubset, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
ScienceSubset_newUMAP <- FindVariableFeatures(ScienceSubset_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
ScienceSubset_newUMAP <- ScaleData(ScienceSubset_newUMAP, features = rownames(ScienceSubset_newUMAP))

# 4️⃣ Run PCA
ScienceSubset_newUMAP <- RunPCA(ScienceSubset_newUMAP, features = VariableFeatures(ScienceSubset_newUMAP))

# 5️⃣ Find neighbors
ScienceSubset_newUMAP <- FindNeighbors(ScienceSubset_newUMAP, dims = 1:50)

# retain previous cluster identities & order
ScienceSubset_newUMAP$clusters_refined <- factor(
  ScienceSubset_newUMAP$clusters_refined,
  levels = intersect(celltype_order, unique(ScienceSubset_newUMAP$clusters_refined))
)
Idents(ScienceSubset_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
ScienceSubset_newUMAP <- RunUMAP(ScienceSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP with color_palette.R colors
p2 <- DimPlot(ScienceSubset_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "ScienceUMAP_newClusters_5k_newUMAP50v2.tiff"), 
          width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

# remove scale.data for all assays to reduce file size ----
# Get assay names
assay_names <- names(ScienceSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  ScienceSubset_newUMAP[[assay]]@scale.data <- matrix()
}

# export new Seurat object as .rds ----
saveRDS(ScienceSubset_newUMAP, file = file.path(RDS_path, "Science_newClusters_newUMAPv2_5k.rds"))

# Reverse the celltype order for the DotPlot so the y-axis renders top-to-bottom
Idents(ScienceSubset_newUMAP) <- factor(
  ScienceSubset_newUMAP$clusters_refined, 
  levels = rev(intersect(celltype_order, unique(ScienceSubset_newUMAP$clusters_refined)))
)

p3 <- DotPlot(
  ScienceSubset_newUMAP,
  features = markers
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

CairoTIFF(filename = file.path(plot_path, "ScienceDotPlot_newClusters_5k_newUMAP50v2_markers.tiff"), 
          width = 10, height = 6, units = "in", res = 600)
print(p3)
dev.off()

# 2. Save Dotplot as PDF
ggplot2::ggsave(
  filename = file.path(plot_path, "ScienceDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"), 
  plot = p3, 
  width = 10, 
  height = 6,
  device = grDevices::cairo_pdf
) 
