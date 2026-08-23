# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
#library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(Cairo)

# Source the color palette ----
source(here("scripts", "color_palette.R"))

# Define and create output directories ----
plot_path <- here("outputs", "SeppPlots")
RDS_path <- here("outputs", "SeppRDS")

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_path, recursive = TRUE, showWarnings = FALSE)

#Load the FC dataset----
# Using the RDS_path to read the output generated from the previous script
Sepp = readRDS(file.path(RDS_path, "Sepp_newClusters_newUMAPv1.rds"))

# Update this path if xenium_5k_genes.rds is stored elsewhere relative to your project root
xenium_genes <- readRDS(here("inputs", "xenium_5k_genes.rds"))

genes_present <- intersect(xenium_genes, rownames(Sepp))

SeppSubset <- subset(
  Sepp,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

# Apply celltype order from color_palette.R as factor levels
SeppSubset$clusters_refined <- factor(
  SeppSubset$clusters_refined,
  levels = intersect(celltype_order, unique(SeppSubset$clusters_refined))
)
Idents(SeppSubset) <- "clusters_refined"

#check UMAP
p <- DimPlot(SeppSubset, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_FC_newClusters_5k.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p)
dev.off()

#export new Seurat object as .rds----
saveRDS(SeppSubset, file = file.path(RDS_path, "Sepp_newClusters_newUMAPv1_5k.rds"))


##redo PCA & UMAP----
# 1️⃣ Normalize data
SeppSubset_newUMAP <- NormalizeData(SeppSubset, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
SeppSubset_newUMAP <- FindVariableFeatures(SeppSubset_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
SeppSubset_newUMAP <- ScaleData(SeppSubset_newUMAP, features = rownames(SeppSubset_newUMAP))

# 4️⃣ Run PCA
SeppSubset_newUMAP <- RunPCA(SeppSubset_newUMAP, features = VariableFeatures(SeppSubset_newUMAP))

# 5️⃣ Find neighbors
SeppSubset_newUMAP <- FindNeighbors(SeppSubset_newUMAP, dims = 1:50)

# retain previous cluster identities and factor order
SeppSubset_newUMAP$clusters_refined <- factor(
  SeppSubset_newUMAP$clusters_refined,
  levels = intersect(celltype_order, unique(SeppSubset_newUMAP$clusters_refined))
)
Idents(SeppSubset_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
SeppSubset_newUMAP <- RunUMAP(SeppSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p2 <- DimPlot(SeppSubset_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_FC_newClusters_newUMAP50v2_5k.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(SeppSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  SeppSubset_newUMAP[[assay]]@scale.data <- matrix()
}

#export new Seurat object as .rds----
saveRDS(SeppSubset_newUMAP, file = file.path(RDS_path, "Sepp_newClusters_newUMAPv2_5k.rds"))


# 1. Ensure DotPlot y-axis is ordered according to your celltype_order (rev for bottom-up plotting)
Idents(SeppSubset_newUMAP) <- factor(
  SeppSubset_newUMAP$clusters_refined, 
  levels = rev(intersect(celltype_order, unique(SeppSubset_newUMAP$clusters_refined)))
)

# 2. Use the 'markers' list sourced directly from color_palette.R
# unlist() converts the named list of genes into a flat character vector for DotPlot
p3 <- DotPlot(
  SeppSubset_newUMAP,
  features = markers
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

CairoTIFF(filename = file.path(plot_path, "SeppDotPlot_FC_newClusters_5k_newUMAP50v2_markers.tiff"), width = 10, height = 6, units = "in", res = 600)
print(p3)
dev.off()

ggplot2::ggsave(filename = file.path(plot_path, "SeppDotPlot_FC_newClusters_5k_newUMAP50v2_markers.pdf"), plot = p3, device = grDevices::cairo_pdf, width = 10, height = 6)
