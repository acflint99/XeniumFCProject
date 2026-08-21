# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(Cairo) # Required for CairoTIFF

# Define and create output directories ----
plot_path <- here("outputs", "AldingerPlots")
RDS_path <- here("outputs", "AldingerRDS")

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_path, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# Source the color palette and aesthetics
# =====================================================
# Adjust the path inside here() if your file is in a subfolder (e.g., here("scripts", "color_palette.R"))
source(here("scripts", "color_palette.R"))

# =====================================================
# Load and Subset the FC dataset
# =====================================================
# Use here() to dynamically build paths relative to the project root
Aldinger = readRDS(here("outputs", "AldingerRDS", "Aldinger_newClusters_newUMAPv1.rds"))
xenium_genes <- readRDS(here("inputs", "xenium_5k_genes.rds"))

genes_present <- intersect(xenium_genes, rownames(Aldinger))

AldingerSubset <- subset(
  Aldinger,
  features = genes_present
)

# number of genes actually present in your object
cat("Number of overlapping genes:", length(genes_present), "\n") 

# export new Seurat object as .rds
saveRDS(AldingerSubset, file = file.path(RDS_path, "Aldinger_newClusters_newUMAPv1_5k.rds"))

# =====================================================
# Redo PCA & UMAP
# =====================================================
# 1️⃣ Normalize data
AldingerSubset_newUMAP <- NormalizeData(AldingerSubset, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
AldingerSubset_newUMAP <- FindVariableFeatures(AldingerSubset_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
AldingerSubset_newUMAP <- ScaleData(AldingerSubset_newUMAP, features = rownames(AldingerSubset_newUMAP))

# 4️⃣ Run PCA
AldingerSubset_newUMAP <- RunPCA(AldingerSubset_newUMAP, features = VariableFeatures(AldingerSubset_newUMAP))

# 5️⃣ Find neighbors
AldingerSubset_newUMAP <- FindNeighbors(AldingerSubset_newUMAP, dims = 1:50)

# 6️⃣ Run UMAP
AldingerSubset_newUMAP <- RunUMAP(AldingerSubset_newUMAP, dims = 1:50)

# =====================================================
# Plot UMAP using imported palette and save as TIFF
# =====================================================
# Set Idents to ensure factor levels match your imported celltype_order
Idents(AldingerSubset_newUMAP) <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = celltype_order
)

# Plot UMAP, mapping colors to your cluster_colors vector
p2 <- DimPlot(AldingerSubset_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors) +
  ggtitle("UMAP of Refined Clusters")

# Save using CairoTIFF
CairoTIFF(filename = file.path(plot_path, "AldingerUMAP_newClusters_5k_newUMAP50v2.tiff"),
          width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

# =====================================================
# Cleanup scale.data to reduce file size
# =====================================================
# Get assay names
assay_names <- names(AldingerSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  AldingerSubset_newUMAP[[assay]]@scale.data <- matrix()
}

# export new Seurat object as .rds
saveRDS(AldingerSubset_newUMAP, file = file.path(RDS_path, "Aldinger_newClusters_newUMAPv2_5k.rds"))

# =====================================================
# Check key cell type markers via DotPlot
# =====================================================
# Reverse the celltype_order so that they display top-to-bottom logically on the y-axis
Idents(AldingerSubset_newUMAP) <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = rev(celltype_order)
)

# Use the 'markers' list directly from color_palette.R
p3 <- DotPlot(
  AldingerSubset_newUMAP,
  features = markers
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

# 1. Save Dotplot using CairoTIFF
CairoTIFF(filename = file.path(plot_path, "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.tiff"),
          width = 10, height = 6, units = "in", res = 600)
print(p3)
dev.off()

# 2. Save Dotplot as PDF
ggsave(
  filename = file.path(plot_path, "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"), 
  plot = p3, 
  width = 10, 
  height = 6,
  device = "pdf"
)