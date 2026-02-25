library(here)
library(Seurat)

# Source the function
source(here("scripts", "XeniumCropCerebellum.R"))

# Load and crop a sample
xenium_cereb <- XeniumCropCerebellum("GZFB5_X_G")

####START FROM HERE####


# ===============================
# Xenium QC: Human Fetal Cerebellum
# ===============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Load Xenium object
# Replace with your path
# xenium.obj <- LoadXenium("path_to_xenium_output",
#                          fov = "fov",
#                          segmentations = "cell",
#                          flip.xy = TRUE)

DefaultAssay(xenium_cereb) <- "Xenium"

# 2. Calculate basic QC metrics
xenium_cereb$nCount_Xenium <- colSums(GetAssayData(xenium_cereb, slot = "counts"))
xenium_cereb$nFeature_Xenium <- colSums(GetAssayData(xenium_cereb, slot = "counts") > 0)

# 3. Visualize distributions
VlnPlot(
  xenium_cereb,
  features = c("nCount_Xenium", "nFeature_Xenium"),
  pt.size = 0.1,
  ncol = 2
)

FeatureScatter(xenium_cereb, "nCount_Xenium", "nFeature_Xenium")

hist(xenium_cereb$nCount_Xenium, breaks = 50, main = "Transcripts per Cell")
hist(xenium_cereb$nFeature_Xenium, breaks = 50, main = "Genes per Cell")

# 4. Suggested filtering (conservative)
# Adjust based on your distributions
xenium_cereb <- subset(
  xenium_cereb,
  subset =
    nCount_Xenium > 30 &
    nFeature_Xenium > 10 &
    nFeature_Xenium < 1200
)

# 5. Spatial QC visualization
SpatialFeaturePlot(xenium_cereb, features = c("nCount_Xenium", "nFeature_Xenium"))

# Optional: flag low-quality cells
xenium_cereb$lowQC <- xenium_cereb$nCount_Xenium < 30 | xenium_cereb$nFeature_Xenium < 10
SpatialDimPlot(xenium_cereb, group.by = "lowQC")

# 6. Save QC'd object
saveRDS(xenium_cereb, file = "xenium_fetalCB_QCed.rds")

cat("QC complete. Final number of cells:", ncol(xenium_cereb), "\n")








## Normalization
## library size based normalization, scale factor - median transcript count
xenium_cereb <- NormalizeData(xenium_cereb, scale.factor=median(xenium_cereb$nCount_Xenium))

## Find top 2k variable genes
xenium_cereb <- FindVariableFeatures(xenium_cereb)

## data scaling
xenium_cereb <- ScaleData(xenium_cereb)

## PCA for dimensionality reduction
xenium_cereb <- RunPCA(xenium_cereb, verbose = FALSE)

## 2D projection by UMAP
xenium_cereb <- RunUMAP(xenium_cereb, dims = 1:30, verbose = FALSE)

## clustering
xenium_cereb <- FindNeighbors(xenium_cereb, dims = 1:30)
xenium_cereb <- FindClusters(xenium_cereb, resolution = .5)

