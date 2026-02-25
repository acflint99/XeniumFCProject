# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)


#Update Seurat object & save to new .rds----
remotes::install_version("Seurat", version = "4.4.0")
library(Seurat)
obj <- readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat.rds")

obj@reductions <- list()   # remove PCA, UMAP, tSNE
obj@graphs <- list()       # remove any graphs (optional, safer)
obj <- UpdateSeuratObject(obj)

DefaultAssay(obj) <- "integrated"

obj <- ScaleData(obj, verbose = FALSE)

obj <- RunPCA(obj, assay = "integrated", npcs = 75, verbose = FALSE)

obj <- RunUMAP(obj, reduction = "pca", dims = 1:75, verbose = FALSE)

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:75, verbose = FALSE)

obj <- FindClusters(obj, resolution = 1.5, verbose = FALSE)

DimPlot(obj, reduction = "umap", group.by = "figure_clusters")
table(obj$figure_clusters)


#add published UMAP dimensions to Seurat Object----
# read the gzipped TSV directly
umap_coords <- read.csv("/data/user/acflint/FC_published/AldingerFC/Seurat_UMAP.coords.tsv.gz.tmp", 
                        sep = "\t", 
                        row.names = 1, 
                        check.names = FALSE)

colnames(umap_coords) <- c("UMAP_1", "UMAP_2")


# 2️⃣ Subset to cells present in UMAP
common_cells <- intersect(colnames(counts_matrix), rownames(umap_coords))
counts_matrix <- counts_matrix[, common_cells]
metadata_df <- metadata_df[common_cells, ]
umap_coords <- umap_coords[common_cells, ]


# Keep only the cells present in both
common.cells <- intersect(colnames(obj), rownames(umap_coords))

# Subset Seurat object temporarily
obj_subset <- subset(obj, cells = common.cells)

# Reorder UMAP coords to match Seurat object
umap_subset <- umap_coords[common.cells, ]

# Create DimReduc object
umap.reduction <- CreateDimReducObject(
  embeddings = as.matrix(umap_subset),
  key = "UMAP_",
  assay = "integrated"
)

# Add to Seurat object
obj_subset[["umap"]] <- umap.reduction

# Quick check
DimPlot(obj_subset, reduction = "umap", label = TRUE)


# 5️⃣ Set cluster identities
Idents(obj_subset) <- obj_subset$figure_clusters

# 6️⃣ Plot UMAP
DimPlot(obj_subset, reduction = "umap", group.by = "figure_clusters")

saveRDS(obj_subset, "/data/user/acflint/FC_published/AldingerFC/Aldinger_seurat_updated.rds")


