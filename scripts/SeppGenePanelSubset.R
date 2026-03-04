# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)

#Load the FC dataset----
Sepp = readRDS("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters_newUMAPv1.rds")

xenium_genes <- readRDS("/home/acflint/R/Projects/XeniumFCProject/inputs/xenium_5k_genes.rds")

genes_present <- intersect(xenium_genes, rownames(Sepp))

SeppSubset <- subset(
  Sepp,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

#check UMAP
p <- DimPlot(SeppSubset, reduction = "umap", group.by = "clusters_refined")
p

#export new Seurat object as .rds----
saveRDS(SeppSubset, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters_newUMAPv1_5k.rds")


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

# retain previous cluster identities
Idents(SeppSubset_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
SeppSubset_newUMAP <- RunUMAP(SeppSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p2 <- DimPlot(SeppSubset_newUMAP, reduction = "umap", group.by = "clusters_refined")
p2
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppUMAP_FC_newClusters_newUMAP50v2_5k.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(SeppSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  SeppSubset_newUMAP[[assay]]@scale.data <- matrix()
}

#export new Seurat object as .rds----
saveRDS(SeppSubset_newUMAP, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters_newUMAPv2_5k.rds")


#check key cell type markers----
markers <- c(
  "FOXP2", "CALB1", "DAB1",      # Purkinje
  "PAX2", "GAD1", "GAD2",        # GABA
  "MKI67", "LTBP1", "OTX2",      # RL
  "EOMES",                        # UBC
  "ATOH1", "PAX6", "NEUROD1", "RELN",         # Granule
  "PRDM13", "DLL1", "ASCL1",      # VZ
  "SOX9", "ADCY2",                # Glia
  "PDGFRA", "OLIG1",              # OPC
  "FOXC1", "SLC7A11",             # Meninges
  "CLDN5", "PECAM1",              # Endothelial
  "P2RY12"                        # Immune
)


# markers <- c(
#   "NES", "MEIS2", "PAX5", "CHST8", "TSHZ1"
# )

# 1. Define the desired cell type order
celltype_order <- c(
  "Purkinje", "GABA", "RL", "UBC", "Granule", "VZ", 
  "Glia", "OPC", "Meninges", "Endothelial", "Immune"
)

# 2. Set Idents to your refined clusters
Idents(SeppSubset_newUMAP) <- factor(
  SeppSubset_newUMAP$clusters_refined, 
  levels = rev(celltype_order)
)

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
p3
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppDotPlot_FC_newClusters_5k_newUMAP50v2_markers.pdf", plot = p3, width = 10, height = 6)
