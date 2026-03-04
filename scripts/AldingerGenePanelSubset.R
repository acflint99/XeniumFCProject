# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)

#Load the FC dataset----
Aldinger = readRDS("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters_newUMAPv1.rds")

xenium_genes <- readRDS("/home/acflint/R/Projects/XeniumFCProject/inputs/xenium_5k_genes.rds")

genes_present <- intersect(xenium_genes, rownames(Aldinger))

AldingerSubset <- subset(
  Aldinger,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

#check UMAP
p <- DimPlot(AldingerSubset, reduction = "umap", group.by = "clusters_refined")
p

#export new Seurat object as .rds----
saveRDS(AldingerSubset, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters_newUMAPv1_5k.rds")


##redo PCA & UMAP----
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


# 7️⃣ Run UMAP
AldingerSubset_newUMAP <- RunUMAP(AldingerSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p2 <- DimPlot(AldingerSubset_newUMAP, reduction = "umap", group.by = "clusters_refined")
p2
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/AldingerUMAP_newClusters_5k_newUMAP50v2.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(AldingerSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  AldingerSubset_newUMAP[[assay]]@scale.data <- matrix()
}

#export new Seurat object as .rds----
saveRDS(AldingerSubset_newUMAP, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters_newUMAPv2_5k.rds")

#check key cell type markers----
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
Idents(AldingerSubset_newUMAP) <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = rev(celltype_order)
)

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
p3
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.pdf", plot = p3, width = 10, height = 6)

