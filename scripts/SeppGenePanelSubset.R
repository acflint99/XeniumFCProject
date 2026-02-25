# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)

#Load the FC dataset----
Sepp = readRDS("/data/user/acflint/FC_published/SeppFC/Sepp_FC_filtered_newUMAP.rds")

xenium_genes <- readRDS("/data/user/acflint/FC_published/Xenium/xenium_5k_genes.rds")

genes_present <- intersect(xenium_genes, rownames(Sepp))

SeppSubset <- subset(
  Sepp,
  features = genes_present
)

length(genes_present)  # number of genes actually present in your object

#check UMAP
p <- DimPlot(SeppSubset, reduction = "umap", group.by = "cell_type_refined")
p

#export new Seurat object as .rds----
saveRDS(SeppSubset, file = "/data/user/acflint/FC_published/SeppFC/Sepp_FC_filtered_noNA_5kgenes.rds")


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
Idents(SeppSubset_newUMAP) <- "cell_type_refined"

# 7️⃣ Run UMAP
SeppSubset_newUMAP <- RunUMAP(SeppSubset_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p2 <- DimPlot(SeppSubset_newUMAP, reduction = "umap", group.by = "cell_type_refined")
p2
ggsave("/data/user/acflint/FC_published/SeppFC/SeppUMAP_FC_newclusters_newUMAP50_noNA_5kgenes.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(SeppSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  SeppSubset_newUMAP[[assay]]@scale.data <- matrix()
}

#export new Seurat object as .rds----
saveRDS(SeppSubset_newUMAP, file = "/data/user/acflint/FC_published/SeppFC/Sepp_FC_filtered_noNA_5kgenes_newUMAP.rds")


#check key cell type markers----
markers <- c(
  "EBF2",      # Purkinje
  "ROR2",
  "FOXP2",
  "PAX2",      # GABA
  "GAD1",
  "SOX2",      # NSC / VZ
  "HES1",
  "TFAP2A",    # VZ
  "MEIS2",     # GCP / NSC
  "ATOH1",     # GCP / RL
  "PAX6",
  "MKI67",
  "WNT2B",     # RL
  "INHBB",     
  "LTBP1",
  "OTX2",
  "EOMES",     # UBC
  "NEUROD1",   # GN
  "RELN",
  "SLC1A2",     # Glia
  "SOX9",
  "ADCY2",
  "GATA1",     # RBC
  "PHOX2B",    # Brainstem
  "HOXB4",
  "LMX1B",
  "HTR2C",     # Choroid
  "FOXJ1",     # Ependymal
  "P2RY12",    # Microglia
  "PDGFRA",    # OPC
  "OLIG1",
  "PDGFRB",    # Pericytes
  "RGS5",
  "CSPG4",
  "CLDN5",     # Endothelial
  "PECAM1",
  #"FOXC1",      # Meninges
  "SLC7A11"
)
# markers <- c(
#   "NES", "MEIS2", "PAX5", "CHST8", "TSHZ1"
# )

p3 <- DotPlot(
  SeppSubset_newUMAP,
  features = markers,
  group.by = "cell_type_refined"  # change if needed
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")
p3
ggsave("/data/user/acflint/FC_published/SeppFC/SeppDotPlot_FC_filtered_newclusters_5kgenes_newUMAP50_markers.pdf", plot = p3, width = 10, height = 6)
