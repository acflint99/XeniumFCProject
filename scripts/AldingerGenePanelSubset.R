# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)

#Load the FC dataset----
Aldinger = readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_filtered.rds")

xenium_genes <- readRDS("/data/user/acflint/FC_published/Xenium/xenium_5k_genes.rds")

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
saveRDS(AldingerSubset, file = "/data/user/acflint/FC_published/AldingerFC/Aldinger_filtered_5kgenes.rds")


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
ggsave("/data/user/acflint/FC_published/AldingerFC/AldingerUMAP_newclusters_5kgenes_newUMAP50.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(AldingerSubset_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  AldingerSubset_newUMAP[[assay]]@scale.data <- matrix()
}

#export new Seurat object as .rds----
saveRDS(AldingerSubset_newUMAP, file = "/data/user/acflint/FC_published/AldingerFC/Aldinger_filtered_5kgenes_newUMAP.rds")

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
  "FOXC1",      # Meninges
  "SLC7A11"
)

# markers <- c(
#   "LAMA2", "COL3A1", "COL1A2", "SLC7A11"
# )

p3 <- DotPlot(
  AldingerSubset_newUMAP,
  features = markers,
  group.by = "clusters_refined"  # change if needed
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")
p3
ggsave("/data/user/acflint/FC_published/AldingerFC/AldingerDotPlot_newclusters_5kgenes_newUMAP50_markers.pdf", plot = p3, width = 10, height = 6)

