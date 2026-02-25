# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

# Load the FC dataset
Science = readRDS("/data/user/acflint/FC_published/ScienceBraunFC/Luo/Science_9-15pcw-celltype.rds")

Idents(Science) <- "celltype"


#plot UMAP
p <- DimPlot(Science, reduction = "umap", group.by = "celltype")
p
ggsave("/data/user/acflint/FC_published/ScienceBraunFC/BraunLuoUMAP_Luoclusters.pdf", plot = p, width = 9, height = 6)


#rename clusters with simple/harmonized names----
Science@meta.data$clusters_refined <- Science@meta.data$celltype

Science@meta.data$clusters_refined <- dplyr::recode(
  Science@meta.data$clusters_refined,
  "PC" = "Purkinje",
  "GCP1" = "GCP",
  "SPON1_NSC" = "NSC",
  "GABA.inter" = "GABA",
  "TCP?" = "NSC",
  "Glut.DN" = "GCP",
  "UBC" = "UBC",
  "TNC_NSC" = "NSC",
  "GCP-P" = "GCP",
  "Granule1" = "GN",
  "RLVZ" = "VZ",
  "VZP/eTCP" = "VZ",
  "GABA.LHX9" = "GABA",
  "RLSVZ" = "RL",
  "brainstem" = "Brainstem",
  "GPC2" = "GCP",
  "Immune" = "Immune",
  "Granule2" = "GN",
  "GABA.CN" = "GABA",
  "RBC" = "RBC",
  "Cycling" = "NSC",
  "34_NSC" = "NSC",
  "Endothelial" = "Endothelial",
  "Pericyte" = "Pericytes",
  "OPC" = "OPC"
)

Idents(Science) <- "celltype"

#check markers
markers_GlutDN <- FindMarkers(Science, ident.1 = "Glut.DN", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)%>% dplyr::arrange(desc(avg_log2FC)) %>% head(20) %>% rownames()

# plot UMAP again with new cluster labels----
p1 <- DimPlot(Science, reduction = "umap", group.by = "clusters_refined")
p1
ggsave("/data/user/acflint/FC_published/ScienceBraunFC/BraunLuoUMAP_newclusters.pdf", plot = p1, width = 7, height = 6)

#save Seurat object w/ new cluster labels----
saveRDS(Science, file = "/data/user/acflint/FC_published/ScienceBraunFC/BraunLuo_newClusters.rds")

#redo PCA & UMAP----
# 1️⃣ Normalize data
Science_newUMAP <- NormalizeData(Science, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
Science_newUMAP <- FindVariableFeatures(Science_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
Science_newUMAP <- ScaleData(Science_newUMAP, features = rownames(Science_newUMAP))

# 4️⃣ Run PCA
Science_newUMAP <- RunPCA(Science_newUMAP, features = VariableFeatures(Science_newUMAP))

# 5️⃣ Find neighbors
Science_newUMAP <- FindNeighbors(Science_newUMAP, dims = 1:50)

# retain previous cluster identities
Idents(Science_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
Science_newUMAP <- RunUMAP(Science_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p2 <- DimPlot(Science_newUMAP, reduction = "umap", group.by = "clusters_refined")
p2
ggsave("/data/user/acflint/FC_published/ScienceBraunFC/BraunLuoUMAP_newclusters_newUMAP50.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(Science_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  Science_newUMAP[[assay]]@scale.data <- matrix()
}


#save Seurat object w/ new cluster labels & new UMAP reductions----
saveRDS(Science_newUMAP, file = "/data/user/acflint/FC_published/ScienceBraunFC/BraunLuo_newClusters_newUMAP.rds")




