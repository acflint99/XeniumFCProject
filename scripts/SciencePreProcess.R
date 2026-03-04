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
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/ScienceUMAP_Luoclusters.pdf", plot = p, width = 9, height = 6)

#check markers
#markers_Cycling <- FindMarkers(Science, ident.1 = "Cycling", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)%>% dplyr::arrange(desc(avg_log2FC)) %>% head(20) %>% rownames()

#remove clusters----
Science <- subset(Science, idents = c("RBC", "brainstem", "Cycling"), invert = TRUE)

#rename clusters with simple/harmonized names----
Science@meta.data$clusters_refined <- Science@meta.data$celltype

Science@meta.data$clusters_refined <- dplyr::recode(
  Science@meta.data$clusters_refined,
  "PC" = "Purkinje",
  "GCP1" = "Granule",
  "SPON1_NSC" = "Glia",
  "GABA.inter" = "GABA",
  "TCP?" = "VZ",
  "Glut.DN" = "Granule",
  "UBC" = "UBC",
  "TNC_NSC" = "Glia",
  "GCP-P" = "Granule",
  "Granule1" = "Granule",
  "RLVZ" = "RL",
  "VZP/eTCP" = "VZ",
  "GABA.LHX9" = "GABA",
  "RLSVZ" = "RL",
  "GPC2" = "Granule",
  "Immune" = "Immune",
  "Granule2" = "Granule",
  "GABA.CN" = "GABA",
  "34_NSC" = "Glia",
  "Endothelial" = "Endothelial",
  "Pericyte" = "Meninges",
  "OPC" = "OPC"
)

Idents(Science) <- "celltype"

#check markers
#markers_GlutDN <- FindMarkers(Science, ident.1 = "Glut.DN", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)%>% dplyr::arrange(desc(avg_log2FC)) %>% head(20) %>% rownames()

# plot UMAP again with new cluster labels----
p1 <- DimPlot(Science, reduction = "umap", group.by = "clusters_refined")
p1
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/ScienceUMAP_newclusters.pdf", plot = p1, width = 7, height = 6)

#save Seurat object w/ new cluster labels----
saveRDS(Science, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Science_newClusters.rds")

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
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/ScienceUMAP_newclusters_newUMAP50v1.pdf", plot = p2, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(Science_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  Science_newUMAP[[assay]]@scale.data <- matrix()
}


#save Seurat object w/ new cluster labels & new UMAP reductions----
saveRDS(Science_newUMAP, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Science_newClusters_newUMAPv1.rds")
