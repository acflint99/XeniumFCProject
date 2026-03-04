# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(SingleCellExperiment)


#Load the FC dataset----
Sepp = readRDS("/data/user/acflint/FC_published/SeppFC/Sepp_hum_sce_final.rds")

#convert to Seurat object
SeuratObj <- as.Seurat(Sepp, counts = "umi", data = "umi")

#activate cell_type as identities
Idents(SeuratObj) <- SeuratObj$cell_type

#convert ENSG -> HGNC symbols, set default assay to "RNA", and remove "originalexp" assay ----
assay_name <- "originalexp"

counts_mat <- GetAssayData(SeuratObj, assay = assay_name, layer = "counts")
data_mat   <- GetAssayData(SeuratObj, assay = assay_name, layer = "data")

ensg_ids <- rownames(counts_mat)

symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensg_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

new_names <- ifelse(is.na(symbols), ensg_ids, symbols)
new_names <- make.unique(new_names)

rownames(counts_mat) <- new_names
rownames(data_mat)   <- new_names

new_assay <- CreateAssayObject(counts = counts_mat)
new_assay <- SetAssayData(new_assay, layer = "data", new.data = data_mat)
new_assay@meta.features$ENSEMBL_ID <- ensg_ids

SeuratObj[["RNA"]] <- new_assay
DefaultAssay(SeuratObj) <- "RNA"

# Optional: remove old assay
SeuratObj[["originalexp"]] <- NULL

rm(new_assay)
rm(counts_mat)
rm(data_mat)

# plot UMAP ----
p1 <- DimPlot(SeuratObj, reduction = "umap2d", group.by = "cell_type")+
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))
p1
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppUMAP_origClusters.pdf", plot = p1, width = 7, height = 6)


#subset fetal cerebellum cells----
SeuratObj_FC <- subset(SeuratObj, subset = Stage %in% c("11 wpc", "17 wpc", "20 wpc", "7 wpc", "8 wpc", "9 wpc"))

#plot UMAP again
p2 <- DimPlot(SeuratObj_FC, reduction = "umap2d", group.by = "cell_type")+
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))
p2
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppUMAP_FC_origClusters.pdf", plot = p2, width = 7, height = 6)

# #check proportions of clusters----
# cluster_counts <- table(Idents(SeuratObj_FC))
# 
# cluster_props <- prop.table(cluster_counts)
# 
# df <- as.data.frame(cluster_counts)
# colnames(df) <- c("Cluster", "Count")
# df$Proportion <- df$Count / sum(df$Count)
# 
# ggplot(df, aes(x = Cluster, y = Proportion)) +
#   geom_bar(stat = "identity") +
#   ylab("Proportion of cells") +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )



#get cluster markers for mixed clusters----
#markers_GC_UBC <- FindMarkers(SeuratObj, ident.1 = "GC/UBC", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#remove clusters----
SeuratObj_FC_filtered <- subset(SeuratObj_FC, idents = c("GABA_MB", "isth_N", "isthmic_neuroblast", "glut_DN", "erythroid", "NTZ_mixed", "noradrenergic"), invert = TRUE)


#rename clusters with simple/harmonized names----
SeuratObj_FC_filtered$clusters_refined <- SeuratObj_FC_filtered$cell_type

SeuratObj_FC_filtered$clusters_refined <- dplyr::recode(
  SeuratObj_FC_filtered$clusters_refined,
  "astroglia" = "Glia",
  "GABA_DN" = "GABA",
  "GC" = "Granule",
  "immune" = "Immune",
  "interneuron" = "GABA",
  "meningeal" = "Meninges",
  "mural/endoth" = "Endothelial",
  "NTZ_neuroblast" = "RL",
  "oligo" = "OPC",
  "parabrachial" = "VZ",
  "Purkinje" = "Purkinje",
  "UBC" = "UBC",
  "VZ_neuroblast" = "VZ",
  "GC/UBC" = "RL"
  )

# remove cells in "NA" cluster
SeuratObj_FC_filtered_noNA <- subset(SeuratObj_FC_filtered, subset = !is.na(clusters_refined))

# plot UMAP again with new cluster labels----
p4 <- DimPlot(SeuratObj_FC_filtered_noNA, reduction = "umap2d", group.by = "clusters_refined")
p4
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppUMAP_FC_newClusters.pdf", plot = p4, width = 7, height = 6)

#redo PCA & UMAP----
# 1️⃣ Normalize data
SeuratObj_FC_filtered_noNA_newUMAP <- NormalizeData(SeuratObj_FC_filtered_noNA, normalization.method = "LogNormalize", scale.factor = 10000)

# 2️⃣ Find variable features
SeuratObj_FC_filtered_noNA_newUMAP <- FindVariableFeatures(SeuratObj_FC_filtered_noNA_newUMAP, selection.method = "vst", nfeatures = 2000)

# 3️⃣ Scale data
SeuratObj_FC_filtered_noNA_newUMAP <- ScaleData(SeuratObj_FC_filtered_noNA_newUMAP, features = rownames(SeuratObj_FC_filtered_noNA_newUMAP))

# 4️⃣ Run PCA
SeuratObj_FC_filtered_noNA_newUMAP <- RunPCA(SeuratObj_FC_filtered_noNA_newUMAP, features = VariableFeatures(SeuratObj_FC_filtered_noNA_newUMAP))

# 5️⃣ Find neighbors
SeuratObj_FC_filtered_noNA_newUMAP <- FindNeighbors(SeuratObj_FC_filtered_noNA_newUMAP, dims = 1:50)

# retain previous cluster identities
Idents(SeuratObj_FC_filtered_noNA_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
SeuratObj_FC_filtered_noNA_newUMAP <- RunUMAP(SeuratObj_FC_filtered_noNA_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP
p5 <- DimPlot(SeuratObj_FC_filtered_noNA_newUMAP, reduction = "umap", group.by = "clusters_refined")
p5
ggsave("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellPlots/SeppUMAP_FC_newClusters_newUMAP50v1.pdf", plot = p5, width = 7, height = 6)

#remove scale.data for all assays to reduce file size----
# Get assay names
assay_names <- names(SeuratObj_FC_filtered_noNA_newUMAP@assays)

# Loop over assays and clear scale.data
for (assay in assay_names) {
  SeuratObj_FC_filtered_noNA_newUMAP[[assay]]@scale.data <- matrix()
}


# #check proportions of clusters----
# cluster_counts2 <- table(Idents(SeuratObj_FC_filtered_noNA_newUMAP))
# 
# cluster_props2 <- prop.table(cluster_counts2)
# 
# df2 <- as.data.frame(cluster_counts2)
# colnames(df2) <- c("Cluster", "Count")
# df2$Proportion <- df2$Count / sum(df2$Count)
# 
# ggplot(df2, aes(x = Cluster, y = Proportion)) +
#   geom_bar(stat = "identity") +
#   ylab("Proportion of cells") +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )

#export new Seurat object as .rds----
saveRDS(SeuratObj_FC_filtered_noNA, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters.rds")
saveRDS(SeuratObj_FC_filtered_noNA_newUMAP, file = "/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters_newUMAPv1.rds")

