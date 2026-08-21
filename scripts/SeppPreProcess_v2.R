# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
#library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(here)
library(Cairo)

# Source the color palette ----
# Ensure color_palette.R is saved in your project root, or update this path accordingly.
source(here("scripts", "color_palette.R"))

# Define and create output directories ----
plot_path <- here("outputs", "SeppPlots")
RDS_path <- here("outputs", "SeppRDS")

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_path, recursive = TRUE, showWarnings = FALSE)

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

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_origClusters.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p1)
dev.off()


#subset fetal cerebellum cells----
SeuratObj_FC <- subset(SeuratObj, subset = Stage %in% c("11 wpc", "17 wpc", "20 wpc", "7 wpc", "8 wpc", "9 wpc"))

#plot UMAP again
p2 <- DimPlot(SeuratObj_FC, reduction = "umap2d", group.by = "cell_type")+
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 10))

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_FC_origClusters.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p2)
dev.off()

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

# Apply celltype order from color_palette.R as factor levels
SeuratObj_FC_filtered_noNA$clusters_refined <- factor(
  SeuratObj_FC_filtered_noNA$clusters_refined, 
  levels = intersect(celltype_order, unique(SeuratObj_FC_filtered_noNA$clusters_refined))
)

# plot UMAP again with new cluster labels and color_palette.R colors ----
p4 <- DimPlot(SeuratObj_FC_filtered_noNA, reduction = "umap2d", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_FC_newClusters.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p4)
dev.off()

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

# retain previous cluster identities & order
SeuratObj_FC_filtered_noNA_newUMAP$clusters_refined <- factor(
  SeuratObj_FC_filtered_noNA_newUMAP$clusters_refined,
  levels = intersect(celltype_order, unique(SeuratObj_FC_filtered_noNA_newUMAP$clusters_refined))
)
Idents(SeuratObj_FC_filtered_noNA_newUMAP) <- "clusters_refined"

# 7️⃣ Run UMAP
SeuratObj_FC_filtered_noNA_newUMAP <- RunUMAP(SeuratObj_FC_filtered_noNA_newUMAP, dims = 1:50)

# 8️⃣ Plot UMAP with color_palette.R colors ----
p5 <- DimPlot(SeuratObj_FC_filtered_noNA_newUMAP, reduction = "umap", group.by = "clusters_refined") +
  scale_color_manual(values = cluster_colors)

CairoTIFF(filename = file.path(plot_path, "SeppUMAP_FC_newClusters_newUMAP50v1.tiff"), width = 7, height = 6, units = "in", res = 600)
print(p5)
dev.off()

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
saveRDS(SeuratObj_FC_filtered_noNA, file = file.path(RDS_path, "Sepp_newClusters.rds"))
saveRDS(SeuratObj_FC_filtered_noNA_newUMAP, file = file.path(RDS_path, "Sepp_newClusters_newUMAPv1.rds"))
