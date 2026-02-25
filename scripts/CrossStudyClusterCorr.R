# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

Aldinger <- readRDS("/data/user/acflint/FC_published/AldingerFC/Aldinger_filtered_5kgenes_newUMAP.rds")
Sepp <- readRDS("/data/user/acflint/FC_published/SeppFC/Sepp_FC_filtered_noNA_5kgenes_newUMAP.rds")
Science <- readRDS("/data/user/acflint/FC_published/ScienceBraunFC/Science_filtered_5kgenes_newUMAP.rds")


############################################################
# 1. Compute cluster pseudobulk (normalized data slot)
############################################################

avg1 <- AverageExpression(
  Aldinger,
  group.by = "clusters_refined",
  assays = "RNA",
  slot = "data"
)$RNA

avg2 <- AverageExpression(
  Sepp,
  group.by = "cell_type_refined",
  assays = "RNA",
  slot = "data"
)$RNA

avg3 <- AverageExpression(
  Science,
  group.by = "clusters_refined",
  assays = "originalexp",
  slot = "data"
)$originalexp

############################################################
# 2️⃣ Restrict to shared genes across all datasets
############################################################
common_genes <- Reduce(intersect, list(
  rownames(avg1),
  rownames(avg2),
  rownames(avg3)
))

avg1 <- avg1[common_genes, ]
avg2 <- avg2[common_genes, ]
avg3 <- avg3[common_genes, ]

############################################################
# 3️⃣ Z-score genes within each dataset
############################################################
zscore_matrix <- function(mat){
  mat <- t(scale(t(mat)))
  mat[is.na(mat)] <- 0
  return(mat)
}

avg1_z <- zscore_matrix(avg1)
avg2_z <- zscore_matrix(avg2)
avg3_z <- zscore_matrix(avg3)

############################################################
# 4️⃣ Compute pairwise diagonal correlations with shared clusters
############################################################
# Shared clusters per pair
shared_12 <- intersect(colnames(avg1_z), colnames(avg2_z))
shared_13 <- intersect(colnames(avg1_z), colnames(avg3_z))
shared_23 <- intersect(colnames(avg2_z), colnames(avg3_z))

# Diagonal correlation vectors
diag_similarity_12 <- sapply(shared_12, function(cl) cor(avg1_z[, cl], avg2_z[, cl]))
diag_similarity_13 <- sapply(shared_13, function(cl) cor(avg1_z[, cl], avg3_z[, cl]))
diag_similarity_23 <- sapply(shared_23, function(cl) cor(avg2_z[, cl], avg3_z[, cl]))

# Pairwise summary tables
similarity_table_12 <- data.frame(Cluster=shared_12, Obj1_vs_Obj2=diag_similarity_12)
similarity_table_13 <- data.frame(Cluster=shared_13, Obj1_vs_Obj3=diag_similarity_13)
similarity_table_23 <- data.frame(Cluster=shared_23, Obj2_vs_Obj3=diag_similarity_23)

print(similarity_table_12)
print(similarity_table_13)
print(similarity_table_23)

# ############################################################
# # 5️⃣ Optional: full cross-cluster correlation heatmap-
# #    Works even if cluster sets differ
# ############################################################
# # Combine all clusters
# combined <- cbind(avg1_z, avg2_z, avg3_z)
# 
# # Compute correlation matrix
# # This produces a square matrix: clusters x clusters
# cor_all <- cor(combined)  # cor() computes correlation between columns by default
# 
# # Prepend dataset label to cluster names
# colnames(cor_all) <- c(
#   paste0("Aldinger", colnames(avg1_z)),
#   paste0("Sepp", colnames(avg2_z)),
#   paste0("Science", colnames(avg3_z))
# )
# 
# # Plot heatmap
# pheatmap(
#   cor_all,
#   main="All Clusters Across Datasets",
#   fontsize_row=8,
#   fontsize_col=8,
#   angle_col=45,
#   cluster_rows=TRUE,
#   cluster_cols=TRUE
# )

############################################################
# 6️⃣ Optional: pairwise heatmaps with axes labeled
############################################################
# Extract only relevant clusters
cols_13 <- colnames(avg3_z)
rows_13 <- colnames(avg1_z)

cor_13 <- cor(avg1_z[, intersect(rows_13, rows_13)], avg3_z[, intersect(cols_13, cols_13)])

# Alphabetical order of rows
sorted_rows <- sort(rownames(cor_13))

# Alphabetical order of columns
sorted_cols <- sort(colnames(cor_13))

# Reorder matrix
cor_13_sorted <- cor_13[sorted_rows, sorted_cols]


p13 <- pheatmap(
  cor_13_sorted,
  main="Aldinger vs Science Cluster Correlation",
  fontsize_row=8,
  fontsize_col=8,
  angle_col=45,
  cluster_rows=FALSE,
  cluster_cols=FALSE
)
p13
ggsave("/data/user/acflint/FC_published/AldingervScienceClusterHeatmap.pdf", plot = p13, width = 7, height = 6)


############################################################
# 6️⃣ Optional: pairwise heatmaps with axes labeled
############################################################
# Extract only relevant clusters
cols_23 <- colnames(avg3_z)
rows_23 <- colnames(avg2_z)

cor_23 <- cor(avg2_z[, intersect(rows_23, rows_23)], avg3_z[, intersect(cols_23, cols_23)])

# Alphabetical order of rows
sorted_rows <- sort(rownames(cor_23))

# Alphabetical order of columns
sorted_cols <- sort(colnames(cor_23))

# Reorder matrix
cor_23_sorted <- cor_23[sorted_rows, sorted_cols]


p23 <- pheatmap(
  cor_23_sorted,
  main="Sepp vs Science Cluster Correlation",
  fontsize_row=8,
  fontsize_col=8,
  angle_col=45,
  cluster_rows=FALSE,
  cluster_cols=FALSE
)
p23
ggsave("/data/user/acflint/FC_published/SeppvScienceClusterHeatmap.pdf", plot = p23, width = 7, height = 6)


############################################################
# 6️⃣ Optional: pairwise heatmaps with axes labeled
############################################################
# Extract only relevant clusters
cols_21 <- colnames(avg1_z)
rows_21 <- colnames(avg2_z)

cor_21 <- cor(avg2_z[, intersect(rows_21, rows_21)], avg1_z[, intersect(cols_21, cols_21)])

# Alphabetical order of rows
sorted_rows <- sort(rownames(cor_21))

# Alphabetical order of columns
sorted_cols <- sort(colnames(cor_21))

# Reorder matrix
cor_21_sorted <- cor_21[sorted_rows, sorted_cols]


p21 <- pheatmap(
  cor_21_sorted,
  main="Sepp vs Aldinger Cluster Correlation",
  fontsize_row=8,
  fontsize_col=8,
  angle_col=45,
  cluster_rows=FALSE,
  cluster_cols=FALSE
)
p21
ggsave("/data/user/acflint/FC_published/SeppvAldingerClusterHeatmap.pdf", plot = p21, width = 7, height = 6)
