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

# set active identities----
Idents(Aldinger) <- "clusters_refined"
Idents(Sepp) <- "cell_type_refined"
Idents(Science) <- "clusters_refined"

# compute markers ----
markers_Aldinger <- FindAllMarkers(Aldinger, only.pos = TRUE)
markers_Sepp     <- FindAllMarkers(Sepp, only.pos = TRUE)
markers_Science  <- FindAllMarkers(Science, only.pos = TRUE)

# filter significant markers ----
filter_markers <- function(markers) {
  markers %>%
    filter(p_val_adj < 0.05,
           avg_log2FC > 0.1) #changed from 0.25
}

markers_Aldinger_f <- filter_markers(markers_Aldinger)
markers_Sepp_f     <- filter_markers(markers_Sepp)
markers_Science_f  <- filter_markers(markers_Science)

#create cluster -> gene lists----
make_marker_list <- function(markers) {
  split(markers$gene, markers$cluster)
}

Aldinger_list <- make_marker_list(markers_Aldinger_f)
Sepp_list     <- make_marker_list(markers_Sepp_f)
Science_list  <- make_marker_list(markers_Science_f)

# marker overlap btwn Aldinger vs Sepp ----
overlap_results <- list()

for (a in names(Aldinger_list)) {
  for (s in names(Sepp_list)) {
    
    genes_A <- Aldinger_list[[a]]
    genes_S <- Sepp_list[[s]]
    
    intersect_genes <- intersect(genes_A, genes_S)
    union_genes     <- union(genes_A, genes_S)
    
    jaccard <- length(intersect_genes) / length(union_genes)
    
    overlap_results[[paste0("A",a,"_S",s)]] <- data.frame(
      Aldinger_cluster = a,
      Sepp_cluster = s,
      overlap_n = length(intersect_genes),
      jaccard = jaccard
    )
  }
}

overlap_table_ASe <- bind_rows(overlap_results)

# marker overlap btwn Sepp vs Science ----
overlap_results <- list()

for (a in names(Science_list)) {
  for (s in names(Sepp_list)) {
    
    genes_Sc <- Science_list[[a]]
    genes_Se <- Sepp_list[[s]]
    
    intersect_genes <- intersect(genes_Sc, genes_Se)
    union_genes     <- union(genes_Sc, genes_Se)
    
    jaccard <- length(intersect_genes) / length(union_genes)
    
    overlap_results[[paste0("A",a,"_S",s)]] <- data.frame(
      Science_cluster = a,
      Sepp_cluster = s,
      overlap_n = length(intersect_genes),
      jaccard = jaccard
    )
  }
}

overlap_table_SS <- bind_rows(overlap_results)

# marker overlap btwn Aldinger vs Science ----
overlap_results <- list()

for (a in names(Science_list)) {
  for (s in names(Aldinger_list)) {
    
    genes_Sc <- Science_list[[a]]
    genes_A <- Aldinger_list[[s]]
    
    intersect_genes <- intersect(genes_A, genes_Sc)
    union_genes     <- union(genes_A, genes_Sc)
    
    jaccard <- length(intersect_genes) / length(union_genes)
    
    overlap_results[[paste0("A",a,"_S",s)]] <- data.frame(
      Science_cluster = a,
      Aldinger_cluster = s,
      overlap_n = length(intersect_genes),
      jaccard = jaccard
    )
  }
}

overlap_table_ASc <- bind_rows(overlap_results)

# marker overlap btwn all 3 clusters ----
three_way_overlap <- list()

for (a in names(Aldinger_list)) {
  for (s in names(Sepp_list)) {
    for (sc in names(Science_list)) {
      
      genes_common <- Reduce(intersect, list(
        Aldinger_list[[a]],
        Sepp_list[[s]],
        Science_list[[sc]]
      ))
      
      if (length(genes_common) > 0) {
        three_way_overlap[[paste(a,s,sc,sep="_")]] <- data.frame(
          Aldinger_cluster = a,
          Sepp_cluster = s,
          Science_cluster = sc,
          shared_genes = length(genes_common)
        )
      }
    }
  }
}

overlap_table_3way <- bind_rows(three_way_overlap)


#identify best matching clusters ----
##High Jaccard (>0.2–0.3) suggests biological correspondence
overlap_table_ASe %>%
  arrange(desc(jaccard)) %>%
  head(10)

overlap_table_SS %>%
  arrange(desc(jaccard)) %>%
  head(10)

overlap_table_ASc %>%
  arrange(desc(jaccard)) %>%
  head(10)

overlap_table_3way %>%
  arrange(desc(jaccard)) %>%
  head(10)

#heatmap Ald vs Sepp ----
heatmap_mat_ASe <- overlap_table_ASe %>%
  select(Aldinger_cluster, Sepp_cluster, jaccard) %>%
  pivot_wider(
    names_from = Sepp_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()   # ← convert tibble → data.frame

# Set rownames BEFORE removing first column
rownames(heatmap_mat_ASe) <- heatmap_mat_ASe$Aldinger_cluster

# Remove cluster column
heatmap_mat_ASe$Aldinger_cluster <- NULL

# Convert to numeric matrix
heatmap_mat_ASe <- as.matrix(heatmap_mat_ASe)

# Sort row names alphabetically
heatmap_mat_ASe <- heatmap_mat_ASe[
  sort(rownames(heatmap_mat_ASe)),
  sort(colnames(heatmap_mat_ASe))
]

# Plot
pASe <- pheatmap(heatmap_mat_ASe,
         main = "Jaccard Index: Aldinger vs Sepp Cluster Marker Overlap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45
)
pASe
ggsave("/data/user/acflint/FC_published/AldingervSeppJaccardHeatmap.pdf", plot = pASe, width = 7, height = 6)

#heatmap Ald vs Science ----
heatmap_mat_ASc <- overlap_table_ASc %>%
  select(Aldinger_cluster, Science_cluster, jaccard) %>%
  pivot_wider(
    names_from = Science_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()   # ← convert tibble → data.frame

# Set rownames BEFORE removing first column
rownames(heatmap_mat_ASc) <- heatmap_mat_ASc$Aldinger_cluster

# Remove cluster column
heatmap_mat_ASc$Aldinger_cluster <- NULL

# Convert to numeric matrix
heatmap_mat_ASc <- as.matrix(heatmap_mat_ASc)

# Sort row names alphabetically
heatmap_mat_ASc <- heatmap_mat_ASc[
  sort(rownames(heatmap_mat_ASc)),
  sort(colnames(heatmap_mat_ASc))
]

# Plot
pASc <- pheatmap(heatmap_mat_ASc,
         main = "Jaccard Index: Aldinger vs Science Cluster Marker Overlap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45
)
pASc
ggsave("/data/user/acflint/FC_published/AldingervScienceJaccardHeatmap.pdf", plot = pASc, width = 7, height = 6)

#heatmap Sepp vs Science ----
heatmap_mat_SS <- overlap_table_SS %>%
  select(Sepp_cluster, Science_cluster, jaccard) %>%
  pivot_wider(
    names_from = Science_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()   # ← convert tibble → data.frame

# Set rownames BEFORE removing first column
rownames(heatmap_mat_SS) <- heatmap_mat_SS$Sepp_cluster

# Remove cluster column
heatmap_mat_SS$Sepp_cluster <- NULL

# Convert to numeric matrix
heatmap_mat_SS <- as.matrix(heatmap_mat_SS)

# Sort row names alphabetically
heatmap_mat_SS <- heatmap_mat_SS[
  sort(rownames(heatmap_mat_SS)),
  sort(colnames(heatmap_mat_SS))
]


# Plot
pSS <- pheatmap(heatmap_mat_SS,
                 main = "Jaccard Index: Sepp vs Science Cluster Marker Overlap",
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                 fontsize_row = 10,
                 fontsize_col = 10,
                 angle_col = 45
)
pSS
ggsave("/data/user/acflint/FC_published/SeppvScienceJaccardHeatmap.pdf", plot = pSS, width = 7, height = 6)


