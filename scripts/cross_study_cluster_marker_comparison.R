# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(viridis)
library(here)
col_fun <- colorRamp2(c(0, 0.6), viridis(2))

library(future)

source(here("scripts", "R", "config.R"))
config <- load_pipeline_config()

Aldinger <- readRDS(here(config$inputs$references$aldinger))
Sepp <- readRDS(here(config$inputs$references$sepp))
Science <- readRDS(here(config$inputs$references$science))

plot_dir <- here("outputs", "references", "cross_study", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# set active identities----
Idents(Aldinger) <- "clusters_refined"
Idents(Sepp) <- "clusters_refined"
Idents(Science) <- "clusters_refined"

# Use exactly the CPUs allocated by Slurm.
workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
plan("multisession", workers = workers)

# compute markers ----
markers_Aldinger <- FindAllMarkers(Aldinger, only.pos = TRUE)
markers_Sepp     <- FindAllMarkers(Sepp, only.pos = TRUE)
markers_Science  <- FindAllMarkers(Science, only.pos = TRUE)

plan("sequential")  # reset to default

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

# -----------------------------
# 1️⃣ Aldinger vs Sepp
# -----------------------------
heatmap_mat_ASe <- overlap_table_ASe %>%
  select(Aldinger_cluster, Sepp_cluster, jaccard) %>%
  pivot_wider(
    names_from = Sepp_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(heatmap_mat_ASe) <- heatmap_mat_ASe$Aldinger_cluster
heatmap_mat_ASe$Aldinger_cluster <- NULL
heatmap_mat_ASe <- as.matrix(heatmap_mat_ASe)
heatmap_mat_ASe <- heatmap_mat_ASe[
  sort(rownames(heatmap_mat_ASe)),
  sort(colnames(heatmap_mat_ASe))
]

htASe <- Heatmap(
  heatmap_mat_ASe,
  name = "Jaccard Index",
  row_title = "Aldinger Clusters",
  row_title_side = "right",
  column_title = "Sepp Clusters",
  column_title_side = "bottom",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  col = col_fun
)

# -----------------------------
# 2️⃣ Aldinger vs Science
# -----------------------------
heatmap_mat_ASc <- overlap_table_ASc %>%
  select(Aldinger_cluster, Science_cluster, jaccard) %>%
  pivot_wider(
    names_from = Science_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(heatmap_mat_ASc) <- heatmap_mat_ASc$Aldinger_cluster
heatmap_mat_ASc$Aldinger_cluster <- NULL
heatmap_mat_ASc <- as.matrix(heatmap_mat_ASc)
heatmap_mat_ASc <- heatmap_mat_ASc[
  sort(rownames(heatmap_mat_ASc)),
  sort(colnames(heatmap_mat_ASc))
]

htASc <- Heatmap(
  heatmap_mat_ASc,
  name = "Jaccard Index",
  row_title = "Aldinger Clusters",
  row_title_side = "right",
  column_title = "Science Clusters",
  column_title_side = "bottom",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  col = col_fun
)

# -----------------------------
# 3️⃣ Sepp vs Science
# -----------------------------
heatmap_mat_SS <- overlap_table_SS %>%
  select(Sepp_cluster, Science_cluster, jaccard) %>%
  pivot_wider(
    names_from = Science_cluster,
    values_from = jaccard,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(heatmap_mat_SS) <- heatmap_mat_SS$Sepp_cluster
heatmap_mat_SS$Sepp_cluster <- NULL
heatmap_mat_SS <- as.matrix(heatmap_mat_SS)
heatmap_mat_SS <- heatmap_mat_SS[
  sort(rownames(heatmap_mat_SS)),
  sort(colnames(heatmap_mat_SS))
]

htSS <- Heatmap(
  heatmap_mat_SS,
  name = "Jaccard Index",
  row_title = "Sepp Clusters",
  row_title_side = "right",
  column_title = "Science Clusters",
  column_title_side = "bottom",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  col = col_fun
)

# -----------------------------
# Export all three heatmaps to multi-page PDF
# -----------------------------
pdf(file.path(plot_dir, "PairwiseClusterJaccardHeatmaps.pdf"),
    width = 7, height = 6)

draw(htASe, newpage = TRUE)

draw(htASc, newpage = TRUE)

draw(htSS, newpage = TRUE)

dev.off()
