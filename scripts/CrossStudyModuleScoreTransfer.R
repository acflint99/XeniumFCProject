# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

Aldinger <- readRDS("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Aldinger_newClusters_newUMAPv2_5k.rds")
Sepp <- readRDS("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Sepp_FC_newClusters_newUMAPv2_5k.rds")
Science <- readRDS("/home/acflint/R/Projects/XeniumFCProject/outputs/SingleCellRDS/Science_newClusters_5k_UMAPv2.rds")


# set active identities----
Idents(Aldinger) <- "clusters_refined"
Idents(Sepp) <- "clusters_refined"
Idents(Science) <- "clusters_refined"

# setup----
datasets <- list(
  Aldinger = Aldinger,
  Sepp     = Sepp,
  Science  = Science
)

cluster_col <- "clusters_refined"
top_n <- 50

# compute top markers per dataset----
get_top_markers <- function(seurat_obj, n = 50) {
  
  # assumes Idents() already set correctly outside
  
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  
  if (nrow(markers) == 0) {
    stop("No markers found. Check identities.")
  }
  
  markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = n) %>%
    summarise(genes = list(gene), .groups = "drop") %>%
    deframe()
}

marker_lists <- lapply(datasets, get_top_markers, n = top_n)

# restrict to shared genes----
shared_genes <- Reduce(
  intersect,
  lapply(datasets, rownames)
)

marker_lists <- lapply(marker_lists, function(cluster_list) {
  lapply(cluster_list, function(genes) {
    genes[genes %in% shared_genes]
  })
})

# automated pairwise transfer function
run_transfer <- function(source_name, target_name) {
  
  source_obj  <- datasets[[source_name]]
  target_obj  <- datasets[[target_name]]
  source_genes <- marker_lists[[source_name]]
  
  # Add module scores
  target_obj <- AddModuleScore(
    target_obj,
    features = source_genes,
    name = paste0(source_name, "_")
  )
  
  # Rename module score columns
  score_cols <- grep(paste0("^", source_name, "_"),
                     colnames(target_obj@meta.data),
                     value = TRUE)
  
  colnames(target_obj@meta.data)[
    match(score_cols, colnames(target_obj@meta.data))
  ] <- paste0(source_name, "_", names(source_genes))
  
  # Extract cluster identities directly
  cluster_ids <- Idents(target_obj)
  
  # Build cluster-level averages manually
  score_matrix <- target_obj@meta.data[, grep(paste0("^", source_name, "_"),
                                              colnames(target_obj@meta.data))]
  
  avg_scores <- aggregate(
    score_matrix,
    by = list(cluster = cluster_ids),
    FUN = mean
  )
  
  rownames(avg_scores) <- avg_scores$cluster
  avg_scores$cluster <- NULL
  
  return(as.matrix(avg_scores))
}


# run all 6 pairwise transfers automatically----
dataset_names <- names(datasets)

transfer_results <- list()

for (source in dataset_names) {
  for (target in dataset_names) {
    if (source != target) {
      
      key <- paste(source, "to", target, sep = "_")
      
      transfer_results[[key]] <- run_transfer(source, target)
    }
  }
}

# heatmap----
pheatmap(
  transfer_results$Aldinger_to_Sepp,
  main = "Module Score Transfer: Aldinger → Sepp",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)





