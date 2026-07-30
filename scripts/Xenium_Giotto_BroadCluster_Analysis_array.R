# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
library(Seurat)
library(future)
library(here)

# --- THE PYTHON FIX ---
# Force reticulate to lock onto your Conda environment BEFORE loading Giotto
library(reticulate)
use_python("/home/acflint/.conda/envs/giotto_env/bin/python", required = TRUE)

# Now load Giotto
library(Giotto)

# Memory and Threading Management (Protects the S4 objects from corrupting)
options(future.globals.maxSize = 20000 * 1024^3)
options(giotto.warn_sequential = FALSE)
#plan("multicore", workers = 8)

# 1. Paths and Samples
data_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS/"
output_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_GiottoBroad_RDS/"
dir.create(output_dir, showWarnings = FALSE)

# 2. HPC Instructions
options(giotto.core_number = 7)

my_insts <- createGiottoInstructions(
  save_dir = output_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE
)

# 3. Accept Sample Name from Slurm
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No sample name provided. Please pass a sample name from the bash script.")
}
s <- args[1] # This grabs the exact sample name for this specific job array task

message("--- Starting Sample: ", s, " ---")
  
# Include the specific suffix in the regex pattern
file_to_load <- list.files(
  data_dir, 
  pattern = paste0("^", s, "_Ald_VZ_RL_QC_Subclusters\\.rds$"), 
  full.names = TRUE
)

if (length(file_to_load) == 0) {
  stop("Could not find file for sample: ", s)
} else if (length(file_to_load) > 1) {
  stop("Found multiple files for sample: ", s, " - please check your directory.")
}
  
message(">> Tracer: Loading RDS file...")
seurat_obj <- readRDS(file_to_load)

# 1. FILTER CLUSTERS
clusters_to_exclude <- c("Meninges", "Immune", "Endothelial")
seurat_obj <- subset(seurat_obj, subset = cluster_weighted %in% clusters_to_exclude, invert = TRUE)

# 2. EXTRACT DATA AFTER FILTERING (Crucial to keep dimensions in sync)
message(">> Tracer: Extracting Expression Matrix...")
# Extracting the sparse counts matrix
expr_mat <- as(LayerData(seurat_obj, assay = "Xenium", layer = "counts"), "sparseMatrix")

# Removing genes with < 2 counts to prevent k-means crash
expr_mat <- expr_mat[rowSums(expr_mat > 0) >= 2, ]

message(">> Tracer: Extracting Coordinates...")
coords <- GetTissueCoordinates(seurat_obj)

message(">> Tracer: Formatting Coordinates...")
spat_locs <- data.frame(
  sdimx = coords[[1]], 
  sdimy = coords[[2]], 
  row.names = rownames(coords)
)

# 3. EXTRACT METADATA
message(">> Tracer: Extracting Metadata...")
meta_df <- seurat_obj@meta.data
meta_df$cell_ID <- rownames(meta_df)

  # --- NEW: CALCULATE CLUSTER-SPECIFIC PERCENTAGES ---
  # This provides the 'ground truth' to stop transcripts from leaking into 0-expression clusters
  message(">> Tracer: Calculating Cluster Percentages...")
  
  # Get counts per cluster
  counts <- LayerData(seurat_obj, assay = "Xenium", layer = "counts")
  clusters <- as.character(meta_df$cluster_weighted)
  
  # Calculate % of cells expressing each gene per cluster
  # (Returns a list of vectors, one per cluster)
  cluster_list <- split(rownames(meta_df), clusters)
  pct_expressed_list <- lapply(cluster_list, function(cells) {
    rowSums(counts[, cells] > 0) / length(cells) * 100
  })
  
  # --- THE FACTOR FIX ---
  # Strip Seurat's hidden factor levels so Giotto LR interactions do not crash
  meta_df$cluster_weighted <- as.character(meta_df$cluster_weighted)
  
  message(">> Tracer: Building Giotto Object...")
  g_obj <- createGiottoObject(
    expression = expr_mat,
    spatial_locs = spat_locs,
    cell_metadata = meta_df,
    instructions = my_insts
  )
  
  message(">> Tracer: SUCCESS - Giotto Object Built!")
  
  # 4. Basic Normalization
  median_sc <- median(colSums(expr_mat))
  g_obj <- normalizeGiotto(g_obj, scalefactor = median_sc)
  
  # 5. Spatial Networks
  message(">> Tracer: Creating Delaunay (Contact) Network...")
  g_obj <- createSpatialNetwork(g_obj, name = "Delaunay_network", method = "Delaunay")
  
  message(">> Tracer: Creating kNN (Paracrine) Network (k=8)...")
  # k=8 captures a neighborhood of roughly 1-2 cell layers in Xenium
  g_obj <- createSpatialNetwork(g_obj, name = "kNN_network", method = "kNN", k = 8)
  
  # 6. Spatial Variable Genes
  # --- REFINED FILTERING LOGIC ---
  message(">> Tracer: Identifying genes expressed in >10% of at least one cluster...")
  
  gene_logical_check <- sapply(rownames(expr_mat), function(g) {
    any(sapply(pct_expressed_list, function(p) p[g] > 10))
  })
  
  genes_to_test <- names(which(gene_logical_check))
  
  message(">> Tracer: Testing ", length(genes_to_test), " genes that met the 10% threshold.")
  
  # CHANGE: 'feats' becomes 'subset_feats' for Giotto 4.x
  spat_genes <- binSpect(g_obj, 
                         bin_method = 'kmeans', 
                         subset_feats = genes_to_test) 
  
  write.csv(spat_genes, paste0(output_dir, s, "_spatial_genes.csv"))
  
  # 7. Cell-Cell Proximity
  prox_stats <- cellProximityEnrichment(g_obj, 
                                        cluster_column = "cluster_weighted", 
                                        spatial_network_name = "Delaunay_network",
                                        number_of_simulations = 100)
  saveRDS(prox_stats, paste0(output_dir, s, "_proximity_stats.rds"))
  
  # 8. Ligand-Receptor (Using the CellChat Human Database)
  message(">> Tracer: Downloading Human LR Database...")
  
  temp_rda <- tempfile(fileext = ".rda")
  download.file(
    url = "https://raw.githubusercontent.com/jinworks/CellChat/master/data/CellChatDB.human.rda",
    destfile = temp_rda,
    mode = "wb",
    quiet = TRUE
  )
  load(temp_rda)
  unlink(temp_rda)
  
  LR_data <- data.frame(
    ligand = CellChatDB.human$interaction$ligand,
    receptor = CellChatDB.human$interaction$receptor
  )
  
  message(">> Tracer: Filtering LR Database with Bio-Sanity Veto...")
  valid_genes <- rownames(expr_mat)
  
  # Filter the database for genes on your panel
  LR_data_filtered <- LR_data[LR_data$ligand %in% valid_genes & LR_data$receptor %in% valid_genes, ]
  
  # --- THE VETO FILTER ---
  # Only keep pairs where:
  # 1. The ligand is in > 5% of AT LEAST ONE cluster
  # 2. The receptor is in > 5% of AT LEAST ONE cluster
  # This kills background noise while preserving stage-specific signals
  keep_ligands <- names(which(sapply(valid_genes, function(g) any(sapply(pct_expressed_list, function(p) p[g] > 5)))))
  keep_receptors <- names(which(sapply(valid_genes, function(g) any(sapply(pct_expressed_list, function(p) p[g] > 5)))))
  
  LR_data_filtered <- LR_data_filtered[LR_data_filtered$ligand %in% keep_ligands & 
                                         LR_data_filtered$receptor %in% keep_receptors, ]
  
  # --- IMPROVED STEP 8: Dual-Network Communication ---
  if (nrow(LR_data_filtered) == 0) {
    lr_summary_list <- NULL
  } else {
    # 1. Update Metadata
    my_metadata <- pDataDT(g_obj)
    my_metadata$cluster_weighted <- as.character(my_metadata$cluster_weighted)
    g_obj <- addCellMetadata(g_obj, new_metadata = my_metadata, by_column = TRUE, column_cell_ID = "cell_ID")
    
    # 2. RUN DELAUNAY (Contact-Dependent)
    message(">> Tracer: Running Delaunay Signaling (min_obs=10)...")
    lr_delaunay <- spatCellCellcom(
      gobject = g_obj,
      spatial_network_name = "Delaunay_network",
      cluster_column = "cluster_weighted",
      feat_set_1 = as.character(LR_data_filtered$ligand),
      feat_set_2 = as.character(LR_data_filtered$receptor),
      min_observations = 10, # Higher threshold for reliability
      adjust_method = "fdr"
    )
    
    # 3. RUN kNN (Paracrine/Secreted)
    message(">> Tracer: Running kNN Signaling (min_obs=10)...")
    lr_knn <- spatCellCellcom(
      gobject = g_obj,
      spatial_network_name = "kNN_network",
      cluster_column = "cluster_weighted",
      feat_set_1 = as.character(LR_data_filtered$ligand),
      feat_set_2 = as.character(LR_data_filtered$receptor),
      min_observations = 10,
      adjust_method = "fdr"
    )
    
    # Combine results into a list for saving
    lr_summary_list <- list("delaunay" = lr_delaunay, "knn" = lr_knn)
  }
  
  saveRDS(lr_summary_list, paste0(output_dir, s, "_LR_interactions_dual.rds"))
  
  # 9. Final Save and Cleanup
  saveRDS(g_obj, file = paste0(output_dir, s, "_Giotto.rds"))
  rm(seurat_obj, g_obj, lr_summary_list, expr_mat, spat_genes, prox_stats, LR_data, LR_data_filtered, CellChatDB.human); gc()
  
  message("--- Finished Sample: ", s, " ---")