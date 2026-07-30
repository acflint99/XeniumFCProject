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
output_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Giotto_RDS/"
dir.create(output_dir, showWarnings = FALSE)

# 2. HPC Instructions
options(giotto.core_number = 4)

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
  
  message(">> Tracer: Extracting Expression Matrix...")
  expr_mat <- as(LayerData(seurat_obj, assay = "Xenium", layer = "counts"), "sparseMatrix")
  
  # --- THE KMEANS FIX ---
  # Remove genes expressed in fewer than 2 cells to prevent k-means clustering crash
  expr_mat <- expr_mat[rowSums(expr_mat > 0) >= 2, ]
  
  message(">> Tracer: Extracting Coordinates...")
  coords <- GetTissueCoordinates(seurat_obj)
  
  message(">> Tracer: Formatting Coordinates...")
  spat_locs <- data.frame(
    sdimx = coords[[1]], 
    sdimy = coords[[2]], 
    row.names = rownames(coords)
  )
  
  message(">> Tracer: Extracting Metadata...")
  meta_df <- seurat_obj@meta.data
  meta_df$cell_ID <- rownames(meta_df)
  
  # --- THE FACTOR FIX ---
  # Strip Seurat's hidden factor levels so Giotto LR interactions do not crash
  meta_df$comb_subcluster <- as.character(meta_df$comb_subcluster)
  
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
  
  # 5. Spatial Network
  g_obj <- createSpatialNetwork(g_obj, name = "Delaunay_network", method = "Delaunay")
  
  # 6. Spatial Variable Genes
  spat_genes <- binSpect(g_obj, bin_method = 'kmeans')
  write.csv(spat_genes, paste0(output_dir, s, "_spatial_genes.csv"))
  
  # 7. Cell-Cell Proximity
  prox_stats <- cellProximityEnrichment(g_obj, 
                                        cluster_column = "comb_subcluster", 
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
  
  message(">> Tracer: Filtering LR Database for Xenium panel genes...")
  valid_genes <- rownames(expr_mat)
  LR_data_filtered <- LR_data[LR_data$ligand %in% valid_genes & LR_data$receptor %in% valid_genes, ]
  
  if (nrow(LR_data_filtered) == 0) {
    warning("No Ligand-Receptor pairs from CellChat were found in this specific Xenium slice!")
    lr_summary <- NULL
  } else {
    message(">> Tracer: Running Interactions with ", nrow(LR_data_filtered), " valid pairs (Sequentially)...")
    lr_summary <- spatCellCellcom(
      gobject = g_obj,
      spatial_network_name = "Delaunay_network",
      cluster_column = "comb_subcluster",
      feat_set_1 = as.character(LR_data_filtered$ligand),
      feat_set_2 = as.character(LR_data_filtered$receptor)
    )
  }
  
  saveRDS(lr_summary, paste0(output_dir, s, "_LR_interactions.rds"))
  
  # 9. Final Save and Cleanup
  saveRDS(g_obj, file = paste0(output_dir, s, "_Giotto.rds"))
  rm(seurat_obj, g_obj, lr_summary, expr_mat, spat_genes, prox_stats, LR_data, LR_data_filtered, CellChatDB.human); gc()
  
  message("--- Finished Sample: ", s, " ---")