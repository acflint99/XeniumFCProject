# ----------------------------
# Setup & Libraries
# ----------------------------
library(future)
library(Seurat)
library(spacexr) # RCTD package
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)
library(tidyr)
library(Cairo) 

# Increase the global object size limit
options(future.globals.maxSize = 400 * 1024^3)

# =========================================================
# Sample Mapping Logic (Slurm Array Support)
# =========================================================
sample_list <- c(
  # "GZFB4_X_G",
  # "FB124_X_G",
  # "FB198_X_G",
  # "FB328_1_X_G",
  # "FB330_1_X_G",
  # "FB78_X_G",
  # "GZFB5_X_G",
  # "GZFB_12_X_G_1",
  # "GZFB_12_X_G_2",
  # "GZFB_12_X_G_3",
  # "GZFB_12_X_G_4",
  # "GZFB_12_X_G_5",
  # "GZFB_1_X_G",
  # "GZFB_9_X_G_1",
  "GZFB_9_X_G_2",
  "GZFB_9_X_G_3"
)

# Get the Task ID from Slurm (e.g., --array=1-2)
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

if (is.na(task_id) || task_id < 1 || task_id > length(sample_list)) {
  stop("Error: Task ID is out of bounds or not provided.")
}

current_sample <- sample_list[task_id]
message("### Processing Annotation for Sample [", task_id, "]: ", current_sample, " ###")

# =========================================================
# Annotation Function
# =========================================================
annotate_xenium_from_ref <- function(xenium_obj, sample_name, reference_name = "Aldinger") {
  
  ## ----------------------------
  ## 0. Setup & Paths
  ## ----------------------------
  set.seed(42)
  
  # OPTIMIZATION: Keep sequential to let RCTD handle OpenMP C++ multithreading
  plan("sequential")
  
  source(here("scripts", "color_palette.R")) 
  
  output_root <- "outputs"
  
  plots_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Plots"))
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_Tables"))
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  ## ----------------------------
  ## 1 & 2. Load Pre-Computed Reference & Format Xenium
  ## ----------------------------
  message("Loading Pre-computed Single-Cell Reference for RCTD...")
  ref_path <- here(output_root, "SingleCellRDS", paste0(reference_name, "_RCTD_Reference.rds"))
  
  if (!file.exists(ref_path)) stop("RCTD Reference file not found: ", ref_path)
  rctd_ref <- readRDS(ref_path)
  
  DefaultAssay(xenium_obj) <- "Xenium"
  
  message("Formatting Spatial Query for RCTD...")
  x_counts <- GetAssayData(xenium_obj, assay = "Xenium", layer = "counts")
  nUMI_x   <- setNames(xenium_obj$nCount_Xenium, colnames(xenium_obj))
  
  # STRICT COORDINATE VALIDATION
  coords <- GetTissueCoordinates(xenium_obj)
  if("cell" %in% colnames(coords)) {
    rownames(coords) <- coords$cell
  } else if (!any(rownames(coords) %in% colnames(xenium_obj))) {
    stop("CRITICAL ERROR: Cell barcodes not found in coordinate rownames or 'cell' column.")
  }
  
  if (!all(c("x", "y") %in% colnames(coords))) {
    colnames(coords)[1:2] <- c("x", "y")
  }
  
  coords <- coords[, c("x", "y")] 
  coords <- coords[colnames(xenium_obj), ] # Force 1:1 barcode alignment
  
  if (any(is.na(coords$x))) {
    stop("CRITICAL ERROR: Coordinate matching failed. NAs introduced to spatial coordinates.")
  }
  
  puck <- SpatialRNA(coords, x_counts, nUMI_x)
  
  ## ----------------------------
  ## 3. Run RCTD Label Transfer
  ## ----------------------------
  message("Initializing RCTD Object...")
  # Max cores set to 8 to match the Slurm array requested CPUS
  myRCTD <- create.RCTD(puck, rctd_ref, max_cores = 4, test_mode = FALSE)
  
  message("Executing RCTD Mapping...")
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  
  ## ----------------------------
  ## 4. Extract Results & Integrate with Seurat
  ## ----------------------------
  message("Integrating RCTD results back into Seurat Object...")
  
  # In 'full' mode, RCTD outputs a matrix of weights instead of a results_df table.
  # Convert to a standard matrix to safely extract the highest values.
  weight_matrix <- as.matrix(myRCTD@results$weights)
  
  # Find the cell type with the highest weight for each spatial spot
  top_cell_type <- colnames(weight_matrix)[max.col(weight_matrix, ties.method = "first")]
  
  # Assign barcodes as names so Seurat can map them perfectly
  names(top_cell_type) <- rownames(weight_matrix)
  
  # Safely add the metadata to the Seurat object. 
  # Any spot RCTD filtered out will safely become NA, which Section 7 turns into 'Unknown'
  xenium_obj <- AddMetaData(xenium_obj, metadata = top_cell_type, col.name = "cluster_weighted")
  
  ## ----------------------------
  ## 5. Export Tables
  ## ----------------------------
  # We only need the primary count matrix now (doublet summary removed)
  count_matrix <- xenium_obj@meta.data %>%
    group_by(seurat_clusters, cluster_weighted) %>%
    tally() %>%
    pivot_wider(names_from = cluster_weighted, values_from = n, values_fill = 0) %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(count_matrix, 
            file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")), 
            row.names = FALSE)
  
  ## ----------------------------
  ## 6. Save RDS
  ## ----------------------------
  output_dir <- here(output_root, paste0("Xenium_", reference_name, "ABT_Res1.5_RDS"))
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(xenium_obj, file.path(output_dir, paste0(sample_name, "_", reference_name, "_annotated_RCTD.rds")))
  message("Successfully annotated with RCTD and saved: ", sample_name)
  
  ## ----------------------------
  ## 7. Visualizations
  ## ----------------------------
  xenium_obj$cluster_weighted <- as.character(xenium_obj$cluster_weighted)
  xenium_obj$cluster_weighted[is.na(xenium_obj$cluster_weighted)] <- "Unknown"
  
  # Identify all unique cell types predicted by RCTD
  predicted_types <- unique(xenium_obj$cluster_weighted)
  
  # Dynamically add missing colors (including 'Unknown')
  missing_types <- setdiff(predicted_types, names(cluster_colors))
  if (length(missing_types) > 0) {
    message("Assigning default colors to missing types: ", paste(missing_types, collapse = ", "))
    for (mt in missing_types) {
      cluster_colors[mt] <- "grey80" 
      if (exists("celltype_order") && !mt %in% celltype_order) {
        celltype_order <- c(celltype_order, mt)
      }
    }
  }
  
  xenium_obj$cluster_weighted <- factor(xenium_obj$cluster_weighted, levels = celltype_order)
  
  # --- UMAP Plot ---
  p_umap <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors) +
    ggtitle(paste(sample_name, "RCTD Annotations"))
  
  ggsave(here(plots_dir, paste0(sample_name, "_", reference_name, "_cluster_weighted_UMAP.tif")), 
         p_umap, device = "tiff", type = "cairo", width = 8, height = 6, dpi = 600)
  
  # --- Spatial Plots ---
  p_wei <- ImageDimPlot(xenium_obj, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "RCTD Predicted Cell Types"))
  p_wei$layers[[1]]$aes_params$stroke <- 0 
  
  ggsave(here(plots_dir, paste0(sample_name, "_Spatial_Weighted.tif")), 
         p_wei, device = "tiff", type = "cairo", width = 10, height = 10, dpi = 600)
  
  # --- Spatial Facet Plot ---
  coords_plot <- GetTissueCoordinates(xenium_obj) 
  if("cell" %in% colnames(coords_plot)) rownames(coords_plot) <- coords_plot$cell
  coords_plot <- coords_plot[colnames(xenium_obj), ]
  
  plot_data <- cbind(coords_plot, cluster = factor(as.character(xenium_obj$cluster_weighted), levels = celltype_order))
  
  p_facet <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = 0.1) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white")) +
    ggtitle(paste(sample_name, "RCTD Predicted Cell Types"))
  
  ggsave(here(plots_dir, paste0(sample_name, "_FacetSpatial.tif")), 
         p_facet, device = "tiff", type = "cairo", width = 12, height = 10, dpi = 600)
  
  # --- DotPlot ---
  if (exists("markers")) {
    existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium_obj)))
    existing_markers <- existing_markers[sapply(existing_markers, length) > 0]
    
    if (length(existing_markers) > 0) {
      Idents(xenium_obj) <- factor(xenium_obj$cluster_weighted, levels = rev(celltype_order))
      
      p_dot <- DotPlot(xenium_obj, features = existing_markers, assay = "Xenium") + 
        RotatedAxis() + scale_color_gradient(low = "lightgrey", high = "red") +
        ggtitle(paste(sample_name, "Markers (RCTD Annotated)"))
      
      ggsave(here(plots_dir, paste0(sample_name, "_DotPlot.tif")), 
             p_dot, device = "tiff", type = "cairo", width = 10, height = 6, dpi = 600)
      ggsave(here(plots_dir, paste0(sample_name, "_DotPlot.pdf")), 
             p_dot, device = cairo_pdf, width = 10, height = 6)
    }
  } else {
    warning("Object 'markers' not found in environment. Skipping DotPlot generation.")
  }
  
  ## ----------------------------
  ## 8. Return
  ## ----------------------------
  return(xenium_obj)
}

# =========================================================
# Execution
# =========================================================
input_file <- here("outputs", "Xenium_Res1.5_RDS", paste0(current_sample, "_CB_QC_cluster.rds"))

if (file.exists(input_file)) {
  seu <- readRDS(input_file)
  annotate_xenium_from_ref(seu, current_sample)
} else {
  stop("Input file not found: ", input_file)
}