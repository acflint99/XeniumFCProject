library(here)
library(Seurat)
library(ggplot2)
library(Cairo)

XeniumCropCerebellum <- function(sample_name,
                                 fov = "fov",
                                 segmentations = "cell",
                                 flip_xy = TRUE,
                                 cell_stat_file = "GZFB20_5_cerebellum_cells_stats.csv", #edit for duo samples
                                 output_folder = "outputs") {
  
  # ---- 1. Validate paths FIRST ----
  sample_path <- here("data", "FCXeniumProject", "GZFB_20_X_G_9__7__5") #edit for duo samples
  
  if (!dir.exists(sample_path)) {
    stop("Sample folder does not exist: ", sample_path)
  }
  
  cell_stat_path <- file.path(sample_path, cell_stat_file)
  
  if (!file.exists(cell_stat_path)) {
    stop("Cell stat CSV not found: ", cell_stat_path)
  }
  
  # ---- 2. Read and validate CSV BEFORE loading Xenium ----
  cell_stats <- read.csv(
    cell_stat_path,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"Cell ID" %in% colnames(cell_stats)) {
    stop("Cell stat CSV must contain a 'Cell ID' column.")
  }
  
  cells_to_keep <- cell_stats$`Cell ID`
  
  if (length(cells_to_keep) == 0) {
    stop("No cells found in cell stat CSV.")
  }
  
  message("Cell stat CSV validated. Loading Xenium object...")
  
  # ---- 3. Load heavy Xenium object LAST ----
  xenium_obj <- LoadXenium(
    sample_path,
    fov = fov,
    segmentations = segmentations,
    flip.xy = flip_xy
  )
  
  # Subset to cerebellum cells  
  cells_present <- intersect(cells_to_keep, colnames(xenium_obj))
  
  if (length(cells_present) == 0) {
    stop("None of the Cell IDs match cells in the Xenium object.")
  }
  
  xenium_cereb <- subset(xenium_obj, cells = cells_present)
  
  # ---- 4. Inject cell statistics (like cell_area) into metadata ----
  # Match rows of cell_stats to the exact order of cells in xenium_cereb
  rownames(cell_stats) <- cell_stats$`Cell ID`
  cell_stats_matched <- cell_stats[colnames(xenium_cereb), ]
  
  # Standardize or assign cell area (check column name in your CSV, usually 'cell_area' or 'Area')
  if ("cell_area" %in% colnames(cell_stats_matched)) {
    xenium_cereb$cell_area <- cell_stats_matched$cell_area
  } else if ("Area" %in% colnames(cell_stats_matched)) {
    xenium_cereb$cell_area <- cell_stats_matched$Area
  } else {
    # Fallback search for any column containing 'area' case-insensitively
    area_col <- grep("area", colnames(cell_stats_matched), ignore.case = TRUE, value = TRUE)
    if(length(area_col) > 0) {
      xenium_cereb$cell_area <- cell_stats_matched[[area_col[1]]]
    } else {
      warning("Could not automatically locate an area column in cell stats CSV. Setting dummy values to prevent crashes.")
      xenium_cereb$cell_area <- 20 # fallback dummy value
    }
  }
  
  # ---- 5. Ensure output folder exists ----
  plot_dir <- here(output_folder, "XeniumCropPlots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # ---- 6. Generate and save plot ----
  options(repr.plot.width = 8, repr.plot.height = 8)
  
  p <- ImageFeaturePlot(xenium_cereb, fov = "fov", features = c("nCount_Xenium"), cols = c("black", "white"), max.cutoff = "q95") + scale_y_reverse()
  
  CairoTIFF(file.path(plot_dir, paste0(sample_name, "_CB_nCount_FeatPlot.tif")), 
            width = 8, height = 8, units = "in", res = 600)
  print(p)
  dev.off()
  
  message("Cerebellum plot saved to: ", file.path(plot_dir, paste0(sample_name, "_CB_nCount_FeatPlot.tif")))
  
  # -------------------------------
  # 11. Save cerebellum (CB) cropped object
  # -------------------------------
  rds_dir <- here("outputs", "XeniumRDS")
  if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
  
  output_file <- here("outputs", "XeniumRDS", paste0(sample_name, "_CB.rds"))
  saveRDS(xenium_cereb, file = output_file, compress = FALSE) 
  
  message("Saved RDS file to:", output_file)
  
  return(xenium_cereb)
}