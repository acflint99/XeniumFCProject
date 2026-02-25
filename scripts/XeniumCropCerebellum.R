library(here)
library(Seurat)

XeniumCropCerebellum <- function(sample_name,
                                 fov = "fov",
                                 segmentations = "cell",
                                 flip_xy = TRUE,
                                 cell_stat_file = "cerebellum_cells_stats.csv",
                                 output_folder = "outputs") {
  
  # ---- 1. Validate paths FIRST ----
  sample_path <- here("data", "FCXeniumProject", sample_name)
  
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
#subset to cerebellum cells  
  cells_present <- intersect(cells_to_keep, colnames(xenium_obj))
  
  if (length(cells_present) == 0) {
    stop("None of the Cell IDs match cells in the Xenium object.")
  }
  
  xenium_cereb <- subset(xenium_obj, cells = cells_present)
  
  # ---- 5. Ensure output folder exists ----
  out_dir <- here(output_folder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ---- 6. Generate and save plot ----
  # Set plot size options
  options(repr.plot.width = 8, repr.plot.height = 8)
  
  p <- ImageFeaturePlot(xenium_cereb, fov = "fov", features = c("nCount_Xenium"), cols = c("black", "white"), max.cutoff = "q95") + scale_y_reverse()
  
  # Save as PNG
  ggsave(
    filename = paste0(sample_name, "_crbllm_nCount_FeatPlot.png"),
    plot = p,
    path = out_dir,
    width = 8,
    height = 8,
    units = "in",
    dpi = 300
  )
  
  message("Cerebellum plot saved to: ", file.path(out_dir, paste0(sample_name, "_cerebellum_plot.png")))
  
  
  return(xenium_cereb)
}