# 1. INITIALIZATION
library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(ggh4x)

# Load your new palette and order
source(here("scripts", "color_palette.R"))

generate_spatial_reports <- function(sample_id, input_rds_path, base_plot_dir) {
  
  message(paste("--- Processing Sample:", sample_id, "---"))
  
  # Add this to the very top of your function
  on.exit(while (!is.null(dev.list())) dev.off(), add = TRUE)
  
  # 1. Setup Output Directory
  sample_plot_dir <- file.path(base_plot_dir, sample_id)
  if (!dir.exists(sample_plot_dir)) dir.create(sample_plot_dir, recursive = TRUE)
  
  # 2. Load Data
  obj <- readRDS(input_rds_path)
  obj$PC_subcluster <- factor(obj$PC_subcluster, levels = PC_subcluster_order)
  
  # Export Cluster Statistics
  stats <- obj@meta.data %>%
    group_by(PC_subcluster) %>%
    summarise(cell_count = n()) %>%
    mutate(percentage = (cell_count / sum(cell_count)) * 100)
  
  write.csv(stats, file.path(sample_plot_dir, paste0(sample_id, "_cluster_counts.csv")), row.names = FALSE)
  
  
  # ----------------------------
  # GENERATE & SAVE PLOTS
  # ----------------------------

  # A. Global Spatial
  # 1. Identify which cells have NAs vs assigned clusters
  is_na <- is.na(obj@meta.data$PC_subcluster)
  
  # 2. Reorder the object: NAs first, then the clusters
  # This ensures clusters are drawn over the gray dots
  obj_ordered <- obj[, c(which(is_na), which(!is_na))]
  
  # 3. Plot using the reordered object
  p_global <- ImageDimPlot(obj_ordered, 
                           group.by = "PC_subcluster",
                           cols = PC_palette,
                           size = 0.8,
                           na.value = alpha("gray20", 0.3)) +
    ggtitle(paste(sample_id, "PC Subclusters"))

  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Global_Spatial.png")),
         plot = p_global, width = 10, height = 8, dpi = 300)

  #C. Faceted Spatial
  coords <- GetTissueCoordinates(obj, type = "Xenium")
  plot_data <- cbind(coords, cluster = obj$PC_subcluster)
  design <- "ABCD\nFGHI\nJKL#\nNOP#"
  
  p_facet <- plot_data %>%
    filter(!is.na(cluster)) %>%
    ggplot(aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.1) +
    facet_manual(~cluster, design = design) + 
    scale_color_manual(values = PC_palette) +
    scale_y_reverse() + 
    coord_fixed() + 
    theme_void() +
    theme(
      legend.position = "none", 
      panel.background = element_rect(fill = "black", color = "black"),
      plot.background = element_rect(fill = "black", color = "black"), 
      strip.text = element_text(color = "white"),
      strip.placement = "inside",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Faceted_Clusters.png")), 
         plot = p_facet, width = 12, height = 10, dpi = 300)
  
  # D. Marker DotPlot
  rev_levels <- factor(as.character(obj$PC_subcluster), levels = rev(PC_subcluster_order))
  names(rev_levels) <- Cells(obj)
  obj <- AddMetaData(obj, metadata = rev_levels, col.name = "PC_subcluster_rev")

  Idents(obj) <- "PC_subcluster_rev"
  obj_plot <- subset(obj, idents = PC_subcluster_order)

  p_dot <- DotPlot(obj_plot, features = PC_markers, cols = c("lightgrey", "red"),
                   dot.scale = 6, cluster.idents = FALSE) +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 8, face = "italic"),
          axis.text.y = element_text(size = 10, face = "bold")) +
    ggtitle(paste(sample_id, "PC Markers"))

  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Markers_DotPlot.png")),
         plot = p_dot, width = 12, height = 7, dpi = 300)
  
  # Clean up memory
  message(paste("Finished sample:", sample_id))
  rm(obj, obj_plot, p_global, p_facet, p_dot, plot_data, coords, stats)
  gc()
}