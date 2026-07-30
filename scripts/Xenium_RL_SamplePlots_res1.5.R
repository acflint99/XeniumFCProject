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
  obj$RL_subcluster <- factor(obj$RL_subcluster, levels = rl_subcluster_order)
  
  # Export Cluster Statistics
  stats <- obj@meta.data %>%
    group_by(RL_subcluster) %>%
    summarise(cell_count = n()) %>%
    mutate(percentage = (cell_count / sum(cell_count)) * 100)
  
  write.csv(stats, file.path(sample_plot_dir, paste0(sample_id, "_cluster_counts.csv")), row.names = FALSE)
  
  #----------------------------
  #HELPER: Internal Highlight Function (Using png())
  #----------------------------
  save_highlight <- function(seurat_obj, targets, title_suffix, filename) {
    
    group_name <- paste0("plot_", title_suffix)
    
    seurat_obj[[group_name]] <- ifelse(
      seurat_obj$RL_subcluster %in% targets,
      as.character(seurat_obj$RL_subcluster),
      "Other"
    )
    
    # Explicit factor ordering: Other first (bottom), targets later (top)
    level_order <- c("Other", as.character(targets))
    
    seurat_obj[[group_name]] <- factor(
      seurat_obj[[group_name]][,1],
      levels = level_order
    )
    
    current_cols <- c("Other" = bg_color, rl_palette[targets])
    
    p <- ImageDimPlot(
      seurat_obj,
      group.by = group_name,
      cols = current_cols,
      size =2
    ) +
      ggtitle(paste(sample_id, title_suffix)) +
      theme(
        legend.title = element_blank(),
        legend.position = "right",
        plot.margin = margin(0, 0, 0, 0, "pt")
      )
    
    ggsave(
      filename = file.path(sample_plot_dir, filename),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300,
      bg = "black"
    )
  }
  
  
  # ----------------------------
  # GENERATE & SAVE PLOTS
  # ----------------------------

  # A. Global Spatial
  p_global <- ImageDimPlot(obj, # Removed the subset() here
                           group.by = "RL_subcluster",
                           cols = rl_palette,
                           size = 2,
                           na.value = alpha("gray20", 0.3)) +
    ggtitle(paste(sample_id, "RL Subclusters"))

  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Global_Spatial2.png")),
         plot = p_global, width = 10, height = 8, dpi = 300)

  # B. Group Highlights
  save_highlight(obj, c("RL VZ", "RL SVZ"),
                 "RL_Lineage", paste0(sample_id, "_RL_Spatial.png"))
  save_highlight(obj, c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
                 "Granule_Lineage", paste0(sample_id, "_Granule_Spatial2.png"))
  save_highlight(obj, c("RL VZ", "RL SVZ", "Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
                 "Granule&RL_Lineage", paste0(sample_id, "_Granule&RL_Spatial2.png"))
  save_highlight(obj, c("RL VZ", "RL SVZ", "Intermediate Progenitors", "Immature UBCs", "Mature UBCs"),
                 "UBC_Lineage", paste0(sample_id, "_UBC_Spatial2.png"))

  #C. Faceted Spatial
  coords <- GetTissueCoordinates(obj, type = "Xenium")
  plot_data <- cbind(coords, cluster = obj$RL_subcluster)
  design <- "ABCDE\nFGHIJ\nKLMN#"
  
  p_facet <- plot_data %>%
    filter(!is.na(cluster)) %>%
    ggplot(aes(x = y, y = x, color = cluster)) +
    geom_point(size = 0.3) +
    facet_manual(~cluster, design = design) + 
    scale_color_manual(values = rl_palette) +
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
  
  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Faceted_Clusters0.3.png")), 
         plot = p_facet, width = 15, height = 10, dpi = 300)
  
  # D. Marker DotPlot
  rev_levels <- factor(as.character(obj$RL_subcluster), levels = rev(rl_subcluster_order))
  names(rev_levels) <- Cells(obj)
  obj <- AddMetaData(obj, metadata = rev_levels, col.name = "subcluster_rev")

  Idents(obj) <- "subcluster_rev"
  obj_plot <- subset(obj, idents = rl_subcluster_order)

  p_dot <- DotPlot(obj_plot, features = rl_markers, cols = c("lightgrey", "red"),
                   dot.scale = 6, cluster.idents = FALSE) +
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 8, face = "italic"),
          axis.text.y = element_text(size = 10, face = "bold")) +
    ggtitle(paste(sample_id, "RL Markers"))

  ggsave(file.path(sample_plot_dir, paste0(sample_id, "_Markers_DotPlot.png")),
         plot = p_dot, width = 12, height = 7, dpi = 300)
  
  # Clean up memory
  message(paste("Finished sample:", sample_id))
  rm(obj, obj_plot, p_global, p_facet, p_dot, plot_data, coords, stats)
  gc()
}