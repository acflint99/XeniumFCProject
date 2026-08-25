# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# Load your new palette and order
source(here("scripts", "color_palette.R"))

input_path <- here("outputs", "XenAld_RL_Subclusters_Res1.5_RDS", "Xenium_RL_Res1.5_newSubclusters_4-3-26.rds")

obj <- readRDS(input_path)

# 1. Define a named list of your cluster sets
cluster_sets <- list(
  "Granule_Lineage" = c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
  "UBC_Lineage" = c("RL VZ", "RL SVZ", "Intermediate Progenitors", "Immature UBCs", "Mature UBCs")
)

target_samples <- c("FB328_1_X_G", "GZFB_12_X_G_3", "GZFB5_X_G",
                    "GZFB_1_X_G", "FB330_1_X_G", "FB78_X_G",
                    "GZFB4_X_G", "FB124_X_G") 

# 2. Iterate through the sets
for (set_name in names(cluster_sets)) {
  
  # Get the specific clusters for this iteration
  current_clusters <- cluster_sets[[set_name]]
  
  # Filter and Aggregate
  plot_df <- obj@meta.data %>%
    select(PCW, RL_subcluster, orig.ident) %>%
    filter(orig.ident %in% target_samples) %>%
    filter(RL_subcluster %in% current_clusters) %>%
    mutate(PCW_num = as.numeric(gsub("PCW", "", PCW))) %>%
    group_by(PCW_num, orig.ident, RL_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop") 
  
  # Skip if the dataframe is empty (e.g., clusters not found in target samples)
  if (nrow(plot_df) == 0) next
  
  # Set Factor Levels
  plot_df$PCW_num <- factor(plot_df$PCW_num, levels = sort(unique(plot_df$PCW_num)))
  plot_df$RL_subcluster <- factor(plot_df$RL_subcluster, 
                                  levels = intersect(rl_subcluster_order, current_clusters))
  
  # Plot
  p <- ggplot(plot_df, aes(x = PCW_num, y = cell_count, fill = RL_subcluster)) +
    geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
    scale_fill_manual(values = rl_palette) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste("Cluster Set:", set_name),
      x = "Age (PCW)", 
      y = "Number of Cells", 
      fill = "Cell Type",
      subtitle = paste("Included samples:", length(target_samples), "samples selected")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))
  
  # Save with a dynamic filename based on the set_name
  file_name <- paste0("XenAld_RL_", set_name, "_ClusterCountPlot.tif")
  
  Cairo::CairoTIFF(
    filename = here("outputs", "XenAld_RL_Subclusters_Res1.5_Plots", file_name),
    width = 10,
    height = 6,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
  ggplot2::ggsave(
    here("outputs", "XenAld_RL_Subclusters_Res1.5_Plots", sub("\\.tif$", ".pdf", file_name)),
    p, device = grDevices::cairo_pdf, width = 10, height = 6
  )
  
  message(paste("Successfully saved plot for:", set_name))
}
