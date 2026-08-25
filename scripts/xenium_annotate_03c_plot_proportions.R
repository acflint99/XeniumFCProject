rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(purrr)
  library(dplyr)
  library(ggplot2)
  library(here)
})

# 1. Load Master Palette & Settings
source(here("scripts", "color_palette.R")) 

# SPECIFY CLUSTERS TO INCLUDE: 
# Edit this list to include only the clusters you want to see in the plot
target_clusters <- c("RL", "UBC", "Granule", "Purkinje", "GABA") 

# Create the plot directory. PCW is already stored in the consensus objects.
plot_dir <- here(
  "outputs", "xenium", "annotation", "03_consensus_labels", "plots", "proportions"
)
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# 3. List and Process RDS Files
data_path <- here("outputs", "xenium", "annotation", "03_consensus_labels", "rds")
sample_files <- list.files(
  path = data_path,
  pattern = "_Consensus_annotated\\.rds$",
  full.names = TRUE
)

all_data <- map_dfr(sample_files, function(f) {
  obj <- readRDS(f)
  if (!all(c("PCW", "consensus_label") %in% colnames(obj[[]]))) {
    stop("Consensus object lacks PCW or consensus_label metadata: ", f)
  }
  
  df <- obj@meta.data %>%
    mutate(PCW_label = as.character(PCW)) %>%
    select(PCW_label, consensus_label)
  
  rm(obj); gc()
  return(df)
})

# 4. Filter and Aggregate
plot_df <- all_data %>%
  # NEW: Only keep the clusters you listed in target_clusters
  filter(consensus_label %in% target_clusters) %>%
  
  # Clean PCW for numeric sorting
  mutate(PCW_num = as.numeric(gsub("PCW", "", PCW_label))) %>%
  
  # Calculate proportions based ONLY on the filtered set
  group_by(PCW_num, consensus_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(PCW_num) %>%
  mutate(relative_percent = (count / sum(count)) * 100) %>%
  ungroup()

# 5. Set Factor Levels for Plotting
plot_df$PCW_num <- factor(plot_df$PCW_num, levels = sort(unique(plot_df$PCW_num)))

# This ensures the legend and stack order match your master celltype_order
plot_df$consensus_label <- factor(plot_df$consensus_label, 
                                   levels = intersect(celltype_order, target_clusters))

# 6. Plot with your custom colors
p <- ggplot(plot_df, aes(x = PCW_num, y = relative_percent, fill = consensus_label)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
  scale_fill_manual(values = cluster_colors) + 
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    x = "Age (PCW)", 
    y = "Relative percent (of selected clusters)", 
    fill = "Cell Type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black")
  )

Cairo::CairoTIFF(
  filename = here(plot_dir, "XeniumConsensusABT_res1.5_ClusterPropPlot.tif"),
  width = 10,
  height = 6,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()
ggplot2::ggsave(
  filename = here(plot_dir, "XeniumConsensusABT_res1.5_ClusterPropPlot.pdf"),
  plot = p, device = grDevices::cairo_pdf, width = 10, height = 6
)
