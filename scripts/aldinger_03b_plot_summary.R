# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

# Load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(here)
library(future) # New: for parallelization

source(here("scripts", "color_palette.R")) # Adjust path as necessary

# --- HPC OPTIMIZATION SETUP ---
# With 440GB RAM, we can be aggressive. 
# Using 12-16 workers is the "sweet spot" to avoid overhead while maximizing speed.
plan("sequential")
     
# Set a 50GB memory limit for passing objects between cores
options(future.globals.maxSize = 50 * 1024^3) 

# --- Load the FC dataset ---
AldingerSubset_newUMAP <- readRDS(here("outputs", "SingleCellRDS", "Aldinger_newClusters_newUMAPv2_5k.rds"))


# Lock the cell type order based on your master palette
AldingerSubset_newUMAP$clusters_refined <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = celltype_order
)

# Plot UMAP
# Validate colors first (using your custom function)
validate_palette(unique(AldingerSubset_newUMAP$clusters_refined))

p2 <- DimPlot(AldingerSubset_newUMAP, reduction = "umap", group.by = "clusters_refined", cols = cluster_colors)
p2
Cairo::CairoTIFF(
  filename = here("outputs", "SingleCellPlots", "AldingerUMAP_newClusters_5k_newUMAP50v2.tif"),
  width = 7,
  height = 6,
  units = "in",
  res = 600
)
print(p2)
grDevices::dev.off()


# --- check key cell type markers ---
# Set Idents and Order
Idents(AldingerSubset_newUMAP) <- factor(
  AldingerSubset_newUMAP$clusters_refined, 
  levels = rev(celltype_order)
)

# markers <- c("OTX2", "EOMES", "RELN", "FOXP2", "PAX2", "TNC",
#                                      "OLIG1", "FOXC1", "CLDN5", "P2RY12")

p3 <- DotPlot(
  AldingerSubset_newUMAP,
  features = markers
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("High-Specificity Marker Expression by Cluster")

p3
Cairo::CairoTIFF(
  filename = here("outputs", "SingleCellPlots", "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.tif"),
  width = 14,
  height = 6,
  units = "in",
  res = 600
)
print(p3)
grDevices::dev.off()
ggplot2::ggsave(
  filename = here("outputs", "SingleCellPlots", "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"),
  plot = p3, device = grDevices::cairo_pdf, width = 14, height = 6
)

#### --- Proportion Plot: clusters_refined by PCW --- ####
target_clusters <- c("RL", "UBC", "Granule", "Purkinje", "GABA") 
target_ages <- c("9 PCW", "10 PCW", "11 PCW", "12 PCW", "14 PCW", "17 PCW", "18 PCW", "20 PCW")

# 1. Calculate proportions
prop_data <- AldingerSubset_newUMAP@meta.data %>%
  filter(clusters_refined %in% target_clusters) %>%
  filter(age %in% target_ages) %>%
  group_by(age, clusters_refined) %>%
  tally() %>%
  group_by(age) %>%
  mutate(percent = (n / sum(n) * 100))

# 2. Generate the Stacked Bar Plot
p4 <- ggplot(prop_data, aes(x = age, y = percent, fill = clusters_refined)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
  scale_fill_manual(values = cluster_colors) + 
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(
    x = "Age (PCW)",
    y = "Relative Proportion",
    fill = "Cell Type"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black")
  )

# Show and Save
p4
Cairo::CairoTIFF(
  filename = here("outputs", "SingleCellPlots", "AldingerBarPlot_newClusters_5k_newUMAP50v2_rmPCW16PCW21_prop.tif"),
  width = 8,
  height = 6,
  units = "in",
  res = 600
)
print(p4)
grDevices::dev.off()
ggplot2::ggsave(
  filename = here("outputs", "SingleCellPlots", "AldingerBarPlot_newClusters_5k_newUMAP50v2_rmPCW16PCW21_prop.pdf"),
  plot = p4, device = grDevices::cairo_pdf, width = 8, height = 6
)
