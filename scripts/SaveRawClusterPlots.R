# Clear the environment
rm(list = ls())

xenium_obj <- readRDS(here("outputs", "GZFB_12_X_G_1_CB_QC_cluster.rds"))
sample_name = "GZFB_12_X_G_1"

xenium_obj$seurat_clusters <- factor(Idents(xenium_obj))
clusters <- levels(xenium_obj$seurat_clusters)
n_clusters <- length(clusters)


cluster_colors <- distinctColorPalette(n_clusters)
names(cluster_colors) <- clusters

global_cluster_plot <- ImageDimPlot(
  object = xenium_obj,
  fov = "fov",
  group.by = "seurat_clusters",
  cols = cluster_colors,
  size = 0.75,
) + scale_y_reverse() +
  ggtitle(paste(sample_name, "- Raw Clusters"))

# This line removes the stroke from the actual points in the plot
global_cluster_plot$layers[[1]]$aes_params$stroke <- 0

# Export Global Plot
png_file_global <- file.path(here("outputs"), paste0(sample_name, "_GlobalRawClustersSpatialPlot.png"))
# We use a high res (300-600 DPI) for Xenium to keep the dots sharp
png(png_file_global, width = 10, height = 10, units = "in", res = 300)
print(global_cluster_plot)
dev.off()

coords <- GetTissueCoordinates(xenium_obj)
metadata <- xenium_obj@meta.data
plot_data <- cbind(coords, cluster = as.character(metadata$seurat_clusters))

facet_cluster_plot <- ggplot(plot_data, aes(x = y, y = x, color = cluster)) +
  geom_point(size = 0.2) +
  facet_wrap(~cluster) +
  scale_color_manual(values = cluster_colors) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "black", color = NA),
        plot.background = element_rect(fill = "black", color = NA),
        strip.text = element_text(color = "white", face = "bold", margin = margin(t = 5, b = 5)),
        plot.title = element_text(color = "white", hjust = 0.5, size = 14)) +
  ggtitle(paste(sample_name, "- Raw Clusters"))


# Export Facet Plot
png_file_facet <- file.path(here("outputs"), paste0(sample_name, "_FacetRawClustersSpatialPlot.png"))
png(png_file_facet, width = 12, height = 10, units = "in", res = 300)
print(facet_cluster_plot)
dev.off()

message("Saved spatial plots PNGs")