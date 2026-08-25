# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(harmony)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

# Load your new palette and order
source(here("scripts", "color_palette.R"))

check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 1. PARALLELIZATION 
plan("multisession", workers = 8) 
options(future.globals.maxSize = 200 * 1024^3)

plot_path <- here("outputs", "XenAld_VZ_Subclusters_Res1.5_Plots")
if(!dir.exists(plot_path)) dir.create(plot_path)

table_path <- here("outputs", "XenAld_VZ_Subclusters_Res1.5_Tables")
if(!dir.exists(table_path)) dir.create(table_path, recursive = TRUE)

merged_path <- here("outputs", "XenAld_VZ_postQC_Res1.5_RDS", "Xenium_VZ_postQC_Res1.5_4-2-26.rds")
obj <- readRDS(merged_path)

# 2. SET THE IDENTITY
# Ensure we are looking at the resolution you liked best (e.g., 0.3)
Idents(obj) <- "Xenium_snn_res.0.5"

check_mem("STARTING FINDALLMARKERS")

# # 3. RUN FINDALLMARKERS (The Optimized Way)
# message(Sys.time(), ": Starting Marker Identification...")
# 
# all_markers <- FindAllMarkers(
#   obj,
#   only.pos = TRUE,          # Only look for upregulated genes (standard for cell types)
#   min.pct = 0.25,           # Gene must be in 25% of the cluster
#   logfc.threshold = 0.25,   # Minimum 1.28x fold change
#   test.use = "wilcox",      # Standard fast test
# 
#   # --- THE SPEED TRICK ---
#   # Downsampling to 1000 cells per cluster gives 99% the same results
#   # but runs 5x faster on massive Xenium datasets.
#   max.cells.per.ident = 1000
# )
# 
# check_mem("MARKERS COMPLETE")
# 
# # 4. FILTER & SAVE TOP 10
# top10_markers <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# write.csv(top10_markers,
#           file.path(table_path, "VZ_RawSubcluster_postQC_top10_Markers_Res0.5.csv"),
#           row.names = FALSE)
# 
# message("Marker analysis complete! Results saved to CSV.")


####input markers into Gemini####

# Ensure we are looking at the resolution you liked best (e.g., 0.3)
Idents(obj) <- "Xenium_snn_res.0.5"

# 1. Create a named vector for the mapping
# The names (0, 1, 2...) must match your current cluster IDs exactly
new_labels <- c(
  "0"  = "Early-born PCs",
  "1"  = "Maturing PCs",
  "2"  = "GABA Progenitors",
  "3"  = "BG",
  "4"  = "Cycling Cells",
  "5"  = "RG Progenitors",
  "6"  = "Astrocytes/Ependyma",
  "7"  = "eCN",
  "8"  = "Patterning PCs",
  "9"  = "iCN",
  "10" = "VZPs",
  "11" = "Golgi Cells",
  "12" = "OPCs",
  "13" = "Late-born PCs",
  "14" = "Late-born PCs",
  "15" = "GCPs",
  "16" = "MLIs",
  "17" = "eCN",
  "18" = "Maturing PCs"
  
)

obj <- RenameIdents(obj, new_labels)

# Final metadata assignment
obj$VZ_subcluster <- Idents(obj)

output_path <- here("outputs", "XenAld_VZ_Subclusters_Res1.5_RDS")
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
saveRDS(obj, file.path(output_path, "Xenium_VZ_Res1.5_newSubclusters_4-3-26.rds"), compress = FALSE)

# 3. APPLY STANDARDIZED ORDERING
# We use 'subcluster_order' from your color_palette.R script
obj$VZ_subcluster <- factor(Idents(obj), levels = vz_subcluster_order)
Idents(obj) <- "VZ_subcluster"

DefaultAssay(obj) <- "Xenium"

# Export Merged Cluster QC Statistics
merged_qc_stats <- obj@meta.data %>%
  group_by(VZ_subcluster) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(merged_qc_stats, 
          file.path(table_path, "XenAld_VZ_Subcluster_QC_Summary.csv"), 
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "XenAld_VZ_Subcluster_QC_Violins.pdf")
pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)
  
  p <- VlnPlot(
    obj, 
    features = feat, 
    group.by = "VZ_subcluster", 
    cols = vz_palette, 
    pt.size = 0      # Remove points for speed and smaller file size
  ) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"), # Explicitly style Y
      legend.position = "none"
    ) +
    ylab(feat) +
    ggtitle(paste("Xenium VZ Subcluster QC-", feat))
  
  # Printing the plot inside the loop sends it to the PDF device
  print(p)
}

# 5. Close the device to finalize the file
dev.off()

message("Individual QC plots saved to: ", pdf_path)

# 4. UMAP WITH CONSISTENT COLORS
p <- DimPlot(obj, 
             reduction = "umap_clean", 
             label = TRUE, 
             label.size = 4, 
             repel = TRUE,
             cols = vz_palette) + # Uses your defined palette
  ggtitle("Xenium Merged VZ Subclusters") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_VZ_Subcluster_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()

# 4. UMAP WITH "Other" Clusters removed
# Define the clusters you want to REMOVE
clusters_to_remove <- c("GCPs", "eCN", "Cycling Cells")

# Create a subset of the object
# ! means 'not', %in% checks if the identity is in your removal list
obj_subset <- subset(obj, idents = clusters_to_remove, invert = TRUE)

# Now run your plotting code on 'obj_subset' instead of 'obj'
p2 <- DimPlot(obj_subset, 
              reduction = "umap_clean", 
              label = TRUE, 
              label.size = 4, 
              repel = TRUE,
              cols = vz_palette) + 
  ggtitle("Xenium Merged VZ Subclusters (Filtered)") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_VZ_Subcluster_rmGCP,eCN,Cyc_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p2)
grDevices::dev.off()

rm(obj_subset)

# 5. DOTPLOT WITH REVERSED ORDER
# Create a temporary reversed factor for the Y-axis
# (Reverse the levels of the factor to get VZP at the top)
obj$VZ_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(vz_subcluster_order)) 
Idents(obj) <- "VZ_subcluster_rev"

# 3. Create the Plot
# Adding 'assay = "Xenium"' is the safest way to avoid multi-assay conflicts
p1 <- DotPlot(obj,
              features = vz_markers, 
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged VZ Subcluster Markers")

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_VZ_SubclusterMarker_DotPlot.tif")),
  width = 12,
  height = 8,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "XenAld_VZ_SubclusterMarker_DotPlot.pdf"),
  p1, device = grDevices::cairo_pdf, width = 12, height = 8
)

# ----------------------------------------------------------------              
# 6. GENERATE TOP 5 MARKERS HEATMAP
# ----------------------------------------------------------------
check_mem("STARTING HEATMAP GENERATION")

# 1. Identify Markers (if not already in environment from earlier)
# We use a lower max.cells.per.ident to speed up the heatmap calculation
message("Finding top 5 markers for heatmap...")
heatmap_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 500 # Faster for visualization purposes
)

# 2. Extract top 5 per cluster
top5_markers <- heatmap_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

top5_markers <- rev(top5_markers)

# 3. Scale the data for the specific markers only 
# (Necessary for DoHeatmap to show relative expression)
# Just run the plot on a temporary scaled version:
temp_obj <- ScaleData(obj, features = top5_markers, verbose = FALSE)

# 4. Create Heatmap
# We downsample the plot to 100 cells per group so the labels are readable
# 1. Prepare the plot without the default lines
p4 <- DoHeatmap(
  subset(temp_obj, downsample = 100), 
  features = top5_markers,
  group.by = "VZ_subcluster",
  group.colors = vz_palette,
  size = 4,          
  angle = 45,        
  draw.lines = TRUE,
  raster = FALSE
) + 
  scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
  guides(color = "none") +
  theme(
    axis.text.y = element_text(size = 6, face = "italic"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Top 5 Markers per VZ Subcluster")

# 4. SAVE AS PDF
ggplot2::ggsave(
  file.path(plot_path, paste0("XenAld_VZ_Subcluster_Top5_Heatmap.pdf")), 
  p4, width = 14, height = 12, device = grDevices::cairo_pdf, useDingbats = FALSE
)

# 5. Save the Heatmap
Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("XenAld_VZ_Subcluster_Top5_Heatmap.tif")),
  width = 14,
  height = 12,
  units = "in",
  res = 600
)
print(p4)
grDevices::dev.off()

check_mem("HEATMAP COMPLETE")

plan("sequential")

check_mem("PIPELINE COMPLETE - FILE SAVED")
