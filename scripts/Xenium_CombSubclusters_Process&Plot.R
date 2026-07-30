# Clear the environment
rm(list = ls())

library(Seurat)
library(harmony)
library(ggplot2)
library(here)
library(future)

# Set the number of workers (cores). 
# Be careful not to exceed your RAM capacity, as each worker copies the object.
plan("sequential")

# Load your new palette and order
source(here("scripts", "color_palette.R"))

# Ensure both directories exist
plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_Plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# # Load the merged object created in your previous script
# merged_obj <- readRDS(here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS", "XenAld_VZ&RL_clean_merge_4-6-26.rds"))
# 
# # 1. Identify the 'prefix' from the full cell name
# # This regex looks for the start of the 16-character Xenium barcode (e.g., aaeb...) 
# # and keeps everything before it.
# full_names <- rownames(merged_obj@meta.data)
# sample_ids <- gsub("_[a-z]{8,}.*", "", full_names)
# 
# # 2. Assign and check levels
# merged_obj$orig.ident <- as.factor(sample_ids)
# 
# # 3. Verification
# n_samples <- length(unique(merged_obj$orig.ident))
# message("Total samples found: ", n_samples)
# print(table(merged_obj$orig.ident))
# 
# # 1. Standard Pre-processing
# merged_obj <- NormalizeData(merged_obj)
# merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000)
# merged_obj <- JoinLayers(merged_obj)
# merged_obj <- ScaleData(merged_obj)
# merged_obj <- RunPCA(merged_obj, npcs = 30)
# 
# # Run UMAP on PCA (PRE-HARMONY) - Name it "umap.unintegrated"
# merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
# 
# # 2. Run Harmony
# merged_obj <- RunHarmony(merged_obj, group.by.vars = "orig.ident", assay.type = "Xenium")
# 
# # 3. Run UMAP on Harmony (POST-HARMONY) - Name it "umap.harmony"
# merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
# 
# saveRDS(merged_obj, here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS", "XenAld_VZ&RL_clean_merge_processed_4-7-26.rds"))

# --- PLOTTING ---

merged_obj <- readRDS(here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Clean_RDS", "XenAld_VZ&RL_clean_merge_processed_4-7-26.rds"))

# # Now p1 and p2 will actually look different
# p1 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "orig.ident", raster = TRUE) + 
#   ggtitle("Pre-Harmony (PCA UMAP)")
# 
# p2 <- DimPlot(merged_obj, reduction = "umap.harmony", group.by = "orig.ident", raster = TRUE) + 
#   ggtitle("Post-Harmony (Harmony UMAP)")
# 
# ggsave(filename = file.path(plot_dir, "XenAld_Batch_Comp_UMAP.png"), p1 + p2, width = 16, height = 7)
# 
# # 2. Generate the Plot
# p_orig <- DimPlot(merged_obj, 
#                   reduction = "umap.harmony", 
#                   group.by = "cluster_weighted", # Use your original label column here
#                   label = TRUE, 
#                   label.size = 4,
#                   label.box = TRUE,       # Makes original labels easier to see
#                   raster = TRUE, 
#                   pt.size = 0.5, 
#                   alpha = 0.8,
#                   cols = cluster_colors) +
#   ggtitle("Xenium Aldinger Merged UMAP: Original Cluster Weighted Labels") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#     axis.text = element_blank(),
#     axis.ticks = element_blank()
#   )
# 
# # 3. Save the Plot
# ggsave(filename = file.path(plot_dir, "XenAld_Merged_OrigClusterWeighted_UMAP.png"), 
#        p_orig, width = 12, height = 9, dpi = 300)
# 
# # Filter markers to only those present in the Xenium assay
# existing_markers <- lapply(markers, function(x) intersect(x, rownames(merged_obj)))
# 
# # existing_markers <- c("OTX2", "EOMES", "RELN", "FOXP2", "PAX2", "TNC",
# #                        "OLIG1", "FOXC1", "CLDN5", "P2RY12")
# 
# Idents(merged_obj) <- factor(merged_obj$cluster_weighted, levels = rev(celltype_order))
# 
# p5 <- DotPlot(merged_obj, features = existing_markers, assay = "Xenium") + 
#   RotatedAxis() +
#   scale_color_gradient(low = "lightgrey", high = "red") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Xenium Merged Marker Expression (LogNorm)")
# 
# ggsave(filename = file.path(plot_dir, "XenAld_Merged_cluster_weighted_DotPlot_markers.png"), 
#     p5, width = 14, height = 6, dpi = 300)
# 
# 
# # 1. Define the Purkinje Cell clusters and your markers of interest
# pc_clusters <- c("Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs")
# 
# # Replace these with your actual gene names
# pc_specific_markers <- c("FOXP1", "ITPR1","COL5A1", 
#                          "VSTM2L", "NDNF", "EBF2",
#                          "CALB1","TRPC3", "NEFL", "ETV1",
#                          "PCDH10", "EBF1", "BCL11A", "RORB", "EN1")
# 
# # 2. Subset the object
# # We use cluster_weighted because that is where your labels are stored
# pc_subset <- subset(merged_obj, subset = VZ_subcluster %in% pc_clusters)
# 
# # 3. Clean up the factor levels
# # This ensures the DotPlot only shows these 4 clusters in a specific order
# pc_subset$VZ_subcluster <- factor(pc_subset$VZ_subcluster, levels = rev(pc_clusters))
# Idents(pc_subset) <- "VZ_subcluster"
# 
# # 4. Filter markers for availability in the Xenium panel
# # Just in case some markers in your list aren't in the specific Xenium assay
# existing_pc_markers <- intersect(pc_specific_markers, rownames(pc_subset))
# 
# # 5. Generate the DotPlot
# p_pc_dots <- DotPlot(
#   pc_subset, 
#   features = existing_pc_markers, 
#   assay = "Xenium",
#   cols = c("lightgrey", "red"), # Using blue for a distinct look, or keep "red"
#   dot.scale = 8
# ) + 
#   RotatedAxis() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   ) +
#   ggtitle("Purkinje Cell Lineage Marker Expression")
# 
# # 6. Save the plot
# ggsave(
#   filename = file.path(plot_dir, "XenAld_Purkinje_Specific_DotPlot.png"), 
#   plot = p_pc_dots, 
#   width = 8, 
#   height = 5, 
#   dpi = 300
# )
# 
# # ----------------------------------------------------------------              
# # 6. GENERATE TOP 5 MARKERS HEATMAP
# # ----------------------------------------------------------------
# # 1. Identify Markers (if not already in environment from earlier)
# # We use a lower max.cells.per.ident to speed up the heatmap calculation
# message("Finding top 5 markers for heatmap...")
# heatmap_markers <- FindAllMarkers(
#   merged_obj,
#   only.pos = TRUE,
#   min.pct = 0.25,
#   logfc.threshold = 0.25,
#   max.cells.per.ident = 500 # Faster for visualization purposes
# )
# 
# # 2. Extract top 5 per cluster
# top5_markers <- heatmap_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 5, order_by = avg_log2FC) %>%
#   pull(gene) %>%
#   unique()
# 
# top5_markers <- rev(top5_markers)
# 
# # 3. Scale the data for the specific markers only 
# # (Necessary for DoHeatmap to show relative expression)
# # Just run the plot on a temporary scaled version:
# temp_obj <- ScaleData(merged_obj, features = top5_markers, verbose = FALSE)
# 
# # ----------------------------------------------------------------
# # 3.5 RE-ORDER FACTORS (The Fix)
# # ----------------------------------------------------------------
# # Ensure the metadata column used for grouping is a factor with your specific order
# celltype_order <- c(
#   "RL", "UBC", "Granule",
#   "VZ", "Purkinje", "GABA", "Glia", "OPC", "Meninges",
#   "Endothelial", "Immune"
# )
# 
# # Apply the factor levels to the temporary object
# temp_obj$cluster_weighted <- factor(
#   temp_obj$cluster_weighted, 
#   levels = celltype_order
# )
# 
# # 4. Create Heatmap
# # We downsample the plot to 100 cells per group so the labels are readable
# # 1. Prepare the plot without the default lines
# p4 <- DoHeatmap(
#   subset(temp_obj, downsample = 100), 
#   features = top5_markers,
#   group.by = "cluster_weighted",
#   group.colors = cluster_colors,
#   size = 4,           
#   angle = 45,         
#   draw.lines = TRUE,
#   raster = FALSE
# ) + 
#   scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
#   guides(color = "none") +
#   theme(
#     axis.text.y = element_text(size = 6, face = "italic"),
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   ) +
#   ggtitle("Top 5 Markers per Cluster")
# 
# # 4. SAVE AS PDF
# ggsave(
#   file.path(plot_dir, paste0("XenAld_Merged_Cluster_Top5_Heatmap.pdf")), 
#   p4, width = 14, height = 12, device = "pdf", useDingbats = FALSE
# )
# 
# # 5. Save the Heatmap
# ggsave(
#   file.path(plot_dir, paste0("XenAld_Merged_Cluster_Top5_Heatmap.png")), 
#   p4, width = 14, height = 12, dpi = 300
# )


# ==============================================================================
#  EXPRESSION VALIDATION
# ==============================================================================

# 1. Set Ident to your broad cluster labels
Idents(merged_obj) <- "cluster_weighted"

# 2. Define the genes of interest
kit_genes <- c("APP", "TNFRSF21", "CADM3", "CADM4", "NECTIN3", "CNTN2", "L1CAM", "CXCL12","CXCR4", "EFNB2", "EPHA4", "GJA1", "NCAM1",
               "L1CAM", "NRXN2", "CLSTN1", "ADGRL1", "DAG1", "LRRTM1", "RELN", "VLDLR", "CD99", "NTF3", "NTRK2")

# 3. Create the Violin Plot
# We use stack = TRUE and flip = TRUE to make a nice comparison of the two genes
p_kit_vln <- VlnPlot(
  merged_obj, 
  features = kit_genes, 
  pt.size = 0,           # Set to 0 to avoid overplotting dots in Xenium data
  ncol = 3, 
  cols = cluster_colors, # Uses your pre-defined palette
  idents = c("Granule", "Purkinje", "Glia",  "GABA", "OPC", "UBC") # Focus on the key players from your heatmap
) & theme(
  plot.title = element_text(face = "bold", size = 14),
  axis.title.x = element_blank()
)

# 4. Save the plot
ggsave(
  filename = file.path(plot_dir, "XenAld_Signaling_Genes_dual_VlnPlot.png"), 
  plot = p_kit_vln, 
  width = 20, 
  height = 20, 
  dpi = 300
)
