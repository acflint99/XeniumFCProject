# Clear the environment
rm(list = ls())

# 1. INITIALIZATION & ENVIRONMENT
source("renv/activate.R")
library(here)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)

check_mem <- function(step_label) {
  # gc() triggers garbage collection and returns a memory report
  m <- gc(full = TRUE)
  # sum(m[,2]) gives the memory in MB currently used by R
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

# 1. PARALLELIZATION (Reducing to 20 workers to stay safe on RAM)
# Even with 440GB, 40 workers * 50GB object = Disaster. 
# 20 workers is the "Sweet Spot" for memory overhead.
plan("multicore", workers = 20) 
options(future.globals.maxSize = 100 * 1024^3) # 100GB limit

merged_path <- here("outputs", "XenAld_VZ_RDS", "Xenium_VZ_Integrated_31226.rds")
obj <- readRDS(merged_path)

# This merges the 15 separate sample layers into one unified matrix
obj <- JoinLayers(obj)

# 2. SET THE IDENTITY
# Ensure we are looking at the resolution you liked best (e.g., 0.3)
Idents(obj) <- "Xenium_snn_res.0.5"

check_mem("STARTING FINDALLMARKERS")

# 3. RUN FINDALLMARKERS (The Optimized Way)
message(Sys.time(), ": Starting Marker Identification...")

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

check_mem("MARKERS COMPLETE")

# # 4. FILTER & SAVE TOP 10
# top10_markers <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# write.csv(top10_markers,
#           here("outputs", "XenAld_VZ_Tables", "VZ_Integrated_top10_Markers_Res0.5.csv"),
#           row.names = FALSE)
# 
# message("Marker analysis complete! Results saved to CSV.")
# 
# markers <- c("SLC17A7", "SLC17A6", "PVALB", "GAD1", "SST")
# 
# DotPlot(obj, 
#         features = markers, 
#         cols = c("lightgrey", "red"), # Blue to Red is also a popular choice
#         dot.scale = 6, 
#         cluster.idents = FALSE) + # Keep them in numeric order 0, 1, 2...
#   RotatedAxis() + # Tilts the gene names for readability
#   theme(
#     axis.text.x = element_text(size = 8, face = "italic"),
#     axis.text.y = element_text(size = 10, face = "bold")
#   ) +
#   ggtitle("Top Markers Across VZ Subpopulations")


# 1. Create a named vector for the mapping
# The names (0, 1, 2...) must match your current cluster IDs exactly
new_labels <- c(
  "0" = "Early PCs",
  "1" = "GABA Progenitors",
  "2" = "Migrating PCs",
  "3" = "BG",
  "4" = "Astrocytes",
  "5" = "Differentiated PCs",
  "6" = "Prolif RG",
  "7" = "RG Progenitors",
  "8" = "eCN",
  "9" = "MLI",
  "10" = "OPCs",
  "11" = "Golgi Cells",
  "12" = "VZP",
  "13" = "Mature PCs",
  "14" = "GABA Interneurons",
  "15" = "GCP"
)

# 2. Update the Idents (active labels) of the object
# Assuming 'xe_obj' is your Xenium Seurat object
Idents(obj) <- "Xenium_snn_res.0.5"

obj <- RenameIdents(obj, new_labels)

# 3. Save these biological names into the metadata for persistent use
obj$VZ_subcluster <- Idents(obj)

# 4. Verify the change
table(obj$VZ_subcluster)

p <- DimPlot(obj, 
             reduction = "umap_harmony", 
             label = TRUE, 
             label.size = 4, 
             repel = TRUE) + 
  ggtitle("Xenium Merged VZ Subclusters") +
  theme(legend.text = element_text(size = 8))

ggsave(here("outputs", "XenAld_VZ_Plots", "XenAld_VZ_Subcluster_UMAP.png"), 
       p, width = 12, height = 9, dpi = 300)

seurat_levels <- c("VZP",
                   "Migrating PCs",
                   "Early PCs",
                   "Differentiated PCs",
                   "Mature PCs",
                   "GABA Progenitors",
                   "Golgi Cells",
                   "GABA Interneurons",
                   "eCN",
                   "MLI",
                   "RG Progenitors",
                   "Prolif RG",
                   "BG",
                   "Astrocytes",
                   "GCP",
                   "OPCs")

# 2. Relevel the active identity
levels(obj) <- seurat_levels

# Reverse the current levels of the active identity
levels(obj) <- rev(levels(obj))

markers <- c("ASCL1",
             "BCL11A",
             "KITLG",
             "EBF1",
             "FOXP2",
             "EBF2",
             "CALB1",
             "NEUROG2",
             "THBS1",
             "SP9",
             "SLC17A6",
             "SOX14",
             "SOX2", #could not find an ideal marker for RG Progenitors
             "TOP2A",
             "TNC",
             "AQP4",
             "EOMES",
             "OLIG1")

p1 <- DotPlot(obj,
        features = markers,
        cols = c("lightgrey", "red"), # Blue to Red is also a popular choice
        dot.scale = 6,
        cluster.idents = FALSE) + # Keep them in numeric order 0, 1, 2...
  RotatedAxis() + # Tilts the gene names for readability
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged VZ Subcluster Markers")

ggsave(here("outputs", "XenAld_VZ_Plots", "XenAld_VZ_SubclusterMarker_DotPlot.png"), 
       p1, width = 12, height = 8, dpi = 300)

# 9. SAVE FINAL RESULT
output_path <- here("outputs", "XenAld_VZ_RDS", "Xenium_VZ_Integrated_newSubclusters.rds")
saveRDS(obj, output_path, compress = FALSE)
check_mem("PIPELINE COMPLETE - FILE SAVED")
