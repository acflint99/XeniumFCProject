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

# Load your new palette and order
source(here("scripts", "color_palette.R"))

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
plan("multisession", workers = 20) 
options(future.globals.maxSize = 100 * 1024^3) # 100GB limit

merged_path <- here("outputs", "XenAld_RL_Integrated_RDS", "Xenium_RL_Integrated_31326.rds")
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
# 
# check_mem("MARKERS COMPLETE")
# 
# # 4. FILTER & SAVE TOP 10
# top10_markers <- all_markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
# 
# write.csv(top10_markers,
#           here("outputs", "XenAld_RL_Integrated_Tables", "RL_Integrated_top10_Markers_Res0.5.csv"),
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
#   ggtitle("Top Markers Across RL Subpopulations")


# 1. Create a named vector for the mapping
# The names (0, 1, 2...) must match your current cluster IDs exactly
new_labels <- c(
  "0" = "Cycling GCPs/EGL",
  "1" = "Differentiating GCPs",
  "2" = "Migrating GCs",
  "3" = "RLSVZ",
  "4" = "Mature GCs/IGL",
  "5" = "Mature UBCs",
  "6" = "Cycling Cells",
  "7" = "Purkinje Cells",
  "8" = "UBC Progenitors",
  "9" = "Transitioning UBCs",
  "10" = "RL Transition",
  "11" = "RG Progenitors",
  "12" = "Mature GCs/IGL",
  "13" = "RLVZ",
  "14" = "Astrocytes"
)

obj <- RenameIdents(obj, new_labels)

# 3. APPLY STANDARDIZED ORDERING
# We use 'subcluster_order' from your color_palette.R script
obj$RL_subcluster <- factor(Idents(obj), levels = rl_subcluster_order)
Idents(obj) <- "RL_subcluster"

# 4. UMAP WITH CONSISTENT COLORS
p <- DimPlot(obj, 
             reduction = "umap_harmony", 
             label = TRUE,
             label.size = 4,
             repel = TRUE,
             cols = rl_palette) + # Uses your defined palette
  ggtitle("Xenium Merged RL Subclusters") +
  theme(legend.text = element_text(size = 8))

ggsave(here("outputs", "XenAld_RL_Integrated_Plots", "XenAld_RL_Subcluster_UMAP.png"), 
       p, width = 12, height = 9, dpi = 300)

# 5. DOTPLOT WITH REVERSED ORDER
# Create a temporary reversed factor for the Y-axis
obj$RL_subcluster_rev <- factor(obj$RL_subcluster, levels = rev(rl_subcluster_order))
Idents(obj) <- "RL_subcluster_rev"

p1 <- DotPlot(obj,
              features = rl_markers,
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) + 
  RotatedAxis() + 
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged RL Subcluster Markers")

ggsave(here("outputs", "XenAld_RL_Integrated_Plots", "XenAld_RL_SubclusterMarker_DotPlot.png"), 
       p1, width = 12, height = 8, dpi = 300)

# 6. SAVE FINAL RESULT
# Revert Idents to the standard order before saving
Idents(obj) <- "RL_subcluster"

output_path <- here("outputs", "XenAld_RL_Integrated_RDS", "Xenium_RL_Integrated_newSubclusters.rds")
saveRDS(obj, output_path, compress = FALSE)

plan("sequential") # Shuts down the 20 background workers

check_mem("PIPELINE COMPLETE - FILE SAVED")
