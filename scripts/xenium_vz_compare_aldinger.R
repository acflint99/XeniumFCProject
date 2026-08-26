# Clear the environment
rm(list = ls())

# Load required libraries
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(here)

xen_path <- here("outputs", "xenium", "vz", "06_subclusters", "rds", "Xenium_VZ_subclusters_Res1.5.rds")
xenium_merged <- readRDS(xen_path)

sn_path <- here("outputs", "references", "aldinger", "vz", "rds", "Aldinger_VZ_4126.rds")
sn_obj <- readRDS(sn_path)

plot_dir <- here("outputs", "references", "cross_study", "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


Idents(sn_obj) <- "RNA_snn_res.0.8"

# --- 1. Calculate Average Expression ---
# This creates a matrix where columns are clusters and rows are genes
# We use the 'data' slot (normalized counts) for correlation

# Calculate for Xenium (Query)
message("Calculating average expression for Xenium subclusters...")
avg_xenium <- AverageExpression(
  xenium_merged, 
  group.by = "ident", # Ensure your subclusters are set as the active ident
  slot = "data"
)[[1]]

# Calculate for snRNA-seq (Reference)
message("Calculating average expression for snRNA-seq clusters...")
avg_sn <- AverageExpression(
  sn_obj, 
  group.by = "ident", # Replace with your snRNA-seq metadata column name
  slot = "data"
)[[1]]

# 1. Find the intersection of genes again to be 100% sure
common_genes <- intersect(rownames(avg_xenium), rownames(avg_sn))

# 2. Subset both matrices to these common genes AND the same order
avg_xenium_sub <- avg_xenium[common_genes, ]
avg_sn_sub     <- avg_sn[common_genes, ]

# 3. Verify dimensions match (Number of rows should be identical)
message("Xenium genes: ", nrow(avg_xenium_sub))
message("snRNA-seq genes: ", nrow(avg_sn_sub))

# 4. Compute Correlation safely
message("Computing Spearman correlation matrix...")
cor_matrix <- cor(as.matrix(avg_xenium_sub), as.matrix(avg_sn_sub), method = "spearman")

# --- 2. Compute Correlation ---
# We use Spearman correlation as it is more robust to the different 
# dynamic ranges of Xenium vs. snRNA-seq (targeted vs. whole transcriptome)
message("Computing Spearman correlation matrix...")
cor_matrix <- cor(as.matrix(avg_xenium_sub), as.matrix(avg_sn_sub), method = "spearman")

# --- 3. Visualization ---
# Define a color palette (White to Blue)
my_colors <- colorRampPalette(c("red", "#DEEBF7", "#084594"))(100)

# Generate Heatmap
# If your subclusters have a one-to-one match, you'll see a strong diagonal
p1 <- pheatmap(
  mat = cor_matrix,
  color = my_colors,
  main = "Spearman Correlation: Xenium Subclusters vs. Aldinger",
  display_numbers = TRUE,       # Show the correlation coefficient in the box
  number_color = "black",
  fontsize_number = 8,
  cluster_rows = TRUE,          # Group similar Xenium subclusters
  cluster_cols = TRUE,          # Group similar snRNA-seq types
  angle_col = 45,               # Tilt column labels for readability
  border_color = "white"
)

Cairo::CairoTIFF(
  filename = file.path(plot_dir, paste0("XenAld_VZ_CorrPlot.tif")),
  width = 10,
  height = 10,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(file.path(plot_dir, "XenAld_VZ_CorrPlot.pdf"), p1,
                device = grDevices::cairo_pdf, width = 10, height = 10)

# --- 4. Optional: Save the matrix ---
# write.csv(cor_matrix, "Xenium_snRNAseq_correlation_matrix.csv")

# Find markers for every cluster
all_xenium_markers <- FindAllMarkers(
  xenium_merged, 
  only.pos = TRUE,          # Only look for upregulated genes
  min.pct = 0.25,           # Gene must be in at least 25% of cells in the cluster
  logfc.threshold = 0.25,   # Minimum log-fold change
  test.use = "wilcox"       # Standard statistical test
)
c
# View the top 5 markers for each cluster
top5_markers <- all_xenium_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Convert the top5_markers dataframe into a list of gene vectors
# Split genes by the 'cluster' column
xenium_signature_list <- split(top5_markers$gene, top5_markers$cluster)

# Get the names of the clusters
cluster_names <- names(xenium_signature_list)

# Loop through each cluster signature
for (clust in cluster_names) {
  
  # 1. Extract the specific gene vector for this cluster
  current_genes <- list(xenium_signature_list[[clust]])
  
  # 2. Add Module Score to the snRNA-seq object
  # We use a fixed name 'TempScore' so we can overwrite it in each loop
  sn_obj <- AddModuleScore(
    object = sn_obj,
    features = current_genes,
    name = "TempScore"
  )
  
  # 3. Generate the FeaturePlot
  # Seurat always appends '1' to the name you give it
  p <- FeaturePlot(sn_obj, reduction = "umap_uncorrected", features = "TempScore1") +
    ggplot2::ggtitle(paste0("Xenium Cluster ", clust, " Signature in snRNA-seq")) +
    ggplot2::scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  # 4. Save or Print the plot
  print(p) 
  Cairo::CairoTIFF(
    filename = file.path(plot_dir, paste0("XenAld_VZ_", clust, "_Sign_UMAP.tif")),
    width = 10,
    height = 10,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
}

# --- 5. Signature Heatmap (Module Score Comparison) ---
message("Generating Signature Heatmap...")

# 1. Add ALL module scores to the snRNA-seq object at once
# We use the list created earlier: xenium_signature_list
sn_obj <- AddModuleScore(
  object = sn_obj,
  features = xenium_signature_list,
  name = "XeniumSig_"
)

# 2. Map the generic names (XeniumSig_1, XeniumSig_2...) back to your Cluster Names
# Seurat appends the number based on the order of the list
score_cols <- paste0("XeniumSig_", 1:length(xenium_signature_list))
# Create a mapping vector
meta_names <- colnames(sn_obj@meta.data)
colnames(sn_obj@meta.data)[meta_names %in% score_cols] <- names(xenium_signature_list)

# 3. Calculate the Average Module Score per snRNA-seq cluster
# This uses your active identity (RNA_snn_res.0.8)
score_matrix <- sn_obj@meta.data %>%
  # Use the column that contains your snRNA-seq cluster IDs
  group_by(RNA_snn_res.0.8) %>% 
  summarise(across(all_of(names(xenium_signature_list)), mean, na.rm = TRUE)) %>%
  column_to_rownames("RNA_snn_res.0.8") %>%
  as.matrix()

# Transpose so Xenium Signatures are Rows and snRNA-seq Clusters are Columns
score_matrix_t <- t(score_matrix)

# 4. Visualization
# We use a divergent color scale because Module Scores are centered around 0
sig_colors <- colorRampPalette(c("white", "ghostwhite", "darkred"))(100)

p_sig <- pheatmap(
  mat = score_matrix_t,
  color = sig_colors,
  main = "Xenium Cluster Signatures mapped to Aldinger Clusters",
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  angle_col = 45,
  border_color = "white",
  silent = FALSE
)

# Save the plot
Cairo::CairoTIFF(
  filename = file.path(plot_dir, "XenAld_VZ_Signature_Heatmap.tif"),
  width = 12,
  height = 10,
  units = "in",
  res = 600
)
print(p_sig)
grDevices::dev.off()
ggplot2::ggsave(file.path(plot_dir, "XenAld_VZ_Signature_Heatmap.pdf"), p_sig,
                device = grDevices::cairo_pdf, width = 12, height = 10)
