
## NOT OPTIMIZED




# ==============================================================================
# SCRIPT 5: GENE EXPRESSION GRADIENT ACROSS AN AXIS
# ==============================================================================
# Clear the environment
rm(list = ls())
options(bitmapType = "cairo")

library(Seurat)
library(dplyr)
library(ggplot2)

# Point this to where your SEURAT objects live (from your earlier Giotto array script)
data_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS/"

# 1. Define Metadata
meta_df <- data.frame(
  Sample = c("GZFB_12_X_G_5", "GZFB_12_X_G_4", "GZFB_12_X_G_3", "GZFB_12_X_G_2", "GZFB_12_X_G_1"),
  Axis_Pos = c(1, 2, 3, 4, 5),
  Label = c("Medial", "Mid-Medial", "Central", "Mid-Lateral", "Lateral")
)

# 2. Define Your Target Genes
# (e.g., Take the genes from your SVG "lost_spatial_genes" list to see IF they also turned off)
target_genes <- c("GENE_A", "GENE_B", "GENE_C") # Replace with your actual genes of interest

# 3. Extract Expression Safely (One file at a time to prevent OOM crashes)
expr_list <- list()

for (i in 1:nrow(meta_df)) {
  s_name <- meta_df$Sample[i]
  message("Processing Expression for: ", s_name)
  
  # Find the Seurat file
  file_path <- list.files(data_dir, pattern = paste0("^", s_name, "_Ald_VZ_RL_QC_Subclusters\\.rds$"), full.names = TRUE)
  
  if (length(file_path) == 1) {
    # Load object
    seurat_obj <- readRDS(file_path)
    
    # Safely pull the normalized expression data for just our target genes
    # Note: If a gene isn't in a specific slice, FetchData might throw a warning, 
    # so we intersect with available genes first.
    available_genes <- intersect(target_genes, rownames(seurat_obj))
    expr_data <- FetchData(seurat_obj, vars = available_genes, layer = "data")
    
    # Calculate Mean Expression across all cells in this slice
    mean_expr <- colMeans(expr_data)
    
    # Store in a clean dataframe
    temp_df <- data.frame(
      Sample = s_name,
      Gene = names(mean_expr),
      Mean_Expression = as.numeric(mean_expr),
      Axis_Pos = meta_df$Axis_Pos[i],
      Label = meta_df$Label[i]
    )
    
    expr_list[[s_name]] <- temp_df
    
    # CRITICAL: Delete object and trigger garbage collection to free RAM
    rm(seurat_obj, expr_data); gc()
  } else {
    warning("Could not find (or found multiple) files for: ", s_name)
  }
}

# Combine all the summaries into one table
expression_df <- bind_rows(expr_list)

# ==============================================================================
# 4. PLOTTING THE EXPRESSION GRADIENT
# ==============================================================================

ggplot(expression_df, aes(x = Axis_Pos, y = Mean_Expression, color = Gene)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = meta_df$Axis_Pos, labels = meta_df$Label) +
  facet_wrap(~Gene, scales = "free_y") +
  labs(
    title = "Gene Expression Levels (Medial -> Lateral)",
    subtitle = "Average normalized expression across all cells per slice",
    x = "Anatomical Position", 
    y = "Mean Normalized Expression"
  ) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))