# ==============================================================================
# SCRIPT 4: SPATIAL VARIABLE GENES (SVG) ACROSS AN AXIS
# ==============================================================================
# Clear the environment
rm(list = ls())

options(bitmapType = "cairo")

library(dplyr)
library(ggplot2)
library(stringr)
library(here)

output_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Giotto_RDS/"
plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_GiottoBroad_Plots") 
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. Define Spatial Metadata (Using the 5-Slide Example)
meta_df <- data.frame(
  Sample = c("GZFB_12_X_G_5", "GZFB_12_X_G_4", "GZFB_12_X_G_3", "GZFB_12_X_G_2", "GZFB_12_X_G_1"),
  Axis_Pos = c(1, 2, 3, 4, 5),
  Label = c("Medial", "Mid-Medial", "Central", "Mid-Lateral", "Lateral")
)

# 2. Load and Aggregate the SVG CSV files
svg_files <- list.files(output_dir, pattern = "_spatial_genes\\.csv$", full.names = TRUE)
svg_files <- svg_files[basename(svg_files) %in% paste0(meta_df$Sample, "_spatial_genes.csv")]

all_svgs <- bind_rows(lapply(svg_files, function(f) {
  s_name <- str_remove(basename(f), "_spatial_genes\\.csv")
  df <- read.csv(f)
  df$Sample <- s_name
  return(df)
})) %>% left_join(meta_df, by = "Sample")


# ==============================================================================
# ANALYSIS A: FINDING "LINEAR DECAY" (COMPLETE GENES ONLY)
# ==============================================================================

# 1. Filter for completeness: Gene must be present in all 5 Axis Positions
complete_genes_df <- all_svgs %>%
  group_by(feats) %>%
  filter(n() == 5) %>% # Only keep genes appearing exactly 5 times (one per slide)
  ungroup()

# 2. Calculate linear stats on the complete dataset
decay_stats <- complete_genes_df %>%
  mutate(log10_pval = -log10(adj.p.value + 1e-300)) %>%
  group_by(feats) %>%
  summarize(
    slope = coef(lm(log10_pval ~ Axis_Pos))[2],
    correlation = cor(Axis_Pos, log10_pval, method = "pearson"),
    max_sig = max(log10_pval),
    # Ensure it starts significant (Medial) and ends less significant (Lateral)
    medial_sig = log10_pval[Axis_Pos == 1],
    lateral_sig = log10_pval[Axis_Pos == 5],
    .groups = 'drop'
  )

# 3. Filter for Linear Decay
# We want: 
# - Negative correlation (dropping significance)
# - Medial p-value < 0.05 (significant start)
# - High linearity (correlation < -0.9 is a very strict, clean linear drop)
linear_decay_genes <- decay_stats %>%
  filter(correlation < -0.8, medial_sig > -log10(0.05)) %>%
  arrange(correlation)

print(paste("Genes present in all 5 samples:", length(unique(complete_genes_df$feats))))
print(paste("Of those, genes with linear decay:", nrow(linear_decay_genes)))


# ==============================================================================
# ANALYSIS B: PLOTTING THE LINEAR DECAY
# ==============================================================================

# Select the top 6 genes with the strongest linear correlation
top_linear_genes <- linear_decay_genes %>% head(9) %>% pull(feats)

plot_svg_data <- complete_genes_df %>% 
  filter(feats %in% top_linear_genes) %>%
  mutate(log10_pval = -log10(adj.p.value + 1e-300))

p_svg <- ggplot(plot_svg_data, aes(x = Axis_Pos, y = log10_pval, color = feats)) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "gray") + # Add trend line
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
  scale_x_continuous(breaks = meta_df$Axis_Pos, labels = meta_df$Label) +
  facet_wrap(~feats, scales = "free_y") +
  labs(
    title = "Linear Decay of Spatial Organization",
    subtitle = "Dotted gray line indicates linear fit; Red dashed line is p = 0.05",
    x = "Anatomical Position", 
    y = "Spatial Significance (-log10 P-Value)"
  ) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Save the plot as TIFF
ggsave(file.path(plot_dir, "Global_LinearSVG_LinePlot.tif"), 
       plot = p_svg, device = "tiff", width = 12, height = 12, dpi = 600, compression = "lzw")