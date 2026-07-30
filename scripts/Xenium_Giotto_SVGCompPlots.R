# ==============================================================================
# SCRIPT 4: CONSENSUS SPATIAL VARIABLE GENES (2 BIOLOGICAL REPLICATES)
# ==============================================================================
# Clear the environment
rm(list = ls())
options(bitmapType = "cairo")

library(dplyr)
library(ggplot2)
library(stringr)

output_dir <- "~/R/Projects/XeniumFCProject/outputs/Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Giotto_RDS/"
plot_dir <- here("outputs", "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_GiottoBroad_Plots") 
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# ==============================================================================
# 1. DEFINE METADATA AND LOAD DATA
# ==============================================================================

# --- Dataset 1: The 5-Slide Sample ---
meta_df_5 <- data.frame(
  Sample = c("GZFB_12_X_G_5", "GZFB_12_X_G_4", "GZFB_12_X_G_3", "GZFB_12_X_G_2", "GZFB_12_X_G_1"),
  Axis_Pos = c(1, 2, 3, 4, 5),
  Label = c("Medial", "Mid-Med", "Central", "Mid-Lat", "Lateral"),
  Dataset = "Sample A (5-Slide)"
)

# --- Dataset 2: The 3-Slide Sample (REPLACE NAMES HERE) ---
meta_df_3 <- data.frame(
  Sample = c("GZFB_9_X_G_3", "GZFB_9_X_G_2", "GZFB_9_X_G_1"),
  # Notice we use 1, 3, 5 so they physically align with the 5-slide sample on the graph!
  Axis_Pos = c(1, 3, 5), 
  Label = c("Medial", "Intermediate", "Lateral"),
  Dataset = "Sample B (3-Slide)"
)

# Function to load and aggregate SVGs
load_svgs <- function(meta_df) {
  svg_files <- list.files(output_dir, pattern = "_spatial_genes\\.csv$", full.names = TRUE)
  svg_files <- svg_files[basename(svg_files) %in% paste0(meta_df$Sample, "_spatial_genes.csv")]
  
  bind_rows(lapply(svg_files, function(f) {
    s_name <- str_remove(basename(f), "_spatial_genes\\.csv")
    df <- read.csv(f)
    df$Sample <- s_name
    return(df)
  })) %>% left_join(meta_df, by = "Sample")
}

# Load both datasets
all_svgs_5 <- load_svgs(meta_df_5)
all_svgs_3 <- load_svgs(meta_df_3)

# Combine them into one giant master table for plotting later
all_svgs_master <- bind_rows(all_svgs_5, all_svgs_3)


# ==============================================================================
# 2. FINDING THE CONSENSUS "GRADIENT" GENES
# ==============================================================================

# --- Process Dataset 1 (5-Slide) ---
sig_svgs_5 <- all_svgs_5 %>% filter(adj.p.value < 0.05)
medial_genes_5 <- sig_svgs_5 %>% filter(Axis_Pos == 1) %>% pull(feats)
lateral_genes_5 <- sig_svgs_5 %>% filter(Axis_Pos == 5) %>% pull(feats)
lost_genes_5 <- setdiff(medial_genes_5, lateral_genes_5)

# --- Process Dataset 2 (3-Slide) ---
sig_svgs_3 <- all_svgs_3 %>% filter(adj.p.value < 0.05)
medial_genes_3 <- sig_svgs_3 %>% filter(Axis_Pos == 1) %>% pull(feats)
lateral_genes_3 <- sig_svgs_3 %>% filter(Axis_Pos == 5) %>% pull(feats) # Axis_Pos 5 is the Lateral slice here too!
lost_genes_3 <- setdiff(medial_genes_3, lateral_genes_3)

# --- The Magic Intersection ---
consensus_lost_genes <- intersect(lost_genes_5, lost_genes_3)

print(paste("Genes lost in 5-slide sample:", length(lost_genes_5)))
print(paste("Genes lost in 3-slide sample:", length(lost_genes_3)))
print(paste("CONSENSUS genes lost in BOTH samples:", length(consensus_lost_genes)))


# ==============================================================================
# 3. PLOTTING THE REPRODUCIBLE SPATIAL DECAY
# ==============================================================================

# Take the top 6 consensus genes (ranked by how significant they were in Sample A's Medial slice)
top_consensus_genes <- sig_svgs_5 %>% 
  filter(feats %in% consensus_lost_genes & Axis_Pos == 1) %>%
  arrange(adj.p.value) %>%
  head(9) %>%
  pull(feats)

# Filter the combined master table
plot_svg_data <- all_svgs_master %>% 
  filter(feats %in% top_consensus_genes) %>%
  mutate(log10_pval = -log10(adj.p.value + 1e-300))

# Plotting
p_spatial_decay <- ggplot(plot_svg_data, aes(x = Axis_Pos, y = log10_pval, color = Dataset, shape = Dataset)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + 
  # Set the x-axis to represent the shared anatomical space
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), 
                     labels = c("Medial", "Mid-Med", "Central", "Mid-Lat", "Lateral")) +
  facet_wrap(~feats, scales = "free_y") +
  labs(
    title = "Consensus Decay of Spatial Gene Organization",
    subtitle = "Comparing spatial structure loss Medial -> Lateral across two biological replicates",
    x = "Anatomical Position", 
    y = "Spatial Significance (-log10 P-Value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("Sample A (5-Slide)" = "steelblue", "Sample B (3-Slide)" = "darkorange"))

# 5. Save the plot
ggsave(file.path(plot_dir, "GlobalComp_SVG_LinePlot.png"), 
       plot = p_spatial_decay, width = 16, height = 12, dpi = 300)