# 1. INITIALIZATION & ENVIRONMENT
rm(list = ls())
options(bitmapType = "cairo")

source("renv/activate.R")
library(here)
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr) # For easier list handling

# Load your new palette and order
source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()

# Define paths
data_dir <- here("outputs", "xenium", "vz_rl", "01_combined_labels", "rds")
plot_out_dir <- here("outputs", "xenium", "vz_rl", "nonspatial", "05_plots", "cluster_counts")
if(!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE)

# 1. Define a named list of your cluster sets
cluster_sets <- list(
  "Granule_Lineage" = c(
    "Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"
  ),
  "UBC_Lineage" = c("RL VZ", "RL SVZ", "Immature UBCs", "Mature UBCs"),
  "eCN_Lineage" = "eCN",
  "GABA_Lineage" = c(
    "VZPs", "Cycling VZPs", "GABA Progenitors", "Golgi Cells", "MLIs", "iCN"
  ),
  "Purkinje_Lineage" = c(
    "VZPs", "Cycling VZPs", "Maturing PCs", "Early-born PCs",
    "Late-born PCs", "Patterning PCs"
  ),
  "Glial_Lineage" = c(
    "VZPs", "Cycling VZPs", "RG Progenitors", "Cycling Glia",
    "BG", "BG Progenitors", "Astrocytes/Ependyma", "Pre-OPCs", "OPCs"
  ),
  "Other_Lineage" = c(
    "Cycling Cells", "Meninges", "Microglia", "Vascular", "Immune"
  )
)

target_samples <- load_sample_manifest(config)$sample_id

# 2. Identify and Load Files
# Get all .rds files in the directory
all_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

# Filter files that match your target_samples list
files_to_process <- all_files[grepl(paste(target_samples, collapse="|"), basename(all_files))]

# 3. Extract Counts from each file
all_counts_list <- map(files_to_process, function(f) {
  message("Processing: ", basename(f))
  temp_obj <- readRDS(f)
  
  # Extract metadata
  df <- temp_obj@meta.data %>%
    select(PCW, comb_subcluster, orig.ident) %>%
    mutate(PCW_num = as.numeric(gsub("PCW", "", PCW))) %>%
    group_by(PCW_num, orig.ident, comb_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  rm(temp_obj); gc() # Clean up memory immediately
  return(df)
})

# Combine all small dataframes into one master plotting dataframe
master_plot_df <- bind_rows(all_counts_list)

# 4. Iterate through the Cluster Sets to Plot
for (set_name in names(cluster_sets)) {
  
  current_clusters <- cluster_sets[[set_name]]
  
  # Filter the master dataframe for the current lineage
  plot_df <- master_plot_df %>%
    filter(comb_subcluster %in% current_clusters) %>%
    # Pool counts across samples at each age before normalizing.
    group_by(PCW_num, comb_subcluster) %>%
    summarise(cell_count = sum(cell_count), .groups = "drop") %>%
    group_by(PCW_num) %>%
    mutate(percent_abundance = 100 * cell_count / sum(cell_count)) %>%
    ungroup()
  
  if (nrow(plot_df) == 0) {
    stop("Cannot calculate percent abundance; no lineage cells for: ", set_name)
  }
  
  # Factor leveling for consistent visualization
  plot_df$PCW_num <- factor(plot_df$PCW_num, levels = sort(unique(plot_df$PCW_num)))
  plot_df$comb_subcluster <- factor(plot_df$comb_subcluster, 
                                  levels = intersect(master_subcluster_order, current_clusters))
  
  # Plot
  p <- ggplot(plot_df, aes(x = PCW_num, y = percent_abundance, fill = comb_subcluster)) +
    geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
    scale_fill_manual(values = subcluster_palette) + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste("Cluster Set:", set_name),
      subtitle = "100% = plotted lineage cells pooled across samples at each age",
      x = "Age (PCW)", 
      y = "Abundance (% of lineage cells)", 
      fill = "Cell Type"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(), 
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )
  
  # Save
  file_name <- paste0("Xenium_CombVZ&RL_", set_name, "_ClusterCountPlot.tif")
  Cairo::CairoTIFF(
    filename = file.path(plot_out_dir, file_name),
    width = 10,
    height = 6,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
  ggplot2::ggsave(
    file.path(plot_out_dir, sub("\\.tif$", ".pdf", file_name)),
    p, device = grDevices::cairo_pdf, width = 10, height = 6
  )
  
  message("Saved plot: ", file_name)
}
