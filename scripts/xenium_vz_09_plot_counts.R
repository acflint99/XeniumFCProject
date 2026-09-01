#!/usr/bin/env Rscript

# Plot the intentional focused VZ developmental cohort.
rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_09_plot_counts.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

input_path <- here("outputs", "xenium", "vz", "06_subclusters", "rds", "Xenium_VZ_subclusters_Res1.5.rds")
plot_dir <- here("outputs", "xenium", "vz", "09_cluster_counts", "plots")

# 1. Define a named list of your cluster sets
cluster_sets <- list(
  "GABA_Lineage" = c("VZPs", "GABA Progenitors", "Golgi Cells", "MLIs", "iCN"),
  "Purkinje_Lineage" = c("VZPs", "Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs"),
  "Glial_Lineage"  = c("VZPs", "RG Progenitors", "BG", "Astrocytes/Ependyma")     
)

# Intentional focused cohort for the developmental lineage count comparison.
# This is not the complete 34-sample manifest by design.
target_samples <- c("FB328_1_X_G", "GZFB_12_X_G_3", "GZFB5_X_G",
                    "GZFB_1_X_G", "FB330_1_X_G", "FB78_X_G",
                    "GZFB4_X_G", "FB124_X_G") 

plot_stems <- paste0("XenAld_VZ_", names(cluster_sets), "_ClusterCountPlot")
expected_outputs <- c(
  file.path(plot_dir, paste0(plot_stems, ".tif")),
  file.path(plot_dir, paste0(plot_stems, ".pdf"))
)

if (dry_run) {
  compact_dry_run(
    paste0("VZ count plots [", length(target_samples), "-sample focused cohort]"),
    inputs = input_path,
    outputs = expected_outputs
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("VZ subcluster input not found: ", input_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ count plots:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
source(here("scripts", "color_palette.R"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(input_path)
required_metadata <- c("PCW", "VZ_subcluster", "orig.ident")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("VZ subcluster input lacks metadata: ", paste(missing_metadata, collapse = ", "))
missing_target_samples <- setdiff(target_samples, unique(as.character(obj$orig.ident)))
if (length(missing_target_samples)) {
  stop("VZ count cohort is missing sample(s): ", paste(missing_target_samples, collapse = ", "))
}

# 2. Iterate through the sets
for (set_name in names(cluster_sets)) {
  
  # Get the specific clusters for this iteration
  current_clusters <- cluster_sets[[set_name]]
  
  # Filter and Aggregate
  plot_df <- obj@meta.data %>%
    select(PCW, VZ_subcluster, orig.ident) %>%
    filter(orig.ident %in% target_samples) %>%
    filter(VZ_subcluster %in% current_clusters) %>%
    mutate(PCW_num = as.numeric(gsub("PCW", "", PCW))) %>%
    group_by(PCW_num, orig.ident, VZ_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop") 
  
  # Skip if the dataframe is empty (e.g., clusters not found in target samples)
  if (nrow(plot_df) == 0) next
  
  # Set Factor Levels
  plot_df$PCW_num <- factor(plot_df$PCW_num, levels = sort(unique(plot_df$PCW_num)))
  plot_df$VZ_subcluster <- factor(plot_df$VZ_subcluster, 
                                  levels = intersect(vz_subcluster_order, current_clusters))
  
  # Plot
  p <- ggplot(plot_df, aes(x = PCW_num, y = cell_count, fill = VZ_subcluster)) +
    geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
    scale_fill_manual(values = vz_palette) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste("Cluster Set:", set_name),
      x = "Age (PCW)", 
      y = "Number of Cells", 
      fill = "Cell Type",
      subtitle = paste("Included samples:", length(target_samples), "samples selected")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))
  
  # Save with a dynamic filename based on the set_name
  file_name <- paste0("XenAld_VZ_", set_name, "_ClusterCountPlot.tif")
  
  Cairo::CairoTIFF(
    filename = file.path(plot_dir, file_name),
    width = 10,
    height = 6,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()
  ggplot2::ggsave(
    file.path(plot_dir, sub("\\.tif$", ".pdf", file_name)),
    p, device = grDevices::cairo_pdf, width = 10, height = 6
  )
  
  message(paste("Successfully saved plot for:", set_name))
}
