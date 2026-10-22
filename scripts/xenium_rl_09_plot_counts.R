#!/usr/bin/env Rscript

# Plot the intentional focused RL developmental cohort.
rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_rl_09_plot_counts.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

input_path <- here("outputs", "xenium", "rl", "06_subclusters", "rds", "Xenium_RL_subclusters.rds")
plot_dir <- here("outputs", "xenium", "rl", "09_cluster_counts", "plots")

# 1. Define a named list of your cluster sets
cluster_sets <- list(
  "Granule_Lineage" = c("Prolif GCPs", "Maturing GCPs", "Differentiating GCs", "Migrating GCs", "Mature GCs"),
  "UBC_Lineage" = c("RL SVZ", "Immature UBCs", "Mature UBCs")
)

# Intentional one-sample-per-PCW cohort for the developmental comparison.
# This is not the complete 34-sample manifest by design.
target_sample_pcw <- c(
  "GZFB_12_X_G_3" = 10,
  "GZFB5_X_G" = 11,
  "FB328_1_X_G" = 12,
  "GZFB_1_X_G" = 12.5,
  "FB330_1_X_G" = 13,
  "FB78_X_G" = 14,
  "GZFB4_X_G" = 16,
  "FB124_X_G" = 17,
  "GZFB_22_X_G_4" = 19
)
target_samples <- names(target_sample_pcw)
target_pcw_levels <- unname(target_sample_pcw)

plot_stems <- paste0("Xenium_RL_", names(cluster_sets), "_ClusterCountPlot")
expected_outputs <- c(
  file.path(plot_dir, paste0(plot_stems, ".tif")),
  file.path(plot_dir, paste0(plot_stems, ".pdf"))
)

if (dry_run) {
  compact_dry_run(
    paste0("RL percent-abundance plots [", length(target_samples), "-sample focused cohort]"),
    inputs = input_path,
    outputs = expected_outputs
  )
  quit(save = "no", status = 0L)
}

if (!file.exists(input_path)) stop("RL subcluster input not found: ", input_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL percent-abundance plots:\n- ",
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
required_metadata <- c("PCW", "RL_subcluster", "orig.ident")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("RL subcluster input lacks metadata: ", paste(missing_metadata, collapse = ", "))
missing_target_samples <- setdiff(target_samples, unique(as.character(obj$orig.ident)))
if (length(missing_target_samples)) {
  stop("RL count cohort is missing sample(s): ", paste(missing_target_samples, collapse = ", "))
}
target_metadata <- unique(data.frame(
  sample_id = as.character(obj$orig.ident[obj$orig.ident %in% target_samples]),
  PCW = as.character(obj$PCW[obj$orig.ident %in% target_samples]),
  stringsAsFactors = FALSE
))
pcw_rows_per_sample <- table(target_metadata$sample_id)
invalid_pcw_samples <- target_samples[pcw_rows_per_sample[target_samples] != 1L]
if (length(invalid_pcw_samples)) {
  stop(
    "Each RL count-cohort sample must have exactly one PCW value; failures: ",
    paste(invalid_pcw_samples, collapse = ", ")
  )
}
observed_pcw <- suppressWarnings(as.numeric(gsub(
  "[^0-9.]", "", target_metadata$PCW[match(target_samples, target_metadata$sample_id)]
)))
pcw_mismatches <- target_samples[
  !is.finite(observed_pcw) | observed_pcw != target_sample_pcw[target_samples]
]
if (length(pcw_mismatches)) {
  stop(
    "RL count-cohort PCW metadata do not match the approved sample-to-PCW map: ",
    paste(pcw_mismatches, collapse = ", ")
  )
}

# 2. Iterate through the sets
for (set_name in names(cluster_sets)) {
  
  # Get the specific clusters for this iteration
  current_clusters <- cluster_sets[[set_name]]
  
  # Normalize within the plotted lineage in each sample (one sample per PCW).
  plot_df <- obj@meta.data %>%
    select(PCW, RL_subcluster, orig.ident) %>%
    filter(orig.ident %in% target_samples) %>%
    filter(RL_subcluster %in% current_clusters) %>%
    mutate(PCW_num = as.numeric(gsub("PCW", "", PCW))) %>%
    group_by(PCW_num, orig.ident, RL_subcluster) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(PCW_num, orig.ident) %>%
    mutate(percent_abundance = 100 * cell_count / sum(cell_count)) %>%
    ungroup() 
  
  # An absent lineage has no defined percentage; do not silently omit outputs.
  missing_lineage_samples <- setdiff(target_samples, as.character(plot_df$orig.ident))
  if (length(missing_lineage_samples)) {
    stop("Cannot calculate percent abundance for ", set_name,
         "; no lineage cells in sample(s): ", paste(missing_lineage_samples, collapse = ", "))
  }
  
  # Set Factor Levels
  plot_df$PCW_num <- factor(plot_df$PCW_num, levels = target_pcw_levels)
  plot_df$RL_subcluster <- factor(plot_df$RL_subcluster, 
                                  levels = intersect(rl_subcluster_order, current_clusters))
  
  # Plot
  p <- ggplot(plot_df, aes(x = PCW_num, y = percent_abundance, fill = RL_subcluster)) +
    geom_bar(stat = "identity", color = "black", width = 0.8, linewidth = 0.2) +
    scale_fill_manual(values = rl_palette) + 
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste("Cluster Set:", set_name),
      x = "Age (PCW)", 
      y = "Abundance (% of lineage cells)", 
      fill = "Cell Type",
      subtitle = paste(length(target_samples), "selected samples; 100% = plotted lineage within each sample")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))
  
  # Save with a dynamic filename based on the set_name
  file_name <- paste0("Xenium_RL_", set_name, "_ClusterCountPlot.tif")
  
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
