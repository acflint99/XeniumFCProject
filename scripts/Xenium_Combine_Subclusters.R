#!/usr/bin/env Rscript

# Combine VZ and RL subcluster labels with the broad consensus labels for all
# manifest-defined Xenium samples. Label priority is VZ -> RL -> consensus.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/Xenium_Combine_Subclusters.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "Xenium_AldingerABT_VZ&RLsubclusters_QC_Res1.5_RDS")
rds_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_RDS")
plot_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Plots")
table_dir <- file.path(output_root, "Xenium_AldingerABT_combVZ&RLsubcluster_Res1.5_Tables")

input_names <- paste0(sample_ids, "_Ald_VZ_RL_QC_Subclusters.rds")
input_paths <- file.path(input_dir, input_names)
rds_paths <- file.path(rds_dir, input_names)
plot_suffixes <- c(
  "_combSubcluster_GlobalSpatial0.6.tif",
  "_combSubcluster_GlobalSpatial0.85.tif",
  "_postQC_GlobalSpatial_plot0.85.tif",
  "_postQC_GlobalSpatial_plot0.6.tif"
)
plot_paths <- lapply(sample_ids, function(sample_id) {
  file.path(plot_dir, paste0(sample_id, plot_suffixes))
})
table_paths <- lapply(sample_ids, function(sample_id) {
  file.path(
    table_dir,
    paste0(sample_id, c("_counts_consensus_label.csv", "_counts_VZ&RLsubclusters.csv"))
  )
})
expected_by_sample <- lapply(seq_along(sample_ids), function(i) {
  c(rds_paths[[i]], plot_paths[[i]], table_paths[[i]])
})

observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  cat("Expected combined-subcluster inputs:", length(input_paths), "\n")
  cat("Existing combined-subcluster inputs:", sum(file.exists(input_paths)), "\n")
  cat("Unexpected top-level RDS files:", length(unexpected_names), "\n")
  write.table(
    data.frame(
      sample_id = sample_ids,
      input = input_paths,
      input_exists = file.exists(input_paths),
      expected_outputs = lengths(expected_by_sample),
      existing_outputs = vapply(expected_by_sample, function(x) sum(file.exists(x)), integer(1))
    ),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  if (length(unexpected_names)) {
    cat("Unexpected files:\n- ", paste(unexpected_names, collapse = "\n- "), "\n")
  }
  quit(save = "no", status = 0L)
}

if (length(missing_names)) {
  stop("Cannot combine subclusters; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop(
    "Cannot combine subclusters; unexpected top-level RDS files:\n- ",
    paste(unexpected_names, collapse = "\n- ")
  )
}
expected_outputs <- unlist(expected_by_sample, use.names = FALSE)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing combined-subcluster outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Processing combined subclusters ", i, " of ", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])

  required_metadata <- c("VZ_subcluster", "RL_subcluster", "consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  }
  pcw_values <- unique(as.character(obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

  combined_labels <- dplyr::coalesce(
    as.character(obj$VZ_subcluster),
    as.character(obj$RL_subcluster),
    as.character(obj$consensus_label)
  )
  if (anyNA(combined_labels) || any(!nzchar(combined_labels))) {
    stop(sample_id, " has cells without a VZ, RL, or consensus label.")
  }
  unexpected_labels <- setdiff(unique(combined_labels), master_subcluster_order)
  if (length(unexpected_labels)) {
    stop(
      sample_id, " has combined labels absent from master_subcluster_order: ",
      paste(unexpected_labels, collapse = ", ")
    )
  }
  present_levels <- master_subcluster_order[master_subcluster_order %in% unique(combined_labels)]
  obj$comb_subcluster <- factor(combined_labels, levels = present_levels)
  Idents(obj) <- "comb_subcluster"
  plot_colors <- subcluster_palette[present_levels]
  if (anyNA(plot_colors)) {
    stop(sample_id, " has combined labels without colors in subcluster_palette.")
  }

  consensus_counts <- as.data.frame(table(obj$consensus_label), stringsAsFactors = FALSE)
  colnames(consensus_counts) <- c("Consensus_Label", "Cell_Count")
  combined_counts <- as.data.frame(table(obj$comb_subcluster), stringsAsFactors = FALSE)
  colnames(combined_counts) <- c("Combined_Subcluster", "Cell_Count")
  write.csv(consensus_counts, table_paths[[i]][[1]], row.names = FALSE)
  write.csv(combined_counts, table_paths[[i]][[2]], row.names = FALSE)

  p_comb_06 <- ImageDimPlot(
    obj, group.by = "comb_subcluster", cols = plot_colors, size = 0.6
  ) + ggtitle(paste(sample_id, "Combined Subclusters"))
  Cairo::CairoTIFF(
    filename = plot_paths[[i]][[1]], width = 14, height = 10,
    units = "in", res = 600, bg = "black"
  )
  print(p_comb_06)
  grDevices::dev.off()

  p_comb_085 <- ImageDimPlot(
    obj, group.by = "comb_subcluster", cols = plot_colors, size = 0.85
  ) + ggtitle(paste(sample_id, "Combined Subclusters"))
  Cairo::CairoTIFF(
    filename = plot_paths[[i]][[2]], width = 14, height = 10,
    units = "in", res = 600, bg = "black"
  )
  print(p_comb_085)
  grDevices::dev.off()

  consensus_levels <- levels(droplevels(factor(obj$consensus_label)))
  consensus_colors <- cluster_colors[consensus_levels]
  if (anyNA(consensus_colors)) {
    stop(sample_id, " has consensus labels without colors in cluster_colors.")
  }
  p_consensus_085 <- ImageDimPlot(
    obj, group.by = "consensus_label", size = 0.85, cols = consensus_colors
  ) + ggtitle(paste(sample_id, "-Post QC (Consensus)"))
  Cairo::CairoTIFF(
    filename = plot_paths[[i]][[3]], width = 12, height = 10,
    units = "in", res = 600, bg = "black"
  )
  print(p_consensus_085)
  grDevices::dev.off()

  p_consensus_06 <- ImageDimPlot(
    obj, group.by = "consensus_label", size = 0.6, cols = consensus_colors
  ) + ggtitle(paste(sample_id, "-Post QC (Consensus)"))
  Cairo::CairoTIFF(
    filename = plot_paths[[i]][[4]], width = 12, height = 10,
    units = "in", res = 600, bg = "black"
  )
  print(p_consensus_06)
  grDevices::dev.off()

  saveRDS(obj, rds_paths[[i]], compress = FALSE)
  message("Saved ", ncol(obj), " cells for ", sample_id, " (", pcw_values[[1]], ").")
  rm(obj, p_comb_06, p_comb_085, p_consensus_06, p_consensus_085)
  gc()
}

message("Combined VZ/RL/consensus labels for all ", length(sample_ids), " configured samples.")
