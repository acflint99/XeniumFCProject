#!/usr/bin/env Rscript

# Add user-approved fine-derived broad labels to completed combined sample
# objects. Original combined objects and consensus plots are never modified.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})
here::i_am("scripts/xenium_vz_rl_02_correct_broad_labels.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))
source(here("scripts", "R", "corrected_broad_labels.R"))
mapping_md5 <- unname(tools::md5sum(here("scripts", "R", "corrected_broad_labels.R")))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_02_correct_broad_labels.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
source_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
stage_dir <- file.path(output_root, "xenium", "vz_rl", "02_corrected_broad")
rds_dir <- file.path(stage_dir, "rds")
plot_dir <- file.path(stage_dir, "plots")
table_dir <- file.path(stage_dir, "tables")
input_names <- paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds")
input_paths <- file.path(source_dir, input_names)
output_paths <- file.path(rds_dir, paste0(sample_ids, "_VZ_RL_QC_CorrectedBroad.rds"))
plot_paths <- lapply(sample_ids, function(id) file.path(plot_dir, paste0(
  id, c("_CorrectedBroad_GlobalSpatial_pt0.6.tif",
        "_CorrectedBroad_GlobalSpatial_pt0.85.tif")
)))
count_paths <- file.path(table_dir, paste0(sample_ids, "_counts_corrected_broad_label.csv"))
transition_paths <- file.path(table_dir, paste0(sample_ids, "_consensus_to_corrected_broad.csv"))
manifest_path <- file.path(table_dir, "Xenium_VZRL_corrected_broad_manifest.csv")
expected_outputs <- c(output_paths, unlist(plot_paths, use.names = FALSE),
                      count_paths, transition_paths, manifest_path)
observed_names <- if (dir.exists(source_dir)) {
  list.files(source_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  compact_dry_run(
    "Fine-derived corrected broad labels",
    inputs = input_paths, outputs = expected_outputs,
    checks = c(
      sample_manifest_unique = !anyDuplicated(sample_ids),
      no_unexpected_combined_rds = !length(unexpected_names)
    )
  )
  if (length(unexpected_names)) {
    cat("  Unexpected combined RDS: ", paste(unexpected_names, collapse = "; "), "\n")
  }
  quit(save = "no", status = 0L)
}
if (anyDuplicated(sample_ids)) stop("Sample manifest contains duplicate sample IDs.")
if (any(!file.exists(input_paths))) {
  stop("Missing combined RDS inputs: ", paste(input_paths[!file.exists(input_paths)], collapse = "; "))
}
if (length(unexpected_names)) {
  stop("Unexpected combined RDS inputs: ", paste(unexpected_names, collapse = "; "))
}
existing <- expected_outputs[file.exists(expected_outputs)]
if (length(existing) && !overwrite) {
  stop("Refusing to replace corrected-broad outputs: ",
       paste(existing, collapse = "; "), ". Review and archive them before rerunning.")
}
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  sample_id = sample_ids, input_path = input_paths, output_path = output_paths,
  input_size = numeric(length(sample_ids)), input_mtime = numeric(length(sample_ids)),
  output_size = numeric(length(sample_ids)), cells = integer(length(sample_ids)),
  mapping_md5 = rep(mapping_md5, length(sample_ids)),
  corrected_unknown_cells = integer(length(sample_ids)),
  changed_broad_cells = integer(length(sample_ids)),
  stringsAsFactors = FALSE
)
for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Correcting broad labels ", i, "/", length(sample_ids), ": ", sample_id)
  obj <- readRDS(input_paths[[i]])
  if (ncol(obj) == 0L || anyDuplicated(colnames(obj))) {
    stop(sample_id, " has zero cells or duplicate cell IDs.")
  }
  required <- c("consensus_label", "comb_subcluster", "regional_assignment", "PCW")
  absent <- setdiff(required, colnames(obj[[]]))
  if (length(absent)) stop(sample_id, " lacks metadata: ", paste(absent, collapse = ", "))
  original_ids <- colnames(obj)
  original_consensus <- as.character(obj$consensus_label)
  if (anyNA(original_consensus) || any(!nzchar(original_consensus))) {
    stop(sample_id, " has missing or blank original consensus labels.")
  }
  fine_labels <- as.character(obj$comb_subcluster)
  corrected <- map_corrected_broad_labels(fine_labels)
  missing_order <- setdiff(unique(corrected), celltype_order)
  if (length(missing_order)) stop("Corrected broad labels lack plot order: ", paste(missing_order, collapse = ", "))
  present <- intersect(celltype_order, unique(corrected))
  obj$corrected_broad_label <- factor(corrected, levels = present)
  if (!identical(colnames(obj), original_ids) ||
      !identical(as.character(obj$consensus_label), original_consensus) ||
      anyNA(obj$corrected_broad_label)) {
    stop(sample_id, " changed cell IDs, original consensus labels, or corrected labels.")
  }
  colors <- resolve_cluster_colors(present)
  if (anyNA(colors)) stop(sample_id, " corrected broad plot has missing colors.")
  p_06 <- ImageDimPlot(
    obj, group.by = "corrected_broad_label", cols = colors, size = 0.6
  ) + ggtitle(paste(sample_id, "Corrected Broad Labels"))
  Cairo::CairoTIFF(plot_paths[[i]][[1]], width = 12, height = 10,
                   units = "in", res = 600, bg = "black")
  print(p_06)
  grDevices::dev.off()
  p_085 <- ImageDimPlot(
    obj, group.by = "corrected_broad_label", cols = colors, size = 0.85
  ) + ggtitle(paste(sample_id, "Corrected Broad Labels"))
  Cairo::CairoTIFF(plot_paths[[i]][[2]], width = 12, height = 10,
                   units = "in", res = 600, bg = "black")
  print(p_085)
  grDevices::dev.off()

  counts <- as.data.frame(table(obj$corrected_broad_label), stringsAsFactors = FALSE)
  names(counts) <- c("Corrected_Broad_Label", "Cell_Count")
  counts <- counts[counts$Cell_Count > 0L, , drop = FALSE]
  write.csv(counts, count_paths[[i]], row.names = FALSE)
  transitions <- as.data.frame(
    table(original_consensus, fine_labels, corrected), stringsAsFactors = FALSE
  )
  names(transitions) <- c("Original_Consensus_Label", "Combined_Subcluster",
                          "Corrected_Broad_Label", "Cell_Count")
  transitions <- transitions[transitions$Cell_Count > 0L, , drop = FALSE]
  write.csv(transitions, transition_paths[[i]], row.names = FALSE)
  if (sum(counts$Cell_Count) != ncol(obj) ||
      sum(transitions$Cell_Count) != ncol(obj)) {
    stop(sample_id, " corrected broad count tables do not cover all cells.")
  }
  saveRDS(obj, output_paths[[i]], compress = FALSE)
  input_info <- file.info(input_paths[[i]])
  output_info <- file.info(output_paths[[i]])
  manifest$input_size[[i]] <- as.numeric(input_info$size)
  manifest$input_mtime[[i]] <- as.numeric(input_info$mtime)
  manifest$output_size[[i]] <- as.numeric(output_info$size)
  manifest$cells[[i]] <- ncol(obj)
  manifest$corrected_unknown_cells[[i]] <- sum(corrected == "Unknown")
  manifest$changed_broad_cells[[i]] <- sum(corrected != original_consensus)
  rm(obj, p_06, p_085)
  invisible(gc())
}
if (nrow(manifest) != length(sample_ids) ||
    !identical(manifest$sample_id, sample_ids) ||
    any(manifest$cells <= 0L) || any(manifest$output_size <= 0)) {
  stop("Corrected-broad manifest is incomplete.")
}
write.csv(manifest, manifest_path, row.names = FALSE)
message("Saved ", length(sample_ids), " corrected-broad sample objects and plots: ", stage_dir)
