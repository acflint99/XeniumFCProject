#!/usr/bin/env Rscript

# Join independently reviewed VZ and RL subcluster labels to the approved
# whole-tissue consensus objects. Shared Unknown cells receive only explicitly
# reviewed labels from the VZ/RL annotation-pair decision manifest; the
# original consensus_label remains unchanged as provenance.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_vz_rl_unknown_relabel.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_rl_01_combine_labels.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
base_input_dir <- production_paths$consensus_rds_dir
vz_input_dir <- file.path(output_root, "xenium", "vz", "07_mapped", "rds")
rl_input_dir <- file.path(output_root, "xenium", "rl", "07_mapped", "rds")
rds_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
plot_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "plots", "combined")
table_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "tables")
review_dir <- file.path(
  output_root, "xenium", "vz_rl", "00_unknown_annotation_review", "tables"
)
comparison_pair_path <- file.path(
  review_dir, "Xenium_VZRL_Unknown_VZ_RL_Pairs.csv"
)
comparison_manifest_path <- file.path(
  review_dir, "Xenium_VZRL_Unknown_Annotation_Comparison_manifest.csv"
)
reviewed_relabel_path <- resolve_config_path(
  config$manifests$reviewed_unknown_relabels, config
)

base_input_names <- paste0(sample_ids, "_Consensus_annotated.rds")
vz_input_names <- paste0(sample_ids, "_VZ_QC_Subclusters.rds")
rl_input_names <- paste0(sample_ids, "_RL_QC_Subclusters.rds")
output_names <- paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds")
base_input_paths <- file.path(base_input_dir, base_input_names)
vz_input_paths <- file.path(vz_input_dir, vz_input_names)
rl_input_paths <- file.path(rl_input_dir, rl_input_names)
rds_paths <- file.path(rds_dir, output_names)
plot_suffixes <- c(
  "_combSubcluster_GlobalSpatial_pt0.6.tif",
  "_combSubcluster_GlobalSpatial_pt0.85.tif",
  "_postQC_GlobalSpatial_pt0.85.tif",
  "_postQC_GlobalSpatial_pt0.6.tif"
)
plot_paths <- lapply(sample_ids, function(sample_id) {
  file.path(plot_dir, paste0(sample_id, plot_suffixes))
})
table_paths <- lapply(sample_ids, function(sample_id) {
  file.path(
    table_dir,
    paste0(
      sample_id,
      c(
        "_counts_consensus_label.csv", "_counts_VZ&RLsubclusters.csv",
        "_counts_regional_assignment.csv", "_reviewed_unknown_relabels.csv"
      )
    )
  )
})
expected_by_sample <- lapply(seq_along(sample_ids), function(i) {
  c(rds_paths[[i]], plot_paths[[i]], table_paths[[i]])
})

observed_vz_names <- if (dir.exists(vz_input_dir)) {
  list.files(vz_input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
observed_rl_names <- if (dir.exists(rl_input_dir)) {
  list.files(rl_input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
all_input_paths <- c(
  base_input_paths, vz_input_paths, rl_input_paths,
  comparison_pair_path, comparison_manifest_path, reviewed_relabel_path
)
missing_paths <- all_input_paths[!file.exists(all_input_paths)]
unexpected_vz_names <- setdiff(observed_vz_names, c(vz_input_names, "Xenium_VZ_Mapping_manifest.csv"))
unexpected_rl_names <- setdiff(observed_rl_names, c(rl_input_names, "Xenium_RL_Mapping_manifest.csv"))

expected_comparison_inputs <- rbind(
  data.frame(
    sample_id = sample_ids, input_type = "consensus",
    input_path = base_input_paths, stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = sample_ids, input_type = "VZ",
    input_path = vz_input_paths, stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = sample_ids, input_type = "RL",
    input_path = rl_input_paths, stringsAsFactors = FALSE
  )
)

review_validation <- tryCatch({
  required_review_files <- c(comparison_pair_path, comparison_manifest_path)
  missing_review_files <- required_review_files[!file.exists(required_review_files)]
  if (length(missing_review_files)) {
    stop(
      "Missing Stage-00 comparison file(s): ",
      paste(missing_review_files, collapse = "; ")
    )
  }
  comparison_pairs <- read.csv(
    comparison_pair_path, stringsAsFactors = FALSE, na.strings = character()
  )
  comparison_manifest <- read.csv(
    comparison_manifest_path, stringsAsFactors = FALSE, na.strings = character()
  )
  reviewed_relabels <- load_reviewed_unknown_relabels(config)
  allowed_final_labels <- setdiff(
    unique(c(celltype_order, master_subcluster_order)), "Unknown"
  )
  validate_unknown_relabel_review(
    reviewed_relabels, sample_ids, allowed_final_labels,
    require_complete_labels = TRUE
  )
  validate_unknown_comparison_provenance(
    comparison_manifest, sample_ids, expected_comparison_inputs
  )
  comparison_keys <- unknown_relabel_key(comparison_pairs)
  reviewed_keys <- unknown_relabel_key(reviewed_relabels)
  if (anyDuplicated(comparison_keys)) {
    stop("Unknown annotation-pair comparison contains duplicate keys.")
  }
  if (!setequal(comparison_keys, reviewed_keys)) {
    stop(
      "Reviewed Unknown decisions do not exactly match the current comparison ",
      "template. Missing decisions: ", length(setdiff(comparison_keys, reviewed_keys)),
      "; unused decisions: ", length(setdiff(reviewed_keys, comparison_keys)), "."
    )
  }
  missing_colors <- setdiff(unique(reviewed_relabels$final_label), names(subcluster_palette))
  if (length(missing_colors)) {
    stop(
      "Reviewed Unknown final labels lack plotting colors: ",
      paste(missing_colors, collapse = ", ")
    )
  }
  list(
    valid = TRUE,
    reviewed_relabels = reviewed_relabels,
    reviewed_keys = reviewed_keys,
    allowed_final_labels = allowed_final_labels
  )
}, error = function(error) list(valid = FALSE, error = conditionMessage(error)))

if (dry_run) {
  compact_dry_run(
    "Combined regional labels",
    inputs = all_input_paths,
    outputs = unlist(expected_by_sample, use.names = FALSE),
    checks = c(
      no_unexpected_vz_inputs = !length(unexpected_vz_names),
      no_unexpected_rl_inputs = !length(unexpected_rl_names),
      reviewed_unknown_relabels = isTRUE(review_validation$valid)
    )
  )
  if (length(unexpected_vz_names)) cat(
    "  Unexpected VZ inputs: ", paste(unexpected_vz_names, collapse = "; "), "\n"
  )
  if (length(unexpected_rl_names)) cat(
    "  Unexpected RL inputs: ", paste(unexpected_rl_names, collapse = "; "), "\n"
  )
  if (!isTRUE(review_validation$valid)) {
    cat("  Unknown-relabel review: ", review_validation$error, "\n")
  }
  quit(save = "no", status = 0L)
}

if (length(missing_paths)) {
  stop("Cannot combine subclusters; missing files:\n- ", paste(missing_paths, collapse = "\n- "))
}
if (length(unexpected_vz_names) || length(unexpected_rl_names)) {
  stop(
    "Cannot combine subclusters; unexpected top-level regional RDS files:\n- ",
    paste(c(unexpected_vz_names, unexpected_rl_names), collapse = "\n- ")
  )
}
if (!isTRUE(review_validation$valid)) {
  stop("Cannot combine subclusters; Unknown-relabel review is not ready: ", review_validation$error)
}
reviewed_relabels <- review_validation$reviewed_relabels
reviewed_keys <- review_validation$reviewed_keys
reviewed_final_labels <- setNames(
  as.character(reviewed_relabels$final_label), reviewed_keys
)
reviewed_notes <- setNames(as.character(reviewed_relabels$review_notes), reviewed_keys)
combined_label_order <- unique(c(master_subcluster_order, celltype_order))
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
  base_obj <- readRDS(base_input_paths[[i]])
  vz_obj <- readRDS(vz_input_paths[[i]])
  rl_obj <- readRDS(rl_input_paths[[i]])

  required_base_metadata <- c("consensus_label", "PCW")
  required_vz_metadata <- c("consensus_label", "VZ_subcluster", "VZ_subset_reason")
  required_rl_metadata <- c("consensus_label", "RL_subcluster", "RL_subset_reason")
  missing_base <- setdiff(required_base_metadata, colnames(base_obj[[]]))
  missing_vz <- setdiff(required_vz_metadata, colnames(vz_obj[[]]))
  missing_rl <- setdiff(required_rl_metadata, colnames(rl_obj[[]]))
  if (length(c(missing_base, missing_vz, missing_rl))) {
    stop(
      sample_id, " lacks required metadata (base/VZ/RL): ",
      paste(c(missing_base, missing_vz, missing_rl), collapse = ", ")
    )
  }
  if (anyDuplicated(colnames(base_obj)) || anyDuplicated(colnames(vz_obj)) ||
      anyDuplicated(colnames(rl_obj))) {
    stop(sample_id, " has duplicate cell IDs in a base or regional object.")
  }
  if (!all(colnames(vz_obj) %in% colnames(base_obj)) ||
      !all(colnames(rl_obj) %in% colnames(base_obj))) {
    stop(sample_id, " regional object contains cell IDs absent from the consensus object.")
  }
  base_consensus <- setNames(as.character(base_obj$consensus_label), colnames(base_obj))
  vz_consensus <- as.character(vz_obj$consensus_label)
  rl_consensus <- as.character(rl_obj$consensus_label)
  expected_vz_consensus <- unname(base_consensus[colnames(vz_obj)])
  expected_rl_consensus <- unname(base_consensus[colnames(rl_obj)])
  if (anyNA(c(vz_consensus, rl_consensus, expected_vz_consensus, expected_rl_consensus)) ||
      !identical(unname(vz_consensus), expected_vz_consensus) ||
      !identical(unname(rl_consensus), expected_rl_consensus)) {
    stop(sample_id, " consensus labels changed between the base and regional objects.")
  }
  pcw_values <- unique(as.character(base_obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

  vz_labels <- setNames(as.character(vz_obj$VZ_subcluster), colnames(vz_obj))
  rl_labels <- setNames(as.character(rl_obj$RL_subcluster), colnames(rl_obj))
  vz_reasons <- setNames(as.character(vz_obj$VZ_subset_reason), colnames(vz_obj))
  rl_reasons <- setNames(as.character(rl_obj$RL_subset_reason), colnames(rl_obj))
  base_cells <- colnames(base_obj)
  is_unknown <- is_unknown_consensus_label(base_consensus)
  is_vz_target <- base_consensus %in% unlist(
    config$regional_subsets$vz_broad_labels, use.names = FALSE
  ) | is_unknown
  is_rl_target <- base_consensus %in% unlist(
    config$regional_subsets$rl_broad_labels, use.names = FALSE
  ) | is_unknown
  vz_kept <- !is.na(vz_labels[base_cells]) & nzchar(vz_labels[base_cells])
  rl_kept <- !is.na(rl_labels[base_cells]) & nzchar(rl_labels[base_cells])
  keep_cells <- (!is_vz_target | vz_kept) & (!is_rl_target | rl_kept)
  obj <- subset(base_obj, cells = base_cells[keep_cells])
  retained_cells <- colnames(obj)
  obj$VZ_subcluster <- unname(vz_labels[retained_cells])
  obj$RL_subcluster <- unname(rl_labels[retained_cells])
  obj$VZ_subset_reason <- unname(vz_reasons[retained_cells])
  obj$RL_subset_reason <- unname(rl_reasons[retained_cells])

  has_vz <- !is.na(obj$VZ_subcluster) & nzchar(obj$VZ_subcluster)
  has_rl <- !is.na(obj$RL_subcluster) & nzchar(obj$RL_subcluster)
  retained_unknown <- is_unknown_consensus_label(obj$consensus_label)
  expected_vz <- as.character(obj$consensus_label) %in% unlist(
    config$regional_subsets$vz_broad_labels, use.names = FALSE
  ) | retained_unknown
  expected_rl <- as.character(obj$consensus_label) %in% unlist(
    config$regional_subsets$rl_broad_labels, use.names = FALSE
  ) | retained_unknown
  if (!all(has_vz == expected_vz) || !all(has_rl == expected_rl)) {
    stop(sample_id, " retained regional targets do not have one mapped branch label each.")
  }
  invalid_vz_reason <- has_vz & (
    is.na(obj$VZ_subset_reason) |
      (obj$VZ_subset_reason == "shared_unknown") != retained_unknown
  )
  invalid_rl_reason <- has_rl & (
    is.na(obj$RL_subset_reason) |
      (obj$RL_subset_reason == "shared_unknown") != retained_unknown
  )
  if (any(invalid_vz_reason) || any(invalid_rl_reason)) {
    stop(sample_id, " regional subset provenance does not match consensus Unknown identities.")
  }
  overlap <- has_vz & has_rl
  if (any(overlap & !is_unknown_consensus_label(obj$consensus_label))) {
    stop(sample_id, " has a VZ/RL label overlap outside a consensus Unknown identity.")
  }
  if (any(overlap & (
    is.na(obj$VZ_subset_reason) | is.na(obj$RL_subset_reason) |
      obj$VZ_subset_reason != "shared_unknown" |
      obj$RL_subset_reason != "shared_unknown"
  ))) {
    stop(sample_id, " has overlapping regional labels without shared-Unknown provenance.")
  }
  combined_labels <- as.character(obj$consensus_label)
  combined_labels[has_vz & !has_rl] <- as.character(obj$VZ_subcluster[has_vz & !has_rl])
  combined_labels[has_rl & !has_vz] <- as.character(obj$RL_subcluster[has_rl & !has_vz])
  unknown_relabel <- rep(NA_character_, ncol(obj))
  unknown_relabel_key_value <- rep(NA_character_, ncol(obj))
  if (any(overlap)) {
    overlap_data <- data.frame(
      sample_id = sample_id,
      consensus_label = as.character(obj$consensus_label[overlap]),
      VZ_subcluster = as.character(obj$VZ_subcluster[overlap]),
      RL_subcluster = as.character(obj$RL_subcluster[overlap]),
      stringsAsFactors = FALSE
    )
    overlap_keys <- unknown_relabel_key(overlap_data)
    if (any(!overlap_keys %in% reviewed_keys)) {
      stop(sample_id, " has shared Unknown annotation pairs without reviewed decisions.")
    }
    unknown_relabel[overlap] <- unname(reviewed_final_labels[overlap_keys])
    unknown_relabel_key_value[overlap] <- paste0(
      "sample_id=", overlap_data$sample_id,
      ";consensus_label=", overlap_data$consensus_label,
      ";VZ_subcluster=", overlap_data$VZ_subcluster,
      ";RL_subcluster=", overlap_data$RL_subcluster
    )
    combined_labels[overlap] <- unknown_relabel[overlap]
  }
  obj$unknown_relabel <- unknown_relabel
  obj$unknown_relabel_key <- unknown_relabel_key_value
  obj$unknown_relabel_status <- factor(
    ifelse(overlap, "reviewed_annotation_pair", "not_applicable"),
    levels = c("not_applicable", "reviewed_annotation_pair")
  )
  regional_assignment <- rep("unassigned_consensus", ncol(obj))
  regional_assignment[has_vz & !has_rl] <- "VZ"
  regional_assignment[has_rl & !has_vz] <- "RL"
  regional_assignment[overlap] <- "shared_unknown_overlap"
  obj$regional_assignment <- factor(
    regional_assignment,
    levels = c("VZ", "RL", "shared_unknown_overlap", "unassigned_consensus")
  )
  if (anyNA(combined_labels) || any(!nzchar(combined_labels))) {
    stop(sample_id, " has cells without a VZ, RL, or consensus label.")
  }
  unexpected_labels <- setdiff(
    unique(combined_labels)[!is_unknown_consensus_label(unique(combined_labels))],
    combined_label_order
  )
  if (length(unexpected_labels)) {
    stop(
      sample_id, " has combined labels absent from master_subcluster_order: ",
      paste(unexpected_labels, collapse = ", ")
    )
  }
  present_levels <- order_consensus_levels(combined_labels, combined_label_order)
  obj$comb_subcluster <- factor(combined_labels, levels = present_levels)
  Idents(obj) <- "comb_subcluster"
  plot_colors <- subcluster_palette[present_levels]
  unknown_present <- present_levels[is_unknown_consensus_label(present_levels)]
  plot_colors[unknown_present] <- resolve_cluster_colors(unknown_present)
  if (anyNA(plot_colors)) {
    stop(sample_id, " has combined labels without colors in subcluster_palette.")
  }

  consensus_counts <- as.data.frame(table(obj$consensus_label), stringsAsFactors = FALSE)
  colnames(consensus_counts) <- c("Consensus_Label", "Cell_Count")
  combined_counts <- as.data.frame(table(obj$comb_subcluster), stringsAsFactors = FALSE)
  colnames(combined_counts) <- c("Combined_Subcluster", "Cell_Count")
  assignment_counts <- as.data.frame(table(obj$regional_assignment), stringsAsFactors = FALSE)
  colnames(assignment_counts) <- c("Regional_Assignment", "Cell_Count")
  assignment_counts$Base_Input_Cells <- ncol(base_obj)
  assignment_counts$Retained_Cells <- ncol(obj)
  assignment_counts$Removed_By_VZ <- sum(is_vz_target & !vz_kept)
  assignment_counts$Removed_By_RL <- sum(is_rl_target & !rl_kept)
  write.csv(consensus_counts, table_paths[[i]][[1]], row.names = FALSE)
  write.csv(combined_counts, table_paths[[i]][[2]], row.names = FALSE)
  write.csv(assignment_counts, table_paths[[i]][[3]], row.names = FALSE)
  if (any(overlap)) {
    overlap_key_data <- data.frame(
      sample_id = sample_id,
      consensus_label = as.character(obj$consensus_label[overlap]),
      VZ_subcluster = as.character(obj$VZ_subcluster[overlap]),
      RL_subcluster = as.character(obj$RL_subcluster[overlap]),
      stringsAsFactors = FALSE
    )
    overlap_keys <- unknown_relabel_key(overlap_key_data)
    relabel_counts <- aggregate(
      list(Cell_Count = rep(1L, sum(overlap))),
      by = data.frame(
        Original_Consensus_Label = as.character(obj$consensus_label[overlap]),
        VZ_Subcluster = as.character(obj$VZ_subcluster[overlap]),
        RL_Subcluster = as.character(obj$RL_subcluster[overlap]),
        Reviewed_Final_Label = as.character(obj$unknown_relabel[overlap]),
        Review_Notes = unname(reviewed_notes[overlap_keys]),
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
  } else {
    relabel_counts <- data.frame(
      Original_Consensus_Label = character(), VZ_Subcluster = character(),
      RL_Subcluster = character(), Reviewed_Final_Label = character(),
      Review_Notes = character(), Cell_Count = integer(),
      stringsAsFactors = FALSE
    )
  }
  write.csv(relabel_counts, table_paths[[i]][[4]], row.names = FALSE)

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

  consensus_levels <- order_consensus_levels(obj$consensus_label, celltype_order)
  consensus_colors <- resolve_cluster_colors(consensus_levels)
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
  rm(base_obj, vz_obj, rl_obj, obj, p_comb_06, p_comb_085, p_consensus_06, p_consensus_085)
  gc()
}

message("Combined VZ/RL/consensus labels for all ", length(sample_ids), " configured samples.")
