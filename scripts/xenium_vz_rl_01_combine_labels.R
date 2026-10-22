#!/usr/bin/env Rscript

# Join independently reviewed VZ, RL, and residual subcluster labels to the
# reviewed whole-tissue inputs. Each Unknown population can enter at most one
# branch; residual labels may remain Unknown when no defensible winner exists.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
here::i_am("scripts/xenium_vz_rl_01_combine_labels.R")

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))
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
residual_qc_state <- load_residual_qc_state(config)

output_root <- here(config$project$outputs_dir)
vz_input_dir <- file.path(output_root, "xenium", "vz", "07_mapped", "rds")
rl_input_dir <- file.path(output_root, "xenium", "rl", "07_mapped", "rds")
residual_input_dir <- file.path(output_root, "xenium", "residual", "06_mapped", "rds")
rds_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "rds")
plot_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "plots")
table_dir <- file.path(output_root, "xenium", "vz_rl", "01_combined_labels", "tables")
route_paths <- unknown_reference_route_paths(config)
routing_state <- tryCatch({
  input_manifest <- load_unknown_reviewed_input_manifest(config)
  routes <- load_unknown_reference_routes(config)
  list(valid = TRUE, input_manifest = input_manifest, routes = routes)
}, error = function(error) list(valid = FALSE, error = conditionMessage(error)))
base_input_paths <- if (isTRUE(routing_state$valid)) {
  as.character(routing_state$input_manifest$mapping_input_path)
} else {
  file.path(
    production_paths$consensus_rds_dir,
    paste0(sample_ids, "_Consensus_annotated.rds")
  )
}
vz_input_names <- paste0(sample_ids, "_VZ_QC_Subclusters.rds")
rl_input_names <- paste0(sample_ids, "_RL_QC_Subclusters.rds")
residual_input_names <- paste0(sample_ids, "_Residual_Subclusters.rds")
output_names <- paste0(sample_ids, "_VZ_RL_QC_Subclusters.rds")
vz_input_paths <- file.path(vz_input_dir, vz_input_names)
rl_input_paths <- file.path(rl_input_dir, rl_input_names)
residual_input_paths <- file.path(residual_input_dir, residual_input_names)
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
        "_counts_regional_assignment.csv", "_unknown_reference_routing.csv"
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
observed_residual_names <- if (dir.exists(residual_input_dir)) {
  list.files(residual_input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
all_input_paths <- c(
  base_input_paths, vz_input_paths, rl_input_paths, residual_input_paths,
  route_paths$input_manifest, route_paths$route_manifest
)
missing_paths <- all_input_paths[!file.exists(all_input_paths)]
unexpected_vz_names <- setdiff(observed_vz_names, c(vz_input_names, "Xenium_VZ_Mapping_manifest.csv"))
unexpected_rl_names <- setdiff(observed_rl_names, c(rl_input_names, "Xenium_RL_Mapping_manifest.csv"))
unexpected_residual_names <- setdiff(observed_residual_names, residual_input_names)

if (dry_run) {
  compact_dry_run(
    "Combined regional labels",
    inputs = all_input_paths,
    outputs = unlist(expected_by_sample, use.names = FALSE),
    checks = c(
      no_unexpected_vz_inputs = !length(unexpected_vz_names),
      no_unexpected_rl_inputs = !length(unexpected_rl_names),
      no_unexpected_residual_inputs = !length(unexpected_residual_names),
      reviewed_reference_routing = isTRUE(routing_state$valid)
    )
  )
  if (length(unexpected_vz_names)) cat(
    "  Unexpected VZ inputs: ", paste(unexpected_vz_names, collapse = "; "), "\n"
  )
  if (length(unexpected_rl_names)) cat(
    "  Unexpected RL inputs: ", paste(unexpected_rl_names, collapse = "; "), "\n"
  )
  if (length(unexpected_residual_names)) cat(
    "  Unexpected residual inputs: ",
    paste(unexpected_residual_names, collapse = "; "), "\n"
  )
  if (!isTRUE(routing_state$valid)) {
    cat("  Unknown routing: ", routing_state$error, "\n")
  }
  quit(save = "no", status = 0L)
}

if (length(missing_paths)) {
  stop("Cannot combine subclusters; missing files:\n- ", paste(missing_paths, collapse = "\n- "))
}
if (length(unexpected_vz_names) || length(unexpected_rl_names) ||
    length(unexpected_residual_names)) {
  stop(
    "Cannot combine subclusters; unexpected top-level regional RDS files:\n- ",
    paste(
      c(unexpected_vz_names, unexpected_rl_names, unexpected_residual_names),
      collapse = "\n- "
    )
  )
}
if (!isTRUE(routing_state$valid)) {
  stop("Cannot combine subclusters; Unknown routing is not ready: ", routing_state$error)
}
routes <- routing_state$routes
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
  residual_obj <- readRDS(residual_input_paths[[i]])

  required_base_metadata <- c("consensus_label", "PCW")
  required_vz_metadata <- c("consensus_label", "VZ_subcluster", "VZ_subset_reason")
  required_rl_metadata <- c("consensus_label", "RL_subcluster", "RL_subset_reason")
  required_residual_metadata <- c(
    "consensus_label", "residual_subcluster", "residual_reference_route",
    "residual_subset_reason", "residual_annotation_confidence",
    "residual_annotation_status"
  )
  missing_base <- setdiff(required_base_metadata, colnames(base_obj[[]]))
  missing_vz <- setdiff(required_vz_metadata, colnames(vz_obj[[]]))
  missing_rl <- setdiff(required_rl_metadata, colnames(rl_obj[[]]))
  missing_residual <- setdiff(required_residual_metadata, colnames(residual_obj[[]]))
  if (length(c(missing_base, missing_vz, missing_rl, missing_residual))) {
    stop(
      sample_id, " lacks required metadata (base/VZ/RL/residual): ",
      paste(c(missing_base, missing_vz, missing_rl, missing_residual), collapse = ", ")
    )
  }
  if (anyDuplicated(colnames(base_obj)) || anyDuplicated(colnames(vz_obj)) ||
      anyDuplicated(colnames(rl_obj)) || anyDuplicated(colnames(residual_obj))) {
    stop(sample_id, " has duplicate cell IDs in a base or regional object.")
  }
  if (!all(colnames(vz_obj) %in% colnames(base_obj)) ||
      !all(colnames(rl_obj) %in% colnames(base_obj)) ||
      !identical(colnames(residual_obj), colnames(base_obj))) {
    stop(sample_id, " regional object contains cell IDs absent from the consensus object.")
  }
  base_consensus <- setNames(as.character(base_obj$consensus_label), colnames(base_obj))
  vz_consensus <- as.character(vz_obj$consensus_label)
  rl_consensus <- as.character(rl_obj$consensus_label)
  residual_consensus <- as.character(residual_obj$consensus_label)
  expected_vz_consensus <- unname(base_consensus[colnames(vz_obj)])
  expected_rl_consensus <- unname(base_consensus[colnames(rl_obj)])
  expected_residual_consensus <- unname(base_consensus[colnames(residual_obj)])
  if (anyNA(c(
        vz_consensus, rl_consensus, residual_consensus,
        expected_vz_consensus, expected_rl_consensus, expected_residual_consensus
      )) ||
      !identical(unname(vz_consensus), expected_vz_consensus) ||
      !identical(unname(rl_consensus), expected_rl_consensus) ||
      !identical(unname(residual_consensus), expected_residual_consensus)) {
    stop(sample_id, " consensus labels changed between the base and regional objects.")
  }
  pcw_values <- unique(as.character(base_obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

  vz_labels <- setNames(as.character(vz_obj$VZ_subcluster), colnames(vz_obj))
  rl_labels <- setNames(as.character(rl_obj$RL_subcluster), colnames(rl_obj))
  vz_reasons <- setNames(as.character(vz_obj$VZ_subset_reason), colnames(vz_obj))
  rl_reasons <- setNames(as.character(rl_obj$RL_subset_reason), colnames(rl_obj))
  residual_labels <- setNames(
    as.character(residual_obj$residual_subcluster), colnames(residual_obj)
  )
  residual_routes <- setNames(
    as.character(residual_obj$residual_reference_route), colnames(residual_obj)
  )
  residual_reasons <- setNames(
    as.character(residual_obj$residual_subset_reason), colnames(residual_obj)
  )
  residual_confidence <- setNames(
    as.character(residual_obj$residual_annotation_confidence), colnames(residual_obj)
  )
  residual_status <- setNames(
    as.character(residual_obj$residual_annotation_status), colnames(residual_obj)
  )
  base_cells <- colnames(base_obj)
  removed_residual <- residual_qc_exclusions_for_sample(
    base_cells, base_obj$consensus_label, sample_id, residual_qc_state)
  if (residual_qc_state$active && (
      !"residual_qc_status" %in% names(residual_obj[[]]) ||
      !identical(residual_obj@misc$residual_qc_decision_md5, residual_qc_state$decision_md5) ||
      anyNA(residual_obj$residual_qc_status) ||
      !identical(as.character(residual_obj$residual_qc_status) == "removed", removed_residual))) {
    stop(sample_id, " residual mapping lacks the current QC exclusion provenance.")
  }
  is_unknown <- is_unknown_consensus_label(base_consensus)
  sample_routes <- routes[routes$sample_id == sample_id, , drop = FALSE]
  sample_route_lookup <- setNames(
    as.character(sample_routes$reference_route),
    as.character(sample_routes$consensus_label)
  )
  sample_route_support_lookup <- setNames(
    as.integer(sample_routes$reference_route_support_n),
    as.character(sample_routes$consensus_label)
  )
  sample_route_status_lookup <- setNames(
    as.character(sample_routes$reference_route_status),
    as.character(sample_routes$consensus_label)
  )
  unknown_routes <- unname(sample_route_lookup[base_consensus])
  if (any(is_unknown & is.na(unknown_routes))) {
    stop(sample_id, " has a retained Unknown population without a reference route.")
  }
  is_vz_target <- base_consensus %in% unlist(
    config$regional_subsets$vz_broad_labels, use.names = FALSE
  ) | (is_unknown & !is.na(unknown_routes) & unknown_routes == "route_vz")
  is_rl_target <- base_consensus %in% unlist(
    config$regional_subsets$rl_broad_labels, use.names = FALSE
  ) | (is_unknown & !is.na(unknown_routes) & unknown_routes == "route_rl")
  is_residual_target <- is_unknown & !is.na(unknown_routes) &
    unknown_routes %in% c("route_other", "route_unresolved")
  vz_kept <- !is.na(vz_labels[base_cells]) & nzchar(vz_labels[base_cells])
  rl_kept <- !is.na(rl_labels[base_cells]) & nzchar(rl_labels[base_cells])
  residual_kept <- !is.na(residual_labels[base_cells]) &
    nzchar(residual_labels[base_cells])
  if (any(removed_residual & !is_residual_target) ||
      !identical(unname(residual_kept), unname(is_residual_target & !removed_residual))) {
    stop(sample_id, " residual mapping differs from targets minus explicitly approved QC exclusions.")
  }
  keep_cells <- (!is_vz_target | vz_kept) & (!is_rl_target | rl_kept) & !removed_residual
  obj <- subset(base_obj, cells = base_cells[keep_cells])
  obj@misc$residual_qc_removal <- list(decision_md5 = residual_qc_state$decision_md5,
    removed_cells = sum(removed_residual), source = residual_qc_state$paths$removed)
  retained_cells <- colnames(obj)
  obj$VZ_subcluster <- unname(vz_labels[retained_cells])
  obj$RL_subcluster <- unname(rl_labels[retained_cells])
  obj$VZ_subset_reason <- unname(vz_reasons[retained_cells])
  obj$RL_subset_reason <- unname(rl_reasons[retained_cells])
  obj$residual_subcluster <- unname(residual_labels[retained_cells])
  obj$residual_reference_route <- unname(residual_routes[retained_cells])
  obj$residual_subset_reason <- unname(residual_reasons[retained_cells])
  obj$residual_annotation_confidence <- unname(residual_confidence[retained_cells])
  obj$residual_annotation_status <- unname(residual_status[retained_cells])

  has_vz <- !is.na(obj$VZ_subcluster) & nzchar(obj$VZ_subcluster)
  has_rl <- !is.na(obj$RL_subcluster) & nzchar(obj$RL_subcluster)
  has_residual <- !is.na(obj$residual_subcluster) & nzchar(obj$residual_subcluster)
  retained_unknown <- is_unknown_consensus_label(obj$consensus_label)
  retained_routes <- unname(sample_route_lookup[as.character(obj$consensus_label)])
  expected_vz <- as.character(obj$consensus_label) %in% unlist(
    config$regional_subsets$vz_broad_labels, use.names = FALSE
  ) | (retained_unknown & !is.na(retained_routes) & retained_routes == "route_vz")
  expected_rl <- as.character(obj$consensus_label) %in% unlist(
    config$regional_subsets$rl_broad_labels, use.names = FALSE
  ) | (retained_unknown & !is.na(retained_routes) & retained_routes == "route_rl")
  expected_residual <- retained_unknown & !is.na(retained_routes) &
    retained_routes %in% c("route_other", "route_unresolved")
  if (!all(has_vz == expected_vz) || !all(has_rl == expected_rl) ||
      !all(has_residual == expected_residual)) {
    stop(sample_id, " retained regional targets do not have one mapped branch label each.")
  }
  invalid_vz_reason <- has_vz & (
    is.na(obj$VZ_subset_reason) |
      (obj$VZ_subset_reason == "reference_route_vz") != retained_unknown
  )
  invalid_rl_reason <- has_rl & (
    is.na(obj$RL_subset_reason) |
      (obj$RL_subset_reason == "reference_route_rl") != retained_unknown
  )
  if (any(invalid_vz_reason) || any(invalid_rl_reason)) {
    stop(sample_id, " regional subset provenance does not match consensus Unknown identities.")
  }
  invalid_residual <- has_residual & (
    is.na(obj$residual_reference_route) |
      !obj$residual_reference_route %in% c("route_other", "route_unresolved") |
      obj$residual_reference_route != retained_routes |
      is.na(obj$residual_subset_reason) |
      obj$residual_subset_reason != "retained_high_quality_other_or_unresolved" |
      is.na(obj$residual_annotation_status) |
      !obj$residual_annotation_status %in% c("reviewed_label", "retained_unknown")
  )
  if (any(invalid_residual)) {
    stop(sample_id, " residual metadata does not match reviewed route and label provenance.")
  }
  branch_count <- as.integer(has_vz) + as.integer(has_rl) + as.integer(has_residual)
  if (any(branch_count > 1L)) {
    stop(sample_id, " has cells assigned to more than one regional/residual branch.")
  }
  combined_labels <- as.character(obj$consensus_label)
  combined_labels[has_vz] <- as.character(obj$VZ_subcluster[has_vz])
  combined_labels[has_rl] <- as.character(obj$RL_subcluster[has_rl])
  combined_labels[has_residual] <- as.character(obj$residual_subcluster[has_residual])
  obj$unknown_reference_route <- ifelse(
    retained_unknown, retained_routes, NA_character_
  )
  obj$unknown_reference_route_support_n <- ifelse(
    retained_unknown,
    unname(sample_route_support_lookup[as.character(obj$consensus_label)]),
    NA_integer_
  )
  obj$unknown_reference_route_status <- ifelse(
    retained_unknown,
    unname(sample_route_status_lookup[as.character(obj$consensus_label)]),
    NA_character_
  )
  obj$unknown_relabel_status <- factor(
    ifelse(
      has_residual & obj$residual_annotation_status == "reviewed_label",
      "reviewed_residual_label",
      ifelse(
        has_residual & obj$residual_annotation_status == "retained_unknown",
        "retained_unknown", "not_applicable"
      )
    ),
    levels = c("not_applicable", "reviewed_residual_label", "retained_unknown")
  )
  regional_assignment <- rep("unassigned_consensus", ncol(obj))
  regional_assignment[has_vz] <- "VZ"
  regional_assignment[has_rl] <- "RL"
  regional_assignment[has_residual] <- "Residual"
  obj$regional_assignment <- factor(
    regional_assignment,
    levels = c("VZ", "RL", "Residual", "unassigned_consensus")
  )
  if (anyNA(combined_labels) || any(!nzchar(combined_labels))) {
    stop(sample_id, " has cells without a VZ, RL, residual, or consensus label.")
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
  assignment_counts$Removed_By_Residual <- sum(is_residual_target & !residual_kept)
  write.csv(consensus_counts, table_paths[[i]][[1]], row.names = FALSE)
  write.csv(combined_counts, table_paths[[i]][[2]], row.names = FALSE)
  write.csv(assignment_counts, table_paths[[i]][[3]], row.names = FALSE)
  unknown_index <- is_unknown_consensus_label(obj$consensus_label)
  if (any(unknown_index)) {
    route_counts <- aggregate(
      list(Cell_Count = rep(1L, sum(unknown_index))),
      by = data.frame(
        Original_Consensus_Label = as.character(obj$consensus_label[unknown_index]),
        Reference_Route = as.character(obj$unknown_reference_route[unknown_index]),
        Final_Combined_Label = as.character(obj$comb_subcluster[unknown_index]),
        Regional_Assignment = as.character(obj$regional_assignment[unknown_index]),
        stringsAsFactors = FALSE
      ),
      FUN = sum
    )
  } else {
    route_counts <- data.frame(
      Original_Consensus_Label = character(), Reference_Route = character(),
      Final_Combined_Label = character(), Regional_Assignment = character(),
      Cell_Count = integer(), stringsAsFactors = FALSE
    )
  }
  write.csv(route_counts, table_paths[[i]][[4]], row.names = FALSE)

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
  rm(
    base_obj, vz_obj, rl_obj, residual_obj, obj,
    p_comb_06, p_comb_085, p_consensus_06, p_consensus_085
  )
  gc()
}

message(
  "Combined VZ/RL/residual/consensus labels for all ",
  length(sample_ids), " configured samples."
)
