#!/usr/bin/env Rscript

# Review post-QC VZ clusters or explicitly apply the reviewed labels.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "subcluster_review_markers.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
resolution_args <- args[grepl("^--resolution=", args)]
if (length(resolution_args) > 1L) stop("Specify --resolution only once.")
subcluster_resolution <- if (length(resolution_args)) {
  suppressWarnings(as.numeric(sub("^--resolution=", "", resolution_args)))
} else {
  0.5
}
supported_review_resolutions <- c(0.3, 0.5, 0.8)
if (!is.finite(subcluster_resolution) ||
    !subcluster_resolution %in% supported_review_resolutions) {
  stop(
    "--resolution must be one of: ",
    paste(supported_review_resolutions, collapse = ", ")
  )
}
subcluster_resolution_label <- format(
  subcluster_resolution, scientific = FALSE, trim = TRUE
)
cluster_column <- paste0("Xenium_snn_res.", subcluster_resolution_label)
args <- args[!grepl("^--resolution=", args)]
valid_options <- c(
  "--review-only", "--apply-reviewed-labels", "--plots-only", "--umaps-only",
  "--dry-run", "--overwrite"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_vz_06_annotate_subclusters.R ",
    "(--review-only|--apply-reviewed-labels|--plots-only|--umaps-only) ",
    "[--resolution=0.3|0.5|0.8] [--dry-run|--overwrite]"
  )
}
review_only <- "--review-only" %in% args
apply_reviewed_labels <- "--apply-reviewed-labels" %in% args
plots_only <- "--plots-only" %in% args
umaps_only <- "--umaps-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (sum(c(review_only, apply_reviewed_labels, plots_only, umaps_only)) != 1L) {
  stop(
    "Choose exactly one mode: --review-only, --apply-reviewed-labels, or ",
    "--plots-only/--umaps-only."
  )
}
if ((plots_only || umaps_only) && !identical(subcluster_resolution_label, "0.8")) {
  stop("--plots-only and --umaps-only require --resolution=0.8 for the approved production object.")
}

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "vz", "06_subclusters", "plots")
table_path <- file.path(output_root, "xenium", "vz", "06_subclusters", "tables")
rds_dir <- file.path(output_root, "xenium", "vz", "06_subclusters", "rds")
merged_path <- file.path(output_root, "xenium", "vz", "05_post_qc", "rds", "Xenium_VZ_postQC.rds")
output_path <- file.path(rds_dir, "Xenium_VZ_subclusters.rds")
label_manifest_path <- file.path(table_path, "Xenium_VZ_subcluster_labels.csv")
review_top10_marker_path <- file.path(
  table_path,
  paste0(
    "VZ_RawSubcluster_postQC_top10_Markers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_top50_marker_path <- file.path(
  table_path,
  paste0(
    "VZ_RawSubcluster_postQC_top50_Markers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_top50_antimarker_path <- file.path(
  table_path,
  paste0(
    "VZ_RawSubcluster_postQC_top50_AntiMarkers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_manifest_path <- file.path(
  table_path,
  paste0(
    "Xenium_VZ_subcluster_review_manifest_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_umap_path <- file.path(
  plot_path,
  paste0(
    "Xenium_VZ_RawSubcluster_postQC_UMAP_Res",
    subcluster_resolution_label,
    ".tif"
  )
)
review_focused_marker_availability_path <- file.path(
  table_path,
  paste0(
    "VZ_RawSubcluster_FocusedMarker_Availability_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_focused_dotplot_tiff_path <- file.path(
  plot_path,
  paste0(
    "Xenium_VZ_RawSubcluster_FocusedMarker_DotPlot_Res",
    subcluster_resolution_label,
    ".tif"
  )
)
review_focused_dotplot_pdf_path <- file.path(
  plot_path,
  paste0(
    "Xenium_VZ_RawSubcluster_FocusedMarker_DotPlot_Res",
    subcluster_resolution_label,
    ".pdf"
  )
)
replicate_marker_dir <- file.path(
  table_path, "replicate_marker_evidence"
)
replicate_marker_paths <- file.path(replicate_marker_dir, c(
  "Xenium_VZ_ReplicateAware_Marker_Evidence.csv",
  "Xenium_VZ_Positive_Marker_Candidates.csv",
  "Xenium_VZ_AntiMarker_Candidates.csv",
  "Xenium_VZ_Marker_Sample_Coverage.csv",
  "Xenium_VZ_ReplicateAware_Marker_Provenance.csv"
))
production_replicate_evidence_required <- identical(
  subcluster_resolution_label, "0.8"
)
review_outputs <- c(
  review_top10_marker_path, review_top50_marker_path,
  review_top50_antimarker_path,
  review_manifest_path, review_umap_path,
  review_focused_marker_availability_path,
  review_focused_dotplot_tiff_path, review_focused_dotplot_pdf_path
)
expected_outputs <- c(
  output_path,
  label_manifest_path,
  file.path(table_path, "Xenium_VZ_Subcluster_QC_Summary.csv"),
  file.path(plot_path, "Xenium_VZ_Subcluster_QC_Violins.pdf"),
  file.path(plot_path, "Xenium_VZ_Subcluster_UMAP.tif"),
  file.path(plot_path, "Xenium_VZ_Lineage_Subcluster_UMAP.tif"),
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_VZ_Lineage_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_VZ_Lineage_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_VZ_Other_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_VZ_Other_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_VZ_Subcluster_Top5_Heatmap.tif"),
  file.path(plot_path, "Xenium_VZ_Subcluster_Top5_Heatmap.pdf")
)
plot_refresh_outputs <- setdiff(expected_outputs, c(output_path, label_manifest_path))
umap_refresh_outputs <- c(
  file.path(plot_path, "Xenium_VZ_Subcluster_UMAP.tif"),
  file.path(plot_path, "Xenium_VZ_Lineage_Subcluster_UMAP.tif")
)

file_is_nonempty <- function(path) {
  info <- suppressWarnings(file.info(path))
  file.exists(path) && !is.na(info$size) && info$size > 0
}

file_set_signature <- function(paths) {
  if (!all(vapply(paths, file_is_nonempty, logical(1)))) return(NA_character_)
  info <- file.info(paths)
  paste(
    paste(
      normalizePath(paths, mustWork = TRUE),
      as.numeric(info$size), as.numeric(info$mtime), sep = "::"
    ),
    collapse = "||"
  )
}

replicate_marker_provenance_matches_input <- function() {
  provenance_path <- tail(replicate_marker_paths, 1L)
  if (!file.exists(merged_path) || !file.exists(provenance_path)) return(FALSE)
  provenance <- tryCatch(
    read.csv(provenance_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  required <- c(
    "input_rds", "input_size", "input_mtime", "cluster_column",
    "inference_unit", "cell_level_p_values_used_for_inference"
  )
  if (is.null(provenance) || nrow(provenance) != 1L ||
      !all(required %in% names(provenance))) return(FALSE)
  input_info <- file.info(merged_path)
  identical(
    as.character(provenance$input_rds),
    normalizePath(merged_path, mustWork = TRUE)
  ) &&
    isTRUE(all.equal(as.numeric(provenance$input_size), as.numeric(input_info$size))) &&
    isTRUE(all.equal(as.numeric(provenance$input_mtime), as.numeric(input_info$mtime))) &&
    identical(as.character(provenance$cluster_column), cluster_column) &&
    identical(as.character(provenance$inference_unit), "biological sample") &&
    identical(as.logical(provenance$cell_level_p_values_used_for_inference), FALSE)
}

review_manifest_matches_input_file <- function() {
  if (!file.exists(merged_path) || !file.exists(review_manifest_path)) return(FALSE)
  review <- tryCatch(
    read.csv(review_manifest_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(review) || nrow(review) != 1L ||
      !all(c("input_size", "input_mtime", "cluster_column") %in% names(review))) {
    return(FALSE)
  }
  input_info <- file.info(merged_path)
  isTRUE(all.equal(as.numeric(review$input_size), as.numeric(input_info$size))) &&
    isTRUE(all.equal(as.numeric(review$input_mtime), as.numeric(input_info$mtime))) &&
    identical(as.character(review$cluster_column), cluster_column)
}

if (dry_run) {
  mode_label <- if (review_only) {
    "VZ post-QC subcluster review"
  } else if (apply_reviewed_labels) {
    "VZ reviewed-label application"
  } else if (plots_only) {
    "VZ plot refresh"
  } else {
    "VZ UMAP refresh"
  }
  dry_inputs <- if (review_only) {
    c(
      merged_path,
      if (production_replicate_evidence_required) replicate_marker_paths else character()
    )
  } else if (apply_reviewed_labels) {
    c(
      merged_path, review_outputs,
      if (production_replicate_evidence_required) replicate_marker_paths else character()
    )
  } else {
    output_path
  }
  dry_checks <- if (review_only && production_replicate_evidence_required) {
    c(
      replicate_marker_outputs_nonempty = all(vapply(
        replicate_marker_paths, file_is_nonempty, logical(1)
      )),
      replicate_marker_provenance_matches_input =
        replicate_marker_provenance_matches_input()
    )
  } else if (apply_reviewed_labels) {
    c(
      review_outputs_nonempty = all(vapply(review_outputs, file_is_nonempty, logical(1))),
      review_manifest_matches_input = review_manifest_matches_input_file(),
      replicate_marker_outputs_nonempty = !production_replicate_evidence_required ||
        all(vapply(replicate_marker_paths, file_is_nonempty, logical(1))),
      replicate_marker_provenance_matches_input = !production_replicate_evidence_required ||
        replicate_marker_provenance_matches_input()
    )
  } else logical()
  compact_dry_run(
    mode_label,
    inputs = dry_inputs,
    outputs = if (review_only) {
      review_outputs
    } else if (apply_reviewed_labels) {
      expected_outputs
    } else if (plots_only) {
      plot_refresh_outputs
    } else {
      umap_refresh_outputs
    },
    checks = dry_checks
  )
  quit(
    save = "no",
    status = if (all(file.exists(dry_inputs)) && all(dry_checks)) 0L else 1L
  )
}

if (plots_only || umaps_only) {
  if (!file.exists(output_path)) stop("Annotated VZ input not found: ", output_path)
} else if (!file.exists(merged_path)) {
  stop("VZ post-QC input not found: ", merged_path)
}
if (apply_reviewed_labels) {
  missing_review_outputs <- review_outputs[!file.exists(review_outputs)]
  if (length(missing_review_outputs)) {
    stop(
      "Run --review-only and inspect its outputs before applying VZ labels. Missing:\n- ",
      paste(missing_review_outputs, collapse = "\n- ")
    )
  }
  empty_review_outputs <- review_outputs[
    !vapply(review_outputs, file_is_nonempty, logical(1))
  ]
  if (length(empty_review_outputs)) {
    stop("VZ review outputs are empty:\n- ", paste(empty_review_outputs, collapse = "\n- "))
  }
}
if ((review_only || apply_reviewed_labels) &&
    production_replicate_evidence_required) {
  missing_replicate_outputs <- replicate_marker_paths[
    !vapply(replicate_marker_paths, file_is_nonempty, logical(1))
  ]
  if (length(missing_replicate_outputs)) {
    stop(
      "Run xenium_vz_06a_build_replicate_marker_evidence.R before the ",
      "production resolution-0.8 review. Missing or empty:\n- ",
      paste(missing_replicate_outputs, collapse = "\n- ")
    )
  }
  if (!replicate_marker_provenance_matches_input()) {
    stop(
      "VZ replicate-aware marker evidence does not match the current ",
      "post-QC object and resolution-0.8 cluster column; rerun Stage 06a."
    )
  }
}
mode_outputs <- if (review_only) {
  review_outputs
} else if (umaps_only) {
  umap_refresh_outputs
} else if (plots_only) {
  plot_refresh_outputs
} else {
  expected_outputs
}
existing_outputs <- mode_outputs[file.exists(mode_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ ",
    if (review_only) "review" else if (umaps_only) "UMAP" else if (plots_only) "plot" else "subcluster",
    " outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(patchwork)
})
source(here("scripts", "color_palette.R"))
vz_review_marker_sets <- if (!plots_only && !umaps_only) {
  add_palette_markers_to_review(
    focused_marker_sets = vz_focused_review_markers,
    palette_marker_sets = vz_markers,
    palette_object_name = "vz_markers"
  )
} else NULL

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
}
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(config$runtime$random_seed)

with_future_workers <- function(expr) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  if (workers > 1L) plan(multisession, workers = workers)
  force(expr)
}

check_mem <- function(step_label) {
  m <- gc(full = TRUE)
  message(paste0("\n[", Sys.time(), "] --- ", step_label, " ---"))
  message("Memory in use: ", round(sum(m[, 2]), 1), " MB\n")
}

dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

vz_umap_label_locations <- function(object) {
  locations <- as.data.frame(
    Embeddings(object, reduction = "umap_clean")[, 1:2, drop = FALSE]
  )
  colnames(locations) <- c("x", "y")
  locations$VZ_subcluster <- as.character(Idents(object))
  locations %>%
    group_by(VZ_subcluster) %>%
    summarise(x = median(x), y = median(y), .groups = "drop")
}

add_vz_umap_labels <- function(plot, label_locations) {
  plot +
    ggrepel::geom_label_repel(
      data = label_locations,
      aes(x = x, y = y, label = VZ_subcluster),
      inherit.aes = FALSE,
      size = 4,
      fill = "white",
      color = "black",
      label.padding = grid::unit(0.15, "lines"),
      label.r = grid::unit(0.12, "lines"),
      linewidth = 0.25,
      box.padding = 0.35,
      point.padding = 0.1,
      max.overlaps = Inf,
      segment.color = NA,
      seed = config$runtime$random_seed
    )
}

plot_vz_umaps <- function(object) {
  full_label_locations <- vz_umap_label_locations(object)
  p <- DimPlot(
    object,
    reduction = "umap_clean",
    label = FALSE,
    cols = vz_palette
  ) +
    ggtitle("Xenium Merged VZ Subclusters") +
    theme(legend.text = element_text(size = 8))
  p <- add_vz_umap_labels(p, full_label_locations)

  Cairo::CairoTIFF(
    filename = file.path(plot_path, "Xenium_VZ_Subcluster_UMAP.tif"),
    width = 12,
    height = 9,
    units = "in",
    res = 600
  )
  print(p)
  grDevices::dev.off()

  clusters_to_remove <- intersect(
    c(
      "RL VZ", "Immature UBCs", "Mature UBCs",
      "eCN", "Microglia", "Cycling Cells"
    ),
    unique(as.character(Idents(object)))
  )
  object_subset <- if (length(clusters_to_remove)) {
    subset(object, idents = clusters_to_remove, invert = TRUE)
  } else {
    object
  }
  filtered_label_locations <- vz_umap_label_locations(object_subset)
  p2 <- DimPlot(
    object_subset,
    reduction = "umap_clean",
    label = FALSE,
    cols = vz_palette
  ) +
    ggtitle("Xenium VZ Lineage Subclusters") +
    theme(legend.text = element_text(size = 8))
  p2 <- add_vz_umap_labels(p2, filtered_label_locations)

  Cairo::CairoTIFF(
    filename = file.path(plot_path, "Xenium_VZ_Lineage_Subcluster_UMAP.tif"),
    width = 12,
    height = 9,
    units = "in",
    res = 600
  )
  print(p2)
  grDevices::dev.off()
  invisible(NULL)
}

if (umaps_only) {
  message("Loading the existing annotated VZ object for UMAP refresh...")
  obj <- readRDS(output_path)
  if (!"VZ_subcluster" %in% colnames(obj[[]])) {
    stop("Annotated VZ input lacks VZ_subcluster metadata.")
  }
  if (!"umap_clean" %in% Reductions(obj)) {
    stop("Annotated VZ input lacks the umap_clean reduction.")
  }
  if (!"Xenium" %in% Assays(obj)) stop("Annotated VZ input lacks the Xenium assay.")
  if (ncol(obj) == 0L) stop("Annotated VZ input contains zero cells.")

  observed_vz_labels <- unique(as.character(obj$VZ_subcluster))
  observed_vz_labels <- observed_vz_labels[
    !is.na(observed_vz_labels) & nzchar(observed_vz_labels)
  ]
  unsupported_vz_labels <- setdiff(observed_vz_labels, vz_subcluster_order)
  missing_vz_colors <- setdiff(observed_vz_labels, names(vz_palette))
  if (length(unsupported_vz_labels) || length(missing_vz_colors)) {
    stop(
      "Annotated VZ labels are not fully registered. Missing from order: ",
      if (length(unsupported_vz_labels)) {
        paste(unsupported_vz_labels, collapse = ", ")
      } else "none",
      "; missing from palette: ",
      if (length(missing_vz_colors)) {
        paste(missing_vz_colors, collapse = ", ")
      } else "none"
    )
  }
  active_vz_order <- vz_subcluster_order[vz_subcluster_order %in% observed_vz_labels]
  obj$VZ_subcluster <- factor(
    as.character(obj$VZ_subcluster),
    levels = active_vz_order
  )
  if (anyNA(obj$VZ_subcluster)) {
    stop("VZ_subcluster factor conversion introduced missing labels during UMAP refresh.")
  }
  Idents(obj) <- "VZ_subcluster"
  plot_vz_umaps(obj)
  message("VZ UMAP refresh complete; the annotated object was not modified.")
  quit(save = "no", status = 0L)
}

if (plots_only) {
  message("Loading the existing annotated VZ object for plot refresh...")
  obj <- readRDS(output_path)
  required_metadata <- c(
    "orig.ident", "consensus_label", "PCW", "VZ_subcluster",
    "nCount_Xenium", "nFeature_Xenium"
  )
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop("Annotated VZ input lacks metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!"umap_clean" %in% Reductions(obj)) stop("Annotated VZ input lacks the umap_clean reduction.")
  if (!"Xenium" %in% Assays(obj)) stop("Annotated VZ input lacks the Xenium assay.")
  if (ncol(obj) == 0L) stop("Annotated VZ input contains zero cells.")

  observed_vz_labels <- unique(as.character(obj$VZ_subcluster))
  observed_vz_labels <- observed_vz_labels[
    !is.na(observed_vz_labels) & nzchar(observed_vz_labels)
  ]
  unsupported_vz_labels <- setdiff(observed_vz_labels, vz_subcluster_order)
  missing_vz_colors <- setdiff(observed_vz_labels, names(vz_palette))
  if (length(unsupported_vz_labels) || length(missing_vz_colors)) {
    stop(
      "Annotated VZ labels are not fully registered. Missing from order: ",
      if (length(unsupported_vz_labels)) {
        paste(unsupported_vz_labels, collapse = ", ")
      } else "none",
      "; missing from palette: ",
      if (length(missing_vz_colors)) {
        paste(missing_vz_colors, collapse = ", ")
      } else "none"
    )
  }
  active_vz_order <- vz_subcluster_order[vz_subcluster_order %in% observed_vz_labels]
  obj$VZ_subcluster <- factor(as.character(obj$VZ_subcluster), levels = active_vz_order)
  if (anyNA(obj$VZ_subcluster)) {
    stop("VZ_subcluster factor conversion introduced missing labels during plot refresh.")
  }
  Idents(obj) <- "VZ_subcluster"
} else {
obj <- readRDS(merged_path)
required_metadata <- c(
  "orig.ident", "consensus_label", "regional_consensus_label",
  "regional_consensus_label_pooled", "regional_unknown_label",
  "PCW", cluster_column
)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("VZ post-QC input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
  stop("VZ post-QC input has invalid sample-qualified Unknown labels.")
}
if (!"umap_clean" %in% Reductions(obj)) stop("VZ post-QC input lacks the umap_clean reduction.")
if (!"Xenium" %in% Assays(obj)) stop("VZ post-QC input lacks the Xenium assay.")
if (ncol(obj) == 0L) stop("VZ post-QC input contains zero cells.")

Idents(obj) <- cluster_column
current_clusters <- levels(Idents(obj))

if (review_only) {
  DefaultAssay(obj) <- "Xenium"
  check_mem("STARTING POST-QC SUBCLUSTER REVIEW")
  message(Sys.time(), ": identifying raw VZ cluster markers for review...")
  all_markers <- with_future_workers(
    FindAllMarkers(
      obj,
      only.pos = FALSE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident = 1000,
      random.seed = config$runtime$random_seed
    )
  )
  if (!nrow(all_markers)) stop("No VZ post-QC markers were returned for review.")
  positive_markers <- all_markers %>% filter(avg_log2FC > 0)
  anti_markers <- all_markers %>% filter(avg_log2FC < 0)
  if (!nrow(positive_markers)) stop("No positive VZ post-QC markers were returned for review.")
  if (!nrow(anti_markers)) stop("No VZ post-QC anti-markers were returned for review.")
  top10_markers <- positive_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  top50_markers <- positive_markers %>%
    group_by(cluster) %>%
    slice_max(n = 50, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  write.csv(top10_markers, review_top10_marker_path, row.names = FALSE)
  write.csv(top50_markers, review_top50_marker_path, row.names = FALSE)
  top50_antimarkers <- anti_markers %>%
    group_by(cluster) %>%
    slice_min(n = 50, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  write.csv(top50_antimarkers, review_top50_antimarker_path, row.names = FALSE)

  p_review <- DimPlot(
    obj,
    reduction = "umap_clean",
    group.by = cluster_column,
    label = TRUE,
    label.size = 4,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle(
      paste0(
        "VZ post-QC raw clusters - resolution ",
        subcluster_resolution_label
      )
    ) +
    theme_classic()
  Cairo::CairoTIFF(
    filename = review_umap_path,
    width = 12,
    height = 9,
    units = "in",
    res = 600
  )
  print(p_review)
  grDevices::dev.off()

  message(Sys.time(), ": generating focused VZ lineage-marker DotPlot...")
  write_focused_marker_dotplot(
    object = obj,
    marker_sets = vz_review_marker_sets,
    assay = "Xenium",
    group_by = cluster_column,
    branch = "VZ",
    resolution_label = subcluster_resolution_label,
    availability_path = review_focused_marker_availability_path,
    tiff_path = review_focused_dotplot_tiff_path,
    pdf_path = review_focused_dotplot_pdf_path
  )

  input_info <- file.info(merged_path)
  write.csv(
    data.frame(
      input_rds = normalizePath(merged_path, mustWork = TRUE),
      input_size = as.numeric(input_info$size),
      input_mtime = as.numeric(input_info$mtime),
      cells = ncol(obj),
      cluster_column = cluster_column,
      clusters = paste(current_clusters, collapse = "|"),
      random_seed = config$runtime$random_seed,
      min_pct = 0.25,
      logfc_threshold = 0.25,
      max_cells_per_ident = 1000,
      marker_table_depths = "10_positive|50_positive|50_anti",
      focused_marker_panel = paste(
        "vz_focused_review_markers",
        "color_palette.R::vz_markers",
        sep = " + "
      ),
      focused_marker_signature = focused_marker_set_signature(
        vz_review_marker_sets
      ),
      replicate_marker_evidence_signature = if (
        production_replicate_evidence_required
      ) file_set_signature(replicate_marker_paths) else "not_required_nonproduction_resolution",
      replicate_marker_inference_unit = if (
        production_replicate_evidence_required
      ) "biological sample" else "not_required_nonproduction_resolution",
      cell_level_marker_use = "descriptive ranking and visual review only",
      stringsAsFactors = FALSE
    ),
    review_manifest_path,
    row.names = FALSE
  )
  message("VZ post-QC review outputs are ready. No subcluster labels were applied.")
  quit(save = "no", status = 0L)
}

review <- read.csv(review_manifest_path, stringsAsFactors = FALSE)
required_review_columns <- c(
  "input_rds", "input_size", "input_mtime", "cells", "cluster_column", "clusters",
  "random_seed", "marker_table_depths", "focused_marker_panel",
  "focused_marker_signature", "replicate_marker_evidence_signature",
  "replicate_marker_inference_unit", "cell_level_marker_use"
)
if (nrow(review) != 1L || !all(required_review_columns %in% names(review))) {
  stop("VZ review manifest is invalid; rerun --review-only.")
}
input_info <- file.info(merged_path)
review_is_current <-
  identical(as.character(review$input_rds), normalizePath(merged_path, mustWork = TRUE)) &&
  isTRUE(all.equal(as.numeric(review$input_size), as.numeric(input_info$size))) &&
  isTRUE(all.equal(as.numeric(review$input_mtime), as.numeric(input_info$mtime))) &&
  identical(as.integer(review$cells), as.integer(ncol(obj))) &&
  identical(as.character(review$cluster_column), cluster_column) &&
  identical(as.character(review$clusters), paste(current_clusters, collapse = "|")) &&
  identical(as.integer(review$random_seed), as.integer(config$runtime$random_seed)) &&
  identical(
    as.character(review$marker_table_depths),
    "10_positive|50_positive|50_anti"
  ) &&
  identical(
    as.character(review$focused_marker_panel),
    "vz_focused_review_markers + color_palette.R::vz_markers"
  ) &&
  identical(
    as.character(review$focused_marker_signature),
    focused_marker_set_signature(vz_review_marker_sets)
  ) &&
  identical(
    as.character(review$replicate_marker_evidence_signature),
    if (production_replicate_evidence_required) {
      file_set_signature(replicate_marker_paths)
    } else "not_required_nonproduction_resolution"
  ) &&
  identical(
    as.character(review$replicate_marker_inference_unit),
    if (production_replicate_evidence_required) {
      "biological sample"
    } else "not_required_nonproduction_resolution"
  ) &&
  identical(
    as.character(review$cell_level_marker_use),
    "descriptive ranking and visual review only"
  )
if (!review_is_current) {
  stop("VZ post-QC input or clusters changed after review; rerun --review-only.")
}

# The resolution-0.8 map was reviewed against the current 30-dimension post-QC
# UMAP, top-50 marker table, focused marker DotPlot, and published developmental
# cerebellar markers. Cell-cycle subdivisions retain their reviewed lineage
# where the lineage signal is supported.
reviewed_label_maps <- list(
  "0.8" = c(
    "0" = "BG",
    "1" = "Maturing PCs",
    "2" = "GABA Progenitors",
    "3" = "Early-born PCs",
    "4" = "RG Progenitors",
    "5" = "Late-born PCs",
    "6" = "Astrocytes/Ependyma",
    "7" = "Cycling Glia",
    "8" = "Patterning PCs",
    "9" = "eCN",
    "10" = "GABA Progenitors",
    "11" = "VZPs",
    "12" = "Cycling VZPs",
    "13" = "iCN",
    "14" = "Maturing PCs",
    "15" = "Golgi Cells",
    "16" = "Pre-OPCs",
    "17" = "OPCs",
    "18" = "RL VZ",
    "19" = "Immature UBCs",
    "20" = "Microglia",
    "21" = "Mature UBCs",
    "22" = "MLIs",
    "23" = "Cycling VZPs",
    "24" = "Cycling Cells",
    "25" = "iCN"
  )
)
new_labels <- reviewed_label_maps[[subcluster_resolution_label]]
if (is.null(new_labels) || !length(new_labels)) {
  stop(
    "No reviewed VZ label map is configured for ", cluster_column, ". ",
    "Run --review-only, inspect the resolution-specific UMAP and top-10/top-50 ",
    "marker tables, then add the approved cluster-to-label mapping to ",
    "reviewed_label_maps."
  )
}

missing_label_mappings <- setdiff(current_clusters, names(new_labels))
unused_label_mappings <- setdiff(names(new_labels), current_clusters)
if (length(missing_label_mappings) || length(unused_label_mappings)) {
  stop(
    "VZ label mapping does not exactly match ", cluster_column, ". Missing: ",
    if (length(missing_label_mappings)) paste(missing_label_mappings, collapse = ", ") else "none",
    "; unused: ",
    if (length(unused_label_mappings)) paste(unused_label_mappings, collapse = ", ") else "none"
  )
}
if (anyNA(new_labels) || any(!nzchar(new_labels))) stop("VZ label mapping contains missing or blank labels.")

reviewed_labels <- unique(unname(new_labels))
labels_missing_from_order <- setdiff(reviewed_labels, vz_subcluster_order)
labels_missing_from_palette <- setdiff(reviewed_labels, names(vz_palette))
if (length(labels_missing_from_order) || length(labels_missing_from_palette)) {
  stop(
    "Reviewed VZ labels are not fully registered. Missing from order: ",
    if (length(labels_missing_from_order)) {
      paste(labels_missing_from_order, collapse = ", ")
    } else "none",
    "; missing from palette: ",
    if (length(labels_missing_from_palette)) {
      paste(labels_missing_from_palette, collapse = ", ")
    } else "none"
  )
}

obj <- RenameIdents(obj, new_labels)

# Final metadata assignment
obj$VZ_subcluster <- Idents(obj)

write.csv(
  data.frame(
    cluster_id = names(new_labels),
    VZ_subcluster = unname(new_labels),
    source_column = cluster_column,
    input_rds = merged_path,
    output_rds = output_path,
    stringsAsFactors = FALSE
  ),
  label_manifest_path,
  row.names = FALSE
)
saveRDS(obj, output_path, compress = FALSE)

# 3. APPLY STANDARDIZED ORDERING
# We use 'subcluster_order' from your color_palette.R script
obj$VZ_subcluster <- factor(Idents(obj), levels = vz_subcluster_order)
if (anyNA(obj$VZ_subcluster)) {
  stop("VZ_subcluster factor conversion introduced missing labels.")
}
Idents(obj) <- "VZ_subcluster"
active_vz_order <- vz_subcluster_order[
  vz_subcluster_order %in% unique(as.character(obj$VZ_subcluster))
]
}

DefaultAssay(obj) <- "Xenium"

# Export Merged Cluster QC Statistics
merged_qc_stats <- obj@meta.data %>%
  group_by(VZ_subcluster) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(merged_qc_stats,
          file.path(table_path, "Xenium_VZ_Subcluster_QC_Summary.csv"),
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "Xenium_VZ_Subcluster_QC_Violins.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)

  p <- VlnPlot(
    obj,
    features = feat,
    group.by = "VZ_subcluster",
    cols = vz_palette,
    pt.size = 0      # Remove points for speed and smaller file size
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"), # Explicitly style Y
      legend.position = "none"
    ) +
    ylab(feat) +
    ggtitle(paste("Xenium VZ Subcluster QC-", feat))

  # Printing the plot inside the loop sends it to the PDF device
  print(p)
}

# 5. Close the device to finalize the file
dev.off()

message("Individual QC plots saved to: ", pdf_path)

# 4. UMAPS WITH CONSISTENT COLORS AND WHITE, LINE-FREE LABEL BOXES
plot_vz_umaps(obj)

# 5. DOTPLOT WITH REVERSED ORDER
# Create a temporary reversed factor for the Y-axis
# (Reverse the levels of the factor to get VZP at the top)
obj$VZ_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(vz_subcluster_order))
Idents(obj) <- "VZ_subcluster_rev"

# The production DotPlots use only the approved compact marker lists. The main
# plot contains their unique union across every VZ subcluster; the split plots
# below show the same markers within the VZ-lineage and Other identity groups.
available_vz_features <- rownames(obj[["Xenium"]])
requested_compact_vz_markers <- unlist(
  vz_all_compact_markers, use.names = FALSE
)
missing_compact_vz_markers <- setdiff(
  requested_compact_vz_markers, available_vz_features
)
if (length(missing_compact_vz_markers)) {
  stop(
    "Compact VZ DotPlot marker genes are absent from the Xenium assay: ",
    paste(missing_compact_vz_markers, collapse = ", ")
  )
}
p1 <- workflow_dotplot(obj,
              features = vz_all_compact_markers,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 12, face = "italic"),
    axis.text.y = element_text(size = 12, face = "plain")
  ) +
  ggtitle("Xenium VZ Subcluster Markers")

dotplot_dimensions <- workflow_dotplot_dimensions(p1)
dotplot_width <- dotplot_dimensions[["width"]]
dotplot_height <- dotplot_dimensions[["height"]]

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_SubclusterMarker_DotPlot.tif")),
  width = dotplot_width,
  height = dotplot_height,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "Xenium_VZ_SubclusterMarker_DotPlot.pdf"),
  p1, device = grDevices::cairo_pdf, width = dotplot_width, height = dotplot_height
)

active_vz_lineage <- intersect(vz_lineage_subclusters, active_vz_order)
active_vz_other <- intersect(vz_other_subclusters, active_vz_order)
if (!length(active_vz_lineage) || !length(active_vz_other)) {
  stop("Both VZ lineage and Other identities are required for split DotPlots.")
}

build_vz_split_dotplot <- function(marker_sets, identities, title) {
  workflow_dotplot(
    obj,
    features = marker_sets,
    assay = "Xenium",
    idents = rev(identities),
    cols = c("lightgrey", "red"),
    dot.scale = 6,
    cluster.idents = FALSE,
    scale = length(identities) > 1L
  ) +
    RotatedAxis() +
    theme(
      axis.text.x = element_text(size = 12, face = "italic"),
      axis.text.y = element_text(size = 12, face = "plain")
    ) +
    ggtitle(title)
}

save_vz_split_dotplot <- function(plot, filename_stem) {
  dimensions <- workflow_dotplot_dimensions(plot)
  width <- dimensions[["width"]]
  height <- dimensions[["height"]]
  Cairo::CairoTIFF(
    filename = file.path(plot_path, paste0(filename_stem, ".tif")),
    width = width,
    height = height,
    units = "in",
    res = 600
  )
  print(plot)
  grDevices::dev.off()
  ggplot2::ggsave(
    file.path(plot_path, paste0(filename_stem, ".pdf")),
    plot, device = grDevices::cairo_pdf, width = width, height = height
  )
}

p_vz_lineage <- build_vz_split_dotplot(
  vz_compact_lineage_markers,
  active_vz_lineage,
  "VZ Lineage Subclusters: Compact Markers"
)
save_vz_split_dotplot(
  p_vz_lineage,
  "Xenium_VZ_Lineage_SubclusterMarker_DotPlot"
)

p_vz_other <- build_vz_split_dotplot(
  vz_other_markers,
  active_vz_other,
  "VZ Other Subclusters: Compact Markers"
)
save_vz_split_dotplot(
  p_vz_other,
  "Xenium_VZ_Other_SubclusterMarker_DotPlot"
)

# ----------------------------------------------------------------
# 6. GENERATE TOP 5 MARKERS HEATMAP
# ----------------------------------------------------------------
check_mem("STARTING HEATMAP GENERATION")

# 1. Identify Markers (if not already in environment from earlier)
# We use a lower max.cells.per.ident to speed up the heatmap calculation
message("Finding top 5 markers for heatmap...")
heatmap_markers <- with_future_workers(
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    random.seed = config$runtime$random_seed
  )
)

# 2. Extract top 5 per cluster
top5_markers <- heatmap_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

top5_markers <- rev(top5_markers)

# 3. Scale the data for the specific markers only
# (Necessary for DoHeatmap to show relative expression)
# Just run the plot on a temporary scaled version:
temp_obj <- ScaleData(obj, features = top5_markers, verbose = FALSE)

# 4. Create Heatmap
# We downsample the plot to 100 cells per group so the labels are readable
# 1. Prepare the plot without the default lines
p4 <- DoHeatmap(
  subset(temp_obj, downsample = 100),
  features = top5_markers,
  group.by = "VZ_subcluster",
  group.colors = vz_palette,
  size = 4,
  angle = 45,
  draw.lines = TRUE,
  raster = FALSE
) +
  scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
  guides(color = "none") +
  theme(
    axis.text.y = element_text(size = 6, face = "italic"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Top 5 Markers per VZ Subcluster")

# 4. SAVE AS PDF
ggplot2::ggsave(
  file.path(plot_path, paste0("Xenium_VZ_Subcluster_Top5_Heatmap.pdf")),
  p4, width = 14, height = 12, device = grDevices::cairo_pdf
)

# 5. Save the Heatmap
Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_VZ_Subcluster_Top5_Heatmap.tif")),
  width = 14,
  height = 12,
  units = "in",
  res = 600
)
print(p4)
grDevices::dev.off()

check_mem("HEATMAP COMPLETE")

plan("sequential")

missing_or_empty_outputs <- mode_outputs[
  !vapply(mode_outputs, file_is_nonempty, logical(1))
]
if (length(missing_or_empty_outputs)) {
  stop(
    "VZ Stage 6 did not create every required nonempty output:\n- ",
    paste(missing_or_empty_outputs, collapse = "\n- ")
  )
}

check_mem("PIPELINE COMPLETE - FILE SAVED")
