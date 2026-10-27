#!/usr/bin/env Rscript

# Review post-QC RL clusters, apply reviewed labels, or refresh production plots.
rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "subcluster_review_markers.R"))
source(here("scripts", "R", "xenium_replicate_marker_evidence.R"))
source(here("scripts", "R", "rl_label_approval.R"))

config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
resolution_args <- args[grepl("^--resolution=", args)]
if (length(resolution_args) > 1L) stop("Specify --resolution only once.")
subcluster_resolution <- if (length(resolution_args)) {
  suppressWarnings(as.numeric(sub("^--resolution=", "", resolution_args)))
} else {
  0.8
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
  "--review-only", "--apply-reviewed-labels", "--plots-only",
  "--dry-run", "--overwrite"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_rl_06_annotate_subclusters.R ",
    "(--review-only|--apply-reviewed-labels|--plots-only) ",
    "[--resolution=0.3|0.5|0.8] [--dry-run|--overwrite]"
  )
}
review_only <- "--review-only" %in% args
apply_reviewed_labels <- "--apply-reviewed-labels" %in% args
plots_only <- "--plots-only" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (sum(c(review_only, apply_reviewed_labels, plots_only)) != 1L) {
  stop(
    "Choose exactly one mode: --review-only, --apply-reviewed-labels, or ",
    "--plots-only."
  )
}
if (plots_only && !identical(subcluster_resolution_label, "0.8")) {
  stop("--plots-only requires --resolution=0.8 for the approved production object.")
}

output_root <- here(config$project$outputs_dir)
plot_path <- file.path(output_root, "xenium", "rl", "06_subclusters", "plots")
table_path <- file.path(output_root, "xenium", "rl", "06_subclusters", "tables")
rds_dir <- file.path(output_root, "xenium", "rl", "06_subclusters", "rds")
merged_path <- file.path(output_root, "xenium", "rl", "05_post_qc", "rds", "Xenium_RL_postQC.rds")
output_path <- file.path(rds_dir, "Xenium_RL_subclusters.rds")
label_manifest_path <- file.path(table_path, "Xenium_RL_subcluster_labels.csv")
review_top10_marker_path <- file.path(
  table_path,
  paste0(
    "RL_RawSubcluster_postQC_top10_Markers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_top50_marker_path <- file.path(
  table_path,
  paste0(
    "RL_RawSubcluster_postQC_top50_Markers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_top50_antimarker_path <- file.path(
  table_path,
  paste0(
    "RL_RawSubcluster_postQC_top50_AntiMarkers_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_manifest_path <- file.path(
  table_path,
  paste0(
    "Xenium_RL_subcluster_review_manifest_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
label_approval_template_path <- file.path(
  table_path, paste0("Xenium_RL_Subcluster_Label_Approval_Template_Res",
                     subcluster_resolution_label, ".csv")
)
label_approval_path <- here("scripts", "reviews", paste0(
  "RL_subcluster_labels_Res", subcluster_resolution_label, ".csv"
))
review_umap_path <- file.path(
  plot_path,
  paste0(
    "Xenium_RL_RawSubcluster_postQC_UMAP_Res",
    subcluster_resolution_label,
    ".tif"
  )
)
review_focused_marker_availability_path <- file.path(
  table_path,
  paste0(
    "RL_RawSubcluster_FocusedMarker_Availability_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
review_focused_dotplot_tiff_path <- file.path(
  plot_path,
  paste0(
    "Xenium_RL_RawSubcluster_FocusedMarker_DotPlot_Res",
    subcluster_resolution_label,
    ".tif"
  )
)
review_focused_dotplot_pdf_path <- file.path(
  plot_path,
  paste0(
    "Xenium_RL_RawSubcluster_FocusedMarker_DotPlot_Res",
    subcluster_resolution_label,
    ".pdf"
  )
)
replicate_marker_dir <- file.path(
  table_path, "replicate_marker_evidence"
)
replicate_marker_paths <- file.path(replicate_marker_dir, c(
  "Xenium_RL_ReplicateAware_Marker_Evidence.csv",
  "Xenium_RL_Positive_Marker_Candidates.csv",
  "Xenium_RL_AntiMarker_Candidates.csv",
  "Xenium_RL_Marker_Sample_Coverage.csv",
  "Xenium_RL_ReplicateAware_Marker_Provenance.csv"
))
production_replicate_evidence_required <- identical(
  subcluster_resolution_label, "0.8"
)
review_outputs <- c(
  review_top10_marker_path, review_top50_marker_path,
  review_top50_antimarker_path,
  review_manifest_path, label_approval_template_path, review_umap_path,
  review_focused_marker_availability_path,
  review_focused_dotplot_tiff_path, review_focused_dotplot_pdf_path
)
expected_outputs <- c(
  output_path,
  label_manifest_path,
  file.path(table_path, "Xenium_RL_Subcluster_QC_Summary.csv"),
  file.path(table_path, "Xenium_RL_Subcluster_Markers.csv"),
  file.path(plot_path, "Xenium_RL_Subcluster_QC_Violins.pdf"),
  file.path(plot_path, "Xenium_RL_Subcluster_UMAP.tif"),
  file.path(plot_path, "Xenium_RL_Subcluster_rmBG,Ep,PC_UMAP.tif"),
  file.path(plot_path, "Xenium_RL_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_RL_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_RL_Lineage_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_RL_Lineage_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_RL_Other_SubclusterMarker_DotPlot.tif"),
  file.path(plot_path, "Xenium_RL_Other_SubclusterMarker_DotPlot.pdf"),
  file.path(plot_path, "Xenium_RL_Subcluster_Top5_Heatmap.tif"),
  file.path(plot_path, "Xenium_RL_Subcluster_Top5_Heatmap.pdf")
)
plot_refresh_outputs <- setdiff(expected_outputs, c(output_path, label_manifest_path))

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

# Production review now includes its replicate-evidence prerequisite.
replicate_review_action <- if (review_only && production_replicate_evidence_required) {
  plan_replicate_marker_review(
    replicate_marker_paths, replicate_marker_provenance_matches_input(), overwrite
  )
} else "not_required"

if (dry_run) {
  approval_valid <- if (apply_reviewed_labels) tryCatch({
    validate_rl_label_approval(label_approval_path, merged_path,
                               review_manifest_path, cluster_column)
    TRUE
  }, error = function(error) {
    cat("  RL label approval: ", conditionMessage(error), "\n", sep = "")
    FALSE
  }) else TRUE
  mode_label <- if (review_only) {
    "RL post-QC subcluster review"
  } else if (apply_reviewed_labels) {
    "RL reviewed-label application"
  } else {
    "RL plot refresh"
  }
  dry_inputs <- if (review_only) {
    c(
      merged_path,
      if (identical(replicate_review_action, "reuse")) replicate_marker_paths else character()
    )
  } else if (apply_reviewed_labels) {
    c(
      merged_path, review_outputs, label_approval_path,
      if (production_replicate_evidence_required) replicate_marker_paths else character()
    )
  } else {
    output_path
  }
  dry_checks <- if (review_only && production_replicate_evidence_required) {
    c(replicate_evidence_can_build_or_reuse = replicate_review_action != "blocked")
  } else if (apply_reviewed_labels) {
    c(
      explicit_source_bound_label_approval = approval_valid,
      review_outputs_nonempty = all(vapply(review_outputs, file_is_nonempty, logical(1))),
      review_manifest_matches_input = review_manifest_matches_input_file(),
      replicate_marker_outputs_nonempty = !production_replicate_evidence_required ||
        all(vapply(replicate_marker_paths, file_is_nonempty, logical(1))),
      replicate_marker_provenance_matches_input = !production_replicate_evidence_required ||
        replicate_marker_provenance_matches_input()
    )
  } else logical()
  if (review_only && production_replicate_evidence_required) {
    cat("  Replicate marker evidence: ", replicate_review_action, "\n", sep = "")
    if (replicate_review_action == "blocked") {
      cat("  Existing evidence is partial or stale. Review it before explicitly using --overwrite.\n")
    }
  }
  compact_dry_run(
    mode_label,
    inputs = dry_inputs,
    outputs = if (review_only) {
      c(review_outputs, if (production_replicate_evidence_required &&
        replicate_review_action != "reuse") replicate_marker_paths else character())
    } else if (apply_reviewed_labels) {
      expected_outputs
    } else {
      plot_refresh_outputs
    },
    checks = dry_checks
  )
  quit(
    save = "no",
    status = if (all(file.exists(dry_inputs)) && all(dry_checks)) 0L else 1L
  )
}

if (plots_only) {
  if (!file.exists(output_path)) stop("Annotated RL input not found: ", output_path)
} else if (!file.exists(merged_path)) {
  stop("RL post-QC input not found: ", merged_path)
}
if (apply_reviewed_labels) {
  # Check the explicit human decision before loading the large source object.
  validate_rl_label_approval(label_approval_path, merged_path,
                             review_manifest_path, cluster_column)
  missing_review_outputs <- review_outputs[!file.exists(review_outputs)]
  if (length(missing_review_outputs)) {
    stop(
      "Run --review-only and inspect its outputs before applying RL labels. Missing:\n- ",
      paste(missing_review_outputs, collapse = "\n- ")
    )
  }
  empty_review_outputs <- review_outputs[
    !vapply(review_outputs, file_is_nonempty, logical(1))
  ]
  if (length(empty_review_outputs)) {
    stop("RL review outputs are empty:\n- ", paste(empty_review_outputs, collapse = "\n- "))
  }
}
if (apply_reviewed_labels && production_replicate_evidence_required) {
  missing_replicate_outputs <- replicate_marker_paths[
    !vapply(replicate_marker_paths, file_is_nonempty, logical(1))
  ]
  if (length(missing_replicate_outputs)) {
    stop(
      "Run Stage 6 with --review-only to generate the required ",
      "resolution-0.8 evidence before applying labels. Missing or empty:\n- ",
      paste(missing_replicate_outputs, collapse = "\n- ")
    )
  }
  if (!replicate_marker_provenance_matches_input()) {
    stop(
      "RL replicate-aware marker evidence does not match the current ",
      "post-QC object and resolution-0.8 cluster column; rerun Stage 6 review."
    )
  }
}
mode_outputs <- if (review_only) {
  review_outputs
} else if (apply_reviewed_labels) {
  expected_outputs
} else {
  plot_refresh_outputs
}
existing_outputs <- mode_outputs[file.exists(mode_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL ",
    if (review_only) "review" else if (plots_only) "plot" else "subcluster",
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

if (review_only && production_replicate_evidence_required) {
  if (replicate_review_action == "blocked") {
    stop("Existing RL replicate evidence is partial or stale. ",
         "Inspect it before using --overwrite; no review outputs were written.")
  }
  if (replicate_review_action %in% c("build", "rebuild")) {
    message("Building RL replicate-aware marker evidence before cluster review...")
    build_replicate_marker_evidence(
      "RL", merged_path, cluster_column, replicate_marker_dir,
      overwrite = overwrite
    )
    # Release the builder's large local object before review loads its own object.
    invisible(gc(full = TRUE))
  } else {
    message("Reusing RL replicate-aware marker evidence matched to the input.")
  }
  if (!all(vapply(replicate_marker_paths, file_is_nonempty, logical(1))) ||
      !replicate_marker_provenance_matches_input()) {
    stop("RL replicate evidence is incomplete or fails input provenance validation.")
  }
}

rl_review_marker_sets <- if (!plots_only) {
  add_palette_markers_to_review(
    focused_marker_sets = rl_focused_review_markers,
    palette_marker_sets = rl_markers,
    palette_object_name = "rl_markers"
  )
} else NULL

workers <- suppressWarnings(as.integer(Sys.getenv("RL_SUBCLUSTER_WORKERS", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("RL_SUBCLUSTER_WORKERS must be a positive integer when set.")
}
allocated_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(allocated_cpus) || allocated_cpus < 1L) stop("Invalid SLURM_CPUS_PER_TASK.")
workers <- min(workers, allocated_cpus)
options(future.globals.maxSize = config$runtime$future_globals_max_gb_default * 1024^3)
plan(sequential)
set.seed(config$runtime$random_seed)

with_future_workers <- function(expr) {
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  if (workers > 1L) plan(multisession, workers = workers) else plan(sequential)
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

if (plots_only) {
  message("Loading the existing annotated RL object for plot refresh...")
  obj <- readRDS(output_path)
  required_metadata <- c(
    "orig.ident", "consensus_label", "PCW", "RL_subcluster",
    "nCount_Xenium", "nFeature_Xenium"
  )
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop("Annotated RL input lacks metadata: ", paste(missing_metadata, collapse = ", "))
  }
  if (!"umap_clean" %in% Reductions(obj)) stop("Annotated RL input lacks the umap_clean reduction.")
  if (!"Xenium" %in% Assays(obj)) stop("Annotated RL input lacks the Xenium assay.")
  if (ncol(obj) == 0L) stop("Annotated RL input contains zero cells.")

  observed_rl_labels <- unique(as.character(obj$RL_subcluster))
  observed_rl_labels <- observed_rl_labels[
    !is.na(observed_rl_labels) & nzchar(observed_rl_labels)
  ]
  unsupported_rl_labels <- setdiff(observed_rl_labels, rl_subcluster_order)
  missing_rl_colors <- setdiff(observed_rl_labels, names(rl_palette))
  if (length(unsupported_rl_labels) || length(missing_rl_colors)) {
    stop(
      "Annotated RL labels are not fully registered. Missing from order: ",
      if (length(unsupported_rl_labels)) {
        paste(unsupported_rl_labels, collapse = ", ")
      } else "none",
      "; missing from palette: ",
      if (length(missing_rl_colors)) {
        paste(missing_rl_colors, collapse = ", ")
      } else "none"
    )
  }
  active_rl_order <- rl_subcluster_order[rl_subcluster_order %in% observed_rl_labels]
  obj$RL_subcluster <- factor(as.character(obj$RL_subcluster), levels = active_rl_order)
  if (anyNA(obj$RL_subcluster)) {
    stop("RL_subcluster factor conversion introduced missing labels during plot refresh.")
  }
  Idents(obj) <- "RL_subcluster"
} else {
obj <- readRDS(merged_path)
required_metadata <- c(
  "orig.ident", "consensus_label", "regional_consensus_label",
  "regional_consensus_label_pooled", "regional_unknown_label",
  "PCW", cluster_column
)
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata)) stop("RL post-QC input lacks metadata: ", paste(missing_metadata, collapse = ", "))
if (!validate_regional_unknown_labels(obj, obj$orig.ident)) {
  stop("RL post-QC input has invalid sample-qualified Unknown labels.")
}
if (!"umap_clean" %in% Reductions(obj)) stop("RL post-QC input lacks the umap_clean reduction.")
if (!"Xenium" %in% Assays(obj)) stop("RL post-QC input lacks the Xenium assay.")
if (ncol(obj) == 0L) stop("RL post-QC input contains zero cells.")

Idents(obj) <- cluster_column
current_clusters <- levels(Idents(obj))

if (review_only) {
  DefaultAssay(obj) <- "Xenium"
  check_mem("STARTING POST-QC SUBCLUSTER REVIEW")
  message(Sys.time(), ": identifying raw RL cluster markers for review...")
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
  if (!nrow(all_markers)) stop("No RL post-QC markers were returned for review.")
  positive_markers <- all_markers %>% filter(avg_log2FC > 0)
  anti_markers <- all_markers %>% filter(avg_log2FC < 0)
  if (!nrow(positive_markers)) stop("No positive RL post-QC markers were returned for review.")
  if (!nrow(anti_markers)) stop("No RL post-QC anti-markers were returned for review.")
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
        "RL post-QC raw clusters - resolution ",
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

  message(Sys.time(), ": generating focused RL lineage-marker DotPlot...")
  write_focused_marker_dotplot(
    object = obj,
    marker_sets = rl_review_marker_sets,
    assay = "Xenium",
    group_by = cluster_column,
    branch = "RL",
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
        "rl_focused_review_markers",
        "color_palette.R::rl_markers",
        sep = " + "
      ),
      focused_marker_signature = focused_marker_set_signature(
        rl_review_marker_sets
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
  write.csv(
    rl_label_approval_template(
      merged_path, review_manifest_path, colnames(obj),
      obj[[cluster_column, drop = TRUE]], cluster_column
    ),
    label_approval_template_path, row.names = FALSE
  )
  message("RL post-QC review outputs are ready. No subcluster labels were applied.")
  message("After reviewing markers, complete the blank approval template and save it to: ",
          label_approval_path)
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
  stop("RL review manifest is invalid; rerun --review-only.")
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
    "rl_focused_review_markers + color_palette.R::rl_markers"
  ) &&
  identical(
    as.character(review$focused_marker_signature),
    focused_marker_set_signature(rl_review_marker_sets)
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
  stop("RL post-QC input or clusters changed after review; rerun --review-only.")
}

# User-approved 2026-10-17 for the 827183-cell, 20-cluster review: marker tables
# and focused DotPlot, with tentative UMAP/marker assignments for sparse 18/19.
# Clusters 5, 6, 10, and 12 retain their reviewed uncertainty; approval selects
# existing names, with cluster 6 explicitly named BG Progenitors. No cell removal.
reviewed_label_maps <- list(
  "0.8" = c(
    "0" = "Differentiating GCs",
    "1" = "Mature UBCs",
    "2" = "Differentiating GCs",
    "3" = "Immature UBCs",
    "4" = "Migrating GCs",
    "5" = "Immature UBCs",
    "6" = "BG Progenitors",
    "7" = "Prolif GCPs",
    "8" = "RL SVZ",
    "9" = "Prolif GCPs",
    "10" = "Differentiating GCs",
    "11" = "Mature GCs",
    "12" = "Maturing PCs",
    "13" = "Migrating GCs",
    "14" = "Prolif GCPs",
    "15" = "Prolif GCPs",
    "16" = "Maturing GCPs",
    "17" = "RL SVZ",
    "18" = "Mature UBCs",
    "19" = "Differentiating GCs"
  )
)
new_labels <- reviewed_label_maps[[subcluster_resolution_label]]
if (is.null(new_labels) || !length(new_labels)) {
  stop(
    "No reviewed RL label map is configured for ", cluster_column, ". ",
    "Run --review-only, inspect the resolution-specific UMAP and top-10/top-50 ",
    "marker tables, then add the approved cluster-to-label mapping to ",
    "reviewed_label_maps."
  )
}

missing_label_mappings <- setdiff(current_clusters, names(new_labels))
unused_label_mappings <- setdiff(names(new_labels), current_clusters)
if (length(missing_label_mappings) || length(unused_label_mappings)) {
  stop(
    "RL label mapping does not exactly match ", cluster_column, ". Missing: ",
    if (length(missing_label_mappings)) paste(missing_label_mappings, collapse = ", ") else "none",
    "; unused: ",
    if (length(unused_label_mappings)) paste(unused_label_mappings, collapse = ", ") else "none"
  )
}
if (anyNA(new_labels) || any(!nzchar(new_labels))) stop("RL label mapping contains missing or blank labels.")
label_approval <- validate_rl_label_approval(
  label_approval_path, merged_path, review_manifest_path, cluster_column,
  labels = new_labels, cell_ids = colnames(obj),
  cluster_ids = obj[[cluster_column, drop = TRUE]]
)
label_approval_md5 <- unname(tools::md5sum(label_approval_path))
reviewed_labels <- unique(unname(new_labels))
labels_missing_from_order <- setdiff(reviewed_labels, rl_subcluster_order)
labels_missing_from_palette <- setdiff(reviewed_labels, names(rl_palette))
if (length(labels_missing_from_order) || length(labels_missing_from_palette)) {
  stop(
    "Reviewed RL labels are not fully registered. Missing from rl_subcluster_order: ",
    if (length(labels_missing_from_order)) {
      paste(labels_missing_from_order, collapse = ", ")
    } else "none",
    "; missing from rl_palette: ",
    if (length(labels_missing_from_palette)) {
      paste(labels_missing_from_palette, collapse = ", ")
    } else "none"
  )
}

obj <- RenameIdents(obj, new_labels)

# Apply and validate the shared ordering before any output object is written.
obj$RL_subcluster <- factor(
  as.character(Idents(obj)),
  levels = rl_subcluster_order
)
if (anyNA(obj$RL_subcluster)) {
  stop("RL_subcluster factor conversion introduced missing labels.")
}
Idents(obj) <- "RL_subcluster"

write.csv(
  data.frame(
    cluster_id = names(new_labels),
    RL_subcluster = unname(new_labels),
    annotation_confidence = ifelse(
      names(new_labels) %in% c("5", "6", "10", "12", "18", "19"),
      "approved_provisional", "marker_supported"
    ),
    source_column = cluster_column,
    review_resolution = subcluster_resolution_label,
    review_basis = label_approval$review_notes[match(names(new_labels), label_approval$cluster_id)],
    label_approval_path = normalizePath(label_approval_path, mustWork = TRUE),
    label_approval_md5 = label_approval_md5,
    source_partition_md5 = unique(label_approval$source_partition_md5),
    input_rds = merged_path,
    output_rds = output_path,
    stringsAsFactors = FALSE
  ),
  label_manifest_path,
  row.names = FALSE
)
saveRDS(obj, output_path, compress = FALSE)
active_rl_order <- rl_subcluster_order[rl_subcluster_order %in% reviewed_labels]
}

DefaultAssay(obj) <- "Xenium"

# Plot only observed approved identities; unused global taxonomy levels create
# ghost heatmap labels and marker facets for populations absent from this object.
obj$RL_subcluster <- factor(as.character(obj$RL_subcluster), levels = active_rl_order)
Idents(obj) <- "RL_subcluster"
active_rl_lineage <- intersect(rl_lineage_subclusters, active_rl_order)
active_rl_other <- intersect(rl_other_subclusters, active_rl_order)
production_rl_other_markers <- rl_other_markers[active_rl_other]
production_rl_markers <- rl_production_marker_panel(active_rl_order)

# Validate the complete production marker panel before writing any refreshed
# table or plot, preventing a missing targeted gene from leaving partial output.
requested_rl_markers <- unique(c(
  unlist(production_rl_markers, use.names = FALSE),
  unlist(rl_compact_lineage_markers, use.names = FALSE),
  unlist(production_rl_other_markers, use.names = FALSE)
))
missing_rl_markers <- setdiff(requested_rl_markers, rownames(obj[["Xenium"]]))
if (length(missing_rl_markers)) {
  stop(
    "Production RL marker genes are absent from the Xenium assay: ",
    paste(missing_rl_markers, collapse = ", ")
  )
}

# Export Merged Cluster QC Statistics
merged_qc_stats <- obj@meta.data %>%
  group_by(RL_subcluster) %>%
  summarise(
    cell_count = n(),
    median_counts = median(nCount_Xenium),
    mean_counts = mean(nCount_Xenium),
    median_features = median(nFeature_Xenium),
    .groups = 'drop'
  )

write.csv(merged_qc_stats,
          file.path(table_path, "Xenium_RL_Subcluster_QC_Summary.csv"),
          row.names = FALSE)

# QC Violin Plot
# 1. Define the QC features you want to plot
qc_features <- c("nCount_Xenium", "nFeature_Xenium")

# 3. Open the PDF device
pdf_path <- file.path(plot_path, "Xenium_RL_Subcluster_QC_Violins.pdf")
grDevices::cairo_pdf(pdf_path, width = 12, height = 8)

# 4. Loop through features and print each to a new page
for (feat in qc_features) {
  message("Plotting: ", feat)
  
  p <- VlnPlot(
    obj,
    features = feat,
    group.by = "RL_subcluster",
    cols = rl_palette,
    pt.size = 0      # Remove points for speed and smaller file size
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(size = 12, face = "bold"), # Explicitly style Y
      legend.position = "none"
    ) +
    ylab(feat) +
    ggtitle(paste("Xenium RL Subcluster QC-", feat))
  
  # Printing the plot inside the loop sends it to the PDF device
  print(p)
}

# 5. Close the device to finalize the file
dev.off()

message("Individual QC plots saved to: ", pdf_path)

# 4. UMAP WITH CONSISTENT COLORS
# Use the same centered white-box label style on both UMAPs.
umap_label_locations <- function(object) {
  locations <- as.data.frame(
    Embeddings(object, reduction = "umap_clean")[, 1:2, drop = FALSE]
  )
  colnames(locations) <- c("x", "y")
  locations$RL_subcluster <- as.character(Idents(object))
  locations %>%
    group_by(RL_subcluster) %>%
    summarise(
      label_x = median(x),
      label_y = median(y),
      .groups = "drop"
    )
}

all_label_locations <- umap_label_locations(obj)
p <- DimPlot(obj,
             reduction = "umap_clean",
             label = FALSE,
             cols = rl_palette) + # Uses your defined palette
  geom_label(
    data = all_label_locations,
    aes(x = label_x, y = label_y, label = RL_subcluster),
    inherit.aes = FALSE,
    size = 4,
    fill = "white",
    color = "black",
    label.padding = grid::unit(0.15, "lines"),
    label.r = grid::unit(0.12, "lines"),
    linewidth = 0.25
  ) +
  ggtitle("Xenium Merged RL Subclusters") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_Subcluster_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p)
grDevices::dev.off()

# 4. UMAP focused on the reviewed RL/UBC/granule lineages. The legacy stable
# filename is retained even though the reviewed off-lineage set is now broader.
off_lineage_clusters <- c(
  "BG", "BG Progenitors", "Maturing PCs", "eCN",
  "GABA Progenitors", "Microglia", "Vascular"
)
clusters_to_remove <- intersect(
  off_lineage_clusters,
  unique(as.character(Idents(obj)))
)

# Create a subset of the object
# ! means 'not', %in% checks if the identity is in your removal list
obj_subset <- subset(obj, idents = clusters_to_remove, invert = TRUE)

filtered_label_locations <- umap_label_locations(obj_subset)

# Now run the plotting code on 'obj_subset' instead of 'obj'.
p2 <- DimPlot(obj_subset,
              reduction = "umap_clean",
              label = FALSE,
              cols = rl_palette) +
  geom_label(
    data = filtered_label_locations,
    aes(x = label_x, y = label_y, label = RL_subcluster),
    inherit.aes = FALSE,
    size = 4,
    fill = "white",
    color = "black",
    label.padding = grid::unit(0.15, "lines"),
    label.r = grid::unit(0.12, "lines"),
    linewidth = 0.25
  ) +
  ggtitle("Xenium Merged RL Subclusters (Filtered)") +
  theme(legend.text = element_text(size = 8))

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_Subcluster_rmBG,Ep,PC_UMAP.tif")),
  width = 12,
  height = 9,
  units = "in",
  res = 600
)
print(p2)
grDevices::dev.off()

rm(obj_subset)

# 5. DOTPLOT WITH REVERSED ORDER
# Create a temporary reversed factor for the Y-axis
# (Reverse the levels of the factor to get VZP at the top)
obj$RL_subcluster_rev <- factor(as.character(Idents(obj)), levels = rev(active_rl_order))
Idents(obj) <- "RL_subcluster_rev"

# 3. Create the Plot
# Adding 'assay = "Xenium"' is the safest way to avoid multi-assay conflicts
p1 <- workflow_dotplot(obj,
              features = production_rl_markers,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 12, face = "italic"),
    axis.text.y = element_text(size = 12, face = "plain")
  ) +
  ggtitle("Xenium Merged RL Subcluster Markers")

dotplot_dimensions <- workflow_dotplot_dimensions(p1)
dotplot_width <- dotplot_dimensions[["width"]]
dotplot_height <- dotplot_dimensions[["height"]]
dotplot_grob <- workflow_dotplot_grob(p1)

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_SubclusterMarker_DotPlot.tif")),
  width = dotplot_width,
  height = dotplot_height,
  units = "in",
  res = 600
)
grid::grid.newpage()
grid::grid.draw(dotplot_grob)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "Xenium_RL_SubclusterMarker_DotPlot.pdf"),
  dotplot_grob, device = grDevices::cairo_pdf, width = dotplot_width, height = dotplot_height
)

active_rl_lineage <- intersect(rl_lineage_subclusters, active_rl_order)
active_rl_other <- intersect(rl_other_subclusters, active_rl_order)
if (!length(active_rl_lineage) || !length(active_rl_other)) {
  stop("Both RL lineage and Other identities are required for split DotPlots.")
}

build_rl_split_dotplot <- function(marker_sets, identities, title) {
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

save_rl_split_dotplot <- function(plot, filename_stem) {
  is_other <- identical(filename_stem, "Xenium_RL_Other_SubclusterMarker_DotPlot")
  if (is_other) plot <- plot + theme(legend.box = "vertical")
  dimensions <- workflow_dotplot_dimensions(
    plot, minimum_width = if (is_other) 6 else 9,
    feature_spacing = if (is_other) 0.22 else 0.28
  )
  width <- dimensions[["width"]]
  height <- dimensions[["height"]] + if (is_other) 0.5 else 0
  display <- workflow_dotplot_grob(plot)
  Cairo::CairoTIFF(
    filename = file.path(plot_path, paste0(filename_stem, ".tif")),
    width = width,
    height = height,
    units = "in",
    res = 600
  )
  grid::grid.newpage()
  grid::grid.draw(display)
  grDevices::dev.off()
  ggplot2::ggsave(
    file.path(plot_path, paste0(filename_stem, ".pdf")),
    display, device = grDevices::cairo_pdf, width = width, height = height
  )
}

p_rl_lineage <- build_rl_split_dotplot(
  rl_compact_lineage_markers,
  active_rl_lineage,
  "RL Lineage Subclusters: Lineage and Granule-State Markers"
)
save_rl_split_dotplot(
  p_rl_lineage,
  "Xenium_RL_Lineage_SubclusterMarker_DotPlot"
)

p_rl_other <- build_rl_split_dotplot(
  production_rl_other_markers,
  active_rl_other,
  "RL Other Subclusters: Markers for Observed Identities"
)
save_rl_split_dotplot(
  p_rl_other,
  "Xenium_RL_Other_SubclusterMarker_DotPlot"
)

# ----------------------------------------------------------------
# 6. GENERATE TOP 5 MARKERS HEATMAP
# ----------------------------------------------------------------
check_mem("STARTING HEATMAP GENERATION")
Idents(obj) <- "RL_subcluster"

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

write.csv(heatmap_markers,
          file.path(table_path, "Xenium_RL_Subcluster_Markers.csv"),
          row.names = FALSE)

# 2. Extract top 5 per cluster
top5_markers <- heatmap_markers %>%
  mutate(cluster = factor(as.character(cluster), levels = active_rl_order)) %>%
  filter(!is.na(cluster)) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

# 3. Scale the data for the specific markers only
# (Necessary for DoHeatmap to show relative expression)
# Just run the plot on a temporary scaled version:
temp_obj <- ScaleData(obj, features = top5_markers, verbose = FALSE)
temp_obj$RL_subcluster <- factor(as.character(temp_obj$RL_subcluster), levels = active_rl_order)
Idents(temp_obj) <- "RL_subcluster"
heatmap_obj <- subset(temp_obj, downsample = 100, seed = config$runtime$random_seed)
heatmap_obj$RL_subcluster <- droplevels(heatmap_obj$RL_subcluster)
Idents(heatmap_obj) <- "RL_subcluster"

# 4. Create Heatmap
# We downsample the plot to 100 cells per group so the labels are readable
# 1. Prepare the plot without the default lines
p4 <- DoHeatmap(
  heatmap_obj,
  features = top5_markers,
  group.by = "RL_subcluster",
  group.colors = unname(rl_palette[levels(heatmap_obj$RL_subcluster)]),
  size = 4,
  angle = 45,
  draw.lines = TRUE,
  raster = FALSE
) +
  scale_fill_viridis_c(option = "viridis", name = "Z-Score", na.value = "white") +
  guides(color = "none") +
  theme(
    axis.text.y = element_text(size = 6, face = "italic"),
    plot.margin = margin(t = 75, r = 12, b = 12, l = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Top 5 Markers per RL Subcluster")

# 4. SAVE AS PDF
ggplot2::ggsave(
  file.path(plot_path, paste0("Xenium_RL_Subcluster_Top5_Heatmap.pdf")),
  p4, width = 14, height = 12, device = grDevices::cairo_pdf
)

# 5. Save the Heatmap
Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_Subcluster_Top5_Heatmap.tif")),
  width = 14,
  height = 12,
  units = "in",
  res = 600
)
print(p4)
grDevices::dev.off()

check_mem("HEATMAP COMPLETE")

plan("sequential")

# Preserve a reviewed visualization choice during subsequent plot-only refreshes.
# Full label application creates a new source object and needs a new comparison.
selected_umap_receipt <- file.path(table_path, "Xenium_RL_Selected_UMAP_manifest.csv")
if (plots_only && file.exists(selected_umap_receipt)) {
  source(here("scripts", "R", "rl_umap_selected.R"))
  render_selected_rl_umap(
    config, output_path,
    file.path(output_root, "xenium", "rl", "06c_umap_comparison"),
    overwrite = overwrite, object = obj
  )
}

check_mem(if (plots_only) "PLOT REFRESH COMPLETE" else "PIPELINE COMPLETE - FILE SAVED")
