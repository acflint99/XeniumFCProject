#!/usr/bin/env Rscript

# Review post-QC RL clusters, apply reviewed labels, or refresh production plots.
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
review_manifest_path <- file.path(
  table_path,
  paste0(
    "Xenium_RL_subcluster_review_manifest_Res",
    subcluster_resolution_label,
    ".csv"
  )
)
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
review_outputs <- c(
  review_top10_marker_path, review_top50_marker_path,
  review_manifest_path, review_umap_path,
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
    "RL post-QC subcluster review"
  } else if (apply_reviewed_labels) {
    "RL reviewed-label application"
  } else {
    "RL plot refresh"
  }
  compact_dry_run(
    mode_label,
    inputs = if (review_only) {
      merged_path
    } else if (apply_reviewed_labels) {
      c(merged_path, review_outputs)
    } else {
      output_path
    },
    outputs = if (review_only) {
      review_outputs
    } else if (apply_reviewed_labels) {
      expected_outputs
    } else {
      plot_refresh_outputs
    },
    checks = if (apply_reviewed_labels) {
      c(
        review_outputs_nonempty = all(vapply(review_outputs, file_is_nonempty, logical(1))),
        review_manifest_matches_input = review_manifest_matches_input_file()
      )
    } else logical()
  )
  quit(save = "no", status = 0L)
}

if (plots_only) {
  if (!file.exists(output_path)) stop("Annotated RL input not found: ", output_path)
} else if (!file.exists(merged_path)) {
  stop("RL post-QC input not found: ", merged_path)
}
if (apply_reviewed_labels) {
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
rl_review_marker_sets <- if (!plots_only) {
  add_palette_markers_to_review(
    focused_marker_sets = rl_focused_review_markers,
    palette_marker_sets = rl_markers,
    palette_object_name = "rl_markers"
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
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident = 1000,
      random.seed = config$runtime$random_seed
    )
  )
  if (!nrow(all_markers)) stop("No RL post-QC markers were returned for review.")
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  top50_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 50, order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup()
  write.csv(top10_markers, review_top10_marker_path, row.names = FALSE)
  write.csv(top50_markers, review_top50_marker_path, row.names = FALSE)

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
    pdf_path = review_focused_dotplot_pdf_path,
    width = 21,
    height = 10
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
      marker_table_depths = "10|50",
      focused_marker_panel = paste(
        "rl_focused_review_markers",
        "color_palette.R::rl_markers",
        sep = " + "
      ),
      focused_marker_signature = focused_marker_set_signature(
        rl_review_marker_sets
      ),
      stringsAsFactors = FALSE
    ),
    review_manifest_path,
    row.names = FALSE
  )
  message("RL post-QC review outputs are ready. No subcluster labels were applied.")
  quit(save = "no", status = 0L)
}

review <- read.csv(review_manifest_path, stringsAsFactors = FALSE)
required_review_columns <- c(
  "input_rds", "input_size", "input_mtime", "cells", "cluster_column", "clusters",
  "random_seed", "marker_table_depths", "focused_marker_panel",
  "focused_marker_signature"
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
  identical(as.character(review$marker_table_depths), "10|50") &&
  identical(
    as.character(review$focused_marker_panel),
    "rl_focused_review_markers + color_palette.R::rl_markers"
  ) &&
  identical(
    as.character(review$focused_marker_signature),
    focused_marker_set_signature(rl_review_marker_sets)
  )
if (!review_is_current) {
  stop("RL post-QC input or clusters changed after review; rerun --review-only.")
}

# No earlier resolution-0.5 map is reused here. The resolution-0.8 map below
# was reviewed against the current 30-dimension post-QC UMAP, top-50 marker
# table, focused RL lineage DotPlot, and RL-VZ diagnostic DotPlot. Cell-cycle
# subdivisions retain their biological lineage label rather than becoming
# separate terminal cell types.
reviewed_label_maps <- list(
  "0.8" = c(
    "0" = "Mature UBCs",
    "1" = "Immature UBCs",
    "2" = "Maturing GCPs",
    "3" = "Migrating GCs",
    "4" = "Differentiating GCs",
    "5" = "Migrating GCs",
    "6" = "BG",
    "7" = "Differentiating GCs",
    "8" = "RL SVZ",
    "9" = "Prolif GCPs",
    "10" = "Migrating GCs",
    "11" = "Mature GCs",
    "12" = "RL SVZ",
    "13" = "Maturing PCs",
    "14" = "Differentiating GCs",
    "15" = "Prolif GCPs",
    "16" = "Maturing GCPs",
    "17" = "eCN",
    "18" = "GABA Progenitors",
    "19" = "Microglia",
    "20" = "Endothelial",
    "21" = "Maturing PCs",
    "22" = "Maturing PCs"
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
    source_column = cluster_column,
    review_resolution = subcluster_resolution_label,
    review_basis = paste(
      "current post-QC UMAP, top-50 markers, focused lineage DotPlot,",
      "and RL-VZ diagnostic DotPlot"
    ),
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

# Validate the complete production marker panel before writing any refreshed
# table or plot, preventing a missing targeted gene from leaving partial output.
requested_rl_markers <- unique(c(
  unlist(rl_markers, use.names = FALSE),
  unlist(rl_compact_lineage_markers, use.names = FALSE),
  unlist(rl_other_markers, use.names = FALSE)
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
  "BG", "Maturing PCs", "eCN",
  "GABA Progenitors", "Microglia", "Endothelial"
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
p1 <- DotPlot(obj,
              features = rl_markers,
              assay = "Xenium",
              cols = c("lightgrey", "red"),
              dot.scale = 6,
              cluster.idents = FALSE) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 8, face = "italic"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) +
  ggtitle("Xenium Merged RL Subcluster Markers")

dotplot_width <- max(
  16,
  0.48 * length(unique(unlist(rl_markers, use.names = FALSE)))
)

Cairo::CairoTIFF(
  filename = file.path(plot_path, paste0("Xenium_RL_SubclusterMarker_DotPlot.tif")),
  width = dotplot_width,
  height = 8,
  units = "in",
  res = 600
)
print(p1)
grDevices::dev.off()
ggplot2::ggsave(
  file.path(plot_path, "Xenium_RL_SubclusterMarker_DotPlot.pdf"),
  p1, device = grDevices::cairo_pdf, width = dotplot_width, height = 8
)

active_rl_lineage <- intersect(rl_lineage_subclusters, active_rl_order)
active_rl_other <- intersect(rl_other_subclusters, active_rl_order)
if (!length(active_rl_lineage) || !length(active_rl_other)) {
  stop("Both RL lineage and Other identities are required for split DotPlots.")
}

build_rl_split_dotplot <- function(marker_sets, identities, title) {
  DotPlot(
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
      axis.text.x = element_text(size = 9, face = "italic"),
      axis.text.y = element_text(size = 10, face = "bold")
    ) +
    ggtitle(title)
}

save_rl_split_dotplot <- function(plot, filename_stem, width, height) {
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

p_rl_lineage <- build_rl_split_dotplot(
  rl_compact_lineage_markers,
  active_rl_lineage,
  "RL Lineage Subclusters: Compact Markers"
)
save_rl_split_dotplot(
  p_rl_lineage,
  "Xenium_RL_Lineage_SubclusterMarker_DotPlot",
  width = max(14, 0.45 * length(unlist(rl_compact_lineage_markers))),
  height = 6
)

p_rl_other <- build_rl_split_dotplot(
  rl_other_markers,
  active_rl_other,
  "RL Other Subclusters: Compact Markers"
)
save_rl_split_dotplot(
  p_rl_other,
  "Xenium_RL_Other_SubclusterMarker_DotPlot",
  width = max(12, 0.5 * length(unlist(rl_other_markers))),
  height = 5
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

# 4. Create Heatmap
# We downsample the plot to 100 cells per group so the labels are readable
# 1. Prepare the plot without the default lines
p4 <- DoHeatmap(
  subset(temp_obj, downsample = 100),
  features = top5_markers,
  group.by = "RL_subcluster",
  group.colors = rl_palette,
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

check_mem(if (plots_only) "PLOT REFRESH COMPLETE" else "PIPELINE COMPLETE - FILE SAVED")
