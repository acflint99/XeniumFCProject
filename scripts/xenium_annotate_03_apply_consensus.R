#!/usr/bin/env Rscript

# Attach the cluster-level consensus labels produced by
# xenium_annotate_02_build_consensus.R to one individual Xenium object. The Aldinger object
# is used as the expression/spatial container; its original label columns are
# retained unchanged.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(readxl)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "color_palette.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
pilot_manifest <- load_resolution2_pilot_manifest(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--dry-run", "--overwrite", "--dotplot-only", "--plots-only", "--pilot-res2",
  "--pilot-res3", "--pilot-res4", "--pilot-res5", "--all-samples-res4",
  "--all-samples-res5",
  "--selected-sample", "--weighted-2of3", "--list"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options) > 0L) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

pilot_res2 <- "--pilot-res2" %in% args
pilot_res3 <- "--pilot-res3" %in% args
pilot_res4 <- "--pilot-res4" %in% args
pilot_res5 <- "--pilot-res5" %in% args
all_samples_res4 <- "--all-samples-res4" %in% args
all_samples_res5 <- "--all-samples-res5" %in% args
all_samples_mode <- all_samples_res4 || all_samples_res5
selected_sample_mode <- "--selected-sample" %in% args
pilot_flags <- c(pilot_res2, pilot_res3, pilot_res4, pilot_res5)
if (sum(pilot_flags) > 1L) {
  stop(
    "Choose only one of --pilot-res2, --pilot-res3, --pilot-res4, or --pilot-res5."
  )
}
pilot_mode <- any(pilot_flags)
if (all_samples_res4 && all_samples_res5) {
  stop("Choose only one of --all-samples-res4 or --all-samples-res5.")
}
if (all_samples_mode && pilot_mode) {
  stop("An all-samples mode cannot be combined with a pilot-resolution option.")
}
if (selected_sample_mode && (!pilot_res5 || all_samples_mode)) {
  stop("--selected-sample currently requires --pilot-res5.")
}
selected_sample_id <- if (selected_sample_mode) {
  trimws(Sys.getenv("SELECTED_SAMPLE_ID", unset = ""))
} else {
  ""
}
if (selected_sample_mode) {
  selected_matches <- which(sample_manifest$sample_id == selected_sample_id)
  if (!nzchar(selected_sample_id) || length(selected_matches) != 1L) {
    stop(
      "SELECTED_SAMPLE_ID must match exactly one config/samples.csv row; found ",
      length(selected_matches), " match(es) for '", selected_sample_id, "'."
    )
  }
}
resolution_mode <- pilot_mode || all_samples_mode
pilot_resolution_tag <- if (all_samples_res5) "5.0" else if (all_samples_res4) "4.0" else if (pilot_res5) "5.0" else if (pilot_res4) "4.0" else if (pilot_res3) "3.0" else "2.0"
pilot_stage <- if (selected_sample_mode) {
  "03h_resolution5_selected_sample"
} else if (all_samples_res4) {
  "03g_resolution4_all_samples"
} else if (all_samples_res5) {
  "03i_resolution5_all_samples"
} else if (pilot_res5) {
  "03e_resolution5_pilot"
} else if (pilot_res4) {
  "03d_resolution4_pilot"
} else if (pilot_res3) {
  "03c_resolution3_pilot"
} else {
  "03b_resolution2_pilot"
}
pilot_cluster_column <- paste0("whole_tissue_cluster_res", pilot_resolution_tag)
pilot_graph_column <- paste0(
  "Xenium_snn_res.",
  if (pilot_res5 || all_samples_res5) "5" else if (pilot_res4 || all_samples_res4) "4" else if (pilot_res3) "3" else "2"
)
facet_point_size <- if (all_samples_mode || selected_sample_mode) {
  0.01
} else if (pilot_res3 || pilot_res4 || pilot_res5) {
  0.03
} else {
  0.1
}
weighted_2of3 <- "--weighted-2of3" %in% args
consensus_method <- if (weighted_2of3) "weighted_2of3" else "legacy_six_vote"
samples <- if (selected_sample_mode) {
  sample_manifest[selected_matches, , drop = FALSE]
} else if (pilot_mode) {
  pilot_manifest
} else {
  sample_manifest
}
mode_description <- if (selected_sample_mode) {
  paste0("resolution-5.0 selected-sample pilot: ", selected_sample_id)
} else if (all_samples_res4) {
  "resolution-4.0 all-sample analysis"
} else if (all_samples_res5) {
  "resolution-5.0 all-sample analysis"
} else if (pilot_mode) {
  paste0("resolution-", pilot_resolution_tag, " pilot")
} else {
  "production"
}
task_map <- data.frame(
  task_id = seq_len(nrow(samples)),
  sample_id = samples$sample_id,
  stringsAsFactors = FALSE
)

if ("--list" %in% args) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept TASK_ID.")
  cat(
    "Mode:", mode_description, "\n"
  )
  cat("Consensus method:", consensus_method, "\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
dotplot_only <- "--dotplot-only" %in% args
plots_only <- "--plots-only" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
if (dotplot_only && plots_only) {
  stop("Choose only one of --dotplot-only or --plots-only.")
}
if ((dotplot_only || plots_only) && overwrite) {
  stop(
    "Plot-only modes already authorize replacing their plot files; ",
    "do not combine them with --overwrite."
  )
}

task_args <- args[!args %in% valid_options]
if (length(task_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_03_apply_consensus.R ",
    "[--pilot-res2|--pilot-res3|--pilot-res4|--pilot-res5|--all-samples-res4|--all-samples-res5] ",
    "[--selected-sample] ",
    "[--weighted-2of3] ",
    "[--dry-run|--overwrite|--dotplot-only|--plots-only] [TASK_ID]"
  )
}

task_value <- if (length(task_args) == 1L) {
  task_args[[1]]
} else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(samples)) {
  stop("TASK_ID must be between 1 and ", nrow(samples), ". Use --list to inspect the mapping.")
}

sample_name <- samples$sample_id[[task_id]]
output_root <- here(config$project$outputs_dir)
annotation_root <- if (all_samples_res4) {
  file.path(output_root, "xenium", "annotation", "resolution4_all_samples")
} else if (all_samples_res5) {
  file.path(output_root, "xenium", "annotation", "resolution5_all_samples")
} else if (pilot_mode) {
  file.path(
    output_root, "xenium", "preprocess", pilot_stage, "annotation"
  )
} else {
  file.path(output_root, "xenium", "annotation")
}
input_path <- file.path(
  annotation_root, "01_label_transfer", "aldinger", "rds",
  paste0(sample_name, "_Aldinger_annotated.rds")
)
comparison_path <- file.path(
  annotation_root,
  if (weighted_2of3) "02_consensus_weighted_2of3" else "02_consensus",
  "tables",
  paste0(sample_name, "_comparison_merged.csv")
)
metadata_path <- resolve_config_path(config$manifests$sample_metadata, config)
output_dir <- file.path(
  annotation_root,
  if (weighted_2of3) "03_consensus_labels_weighted_2of3" else "03_consensus_labels",
  "rds"
)
plot_dir <- file.path(
  annotation_root,
  if (weighted_2of3) "03_consensus_labels_weighted_2of3" else "03_consensus_labels",
  "plots"
)
output_path <- file.path(output_dir, paste0(sample_name, "_Consensus_annotated.rds"))
plot_paths <- file.path(
  plot_dir,
  paste0(
    sample_name,
    c("_Consensus_UMAP.tif", "_Consensus_GlobalSpatial.tif", "_Consensus_FacetSpatial.tif")
  )
)
known_plot_paths <- file.path(
  plot_dir,
  paste0(
    sample_name,
    c(
      "_Consensus_KnownOnly_UMAP.tif",
      "_Consensus_KnownOnly_GlobalSpatial.tif",
      "_Consensus_KnownOnly_FacetSpatial.tif"
    )
  )
)
dotplot_paths <- file.path(
  plot_dir,
  paste0(sample_name, "_Consensus_Marker_DotPlot", c(".tif", ".pdf"))
)
known_dotplot_paths <- file.path(
  plot_dir,
  paste0(sample_name, "_Consensus_KnownOnly_Marker_DotPlot", c(".tif", ".pdf"))
)
marker_table_dir <- file.path(
  annotation_root,
  if (weighted_2of3) "03_consensus_labels_weighted_2of3" else "03_consensus_labels",
  "tables", "unknown_markers"
)
marker_workbook_path <- file.path(
  marker_table_dir, paste0(sample_name, "_Unknown_cluster_top20_markers.xlsx")
)
expected_outputs <- c(
  output_path, plot_paths, known_plot_paths, dotplot_paths, known_dotplot_paths
)

save_consensus_dotplot <- function(
    obj, sample_name, dotplot_paths, title_suffix = "with Unknown clusters") {
  if (!"consensus_label" %in% colnames(obj[[]])) {
    stop("Consensus object lacks the 'consensus_label' metadata column: ", sample_name)
  }

  consensus_values <- as.character(obj$consensus_label)
  if (anyNA(consensus_values) || any(!nzchar(consensus_values))) {
    stop("Consensus object contains blank or missing consensus labels: ", sample_name)
  }
  consensus_levels <- order_consensus_levels(consensus_values, celltype_order)
  obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
  Idents(obj) <- "consensus_label"

  existing_markers <- lapply(markers, function(features) intersect(features, rownames(obj)))
  existing_markers <- existing_markers[lengths(existing_markers) > 0L]
  if (!length(existing_markers)) {
    stop("None of the configured broad-cell markers are present in ", sample_name, ".")
  }

  p_dot <- DotPlot(
    obj,
    features = existing_markers,
    assay = "Xenium",
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot <- standardize_broad_dotplot(p_dot, consensus_levels) +
    ggtitle(paste(sample_name, "Consensus marker expression —", title_suffix))

  Cairo::CairoTIFF(dotplot_paths[[1]], width = 10, height = 6, units = "in", res = 600)
  print(p_dot)
  grDevices::dev.off()
  ggplot2::ggsave(
    filename = dotplot_paths[[2]],
    plot = p_dot,
    device = grDevices::cairo_pdf,
    width = 10,
    height = 6,
    units = "in"
  )
}

save_unknown_marker_workbook <- function(
    obj, sample_name, marker_workbook_path, overwrite = FALSE) {
  labels <- trimws(as.character(obj$consensus_label))
  unknown_levels <- order_consensus_levels(labels, celltype_order)
  unknown_levels <- unknown_levels[is_unknown_consensus_label(unknown_levels)]
  if (!length(unknown_levels)) {
    message(
      "Skipping Unknown marker workbook for ", sample_name,
      ": no Unknown clusters were found."
    )
    return(invisible(FALSE))
  }
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop(
      "The 'writexl' package is required for Unknown-cluster marker workbooks. ",
      "Restore it through the project renv environment before rerunning this stage."
    )
  }
  if (file.exists(marker_workbook_path) && !overwrite) {
    stop(
      "Refusing to overwrite the existing Unknown marker workbook: ",
      marker_workbook_path, "\nUse --overwrite only after reviewing it."
    )
  }

  obj$consensus_label <- factor(
    labels, levels = order_consensus_levels(labels, celltype_order)
  )
  Idents(obj) <- "consensus_label"
  DefaultAssay(obj) <- "Xenium"
  sheets <- lapply(unknown_levels, function(unknown_label) {
    cluster_ids <- unique(as.character(obj$seurat_clusters[labels == unknown_label]))
    if (length(cluster_ids) != 1L || is.na(cluster_ids) || !nzchar(cluster_ids)) {
      stop(
        sample_name, " ", unknown_label,
        " must map to exactly one nonblank seurat_clusters value; found: ",
        paste(cluster_ids, collapse = ", ")
      )
    }
    markers_table <- Seurat::FindMarkers(
      object = obj,
      ident.1 = unknown_label,
      ident.2 = NULL,
      assay = "Xenium",
      test.use = "wilcox",
      only.pos = TRUE,
      min.pct = 0.01,
      logfc.threshold = 0.1,
      verbose = FALSE
    )
    markers_table$gene <- rownames(markers_table)
    fold_column <- intersect(c("avg_log2FC", "avg_logFC"), names(markers_table))
    if (!length(fold_column)) {
      stop("FindMarkers did not return an average log-fold-change column for ", unknown_label, ".")
    }
    adjusted_p <- if ("p_val_adj" %in% names(markers_table)) {
      markers_table$p_val_adj
    } else {
      rep(NA_real_, nrow(markers_table))
    }
    markers_table <- markers_table[
      order(-markers_table[[fold_column[[1]]]], adjusted_p, markers_table$gene),
      ,
      drop = FALSE
    ]
    markers_table <- head(markers_table, 20L)
    markers_table <- data.frame(
      rank = seq_len(nrow(markers_table)),
      gene = markers_table$gene,
      unknown_label = unknown_label,
      seurat_cluster = cluster_ids,
      contrast = paste0(unknown_label, " vs all other cells"),
      markers_table[setdiff(names(markers_table), "gene")],
      row.names = NULL,
      check.names = FALSE
    )
    markers_table
  })
  names(sheets) <- unknown_levels

  dir.create(dirname(marker_workbook_path), recursive = TRUE, showWarnings = FALSE)
  partial_path <- paste0(marker_workbook_path, ".", Sys.getpid(), ".part.xlsx")
  on.exit(if (file.exists(partial_path)) unlink(partial_path), add = TRUE)
  writexl::write_xlsx(sheets, partial_path)
  backup_path <- NULL
  if (file.exists(marker_workbook_path)) {
    backup_path <- paste0(marker_workbook_path, ".backup_", Sys.getpid())
    if (!file.rename(marker_workbook_path, backup_path)) {
      stop("Could not preserve the existing Unknown marker workbook: ", marker_workbook_path)
    }
  }
  if (!file.rename(partial_path, marker_workbook_path)) {
    if (!is.null(backup_path) && file.exists(backup_path)) {
      file.rename(backup_path, marker_workbook_path)
    }
    stop("Could not move the completed Unknown marker workbook into place: ", marker_workbook_path)
  }
  if (!is.null(backup_path) && file.exists(backup_path)) unlink(backup_path)
  message("Saved Unknown-cluster top-20 marker workbook: ", marker_workbook_path)
  invisible(TRUE)
}

save_consensus_plots <- function(
    obj, sample_name, plot_paths, known_plot_paths, dotplot_paths,
    known_dotplot_paths, facet_point_size,
    resolution_mode, pilot_cluster_column, pilot_graph_column,
    pilot_resolution_tag, weighted_2of3) {
  required_metadata <- c("consensus_label")
  if (resolution_mode) {
    required_metadata <- c(
      required_metadata, "seurat_clusters", "whole_tissue_cluster_res1.5",
      pilot_cluster_column, pilot_graph_column
    )
  }
  if (weighted_2of3) required_metadata <- c(required_metadata, "consensus_method")
  missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
  if (length(missing_metadata)) {
    stop(
      "Consensus object lacks required plotting metadata for ", sample_name,
      ": ", paste(missing_metadata, collapse = ", ")
    )
  }
  if (anyDuplicated(Cells(obj))) {
    stop("Consensus object contains duplicate cell IDs: ", sample_name)
  }
  if (!"umap" %in% Reductions(obj)) stop("Consensus object lacks UMAP: ", sample_name)
  if (!"Xenium" %in% Assays(obj)) stop("Consensus object lacks Xenium assay: ", sample_name)
  if (resolution_mode && !identical(
    as.character(obj$seurat_clusters),
    as.character(obj[[pilot_cluster_column]][, 1])
  )) {
    stop(
      "Resolution-", pilot_resolution_tag,
      " consensus object seurat_clusters does not match ", pilot_cluster_column,
      ": ", sample_name
    )
  }
  if (weighted_2of3) {
    object_methods <- unique(trimws(as.character(obj$consensus_method)))
    if (!identical(object_methods, "weighted_2of3")) {
      stop(
        "Consensus object does not contain weighted_2of3 provenance for ",
        sample_name, "; found: ", paste(object_methods, collapse = ", ")
      )
    }
  }

  consensus_values <- trimws(as.character(obj$consensus_label))
  if (anyNA(consensus_values) || any(!nzchar(consensus_values))) {
    stop("Consensus object contains blank or missing labels: ", sample_name)
  }
  consensus_levels <- order_consensus_levels(consensus_values, celltype_order)
  validate_palette(consensus_levels)
  obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
  Idents(obj) <- "consensus_label"
  save_plot_variant <- function(plot_obj, variant_paths, variant_dotplot_paths, title_suffix) {
    variant_values <- trimws(as.character(plot_obj$consensus_label))
    variant_levels <- order_consensus_levels(variant_values, celltype_order)
    plot_obj$consensus_label <- factor(variant_values, levels = variant_levels)
    Idents(plot_obj) <- "consensus_label"
    plot_colors <- resolve_cluster_colors(variant_levels)

    p_umap <- DimPlot(
      plot_obj, reduction = "umap", group.by = "consensus_label",
      label = TRUE, cols = plot_colors
    ) + ggtitle(paste(sample_name, "Consensus labels —", title_suffix))
    Cairo::CairoTIFF(variant_paths[[1]], width = 8, height = 6, units = "in", res = 600)
    print(p_umap)
    grDevices::dev.off()

    p_spatial <- ImageDimPlot(
      plot_obj, group.by = "consensus_label", size = 0.75, cols = plot_colors
    ) + ggtitle(paste(sample_name, "Consensus labels —", title_suffix))
    Cairo::CairoTIFF(
      variant_paths[[2]], width = 8, height = 8, units = "in", res = 600, bg = "black"
    )
    print(p_spatial)
    grDevices::dev.off()

    coords <- Seurat::GetTissueCoordinates(plot_obj, image = "fov", full = FALSE)
    required_coordinate_columns <- c("x", "y", "cell")
    missing_coordinate_columns <- setdiff(required_coordinate_columns, colnames(coords))
    if (length(missing_coordinate_columns)) {
      stop(
        "Spatial coordinates lack required columns for ", sample_name, ": ",
        paste(missing_coordinate_columns, collapse = ", ")
      )
    }
    coordinate_cell_ids <- as.character(coords$cell)
    if (anyDuplicated(coordinate_cell_ids)) {
      stop("Spatial coordinates contain duplicate cell identifiers: ", sample_name)
    }
    coordinate_cells <- match(coordinate_cell_ids, colnames(plot_obj))
    if (anyNA(coordinate_cells)) {
      missing_coordinate_cells <- unique(coordinate_cell_ids[is.na(coordinate_cells)])
      stop(
        "Spatial coordinates contain cell identifiers absent from ", sample_name,
        ": ", paste(head(missing_coordinate_cells, 10L), collapse = ", ")
      )
    }
    plot_data <- cbind(
      coords,
      consensus_label = factor(
        variant_values[coordinate_cells], levels = variant_levels
      )
    )
    if (anyNA(plot_data$consensus_label)) {
      stop("Failed to join consensus labels to spatial coordinates: ", sample_name)
    }
    p_facet <- ggplot(plot_data, aes(x = y, y = x, color = consensus_label)) +
      geom_point(size = facet_point_size) +
      facet_wrap(~consensus_label) +
      scale_color_manual(values = plot_colors, drop = FALSE) +
      coord_fixed() +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        legend.position = "none",
        strip.text = element_text(color = "white")
      )
    Cairo::CairoTIFF(
      variant_paths[[3]], width = 12, height = 8, units = "in", res = 600, bg = "black"
    )
    print(p_facet)
    grDevices::dev.off()

    save_consensus_dotplot(
      plot_obj, sample_name, variant_dotplot_paths, title_suffix = title_suffix
    )
  }

  save_plot_variant(obj, plot_paths, dotplot_paths, "with Unknown clusters")
  known_cells <- Cells(obj)[!is_unknown_consensus_label(obj$consensus_label)]
  if (!length(known_cells)) stop("No known cells remain for known-only plots: ", sample_name)
  known_obj <- subset(obj, cells = known_cells)
  save_plot_variant(
    known_obj, known_plot_paths, known_dotplot_paths, "known cells only"
  )
}

if (dotplot_only) {
  if (dry_run) {
    compact_dry_run(
      paste0(
        "Consensus DotPlot task ", task_id, "/", nrow(samples),
        " [", sample_name, "]"
      ),
      inputs = output_path,
      outputs = c(dotplot_paths, known_dotplot_paths)
    )
    quit(save = "no", status = 0L)
  }

  if (!file.exists(output_path)) stop("Consensus object not found: ", output_path)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(output_path)
  save_consensus_dotplot(obj, sample_name, dotplot_paths)
  known_cells <- Cells(obj)[!is_unknown_consensus_label(obj$consensus_label)]
  if (!length(known_cells)) stop("No known cells remain for known-only DotPlot: ", sample_name)
  known_obj <- subset(obj, cells = known_cells)
  save_consensus_dotplot(
    known_obj, sample_name, known_dotplot_paths, title_suffix = "known cells only"
  )
  message("Regenerated all-cell and known-only consensus DotPlots for ", sample_name, ".")
  quit(save = "no", status = 0L)
}

if (plots_only) {
  if (dry_run) {
    compact_dry_run(
      paste0(
        "Consensus plots task ", task_id, "/", nrow(samples),
        " [", sample_name, "]"
      ),
      inputs = output_path,
      outputs = c(plot_paths, known_plot_paths, dotplot_paths, known_dotplot_paths)
    )
    quit(save = "no", status = if (file.exists(output_path)) 0L else 1L)
  }

  if (!file.exists(output_path)) stop("Consensus object not found: ", output_path)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- readRDS(output_path)
  save_consensus_plots(
    obj = obj,
    sample_name = sample_name,
    plot_paths = plot_paths,
    known_plot_paths = known_plot_paths,
    dotplot_paths = dotplot_paths,
    known_dotplot_paths = known_dotplot_paths,
    facet_point_size = facet_point_size,
    resolution_mode = resolution_mode,
    pilot_cluster_column = pilot_cluster_column,
    pilot_graph_column = pilot_graph_column,
    pilot_resolution_tag = pilot_resolution_tag,
    weighted_2of3 = weighted_2of3
  )
  message("Regenerated all consensus plots without rewriting RDS: ", sample_name)
  quit(save = "no", status = 0L)
}

if (dry_run) {
  input_checks <- data.frame(
    input_type = c(
      "Aldinger annotated RDS",
      "consensus table",
      if (pilot_mode && !selected_sample_mode) "pilot PCW manifest" else "sample metadata"
    ),
    input = c(
      input_path,
      comparison_path,
      if (pilot_mode && !selected_sample_mode) {
        "config/resolution2_pilot_samples.csv"
      } else {
        metadata_path
      }
    ),
    exists = c(
      file.exists(input_path),
      file.exists(comparison_path),
      if (pilot_mode && !selected_sample_mode) TRUE else file.exists(metadata_path)
    ),
    stringsAsFactors = FALSE
  )
  compact_dry_run(
    paste0(
      "Consensus application task ", task_id, "/", nrow(samples),
      " [", sample_name, "]"
    ),
    inputs = input_checks$input,
    outputs = c(expected_outputs, marker_workbook_path),
    checks = setNames(input_checks$exists, input_checks$input_type)
  )
  quit(save = "no", status = if (all(input_checks$exists)) 0L else 1L)
}

if (pilot_mode && !selected_sample_mode) {
  if (!"PCW" %in% names(samples)) {
    stop("Clustering pilot manifest lacks the required PCW column.")
  }
  pcw_value <- as.character(samples$PCW[[task_id]])
} else {
  if (!file.exists(metadata_path)) stop("Sample metadata not found: ", metadata_path)
  sample_map <- readxl::read_excel(metadata_path)
  required_metadata_columns <- c("sample", "PCW")
  if (!all(required_metadata_columns %in% names(sample_map))) {
    stop("Sample metadata must contain: ", paste(required_metadata_columns, collapse = ", "))
  }
  base_sample_name <- sub("_\\d+$", "", sample_name)
  metadata_matches <- which(as.character(sample_map$sample) == base_sample_name)
  if (length(metadata_matches) != 1L) {
    stop(
      "Expected one PCW metadata match for '", base_sample_name,
      "'; found ", length(metadata_matches), "."
    )
  }
  pcw_value <- as.character(sample_map$PCW[[metadata_matches]])
}
if (is.na(pcw_value) || !nzchar(pcw_value)) {
  stop("PCW metadata is blank for '", sample_name, "'.")
}

if (!file.exists(input_path)) stop("Annotated input not found: ", input_path)
if (!file.exists(comparison_path)) stop("Consensus table not found: ", comparison_path)

protected_outputs <- c(expected_outputs, marker_workbook_path)
existing_outputs <- protected_outputs[file.exists(protected_outputs)]
if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Refusing to overwrite existing consensus outputs for ", sample_name,
    ":\n- ", paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing files."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_path)
if (!"seurat_clusters" %in% colnames(obj[[]])) {
  stop("Annotated object lacks the 'seurat_clusters' metadata column: ", input_path)
}
if (anyDuplicated(Cells(obj))) stop("Annotated object contains duplicate cell IDs: ", input_path)
if (resolution_mode) {
  required_pilot_metadata <- c(
    "whole_tissue_cluster_res1.5",
    pilot_cluster_column,
    pilot_graph_column,
    "seurat_clusters"
  )
  missing_pilot_metadata <- setdiff(required_pilot_metadata, colnames(obj[[]]))
  if (length(missing_pilot_metadata)) {
    stop(
      "Resolution-", pilot_resolution_tag, " Aldinger object lacks metadata: ",
      paste(missing_pilot_metadata, collapse = ", ")
    )
  }
  if (!identical(
    as.character(obj$seurat_clusters),
    as.character(obj[[pilot_cluster_column]][, 1])
  )) {
    stop(
      "Resolution-specific Aldinger object seurat_clusters does not exactly match ",
      pilot_cluster_column, "."
    )
  }
}

obj$PCW <- pcw_value

comparison <- read.csv(comparison_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("seurat_clusters", "consensus_label")
if (weighted_2of3) {
  required_columns <- c(
    required_columns, "consensus_support_n", "consensus_nonunknown_n",
    "consensus_status", "consensus_method"
  )
}
if (!all(required_columns %in% names(comparison))) {
  stop("Consensus table must contain: ", paste(required_columns, collapse = ", "))
}
if (weighted_2of3) {
  table_methods <- unique(as.character(comparison$consensus_method))
  if (!identical(table_methods, "weighted_2of3")) {
    stop("Consensus table was not produced by the weighted_2of3 method: ", comparison_path)
  }
}
if (anyDuplicated(as.character(comparison$seurat_clusters))) {
  stop("Consensus table contains duplicate seurat_clusters values: ", comparison_path)
}
table_labels <- trimws(as.character(comparison$consensus_label))
if (anyNA(table_labels) || any(!nzchar(table_labels))) {
  stop("Consensus table contains blank or missing consensus_label values: ", comparison_path)
}
if (any(tolower(table_labels) == "unknown")) {
  stop(
    "Consensus table contains the legacy exact label 'Unknown'. Rerun ",
    "xenium_annotate_02_build_consensus.R so Unknown clusters are numbered."
  )
}
unknown_table_labels <- unique(table_labels[is_unknown_consensus_label(table_labels)])
unknown_label_counts <- table(table_labels[is_unknown_consensus_label(table_labels)])
if (length(unknown_label_counts) && any(unknown_label_counts != 1L)) {
  stop(
    "Each Unknown-## label must map to exactly one consensus-table cluster; duplicates: ",
    paste(names(unknown_label_counts)[unknown_label_counts != 1L], collapse = ", ")
  )
}
if (length(unknown_table_labels) >= 1L && !requireNamespace("writexl", quietly = TRUE)) {
  stop(
    "The 'writexl' package is required because ", sample_name, " contains ",
    length(unknown_table_labels), " Unknown clusters. Restore it through the project ",
    "renv environment before rerunning this stage."
  )
}

cell_clusters <- as.character(obj$seurat_clusters)
cluster_ids <- as.character(comparison$seurat_clusters)
if (anyNA(cluster_ids) || any(!nzchar(cluster_ids))) {
  stop("Consensus table contains blank seurat_clusters values: ", comparison_path)
}
if (resolution_mode && !setequal(unique(cell_clusters), cluster_ids)) {
  stop(
    "Resolution-", pilot_resolution_tag,
    " object and consensus table do not contain identical cluster IDs for ",
    sample_name, "."
  )
}
cluster_lookup <- setNames(trimws(as.character(comparison$consensus_label)), cluster_ids)
unmapped_clusters <- setdiff(unique(cell_clusters), names(cluster_lookup))
if (length(unmapped_clusters) > 0L) {
  stop("No consensus label for cluster(s): ", paste(unmapped_clusters, collapse = ", "))
}

comparison_label_columns <- setdiff(names(comparison), "seurat_clusters")
for (column in comparison_label_columns) {
  column_lookup <- setNames(as.character(comparison[[column]]), cluster_ids)
  obj[[column]] <- unname(column_lookup[cell_clusters])
}

consensus_values <- trimws(as.character(obj$consensus_label))
if (anyNA(consensus_values) || any(!nzchar(consensus_values))) {
  stop("Consensus table produced blank or missing cell-level labels for ", sample_name, ".")
}
consensus_levels <- order_consensus_levels(consensus_values, celltype_order)
validate_palette(consensus_levels)
obj$consensus_label <- factor(consensus_values, levels = consensus_levels)
Idents(obj) <- "consensus_label"

if (!identical(as.character(Idents(obj)), as.character(obj$consensus_label))) {
  stop("Failed to set consensus_label as the active identity.")
}

save_consensus_plots(
  obj = obj,
  sample_name = sample_name,
  plot_paths = plot_paths,
  known_plot_paths = known_plot_paths,
  dotplot_paths = dotplot_paths,
  known_dotplot_paths = known_dotplot_paths,
  facet_point_size = facet_point_size,
  resolution_mode = resolution_mode,
  pilot_cluster_column = pilot_cluster_column,
  pilot_graph_column = pilot_graph_column,
  pilot_resolution_tag = pilot_resolution_tag,
  weighted_2of3 = weighted_2of3
)

save_unknown_marker_workbook(
  obj = obj,
  sample_name = sample_name,
  marker_workbook_path = marker_workbook_path,
  overwrite = overwrite
)

saveRDS(obj, output_path, compress = FALSE)
message(
  "Saved consensus-labelled object with PCW metadata and plots for ",
  sample_name, "."
)
