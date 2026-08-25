#!/usr/bin/env Rscript

# Read-only validation of pipeline outputs through xenium_annotate_03_apply_consensus.R.
# Run from the XeniumFCProject root:
#   Rscript scripts/validate_outputs_through_consensus.R
#
# Reports are written to outputs/validation/through_consensus/. No analysis
# objects or figures are modified.

options(warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(yaml)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Could not determine this script's location.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
project_root <- dirname(script_dir)
config_path <- file.path(project_root, "config", "config.yml")

# Also support a self-contained scripts checkout where config/ is inside the
# scripts directory. The normal HPC layout uses the first path above.
if (!file.exists(config_path)) {
  project_root <- script_dir
  config_path <- file.path(project_root, "config", "config.yml")
}
if (!file.exists(config_path)) stop("Could not find config/config.yml.")

config <- yaml::read_yaml(config_path)
resolve_project_path <- function(path) {
  candidates <- c(
    file.path(project_root, path),
    file.path(dirname(config_path), basename(path))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) normalizePath(existing[[1]], mustWork = TRUE) else candidates[[1]]
}

samples_path <- resolve_project_path(config$manifests$samples)
samples <- read.csv(
  samples_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character()
)
output_root <- file.path(project_root, config$project$outputs_dir)
report_dir <- file.path(output_root, "validation", "through_consensus")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

references <- c("Aldinger", "Sepp", "Science")
required_qc_metadata <- c(
  "nCount_Xenium", "nFeature_Xenium", "cell_area", "percent_control",
  "seurat_clusters"
)
required_consensus_columns <- c(
  "aldinger_cluster_majority", "aldinger_cluster_weighted",
  "sepp_cluster_majority", "sepp_cluster_weighted",
  "science_cluster_majority", "science_cluster_weighted",
  "consensus_label", "PCW"
)

file_is_nonempty <- function(path) {
  info <- suppressWarnings(file.info(path))
  file.exists(path) && !is.na(info$size) && info$size > 0
}

all_files_nonempty <- function(paths) {
  length(paths) > 0L && all(vapply(paths, file_is_nonempty, logical(1)))
}

collapse_values <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(x, collapse = ";") else ""
}

safe_read_csv <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) structure(list(error = conditionMessage(e)), class = "validation_error")
  )
}

append_issue <- function(issues, text) {
  c(issues, text)
}

initial_plot_paths <- function(sample_name) {
  file.path(
    output_root,
    "xenium", "preprocess", "03_clustered", "plots",
    paste0(
      sample_name,
      c(
        "_UMAP.tif", "_RawCluster_UMAP.tif",
        "_GlobalRawClustersSpatialPlot.tif",
        "_FacetRawClustersSpatialPlot.tif"
      )
    )
  )
}

abt_paths <- function(sample_name, reference) {
  reference_dir <- file.path(
    output_root, "xenium", "annotation", "01_label_transfer", tolower(reference)
  )
  table_dir <- file.path(reference_dir, "tables")
  plot_dir <- file.path(reference_dir, "plots")
  rds_dir <- file.path(reference_dir, "rds")
  list(
    majority_table = file.path(
      table_dir, paste0(sample_name, "_", reference, "_majority_vs_weighted.csv")
    ),
    count_table = file.path(
      table_dir, paste0(sample_name, "_", reference, "_prediction_cellcounts.csv")
    ),
    rds = file.path(
      rds_dir, paste0(sample_name, "_", reference, "_annotated.rds")
    ),
    plots = file.path(
      plot_dir,
      paste0(
        sample_name,
        c(
          paste0("_", reference, "_Broad_ClusterWeighted_UMAP.tif"),
          paste0("_", reference, "_Broad_ClusterMajority_UMAP.tif"),
          "_Broad_PredictionScores_Hist.tif",
          "_Broad_PredictionScores_Hist.pdf",
          "_Broad_GlobalSpatial_ClusterWeighted.tif",
          "_Broad_GlobalSpatial_ClusterMajority.tif",
          "_Broad_FacetSpatial_ClusterWeighted.tif",
          "_Broad_FacetSpatial_ClusterMajority.tif",
          "_Broad_Marker_DotPlot_Weighted.tif",
          "_Broad_Marker_DotPlot_Weighted.pdf",
          "_Broad_Marker_DotPlot_Majority.tif",
          "_Broad_Marker_DotPlot_Majority.pdf"
        )
      )
    )
  )
}

consensus_paths <- function(sample_name) {
  plot_dir <- file.path(
    output_root, "xenium", "annotation", "03_consensus_labels", "plots", "samples"
  )
  list(
    table = file.path(
      output_root, "xenium", "annotation", "02_consensus", "tables",
      paste0(sample_name, "_comparison_merged.csv")
    ),
    rds = file.path(
      output_root, "xenium", "annotation", "03_consensus_labels", "rds",
      paste0(sample_name, "_Consensus_annotated.rds")
    ),
    plots = file.path(
      plot_dir,
      paste0(
        sample_name,
        c(
          "_Consensus_UMAP.tif", "_Consensus_GlobalSpatial.tif",
          "_Consensus_FacetSpatial.tif", "_Consensus_Marker_DotPlot.tif",
          "_Consensus_Marker_DotPlot.pdf"
        )
      )
    )
  )
}

metadata_path <- resolve_project_path(config$manifests$sample_metadata)
if (!file.exists(metadata_path)) stop("Sample metadata not found: ", metadata_path)
sample_metadata <- readxl::read_excel(metadata_path)
if (!all(c("sample", "PCW") %in% names(sample_metadata))) {
  stop("Sample metadata must contain columns named 'sample' and 'PCW'.")
}

validate_reference <- function(reference) {
  config_key <- tolower(reference)
  configured_path <- config$inputs$references[[config_key]]
  path <- resolve_project_path(configured_path)
  issues <- character()
  n_cells <- NA_integer_
  n_features <- NA_integer_
  n_labels <- NA_integer_
  readable <- FALSE

  if (!file_is_nonempty(path)) {
    issues <- append_issue(issues, "missing_or_empty_reference_rds")
  } else {
    obj <- tryCatch(readRDS(path), error = function(e) e)
    if (inherits(obj, "error")) {
      issues <- append_issue(issues, paste0("reference_read_error: ", conditionMessage(obj)))
    } else {
      readable <- TRUE
      n_cells <- ncol(obj)
      n_features <- nrow(obj)
      if (!"clusters_refined" %in% colnames(obj[[]])) {
        issues <- append_issue(issues, "missing_clusters_refined")
      } else {
        labels <- as.character(obj$clusters_refined)
        n_labels <- length(unique(labels[!is.na(labels) & nzchar(labels)]))
        if (anyNA(labels) || any(!nzchar(labels))) {
          issues <- append_issue(issues, "blank_or_missing_clusters_refined")
        }
      }
      if (!"pca" %in% names(obj@reductions)) issues <- append_issue(issues, "missing_pca")
      if (!"umap" %in% names(obj@reductions)) issues <- append_issue(issues, "missing_umap")
      rm(obj)
      gc(verbose = FALSE)
    }
  }

  data.frame(
    reference = reference,
    path = path,
    readable = readable,
    n_cells = n_cells,
    n_features = n_features,
    n_labels = n_labels,
    status = if (length(issues)) "FAIL" else "PASS",
    issues = collapse_values(issues),
    stringsAsFactors = FALSE
  )
}

reference_report <- do.call(rbind, lapply(references, validate_reference))
write.csv(
  reference_report,
  file.path(report_dir, "reference_validation.csv"),
  row.names = FALSE,
  na = ""
)

consensus_label_reports <- list()

validate_sample <- function(sample_name) {
  message("Validating ", sample_name, " ...")
  issues <- character()
  clustered_path <- file.path(
    output_root, "xenium", "preprocess", "03_clustered", "rds",
    paste0(sample_name, "_CB_QC_cluster.rds")
  )
  qc_files <- file.path(
    output_root, "xenium", "preprocess", "02_qc", "reports",
    paste0(sample_name, c("_QCplots.pdf", "_QC_thresholds.txt"))
  )
  cons <- consensus_paths(sample_name)

  n_cells <- NA_integer_
  consensus_cells <- NA_integer_
  n_features <- NA_integer_
  n_clusters <- NA_integer_
  n_consensus_labels <- NA_integer_
  unknown_cells <- NA_integer_
  unknown_percent <- NA_real_
  pcw <- ""
  cluster_ids <- character()
  cluster_sizes <- integer()
  clustered_cells <- character()

  if (!file_is_nonempty(clustered_path)) {
    issues <- append_issue(issues, "missing_or_empty_clustered_rds")
  } else {
    clustered <- tryCatch(readRDS(clustered_path), error = function(e) e)
    if (inherits(clustered, "error")) {
      issues <- append_issue(
        issues, paste0("clustered_rds_read_error: ", conditionMessage(clustered))
      )
    } else {
      n_cells <- ncol(clustered)
      n_features <- nrow(clustered)
      clustered_cells <- colnames(clustered)
      if (anyDuplicated(clustered_cells)) issues <- append_issue(issues, "duplicate_cell_ids")

      missing_meta <- setdiff(required_qc_metadata, colnames(clustered[[]]))
      if (length(missing_meta)) {
        issues <- append_issue(
          issues, paste0("missing_clustered_metadata:", paste(missing_meta, collapse = "|"))
        )
      } else {
        cluster_ids <- unique(as.character(clustered$seurat_clusters))
        n_clusters <- length(cluster_ids)
        cluster_sizes <- table(as.character(clustered$seurat_clusters))
        qc_bad <-
          !is.finite(clustered$nCount_Xenium) |
          !is.finite(clustered$nFeature_Xenium) |
          !is.finite(clustered$cell_area) |
          !is.finite(clustered$percent_control) |
          clustered$nCount_Xenium <= 0 |
          clustered$nFeature_Xenium <= 0 |
          clustered$cell_area <= 0 |
          clustered$percent_control < 0
        if (any(qc_bad)) issues <- append_issue(issues, "invalid_qc_values")
      }
      if (is.null(clustered@misc$QC_thresholds)) {
        issues <- append_issue(issues, "missing_QC_thresholds")
      }
      if (!"Xenium" %in% Assays(clustered)) issues <- append_issue(issues, "missing_Xenium_assay")
      if (!"pca" %in% names(clustered@reductions)) issues <- append_issue(issues, "missing_pca")
      if (!"umap" %in% names(clustered@reductions)) issues <- append_issue(issues, "missing_umap")
      if (!"Xenium_snn" %in% names(clustered@graphs)) {
        issues <- append_issue(issues, "missing_Xenium_snn_graph")
      }
      coords <- tryCatch(GetTissueCoordinates(clustered), error = function(e) e)
      if (inherits(coords, "error") || nrow(coords) != n_cells) {
        issues <- append_issue(issues, "missing_or_incomplete_spatial_coordinates")
      }
      rm(clustered, coords)
      gc(verbose = FALSE)
    }
  }

  if (!all_files_nonempty(qc_files)) issues <- append_issue(issues, "incomplete_qc_reports")
  if (!all_files_nonempty(initial_plot_paths(sample_name))) {
    issues <- append_issue(issues, "incomplete_initial_cluster_plots")
  }

  for (reference in references) {
    paths <- abt_paths(sample_name, reference)
    ref_key <- tolower(reference)
    if (!file_is_nonempty(paths$rds)) {
      issues <- append_issue(issues, paste0(ref_key, "_missing_or_empty_annotated_rds"))
    }
    if (!all_files_nonempty(paths$plots)) {
      issues <- append_issue(issues, paste0(ref_key, "_incomplete_plots"))
    }

    majority <- safe_read_csv(paths$majority_table)
    if (inherits(majority, "validation_error")) {
      issues <- append_issue(issues, paste0(ref_key, "_majority_table_unreadable"))
    } else {
      required <- c("seurat_clusters", "cluster_majority", "cluster_weighted")
      if (!all(required %in% names(majority))) {
        issues <- append_issue(issues, paste0(ref_key, "_majority_table_bad_columns"))
      } else {
        table_clusters <- as.character(majority$seurat_clusters)
        if (anyDuplicated(table_clusters)) {
          issues <- append_issue(issues, paste0(ref_key, "_duplicate_cluster_rows"))
        }
        if (length(cluster_ids) && !setequal(cluster_ids, table_clusters)) {
          issues <- append_issue(issues, paste0(ref_key, "_cluster_coverage_mismatch"))
        }
        if (anyNA(majority$cluster_majority) || any(!nzchar(majority$cluster_majority)) ||
            anyNA(majority$cluster_weighted) || any(!nzchar(majority$cluster_weighted))) {
          issues <- append_issue(issues, paste0(ref_key, "_blank_cluster_labels"))
        }
      }
    }

    counts <- safe_read_csv(paths$count_table)
    if (inherits(counts, "validation_error")) {
      issues <- append_issue(issues, paste0(ref_key, "_count_table_unreadable"))
    } else if (!"seurat_clusters" %in% names(counts) || ncol(counts) < 2L) {
      issues <- append_issue(issues, paste0(ref_key, "_count_table_bad_columns"))
    } else if (length(cluster_sizes)) {
      count_values <- suppressWarnings(as.matrix(data.frame(
        lapply(counts[setdiff(names(counts), "seurat_clusters")], as.numeric),
        check.names = FALSE
      )))
      observed <- rowSums(count_values, na.rm = FALSE)
      expected <- unname(cluster_sizes[as.character(counts$seurat_clusters)])
      if (anyNA(observed) || anyNA(expected) || !identical(as.numeric(observed), as.numeric(expected))) {
        issues <- append_issue(issues, paste0(ref_key, "_prediction_counts_mismatch"))
      }
    }
  }

  comparison <- safe_read_csv(cons$table)
  comparison_ok <- TRUE
  if (inherits(comparison, "validation_error")) {
    comparison_ok <- FALSE
    issues <- append_issue(issues, "consensus_table_unreadable")
  } else {
    needed <- c("seurat_clusters", required_consensus_columns[1:7])
    if (!all(needed %in% names(comparison))) {
      comparison_ok <- FALSE
      issues <- append_issue(issues, "consensus_table_bad_columns")
    } else {
      comparison_clusters <- as.character(comparison$seurat_clusters)
      if (anyDuplicated(comparison_clusters)) {
        comparison_ok <- FALSE
        issues <- append_issue(issues, "consensus_table_duplicate_clusters")
      }
      if (length(cluster_ids) && !setequal(cluster_ids, comparison_clusters)) {
        comparison_ok <- FALSE
        issues <- append_issue(issues, "consensus_table_cluster_coverage_mismatch")
      }
      labels <- trimws(as.character(comparison$consensus_label))
      if (anyNA(labels) || any(!nzchar(labels))) {
        comparison_ok <- FALSE
        issues <- append_issue(issues, "consensus_table_blank_labels")
      }
    }
  }

  if (!file_is_nonempty(cons$rds)) {
    issues <- append_issue(issues, "missing_or_empty_consensus_rds")
  } else {
    obj <- tryCatch(readRDS(cons$rds), error = function(e) e)
    if (inherits(obj, "error")) {
      issues <- append_issue(issues, paste0("consensus_rds_read_error: ", conditionMessage(obj)))
    } else {
      consensus_cells <- ncol(obj)
      if (length(clustered_cells) && !identical(colnames(obj), clustered_cells)) {
        issues <- append_issue(issues, "consensus_cell_ids_or_order_mismatch")
      }
      missing_meta <- setdiff(required_consensus_columns, colnames(obj[[]]))
      if (length(missing_meta)) {
        issues <- append_issue(
          issues, paste0("missing_consensus_metadata:", paste(missing_meta, collapse = "|"))
        )
      } else {
        consensus_values <- as.character(obj$consensus_label)
        n_consensus_labels <- length(unique(consensus_values))
        unknown_cells <- sum(is.na(consensus_values) | !nzchar(consensus_values) |
                               tolower(consensus_values) == "unknown")
        unknown_percent <- 100 * unknown_cells / ncol(obj)
        label_counts <- as.data.frame(table(consensus_label = consensus_values, useNA = "ifany"))
        names(label_counts)[[2]] <- "n_cells"
        label_counts$sample_id <- sample_name
        consensus_label_reports[[sample_name]] <<- label_counts[
          label_counts$n_cells > 0,
          c("sample_id", "consensus_label", "n_cells"),
          drop = FALSE
        ]
        pcw <- collapse_values(obj$PCW)
        if (length(unique(as.character(obj$PCW))) != 1L || !nzchar(pcw)) {
          issues <- append_issue(issues, "PCW_not_one_nonblank_value")
        }
        if (!identical(as.character(Idents(obj)), consensus_values)) {
          issues <- append_issue(issues, "active_idents_do_not_match_consensus")
        }
        if (comparison_ok && "seurat_clusters" %in% colnames(obj[[]])) {
          lookup <- setNames(
            trimws(as.character(comparison$consensus_label)),
            as.character(comparison$seurat_clusters)
          )
          expected_labels <- unname(lookup[as.character(obj$seurat_clusters)])
          expected_labels[is.na(expected_labels) | !nzchar(expected_labels)] <- "Unknown"
          if (!identical(consensus_values, expected_labels)) {
            issues <- append_issue(issues, "consensus_labels_do_not_match_table")
          }
        }
      }
      rm(obj)
      gc(verbose = FALSE)
    }
  }

  base_sample <- sub("_\\d+$", "", sample_name)
  metadata_matches <- which(as.character(sample_metadata$sample) == base_sample)
  if (length(metadata_matches) != 1L) {
    issues <- append_issue(issues, "sample_metadata_match_not_unique")
  } else if (nzchar(pcw) && pcw != as.character(sample_metadata$PCW[[metadata_matches]])) {
    issues <- append_issue(issues, "PCW_does_not_match_sample_metadata")
  }

  if (!all_files_nonempty(cons$plots)) issues <- append_issue(issues, "incomplete_consensus_plots")

  data.frame(
    sample_id = sample_name,
    n_cells = n_cells,
    consensus_cells = consensus_cells,
    n_features = n_features,
    n_clusters = n_clusters,
    n_consensus_labels = n_consensus_labels,
    unknown_cells = unknown_cells,
    unknown_percent = round(unknown_percent, 4),
    PCW = pcw,
    status = if (length(issues)) "FAIL" else "PASS",
    issues = collapse_values(issues),
    stringsAsFactors = FALSE
  )
}

sample_report <- do.call(rbind, lapply(samples$sample_id, validate_sample))
write.csv(
  sample_report,
  file.path(report_dir, "sample_validation.csv"),
  row.names = FALSE,
  na = ""
)

label_report <- if (length(consensus_label_reports)) {
  do.call(rbind, consensus_label_reports)
} else {
  data.frame(sample_id = character(), consensus_label = character(), n_cells = integer())
}
write.csv(
  label_report,
  file.path(report_dir, "consensus_label_counts.csv"),
  row.names = FALSE,
  na = ""
)

failed_references <- reference_report$reference[reference_report$status != "PASS"]
failed_samples <- sample_report$sample_id[sample_report$status != "PASS"]
summary_lines <- c(
  paste("Validation time:", format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
  paste("Project root:", project_root),
  paste("Samples expected:", nrow(samples)),
  paste("Samples passed:", sum(sample_report$status == "PASS")),
  paste("Samples failed:", length(failed_samples)),
  paste("References passed:", sum(reference_report$status == "PASS"), "of", nrow(reference_report)),
  paste("Clustered objects validated:", sum(!is.na(sample_report$n_cells))),
  paste("Total clustered cells:", sum(sample_report$n_cells, na.rm = TRUE)),
  paste("Consensus objects validated:", sum(!is.na(sample_report$consensus_cells))),
  paste("Cells in validated consensus objects:", sum(sample_report$consensus_cells, na.rm = TRUE)),
  paste("Unknown cells in validated consensus objects:", sum(sample_report$unknown_cells, na.rm = TRUE)),
  paste("Failed references:", if (length(failed_references)) paste(failed_references, collapse = ", ") else "none"),
  paste("Failed samples:", if (length(failed_samples)) paste(failed_samples, collapse = ", ") else "none")
)
writeLines(summary_lines, file.path(report_dir, "validation_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Reports written to:", report_dir, "\n")

if (length(failed_references) || length(failed_samples)) quit(save = "no", status = 1L)
message("All outputs through consensus labels passed validation.")
