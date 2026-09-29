# Manifest-driven annotation for every biological sample and reference.
# Use --list or --dry-run before launching heavy Seurat work.

rm(list = ls())

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)
sample_list <- sample_manifest$sample_id

pilot_manifest <- load_resolution2_pilot_manifest(config)
required_pilot_columns <- c("task_id", "sample_id", "PCW", "age_group")
if (!all(required_pilot_columns %in% names(pilot_manifest))) {
  stop(
    "Resolution-2.0 pilot manifest must contain: ",
    paste(required_pilot_columns, collapse = ", ")
  )
}
pilot_manifest_counts <- vapply(
  pilot_manifest$sample_id,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(pilot_manifest_counts != 1L)) {
  stop("Every clustering-pilot sample must map to exactly one config/samples.csv row.")
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

args <- commandArgs(trailingOnly = TRUE)

valid_options <- c(
  "--dry-run", "--overwrite", "--pilot-res2", "--pilot-res3",
  "--pilot-res4", "--pilot-res5", "--all-samples-res4", "--all-samples-res5",
  "--selected-sample",
  "--list"
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
pilot_resolution <- if (all_samples_res5) 5.0 else if (all_samples_res4) 4.0 else if (pilot_res5) 5.0 else if (pilot_res4) 4.0 else if (pilot_res3) 3.0 else 2.0
pilot_resolution_tag <- sprintf("%.1f", pilot_resolution)
pilot_stage <- if (selected_sample_mode) {
  "03h_resolution5_selected_sample"
} else if (all_samples_res4) {
  "03g_resolution4_all_samples"
} else if (all_samples_res5) {
  "04_resolution5_clustered"
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
pilot_graph_column <- paste0("Xenium_snn_res.", format(pilot_resolution, trim = TRUE))
facet_point_size <- if (all_samples_mode || selected_sample_mode) {
  0.01
} else if (pilot_res3 || pilot_res4 || pilot_res5) {
  0.03
} else {
  0.1
}
list_requested <- "--list" %in% args
sample_list <- if (selected_sample_mode) {
  selected_sample_id
} else if (pilot_mode) {
  pilot_manifest$sample_id
} else {
  sample_manifest$sample_id
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

legacy_reference_cap <- as.integer(
  config$label_transfer$reference_cells_per_identity
)
production_res5_settings <- config$label_transfer$resolution5_all_samples
production_res5_contract <- all_samples_res5
reference_sampling_mode <- if (all_samples_res5) {
  tolower(trimws(as.character(
    production_res5_settings$reference_sampling
  )))
} else {
  "cap_per_identity"
}
min_top2_mean_score_margin <- if (all_samples_res5) {
  as.numeric(production_res5_settings$min_top2_mean_score_margin)
} else {
  0
}
prediction_score_threshold <- as.numeric(
  config$label_transfer$prediction_score_threshold
)
random_seed <- as.integer(config$runtime$random_seed)
transfer_dimensions <- seq_len(as.integer(config$label_transfer$dimensions))
k_anchor <- as.integer(config$label_transfer$k_anchor)
k_score <- as.integer(config$label_transfer$k_score)
k_weight <- as.integer(config$label_transfer$k_weight)
variable_features <- as.integer(config$label_transfer$variable_features)

if (!reference_sampling_mode %in% c("full", "cap_per_identity")) {
  stop("Unsupported reference sampling mode: ", reference_sampling_mode)
}
if (!is.finite(min_top2_mean_score_margin) ||
    min_top2_mean_score_margin < 0 || min_top2_mean_score_margin >= 1) {
  stop("The minimum top-two mean-score margin must be at least 0 and less than 1.")
}
if (!is.finite(prediction_score_threshold) ||
    prediction_score_threshold <= 0 || prediction_score_threshold >= 1) {
  stop("The prediction-score threshold must be between 0 and 1.")
}
if (anyNA(c(
  random_seed, max(transfer_dimensions), k_anchor, k_score, k_weight,
  variable_features, legacy_reference_cap
))) {
  stop("Label-transfer integer settings must not be missing.")
}

if (list_requested) {
  list_args <- args[!args %in% valid_options]
  if (length(list_args)) stop("--list does not accept REFERENCE or TASK_ID arguments.")
  cat(
    "Mode:", mode_description, "\n"
  )
  cat("References:", paste(names(reference_paths), collapse = ", "), "\n")
  cat("Reference sampling:", reference_sampling_mode, "\n")
  if (reference_sampling_mode == "cap_per_identity") {
    cat("Reference cap per identity:", legacy_reference_cap, "\n")
  }
  if (production_res5_contract) {
    cat(
      "Reference-vote thresholds: winner mean score >=",
      prediction_score_threshold, "; top-two mean-score margin >=",
      min_top2_mean_score_margin, "\n"
    )
  } else {
    cat(
      "Legacy vote rule: summed cell maximum scores; no reference abstention; ",
      "majority-label mean-score threshold =", prediction_score_threshold, "\n"
    )
  }
  cat(
    "Transfer settings: RPCA; dimensions=1:", max(transfer_dimensions),
    "; k.anchor=", k_anchor, "; k.score=", k_score,
    "; k.weight=", k_weight, "; seed=", random_seed, "\n",
    sep = ""
  )
  write.table(
    data.frame(task_id = seq_along(sample_list), sample_id = sample_list),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  quit(save = "no", status = 0L)
}

dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) {
  stop("--dry-run and --overwrite cannot be used together.")
}

job_args <- args[!args %in% valid_options]
if (length(job_args) != 2L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_01_label_transfer_rpca.R ",
    "[--pilot-res2|--pilot-res3|--pilot-res4|--pilot-res5|--all-samples-res4|--all-samples-res5] ",
    "[--selected-sample] ",
    "[--dry-run|--overwrite] REFERENCE TASK_ID"
  )
}

reference_match <- match(tolower(job_args[[1]]), tolower(names(reference_paths)))
if (is.na(reference_match)) {
  stop("REFERENCE must be Aldinger, Sepp, or Science.")
}
reference_name <- names(reference_paths)[reference_match]
reference_path <- unname(reference_paths[[reference_name]])
reference_key <- tolower(reference_name)

task_id <- suppressWarnings(as.integer(job_args[[2]]))
if (is.na(task_id) || task_id < 1L || task_id > length(sample_list)) {
  stop("TASK_ID must be between 1 and ", length(sample_list), ".")
}

current_sample <- sample_list[[task_id]]
if (resolution_mode) {
  pilot_root <- here(
    "outputs", "xenium", "preprocess", pilot_stage
  )
  input_file <- file.path(
    pilot_root, "rds",
    paste0(
      current_sample, "_whole_tissue_Res", pilot_resolution_tag,
      if (all_samples_mode) "" else "_pilot", ".rds"
    )
  )
  annotation_root <- if (all_samples_res4) {
    here("outputs", "xenium", "annotation", "resolution4_all_samples")
  } else if (all_samples_res5) {
    here("outputs", "xenium", "annotation")
  } else {
    file.path(pilot_root, "annotation")
  }
} else {
  input_file <- here(
    "outputs", "xenium", "preprocess", "03_clustered", "rds",
    paste0(current_sample, "_CB_QC_cluster.rds")
  )
  annotation_root <- here(
    "outputs", "xenium", "annotation", "legacy_resolution1_5"
  )
}
annotation_dir <- file.path(annotation_root, "01_label_transfer", reference_key)
plots_dir <- file.path(annotation_dir, "plots")
tables_dir <- file.path(annotation_dir, "tables")
rds_dir <- file.path(annotation_dir, "rds")

expected_outputs <- c(
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_majority_vs_weighted.csv")),
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_prediction_cellcounts.csv")),
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_reference_class_balance.csv")),
  file.path(tables_dir, paste0(current_sample, "_", reference_name, "_transfer_provenance.csv")),
  file.path(plots_dir, paste0(current_sample, "_", reference_name, "_Broad_ClusterWeighted_UMAP.tif")),
  file.path(plots_dir, paste0(current_sample, "_", reference_name, "_Broad_ClusterMajority_UMAP.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_PredictionScores_Hist.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_PredictionScores_Hist.pdf")),
  file.path(plots_dir, paste0(current_sample, "_Broad_GlobalSpatial_ClusterWeighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_GlobalSpatial_ClusterMajority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_FacetSpatial_ClusterWeighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_FacetSpatial_ClusterMajority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Weighted.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Weighted.pdf")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Majority.tif")),
  file.path(plots_dir, paste0(current_sample, "_Broad_Marker_DotPlot_Majority.pdf")),
  file.path(rds_dir, paste0(current_sample, "_", reference_name, "_annotated.rds"))
)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]

if (dry_run) {
  compact_dry_run(
    paste0(
      reference_name, " transfer task ", task_id, "/", length(sample_list),
      " [", current_sample, "]"
    ),
    inputs = c(input_file, reference_path),
    outputs = expected_outputs
  )
  inputs_ready <- file.exists(input_file) && file.exists(reference_path)
  quit(save = "no", status = if (inputs_ready) 0L else 1L)
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (!file.exists(reference_path)) stop("Reference file not found: ", reference_path)

if (length(existing_outputs) > 0L && !overwrite) {
  stop(
    "Refusing to overwrite existing ", reference_name,
    " annotation outputs for ", current_sample, ":\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}
if (length(existing_outputs) > 0L) {
  warning(
    "Overwriting ", length(existing_outputs), " existing ", reference_name,
    " annotation outputs for ", current_sample, "."
  )
}

library(future)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Cairo)

options(future.globals.maxSize = 400 * 1024^3)

message(
  "### Processing ", reference_name, " annotation for sample [",
  task_id, "/", length(sample_list), "]: ", current_sample, " ###"
)

# =========================================================
# Annotation Function
# =========================================================
annotate_xenium_from_ref <- function(xenium_obj,
                                     sample_name,
                                     reference_name,
                                     reference_path,
                                     annotation_dir,
                                     production_res5_contract,
                                     reference_sampling_mode,
                                     reference_cap_per_identity,
                                     prediction_score_threshold,
                                     min_top2_mean_score_margin,
                                     random_seed,
                                     variable_features,
                                     transfer_dimensions,
                                     k_anchor,
                                     k_score,
                                     k_weight) {
  
  ## ----------------------------
  ## 0. Setup & Paths
  ## ----------------------------
  set.seed(random_seed)
  plan("sequential") 
  
  # Load custom palette and ordering from your specific script
  source(here("scripts", "color_palette.R")) 
  
  pred_score_thresh <- prediction_score_threshold
  plots_dir <- file.path(annotation_dir, "plots")
  if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  tables_dir <- file.path(annotation_dir, "tables")
  if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
  
  ## ----------------------------
  ## 1. Load Reference
  ## ----------------------------
  ref_path <- reference_path
  if (!file.exists(ref_path)) stop("Reference file not found: ", ref_path)
  reference <- readRDS(ref_path)
  
  ## ----------------------------
  ## 2. Preparation & Anchor Finding
  ## ----------------------------
  # Set reference assay safely
  if ("RNA" %in% Assays(reference)) {
    DefaultAssay(reference) <- "RNA"
  } else {
    # Fallback to whatever assay exists (e.g., "originalexp")
    DefaultAssay(reference) <- Assays(reference)[1] 
  }
  
  # Set xenium assay safely
  if ("Xenium" %in% Assays(xenium_obj)) {
    DefaultAssay(xenium_obj) <- "Xenium"
  } else {
    DefaultAssay(xenium_obj) <- Assays(xenium_obj)[1]
  }
  if (anyDuplicated(Cells(xenium_obj))) {
    stop("Query contains duplicate cell IDs: ", sample_name)
  }
  if (!identical(rownames(xenium_obj[[]]), Cells(xenium_obj))) {
    stop("Query metadata rows do not exactly align with cell IDs: ", sample_name)
  }
  
  if (anyDuplicated(Cells(reference))) {
    stop("Reference contains duplicate cell IDs: ", reference_name)
  }
  if (!identical(rownames(reference[[]]), Cells(reference))) {
    stop("Reference metadata rows do not exactly align with cell IDs: ", reference_name)
  }
  if (!"clusters_refined" %in% colnames(reference[[]])) {
    stop(reference_name, " reference lacks clusters_refined metadata.")
  }
  reference_labels <- trimws(as.character(reference$clusters_refined))
  if (anyNA(reference_labels) || any(!nzchar(reference_labels))) {
    stop(reference_name, " reference has blank or missing clusters_refined labels.")
  }
  reference$clusters_refined <- reference_labels

  reference <- FindVariableFeatures(
    reference,
    selection.method = "vst",
    nfeatures = variable_features,
    verbose = FALSE
  )
  shared_genes <- intersect(rownames(reference), rownames(xenium_obj))

  Idents(reference) <- "clusters_refined"
  set.seed(random_seed)
  reference_for_transfer <- if (reference_sampling_mode == "full") {
    reference
  } else {
    subset(reference, downsample = reference_cap_per_identity)
  }

  class_counts <- function(object, count_name) {
    counts <- as.data.frame(
      table(reference_class = as.character(object$clusters_refined)),
      stringsAsFactors = FALSE
    )
    names(counts)[[2]] <- count_name
    counts
  }
  reference_balance <- merge(
    class_counts(reference, "full_reference_cells"),
    class_counts(reference_for_transfer, "used_reference_cells"),
    by = "reference_class", all = TRUE, sort = FALSE
  )
  reference_balance$full_reference_cells[
    is.na(reference_balance$full_reference_cells)
  ] <- 0L
  reference_balance$used_reference_cells[
    is.na(reference_balance$used_reference_cells)
  ] <- 0L
  expected_used_cells <- if (reference_sampling_mode == "full") {
    reference_balance$full_reference_cells
  } else {
    pmin(reference_balance$full_reference_cells, reference_cap_per_identity)
  }
  if (!identical(
    as.integer(reference_balance$used_reference_cells),
    as.integer(expected_used_cells)
  )) {
    stop("Reference sampling retained an unexpected number of cells per identity.")
  }
  reference_balance$retained_fraction <- (
    reference_balance$used_reference_cells /
      reference_balance$full_reference_cells
  )
  reference_balance$sample_id <- sample_name
  reference_balance$reference <- reference_name
  reference_balance$reference_sampling_mode <- reference_sampling_mode
  reference_balance$reference_cap_per_identity <- if (
    reference_sampling_mode == "full"
  ) {
    NA_integer_
  } else {
    reference_cap_per_identity
  }
  reference_balance <- reference_balance[c(
    "sample_id", "reference", "reference_sampling_mode",
    "reference_class", "full_reference_cells", "used_reference_cells",
    "retained_fraction", "reference_cap_per_identity"
  )]
  full_reference_cells_total <- ncol(reference)
  used_reference_cells_total <- ncol(reference_for_transfer)

  # Avoid retaining both the full input binding and its scaled copy.
  rm(reference)
  invisible(gc(verbose = FALSE))

  transfer_features <- intersect(
    VariableFeatures(reference_for_transfer), shared_genes
  )
  if (length(transfer_features) <= max(transfer_dimensions)) {
    stop(
      "Only ", length(transfer_features),
      " transfer features are available; more than ",
      max(transfer_dimensions), " are required."
    )
  }

  reference_for_transfer <- ScaleData(
    reference_for_transfer, features = transfer_features, verbose = FALSE
  )
  xenium_obj         <- ScaleData(xenium_obj, features = transfer_features, verbose = FALSE)

  set.seed(random_seed)
  reference_for_transfer <- RunPCA(
    reference_for_transfer, features = transfer_features, verbose = FALSE
  )
  set.seed(random_seed)
  xenium_obj         <- RunPCA(xenium_obj, features = transfer_features, verbose = FALSE)
  
  cat("Shared features for transfer:", length(transfer_features), "\n")
  
  set.seed(random_seed)
  anchors <- FindTransferAnchors(
    reference = reference_for_transfer,
    query = xenium_obj,
    normalization.method = "LogNormalize", 
    reduction = "rpca", 
    features = transfer_features,
    dims = transfer_dimensions,
    k.anchor = k_anchor,
    k.score = k_score,
    approx.pca = TRUE
  )
  
  ## ----------------------------
  ## 3. Transfer & Voting Logic
  ## ----------------------------
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_for_transfer$clusters_refined,
    dims = transfer_dimensions,
    k.weight = k_weight,
    store.weights = FALSE
  )

  predictions <- as.data.frame(predictions, check.names = FALSE)
  if (!identical(rownames(predictions), Cells(xenium_obj))) {
    if (!setequal(rownames(predictions), Cells(xenium_obj))) {
      stop("Transfer predictions do not align one-to-one with query cells.")
    }
    predictions <- predictions[Cells(xenium_obj), , drop = FALSE]
  }

  xenium_obj <- AddMetaData(xenium_obj, predictions)
  xenium_obj$high_conf <- if (production_res5_contract) {
    xenium_obj$prediction.score.max >= pred_score_thresh
  } else {
    xenium_obj$prediction.score.max > pred_score_thresh
  }

  set.seed(random_seed)
  xenium_obj <- RunUMAP(
    xenium_obj, dims = transfer_dimensions, reduction = "pca"
  )
  
  # -------------------
  # Thresholded Majority Voting
  # -------------------
  majority_labels <- xenium_obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      cluster_majority = names(sort(table(predicted.id), decreasing = TRUE))[1],
      mean_max_score = mean(prediction.score.max), # Calculate cluster average confidence
      .groups = "drop"
    ) %>%
    mutate(
      cluster_majority = ifelse(mean_max_score >= pred_score_thresh, cluster_majority, "Unknown")
    )
  majority_labels$seurat_clusters <- as.character(
    majority_labels$seurat_clusters
  )
  
  xenium_obj$cluster_majority <- majority_labels$cluster_majority[match(xenium_obj$seurat_clusters, majority_labels$seurat_clusters)]
  
  cluster_ids <- as.character(xenium_obj$seurat_clusters)
  cluster_levels <- unique(cluster_ids)
  numeric_clusters <- suppressWarnings(as.numeric(cluster_levels))
  cluster_levels <- cluster_levels[order(
    is.na(numeric_clusters), numeric_clusters, cluster_levels
  )]
  if (production_res5_contract) {
    # Resolution-5 production: use the same per-class mean-score rule tested in
    # the sampling/abstention sensitivity analysis, then allow the reference to
    # abstain when either approved threshold fails.
    score_cols <- grep(
      "^prediction\\.score\\.", colnames(xenium_obj@meta.data), value = TRUE
    )
    score_cols <- score_cols[
      !score_cols %in% c("prediction.score.max", "prediction.score.id")
    ]
    if (length(score_cols) < 2L) {
      stop("TransferData returned fewer than two class-score columns.")
    }
    score_matrix <- as.matrix(
      xenium_obj@meta.data[, score_cols, drop = FALSE]
    )
    storage.mode(score_matrix) <- "double"
    class_labels <- sub("^prediction\\.score\\.", "", score_cols)
    if (anyDuplicated(class_labels)) {
      stop("TransferData returned duplicate class-score labels.")
    }
    colnames(score_matrix) <- class_labels
    if (any(!is.finite(score_matrix))) {
      stop("TransferData returned non-finite class scores.")
    }
    predicted_ids <- trimws(as.character(xenium_obj$predicted.id))
    if (anyNA(predicted_ids) || any(!nzchar(predicted_ids))) {
      stop("TransferData returned blank or missing predicted labels.")
    }

    weighted_rows <- lapply(cluster_levels, function(cluster_id) {
      cell_index <- which(cluster_ids == cluster_id)
      mean_scores <- colMeans(score_matrix[cell_index, , drop = FALSE])
      score_order <- order(mean_scores, decreasing = TRUE)
      winner <- names(mean_scores)[score_order[[1]]]
      runner_up <- names(mean_scores)[score_order[[2]]]
      winner_score <- unname(mean_scores[[winner]])
      runner_up_score <- unname(mean_scores[[runner_up]])
      top2_margin <- winner_score - runner_up_score
      failed_score <- winner_score < pred_score_thresh
      failed_margin <- top2_margin < min_top2_mean_score_margin
      abstained <- failed_score || failed_margin
      data.frame(
        seurat_clusters = cluster_id,
        n_cells = length(cell_index),
        cluster_weighted_before_abstention = winner,
        weighted_winner_mean_score = winner_score,
        weighted_runner_up_label = runner_up,
        weighted_runner_up_mean_score = runner_up_score,
        weighted_top2_mean_score_margin = top2_margin,
        weighted_winner_cell_fraction = mean(predicted_ids[cell_index] == winner),
        weighted_failed_score = failed_score,
        weighted_failed_margin = failed_margin,
        weighted_abstained = abstained,
        cluster_weighted = if (abstained) "Unknown" else winner,
        stringsAsFactors = FALSE
      )
    })
    weighted_labels <- do.call(rbind, weighted_rows)
    rownames(weighted_labels) <- NULL
  } else {
    # Preserve the established vote outside resolution-5 production: add each
    # cell's maximum prediction score within its predicted class, then retain
    # the class with the largest sum for the cluster.
    weighted_labels <- xenium_obj@meta.data %>%
      group_by(seurat_clusters, predicted.id) %>%
      summarise(score_sum = sum(prediction.score.max), .groups = "drop") %>%
      group_by(seurat_clusters) %>%
      slice_max(score_sum, n = 1) %>%
      ungroup() %>%
      select(seurat_clusters, cluster_weighted = predicted.id)
    weighted_labels$seurat_clusters <- as.character(
      weighted_labels$seurat_clusters
    )
  }
  if (anyDuplicated(weighted_labels$seurat_clusters) ||
      !setequal(weighted_labels$seurat_clusters, unique(cluster_ids))) {
    stop(
      "Weighted-vote table does not map one-to-one to query clusters; ",
      "inspect tied class-score sums."
    )
  }

  weighted_match <- match(cluster_ids, weighted_labels$seurat_clusters)
  metadata_columns <- if (production_res5_contract) {
    setdiff(names(weighted_labels), c("seurat_clusters", "n_cells"))
  } else {
    "cluster_weighted"
  }
  for (column in metadata_columns) {
    xenium_obj[[column]] <- weighted_labels[[column]][weighted_match]
  }
  
  ## ----------------------------
  ## 4. Export Comparison Table
  ## ----------------------------
  comparison <- if (production_res5_contract) {
    majority_labels %>%
      select(
        seurat_clusters, cluster_majority,
        majority_mean_max_score = mean_max_score
      ) %>%
      left_join(weighted_labels, by = "seurat_clusters") %>%
      mutate(
        reference_sampling_mode = .env$reference_sampling_mode,
        reference_cap_per_identity = if (
          .env$reference_sampling_mode == "full"
        ) NA_integer_ else .env$reference_cap_per_identity,
        prediction_score_threshold = .env$pred_score_thresh,
        min_top2_mean_score_margin = .env$min_top2_mean_score_margin,
        random_seed = .env$random_seed,
        transfer_dimensions = max(.env$transfer_dimensions),
        k_anchor = .env$k_anchor,
        k_score = .env$k_score,
        k_weight = .env$k_weight
      ) %>%
      select(
        seurat_clusters, cluster_majority, cluster_weighted,
        everything()
      )
  } else {
    majority_labels %>%
      select(seurat_clusters, cluster_majority) %>%
      left_join(
        weighted_labels %>%
          select(seurat_clusters, cluster_weighted),
        by = "seurat_clusters"
      )
  }
  comparison <- comparison %>%
    arrange(as.numeric(as.character(seurat_clusters)))

  if (nrow(comparison) != length(cluster_levels) ||
      anyDuplicated(as.character(comparison$seurat_clusters))) {
    stop("Comparison-table join changed cluster cardinality.")
  }
  
  write.csv(comparison, file = here(tables_dir, paste0(sample_name, "_", reference_name, "_majority_vs_weighted.csv")), row.names = FALSE)

  write.csv(
    reference_balance,
    file.path(
      tables_dir,
      paste0(sample_name, "_", reference_name, "_reference_class_balance.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  anchor_count <- tryCatch(
    nrow(methods::slot(anchors, "anchors")),
    error = function(e) NA_integer_
  )
  transfer_provenance <- data.frame(
    field = c(
      "sample_id", "reference", "query_cells", "reference_path",
      "reference_sampling_mode", "reference_cap_per_identity",
      "full_reference_cells", "used_reference_cells", "shared_genes",
      "transfer_features", "prediction_score_threshold",
      "min_top2_mean_score_margin", "dimensions", "reduction",
      "k_anchor", "k_score", "k_weight", "random_seed", "anchor_count",
      "reference_vote_rule"
    ),
    value = c(
      sample_name, reference_name, as.character(ncol(xenium_obj)), ref_path,
      reference_sampling_mode,
      if (reference_sampling_mode == "full") {
        "none"
      } else {
        as.character(reference_cap_per_identity)
      },
      as.character(full_reference_cells_total),
      as.character(used_reference_cells_total),
      as.character(length(shared_genes)),
      as.character(length(transfer_features)),
      as.character(pred_score_thresh),
      as.character(min_top2_mean_score_margin),
      paste(transfer_dimensions, collapse = ","), "rpca",
      as.character(k_anchor), as.character(k_score), as.character(k_weight),
      as.character(random_seed), as.character(anchor_count),
      if (production_res5_contract) {
        "mean class-score winner retained only when score and top-two margin both pass"
      } else {
        "legacy winner maximizes summed cell prediction.score.max within predicted class"
      }
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    transfer_provenance,
    file.path(
      tables_dir,
      paste0(sample_name, "_", reference_name, "_transfer_provenance.csv")
    ),
    row.names = FALSE,
    na = ""
  )
  
  ## ----------------------------
  ## 4c. Export Cluster-by-CellType Contingency Table
  ## ----------------------------
  # Note: Retains `predicted.id` so you can see the raw background noise/variance that caused a cluster to become "Unknown"
  count_matrix <- xenium_obj@meta.data %>%
    count(seurat_clusters, predicted.id) %>%
    pivot_wider(names_from = predicted.id, values_from = n, values_fill = list(n = 0)) %>%
    arrange(as.numeric(as.character(seurat_clusters)))
  
  write.csv(count_matrix, 
            file = here(tables_dir, paste0(sample_name, "_", reference_name, "_prediction_cellcounts.csv")), 
            row.names = FALSE)
  
  message("Contingency matrix exported for: ", sample_name)
  
  ## ----------------------------
  ## 5. Visualizations
  ## ----------------------------
  
  # Safely add "Unknown" to your palette if it isn't already in color_palette.R
  if (!"Unknown" %in% celltype_order) {
    celltype_order <- c(celltype_order, "Unknown")
    if (!"Unknown" %in% names(cluster_colors)) {
      cluster_colors["Unknown"] <- "grey80" 
    }
  }
  
  # Validate palette (exclude "Unknown" from the strict check so it doesn't crash)
  clusters_to_validate <- setdiff(unique(c(xenium_obj$cluster_weighted, xenium_obj$cluster_majority)), "Unknown")
  validate_palette(clusters_to_validate)
  
  xenium_obj$cluster_weighted <- factor(xenium_obj$cluster_weighted, levels = celltype_order)
  xenium_obj$cluster_majority <- factor(xenium_obj$cluster_majority, levels = celltype_order)
  
  # --- UMAPs ---
  p_umap_wei <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_weighted", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterWeighted_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_wei)
  dev.off()
  
  p_umap_maj <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, group.by = "cluster_majority", cols = cluster_colors)
  CairoTIFF(here(plots_dir, paste0(sample_name, "_", reference_name, "_Broad_ClusterMajority_UMAP.tif")), width = 8, height = 6, units = "in", res = 600)
  print(p_umap_maj)
  dev.off()
  
  # --- Histogram of scores ---
  p_hist <- ggplot(xenium_obj@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pred_score_thresh, linetype = "dashed", color = "red") +
    theme_minimal() + labs(title = paste(sample_name, "Scores"), x = "Max Prediction Score")
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.tif")), width = 6, height = 4, units = "in", res = 600)
  print(p_hist)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_PredictionScores_Hist.pdf")), p_hist, device = grDevices::cairo_pdf, width = 6, height = 4)
  
  # --- Spatial Plots ---
  p_wei <- ImageDimPlot(xenium_obj, group.by = "cluster_weighted", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Weighted"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterWeighted.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_wei)
  dev.off() # Fixed missing parentheses
  
  p_maj <- ImageDimPlot(xenium_obj, group.by = "cluster_majority", size = 0.75, cols = cluster_colors) + 
    ggtitle(paste(sample_name, "Majority")) # Fixed title
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_GlobalSpatial_ClusterMajority.tif")), width = 8, height = 8, units = "in", res = 600)
  print(p_maj)
  dev.off()
  
  # --- Spatial Facet Plots ---
  coords <- GetTissueCoordinates(xenium_obj) 
  
  # Weighted Facet
  plot_data_wei <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_weighted), levels = celltype_order))
  p_facet_wei <- ggplot(plot_data_wei, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = facet_point_size) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterWeighted.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_wei)
  dev.off()
  
  # Majority Facet
  plot_data_maj <- cbind(coords, cluster = factor(as.character(xenium_obj$cluster_majority), levels = celltype_order))
  p_facet_maj <- ggplot(plot_data_maj, aes(x = y, y = x, color = cluster)) + 
    geom_point(size = facet_point_size) + facet_wrap(~cluster) +
    scale_color_manual(values = cluster_colors) + coord_fixed() + theme_void() +
    theme(panel.background = element_rect(fill = "black"), plot.background = element_rect(fill = "black"),
          legend.position = "none", strip.text = element_text(color = "white"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_FacetSpatial_ClusterMajority.tif")), width = 12, height = 8, units = "in", res = 600)
  print(p_facet_maj)
  dev.off()
  
  # --- DotPlots ---
  existing_markers <- lapply(markers, function(x) intersect(x, rownames(xenium_obj)))
  
  # Weighted DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_weighted, levels = rev(celltype_order))
  p_dot_wei <- DotPlot(
    xenium_obj,
    features = existing_markers,
    assay = "Xenium",
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot_wei <- standardize_broad_dotplot(p_dot_wei) +
    ggtitle(paste(sample_name, "Markers (Weighted)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_wei)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Weighted.pdf")), p_dot_wei, device = grDevices::cairo_pdf, width = 10, height = 6)
  
  # Majority DotPlot
  Idents(xenium_obj) <- factor(xenium_obj$cluster_majority, levels = rev(celltype_order))
  p_dot_maj <- DotPlot(
    xenium_obj,
    features = existing_markers,
    assay = "Xenium",
    col.min = broad_dotplot_col_min,
    col.max = broad_dotplot_col_max,
    dot.min = broad_dotplot_dot_min / 100,
    dot.scale = broad_dotplot_dot_scale,
    scale.min = broad_dotplot_dot_min,
    scale.max = broad_dotplot_dot_max
  )
  p_dot_maj <- standardize_broad_dotplot(p_dot_maj) +
    ggtitle(paste(sample_name, "Markers (Majority)"))
  
  CairoTIFF(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.tif")), width = 10, height = 6, units = "in", res = 600)
  print(p_dot_maj)
  dev.off()
  ggplot2::ggsave(here(plots_dir, paste0(sample_name, "_Broad_Marker_DotPlot_Majority.pdf")), p_dot_maj, device = grDevices::cairo_pdf, width = 10, height = 6)
  
  ## ----------------------------
  ## 6. Save & Return
  ## ----------------------------
  output_dir <- file.path(annotation_dir, "rds")
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(xenium_obj, file.path(output_dir, paste0(sample_name, "_", reference_name, "_annotated.rds")))
  message("Successfully annotated and saved: ", sample_name)
  
  return(xenium_obj)
}

# =========================================================
# Execution
# =========================================================
seu <- readRDS(input_file)
if (resolution_mode) {
  required_pilot_metadata <- c(
    "whole_tissue_cluster_res1.5",
    pilot_cluster_column,
    pilot_graph_column,
    "seurat_clusters"
  )
  missing_pilot_metadata <- setdiff(required_pilot_metadata, colnames(seu[[]]))
  if (length(missing_pilot_metadata)) {
    stop(
      "Resolution-", pilot_resolution_tag, " input lacks metadata: ",
      paste(missing_pilot_metadata, collapse = ", ")
    )
  }
  candidate_clusters <- as.character(seu[[pilot_cluster_column]][, 1])
  seurat_clusters <- as.character(seu$seurat_clusters)
  if (anyNA(candidate_clusters) || !identical(candidate_clusters, seurat_clusters)) {
    stop(
      "Resolution-specific input seurat_clusters does not exactly match ",
      pilot_cluster_column, " for every cell."
    )
  }
  if (anyDuplicated(Cells(seu))) {
    stop("Resolution-specific input contains duplicate cell IDs.")
  }
}
annotate_xenium_from_ref(
  xenium_obj = seu,
  sample_name = current_sample,
  reference_name = reference_name,
  reference_path = reference_path,
  annotation_dir = annotation_dir,
  production_res5_contract = production_res5_contract,
  reference_sampling_mode = reference_sampling_mode,
  reference_cap_per_identity = legacy_reference_cap,
  prediction_score_threshold = prediction_score_threshold,
  min_top2_mean_score_margin = min_top2_mean_score_margin,
  random_seed = random_seed,
  variable_features = variable_features,
  transfer_dimensions = transfer_dimensions,
  k_anchor = k_anchor,
  k_score = k_score,
  k_weight = k_weight
)
