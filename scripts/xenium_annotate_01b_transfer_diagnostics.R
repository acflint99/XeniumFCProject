#!/usr/bin/env Rscript

# Focused transfer diagnostic for three resolution-5 queries and three
# references. It compares the existing RPCA preprocessing, matched 10k
# LogNormalize RPCA, and matched 10k LogNormalize CCA. Each anchor set is reused
# for k.weight 25 and 50. Outputs are tables only and never replace production.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
})

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
sample_manifest <- load_sample_manifest(config)

requested_samples <- c("GZFB_9_X_G_1", "GZFB_20_X_G_1", "GZFB_22_X_G_3")
references <- c("Aldinger", "Sepp", "Science")
k_weight_values <- c(25L, 50L)
k_anchor <- 10L
k_score <- as.integer(config$label_transfer$k_score)
transfer_dims <- seq_len(as.integer(config$label_transfer$dimensions))
legacy_threshold <- as.numeric(config$label_transfer$prediction_score_threshold)
normalization_scale_factor <- 10000

conditions <- data.frame(
  condition = c("baseline_rpca", "norm10k_rpca", "norm10k_cca"),
  renormalize_10k = c(FALSE, TRUE, TRUE),
  reduction = c("rpca", "rpca", "cca"),
  stringsAsFactors = FALSE
)

sample_counts <- vapply(
  requested_samples,
  function(sample_id) sum(sample_manifest$sample_id == sample_id),
  integer(1)
)
if (any(sample_counts != 1L)) {
  stop(
    "Every diagnostic sample must match exactly one config/samples.csv row: ",
    paste(requested_samples[sample_counts != 1L], collapse = ", ")
  )
}

target_manifest_candidates <- c(
  here("config", "transfer_diagnostic_clusters.csv"),
  here("scripts", "config", "transfer_diagnostic_clusters.csv")
)
target_manifest_existing <- target_manifest_candidates[
  file.exists(target_manifest_candidates)
]
if (!length(target_manifest_existing)) {
  stop(
    "Missing diagnostic target manifest. Checked:\n- ",
    paste(target_manifest_candidates, collapse = "\n- ")
  )
}
target_manifest_path <- normalizePath(target_manifest_existing[[1]], mustWork = TRUE)
targets <- read.csv(target_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_target_columns <- c(
  "sample_id", "seurat_cluster", "diagnostic_question", "selection_basis"
)
if (!all(required_target_columns %in% names(targets))) {
  stop(
    "Diagnostic target manifest must contain: ",
    paste(required_target_columns, collapse = ", ")
  )
}
targets$sample_id <- trimws(as.character(targets$sample_id))
targets$seurat_cluster <- trimws(as.character(targets$seurat_cluster))
if (anyNA(targets$sample_id) || any(!nzchar(targets$sample_id)) ||
    anyNA(targets$seurat_cluster) || any(!nzchar(targets$seurat_cluster))) {
  stop("Diagnostic target manifest contains blank sample or cluster IDs.")
}
if (!setequal(unique(targets$sample_id), requested_samples)) {
  stop("Diagnostic target samples do not exactly match the requested samples.")
}
if (anyDuplicated(targets[c("sample_id", "seurat_cluster")])) {
  stop("Diagnostic target manifest contains duplicate sample/cluster rows.")
}

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)

task_map <- do.call(rbind, lapply(seq_len(nrow(conditions)), function(i) {
  grid <- expand.grid(
    sample_id = requested_samples,
    reference = references,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$condition <- conditions$condition[[i]]
  grid$renormalize_10k <- conditions$renormalize_10k[[i]]
  grid$reduction <- conditions$reduction[[i]]
  grid
}))
task_map$task_id <- seq_len(nrow(task_map))
task_map <- task_map[c(
  "task_id", "condition", "renormalize_10k", "reduction",
  "sample_id", "reference"
)]
if (nrow(task_map) != 27L) stop("Expected exactly 27 diagnostic anchor tasks.")

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--list", "--dry-run", "--overwrite", "--combine")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
list_requested <- "--list" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
combine <- "--combine" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
positional_args <- args[!args %in% valid_options]

if (list_requested) {
  if (length(positional_args)) stop("--list does not accept TASK_ID.")
  cat("Scope: three resolution-5 samples x three references x three anchor conditions\n")
  cat("Fixed: k.anchor=10, k.score=", k_score,
      ", dimensions=1:", max(transfer_dims), ", threshold=", legacy_threshold, "\n", sep = "")
  cat("Each anchor task evaluates k.weight=25 and 50.\n")
  write.table(task_map, row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}

output_root <- here(config$project$outputs_dir)
query_root <- file.path(
  output_root, "xenium", "preprocess", "04_resolution5_clustered", "rds"
)
diagnostic_root <- file.path(
  output_root, "validation", "resolution5_method_development",
  "transfer_diagnostics_k10"
)

query_path_for <- function(sample_id) {
  file.path(query_root, paste0(sample_id, "_whole_tissue_Res5.0.rds"))
}

task_output_paths <- function(condition, sample_id, reference) {
  out_dir <- file.path(diagnostic_root, condition, tolower(reference), "tables")
  prefix <- paste(sample_id, reference, condition, sep = "_")
  c(
    cluster_summary = file.path(out_dir, paste0(prefix, "_cluster_summary.csv")),
    target_cells = file.path(out_dir, paste0(prefix, "_target_cell_scores.csv")),
    anchor_composition = file.path(out_dir, paste0(prefix, "_anchor_composition.csv")),
    anchor_coverage = file.path(out_dir, paste0(prefix, "_anchor_coverage.csv")),
    mapping_summary = file.path(out_dir, paste0(prefix, "_mapping_summary.csv")),
    reference_balance = file.path(out_dir, paste0(prefix, "_reference_class_balance.csv")),
    reference_strata = file.path(out_dir, paste0(prefix, "_reference_strata.csv")),
    feature_coverage = file.path(out_dir, paste0(prefix, "_feature_coverage.csv")),
    provenance = file.path(out_dir, paste0(prefix, "_provenance.csv"))
  )
}

all_task_outputs <- do.call(rbind, lapply(seq_len(nrow(task_map)), function(i) {
  paths <- task_output_paths(
    task_map$condition[[i]], task_map$sample_id[[i]], task_map$reference[[i]]
  )
  data.frame(
    task_id = task_map$task_id[[i]], output_type = names(paths),
    output_path = unname(paths), stringsAsFactors = FALSE
  )
}))

combined_dir <- file.path(diagnostic_root, "combined_tables")
combined_outputs <- c(
  cluster_summary = file.path(combined_dir, "diagnostic_cluster_summary.csv"),
  target_cells_long = file.path(combined_dir, "diagnostic_target_cells_long.csv"),
  target_cells_crossref = file.path(combined_dir, "diagnostic_target_cells_cross_reference.csv"),
  target_cluster_crossref = file.path(combined_dir, "diagnostic_target_cluster_cross_reference_summary.csv"),
  anchor_composition = file.path(combined_dir, "diagnostic_anchor_composition.csv"),
  anchor_coverage = file.path(combined_dir, "diagnostic_anchor_coverage.csv"),
  mapping_summary = file.path(combined_dir, "diagnostic_mapping_summary.csv"),
  reference_balance = file.path(combined_dir, "diagnostic_reference_class_balance.csv"),
  reference_strata = file.path(combined_dir, "diagnostic_reference_strata.csv"),
  feature_coverage = file.path(combined_dir, "diagnostic_feature_coverage.csv"),
  provenance = file.path(combined_dir, "diagnostic_combined_provenance.csv")
)

rbind_fill <- function(tables) {
  all_names <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(x) {
    missing <- setdiff(all_names, names(x))
    for (column in missing) x[[column]] <- NA
    x[all_names]
  })
  do.call(rbind, tables)
}

read_output_type <- function(output_type) {
  paths <- all_task_outputs$output_path[all_task_outputs$output_type == output_type]
  rbind_fill(lapply(paths, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
}

consensus_label <- function(values) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values) & tolower(values) != "unknown"]
  if (!length(values)) return("Unknown")
  counts <- table(values)
  winners <- names(counts)[counts == max(counts)]
  if (max(counts) < 2L || length(winners) != 1L) "Unknown" else winners[[1]]
}

if (combine) {
  if (length(positional_args)) stop("--combine does not accept TASK_ID.")
  required_inputs <- all_task_outputs$output_path
  if (dry_run) {
    ok <- compact_dry_run(
      "Combine 27 transfer-diagnostic tasks",
      inputs = required_inputs,
      outputs = combined_outputs,
      checks = c(expected_outputs = length(required_inputs) == 27L * 9L)
    )
    quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
  }
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Missing diagnostic task outputs:\n- ", paste(missing_inputs, collapse = "\n- "))
  }
  existing_outputs <- combined_outputs[file.exists(combined_outputs)]
  if (length(existing_outputs) && !overwrite) {
    stop(
      "Refusing to overwrite combined diagnostic outputs:\n- ",
      paste(existing_outputs, collapse = "\n- "),
      "\nUse --overwrite only after reviewing them."
    )
  }

  cluster_summary <- read_output_type("cluster_summary")
  target_long <- read_output_type("target_cells")
  key <- c(
    "sample_id", "condition", "k_weight", "seurat_cluster", "cell_id",
    "diagnostic_question", "selection_basis"
  )
  if (anyDuplicated(target_long[c(key, "reference")])) {
    stop("Target-cell outputs duplicate a sample/condition/k.weight/cluster/cell/reference key.")
  }

  standard_cell_columns <- c(
    key, "reference", "predicted_id", "thresholded_predicted_id",
    "prediction_score_max", "top2_cell_score_margin", "normalized_entropy",
    "mapping_score", "score_RL", "score_Granule", "score_Glia", "score_VZ"
  )
  target_long <- target_long[standard_cell_columns]
  reference_wide <- lapply(references, function(reference) {
    x <- target_long[target_long$reference == reference, , drop = FALSE]
    x$reference <- NULL
    names(x)[!names(x) %in% key] <- paste0(
      tolower(reference), "_", names(x)[!names(x) %in% key]
    )
    x
  })
  target_wide <- Reduce(
    function(x, y) merge(x, y, by = key, all = TRUE, sort = FALSE),
    reference_wide
  )
  vote_columns <- paste0(tolower(references), "_thresholded_predicted_id")
  target_wide$cell_consensus_label <- apply(
    target_wide[vote_columns], 1, consensus_label
  )
  target_wide$cell_nonunknown_reference_n <- apply(
    target_wide[vote_columns], 1,
    function(v) sum(!is.na(v) & nzchar(v) & tolower(v) != "unknown")
  )
  target_wide$cell_distinct_nonunknown_label_n <- apply(
    target_wide[vote_columns], 1,
    function(v) {
      v <- v[!is.na(v) & nzchar(v) & tolower(v) != "unknown"]
      length(unique(v))
    }
  )
  target_wide$all_three_references_agree <- apply(
    target_wide[vote_columns], 1,
    function(v) all(!is.na(v)) && length(unique(v)) == 1L
  )

  cluster_crossref <- do.call(rbind, lapply(
    split(
      target_wide,
      list(
        target_wide$sample_id, target_wide$condition,
        target_wide$k_weight, target_wide$seurat_cluster
      ),
      drop = TRUE
    ),
    function(x) {
      consensus_counts <- sort(table(x$cell_consensus_label), decreasing = TRUE)
      data.frame(
        sample_id = x$sample_id[[1]], condition = x$condition[[1]],
        k_weight = x$k_weight[[1]], seurat_cluster = x$seurat_cluster[[1]],
        n_cells = nrow(x),
        modal_cell_consensus_label = names(consensus_counts)[[1]],
        modal_cell_consensus_fraction = max(consensus_counts) / nrow(x),
        all_three_agree_fraction = mean(x$all_three_references_agree),
        unresolved_cell_fraction = mean(x$cell_consensus_label == "Unknown"),
        mean_distinct_nonunknown_labels = mean(x$cell_distinct_nonunknown_label_n),
        stringsAsFactors = FALSE
      )
    }
  ))

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(cluster_summary, combined_outputs[["cluster_summary"]], row.names = FALSE, na = "")
  write.csv(target_long, combined_outputs[["target_cells_long"]], row.names = FALSE, na = "")
  write.csv(target_wide, combined_outputs[["target_cells_crossref"]], row.names = FALSE, na = "")
  write.csv(cluster_crossref, combined_outputs[["target_cluster_crossref"]], row.names = FALSE, na = "")
  for (output_type in c(
    "anchor_composition", "anchor_coverage", "mapping_summary",
    "reference_balance", "reference_strata", "feature_coverage"
  )) {
    write.csv(
      read_output_type(output_type), combined_outputs[[output_type]],
      row.names = FALSE, na = ""
    )
  }
  provenance <- data.frame(
    field = c(
      "purpose", "samples", "references", "conditions", "k_anchor",
      "k_score", "k_weight_values", "threshold", "normalization_scale_factor",
      "target_manifest", "mapping_score_scope", "production_effect"
    ),
    value = c(
      "minimal transfer-source diagnostic", paste(requested_samples, collapse = "|"),
      paste(references, collapse = "|"), paste(conditions$condition, collapse = "|"),
      k_anchor, k_score, paste(k_weight_values, collapse = "|"), legacy_threshold,
      normalization_scale_factor, target_manifest_path, "RPCA conditions only",
      "none; isolated tables only"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(provenance, combined_outputs[["provenance"]], row.names = FALSE, na = "")
  message("Wrote combined transfer diagnostics beneath: ", combined_dir)
  quit(save = "no", status = 0L)
}

if (length(positional_args) > 1L) {
  stop("Use [--dry-run|--overwrite] TASK_ID, --list, or --combine.")
}
task_value <- if (length(positional_args)) positional_args[[1]] else {
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
}
task_id <- suppressWarnings(as.integer(task_value))
if (is.na(task_id) || task_id < 1L || task_id > nrow(task_map)) {
  stop("TASK_ID must be between 1 and ", nrow(task_map), ". Use --list.")
}

task <- task_map[task_map$task_id == task_id, , drop = FALSE]
condition <- task$condition[[1]]
sample_id <- task$sample_id[[1]]
reference_name <- task$reference[[1]]
reduction <- task$reduction[[1]]
renormalize_10k <- isTRUE(task$renormalize_10k[[1]])
query_path <- query_path_for(sample_id)
reference_path <- unname(reference_paths[[reference_name]])
output_paths <- task_output_paths(condition, sample_id, reference_name)

if (dry_run) {
  target_rows <- targets[targets$sample_id == sample_id, , drop = FALSE]
  ok <- compact_dry_run(
    paste0(
      "Transfer diagnostic task ", task_id, "/27 [", condition, " | ",
      sample_id, " | ", reference_name, "]"
    ),
    inputs = c(query_path, reference_path, target_manifest_path),
    outputs = output_paths,
    checks = c(
      sample_matches_manifest = sum(sample_manifest$sample_id == sample_id) == 1L,
      target_clusters_declared = nrow(target_rows) > 0L,
      fixed_k_anchor_10 = k_anchor == 10L,
      fixed_k_score_30 = k_score == 30L
    )
  )
  quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
}

missing_inputs <- c(query_path, reference_path)[!file.exists(c(query_path, reference_path))]
if (length(missing_inputs)) {
  stop("Missing diagnostic input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite diagnostic outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

set.seed(as.integer(config$runtime$random_seed))
query <- readRDS(query_path)
reference <- readRDS(reference_path)
if (anyDuplicated(Cells(query)) || anyDuplicated(Cells(reference))) {
  stop("Query or reference contains duplicate cell IDs.")
}
if (!"clusters_refined" %in% colnames(reference[[]])) {
  stop(reference_name, " reference lacks clusters_refined metadata.")
}
reference_labels <- trimws(as.character(reference$clusters_refined))
if (anyNA(reference_labels) || any(!nzchar(reference_labels))) {
  stop(reference_name, " reference has blank or missing clusters_refined labels.")
}
reference$clusters_refined <- reference_labels

required_query_metadata <- c(
  "whole_tissue_cluster_res5.0", "Xenium_snn_res.5", "seurat_clusters"
)
missing_query_metadata <- setdiff(required_query_metadata, colnames(query[[]]))
if (length(missing_query_metadata)) {
  stop("Query lacks metadata: ", paste(missing_query_metadata, collapse = ", "))
}
cluster_ids <- as.character(query$seurat_clusters)
if (anyNA(cluster_ids) ||
    !identical(cluster_ids, as.character(query$whole_tissue_cluster_res5.0))) {
  stop("Query seurat_clusters does not exactly match resolution 5 metadata.")
}
sample_targets <- targets[targets$sample_id == sample_id, , drop = FALSE]
missing_targets <- setdiff(sample_targets$seurat_cluster, unique(cluster_ids))
if (length(missing_targets)) {
  stop(sample_id, " lacks declared diagnostic clusters: ", paste(missing_targets, collapse = ", "))
}

DefaultAssay(reference) <- if ("RNA" %in% Assays(reference)) "RNA" else Assays(reference)[[1]]
DefaultAssay(query) <- if ("Xenium" %in% Assays(query)) "Xenium" else Assays(query)[[1]]

if (renormalize_10k) {
  reference <- NormalizeData(
    reference, normalization.method = "LogNormalize",
    scale.factor = normalization_scale_factor, verbose = FALSE
  )
  query <- NormalizeData(
    query, normalization.method = "LogNormalize",
    scale.factor = normalization_scale_factor, verbose = FALSE
  )
}

reference <- FindVariableFeatures(
  reference, selection.method = "vst",
  nfeatures = as.integer(config$label_transfer$variable_features), verbose = FALSE
)
shared_genes <- intersect(rownames(reference), rownames(query))
Idents(reference) <- "clusters_refined"
set.seed(as.integer(config$runtime$random_seed))
reference_balanced <- subset(
  reference,
  downsample = as.integer(config$label_transfer$reference_cells_per_identity)
)
transfer_features <- intersect(VariableFeatures(reference_balanced), shared_genes)
if (length(transfer_features) <= max(transfer_dims)) {
  stop("Too few transfer features: ", length(transfer_features))
}
pca_npcs_available <- min(
  50L, length(transfer_features) - 1L,
  ncol(reference_balanced) - 1L, ncol(query) - 1L
)
if (reduction == "rpca" && pca_npcs_available < max(transfer_dims)) {
  stop(
    "RPCA requires at least ", max(transfer_dims),
    " PCs, but only ", pca_npcs_available, " can be computed."
  )
}
pca_npcs <- if (reduction == "rpca") pca_npcs_available else NA_integer_

class_count_table <- function(object, scope) {
  counts <- as.data.frame(table(as.character(object$clusters_refined)), stringsAsFactors = FALSE)
  names(counts) <- c("reference_class", "n_cells")
  counts$scope <- scope
  counts
}
reference_balance <- rbind(
  class_count_table(reference, "full"),
  class_count_table(reference_balanced, "balanced")
)
reference_balance$sample_id <- sample_id
reference_balance$reference <- reference_name
reference_balance$condition <- condition
reference_balance <- reference_balance[c(
  "sample_id", "reference", "condition", "scope", "reference_class", "n_cells"
)]

strata_candidates <- c(
  "figure_clusters", "cell_type", "science_broad_label", "age", "Age",
  "PCW", "Stage", "stage", "donor", "Donor", "sample", "Sample",
  "sample_id", "orig.ident"
)
strata_fields <- intersect(strata_candidates, colnames(reference[[]]))
strata_fields <- strata_fields[vapply(strata_fields, function(field) {
  values <- as.character(reference[[field]][, 1])
  length(unique(values[!is.na(values) & nzchar(values)])) <= 100L
}, logical(1))]

strata_table_for <- function(object, scope) {
  if (!length(strata_fields)) {
    return(data.frame(
      scope = scope, reference_class = NA_character_, stratum_field = NA_character_,
      stratum_value = NA_character_, n_cells = NA_integer_, stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(strata_fields, function(field) {
    md <- object[[]]
    value <- as.character(md[[field]])
    value[is.na(value) | !nzchar(trimws(value))] <- "<missing>"
    tab <- as.data.frame(table(
      reference_class = as.character(md$clusters_refined),
      stratum_value = value
    ), stringsAsFactors = FALSE)
    names(tab)[[3]] <- "n_cells"
    tab <- tab[tab$n_cells > 0L, , drop = FALSE]
    tab$scope <- scope
    tab$stratum_field <- field
    tab[c("scope", "reference_class", "stratum_field", "stratum_value", "n_cells")]
  }))
}
reference_strata <- rbind(
  strata_table_for(reference, "full"),
  strata_table_for(reference_balanced, "balanced")
)
reference_strata$sample_id <- sample_id
reference_strata$reference <- reference_name
reference_strata$condition <- condition
reference_strata <- reference_strata[c(
  "sample_id", "reference", "condition", "scope", "reference_class",
  "stratum_field", "stratum_value", "n_cells"
)]

diagnostic_markers <- list(
  RL = c("MKI67", "LTBP1", "OTX2"),
  Granule = c("ATOH1", "PAX6", "NEUROD1", "RELN"),
  Glia = c("SOX9", "TNC"),
  VZ = c("PRDM13", "ASCL1")
)
feature_coverage <- do.call(rbind, lapply(names(diagnostic_markers), function(class_name) {
  genes <- diagnostic_markers[[class_name]]
  data.frame(
    marker_class = class_name, gene = genes,
    in_reference = genes %in% rownames(reference),
    in_query = genes %in% rownames(query),
    shared = genes %in% shared_genes,
    used_for_transfer = genes %in% transfer_features,
    stringsAsFactors = FALSE
  )
}))
feature_coverage$sample_id <- sample_id
feature_coverage$reference <- reference_name
feature_coverage$condition <- condition
feature_coverage <- feature_coverage[c(
  "sample_id", "reference", "condition", "marker_class", "gene",
  "in_reference", "in_query", "shared", "used_for_transfer"
)]

reference_balanced <- ScaleData(
  reference_balanced, features = transfer_features, verbose = FALSE
)
query <- ScaleData(query, features = transfer_features, verbose = FALSE)
if (reduction == "rpca") {
  reference_balanced <- RunPCA(
    reference_balanced, features = transfer_features, npcs = pca_npcs, verbose = FALSE
  )
  query <- RunPCA(query, features = transfer_features, npcs = pca_npcs, verbose = FALSE)
}

mapping_k <- if (reduction == "rpca") 100L else NULL
anchors <- FindTransferAnchors(
  reference = reference_balanced, query = query,
  normalization.method = "LogNormalize", reduction = reduction,
  features = transfer_features, dims = transfer_dims,
  k.anchor = k_anchor, k.score = k_score, approx.pca = TRUE,
  mapping.score.k = mapping_k
)

anchor_matrix <- as.data.frame(methods::slot(anchors, "anchors"))
required_anchor_columns <- c("cell1", "cell2", "score")
if (!all(required_anchor_columns %in% names(anchor_matrix)) || !nrow(anchor_matrix)) {
  stop("AnchorSet lacks a non-empty cell1/cell2/score matrix.")
}
reference_lookup <- methods::slot(anchors, "reference.cells")
query_lookup <- methods::slot(anchors, "query.cells")
strip_suffix <- function(values, suffix) sub(paste0("_", suffix, "$"), "", values)
if (!all(reference_lookup %in% Cells(reference_balanced))) {
  reference_lookup <- strip_suffix(reference_lookup, "reference")
}
if (!all(query_lookup %in% Cells(query))) query_lookup <- strip_suffix(query_lookup, "query")
if (!all(reference_lookup %in% Cells(reference_balanced)) || !all(query_lookup %in% Cells(query))) {
  stop("Could not align AnchorSet cell lookups to reference and query cells.")
}
reference_index <- as.integer(anchor_matrix$cell1)
query_index <- as.integer(anchor_matrix$cell2)
if (anyNA(reference_index) || anyNA(query_index) ||
    any(reference_index < 1L | reference_index > length(reference_lookup)) ||
    any(query_index < 1L | query_index > length(query_lookup))) {
  stop("AnchorSet contains invalid reference or query cell indices.")
}
anchor_matrix$reference_cell <- reference_lookup[reference_index]
anchor_matrix$query_cell <- query_lookup[query_index]
reference_label_lookup <- setNames(
  as.character(reference_balanced$clusters_refined), Cells(reference_balanced)
)
query_cluster_lookup <- setNames(cluster_ids, Cells(query))
anchor_matrix$reference_class <- reference_label_lookup[anchor_matrix$reference_cell]
anchor_matrix$seurat_cluster <- query_cluster_lookup[anchor_matrix$query_cell]
if (anyNA(anchor_matrix$reference_class) || anyNA(anchor_matrix$seurat_cluster)) {
  stop("Anchor cells did not map cleanly to reference labels and query clusters.")
}

anchor_composition <- do.call(rbind, lapply(
  split(anchor_matrix, list(anchor_matrix$seurat_cluster, anchor_matrix$reference_class), drop = TRUE),
  function(x) data.frame(
    seurat_cluster = x$seurat_cluster[[1]], reference_class = x$reference_class[[1]],
    n_anchors = nrow(x), n_query_anchor_cells = length(unique(x$query_cell)),
    mean_anchor_score = mean(x$score), median_anchor_score = median(x$score),
    stringsAsFactors = FALSE
  )
))
anchor_composition$sample_id <- sample_id
anchor_composition$reference <- reference_name
anchor_composition$condition <- condition
anchor_composition <- anchor_composition[c(
  "sample_id", "reference", "condition", "seurat_cluster", "reference_class",
  "n_anchors", "n_query_anchor_cells", "mean_anchor_score", "median_anchor_score"
)]

anchor_coverage <- do.call(rbind, lapply(split(anchor_matrix, anchor_matrix$seurat_cluster), function(x) {
  cluster_id <- x$seurat_cluster[[1]]
  cluster_n <- sum(cluster_ids == cluster_id)
  data.frame(
    seurat_cluster = cluster_id, n_cells = cluster_n, n_anchors = nrow(x),
    n_query_anchor_cells = length(unique(x$query_cell)),
    query_anchor_cell_fraction = length(unique(x$query_cell)) / cluster_n,
    mean_anchor_score = mean(x$score), median_anchor_score = median(x$score),
    stringsAsFactors = FALSE
  )
}))
anchor_coverage$sample_id <- sample_id
anchor_coverage$reference <- reference_name
anchor_coverage$condition <- condition
anchor_coverage <- anchor_coverage[c(
  "sample_id", "reference", "condition", "seurat_cluster", "n_cells",
  "n_anchors", "n_query_anchor_cells", "query_anchor_cell_fraction",
  "mean_anchor_score", "median_anchor_score"
)]

mapping_error <- ""
mapping_scores <- tryCatch(
  {
    if (reduction != "rpca") stop("MappingScore is not defined for the CCA diagnostic condition.")
    mapping_anchors <- anchors
    mapping_objects <- methods::slot(mapping_anchors, "object.list")
    if (!"pcaproject.l2" %in% Reductions(mapping_objects[[1]])) {
      mapping_objects[[1]] <- L2Dim(
        mapping_objects[[1]], reduction = "pcaproject"
      )
      methods::slot(mapping_anchors, "object.list") <- mapping_objects
    }
    Seurat::MappingScore(
      mapping_anchors, kanchors = min(50L, nrow(anchor_matrix)),
      ndim = max(transfer_dims), verbose = FALSE
    )
  },
  error = function(e) {
    mapping_error <<- conditionMessage(e)
    NULL
  }
)
mapping_by_cell <- setNames(rep(NA_real_, ncol(query)), Cells(query))
if (!is.null(mapping_scores)) {
  score_names <- names(mapping_scores)
  if (is.null(score_names) && length(mapping_scores) == length(query_lookup)) {
    score_names <- query_lookup
  } else if (!is.null(score_names) && !all(score_names %in% Cells(query))) {
    score_names <- strip_suffix(score_names, "query")
  }
  if (length(mapping_scores) != length(score_names) || !all(score_names %in% Cells(query))) {
    stop("MappingScore results could not be aligned to query cells.")
  }
  mapping_by_cell[score_names] <- as.numeric(mapping_scores)
}

focus_classes <- c("RL", "Granule", "Glia", "VZ")
cluster_summaries <- list()
target_cell_tables <- list()

for (k_weight in k_weight_values) {
  predictions <- TransferData(
    anchorset = anchors, refdata = reference_balanced$clusters_refined,
    dims = transfer_dims, k.weight = k_weight,
    weight.reduction = if (reduction == "cca") "cca" else "pcaproject",
    store.weights = FALSE
  )
  predictions <- as.data.frame(predictions, check.names = FALSE)
  if (!identical(rownames(predictions), Cells(query))) {
    if (!setequal(rownames(predictions), Cells(query))) {
      stop("Transfer predictions do not align one-to-one with query cells.")
    }
    predictions <- predictions[Cells(query), , drop = FALSE]
  }
  score_columns <- grep("^prediction\\.score\\.", names(predictions), value = TRUE)
  score_columns <- setdiff(score_columns, c("prediction.score.max", "prediction.score.id"))
  if (length(score_columns) < 2L) {
    stop("TransferData returned fewer than two class-score columns.")
  }
  score_matrix <- as.matrix(predictions[score_columns])
  storage.mode(score_matrix) <- "double"
  class_labels <- sub("^prediction\\.score\\.", "", score_columns)
  colnames(score_matrix) <- class_labels
  if (any(!is.finite(score_matrix))) stop("Predictions contain non-finite class scores.")
  row_totals <- rowSums(score_matrix)
  probability_matrix <- score_matrix / pmax(row_totals, .Machine$double.eps)
  normalized_entropy <- -rowSums(
    ifelse(probability_matrix > 0, probability_matrix * log(probability_matrix), 0)
  ) / log(ncol(probability_matrix))
  sorted_scores <- t(apply(score_matrix, 1, sort, decreasing = TRUE))
  top2_cell_margin <- sorted_scores[, 1] - sorted_scores[, 2]
  predicted_id <- trimws(as.character(predictions$predicted.id))
  prediction_score_max <- as.numeric(predictions$prediction.score.max)
  thresholded_id <- ifelse(
    prediction_score_max >= legacy_threshold, predicted_id, "Unknown"
  )

  cluster_rows <- lapply(unique(cluster_ids), function(cluster_id) {
    index <- which(cluster_ids == cluster_id)
    mean_scores <- colMeans(score_matrix[index, , drop = FALSE])
    score_order <- order(mean_scores, decreasing = TRUE)
    winner <- names(mean_scores)[score_order[[1]]]
    runner_up <- names(mean_scores)[score_order[[2]]]
    majority_counts <- sort(table(predicted_id[index]), decreasing = TRUE)
    majority <- names(majority_counts)[[1]]
    row <- data.frame(
      sample_id = sample_id, reference = reference_name, condition = condition,
      renormalize_10k = renormalize_10k, reduction = reduction,
      k_anchor = k_anchor, k_score = k_score, k_weight = k_weight,
      seurat_cluster = cluster_id, n_cells = length(index),
      cluster_majority = if (mean(prediction_score_max[index]) >= legacy_threshold) majority else "Unknown",
      cluster_weighted = if (mean_scores[[winner]] >= legacy_threshold) winner else "Unknown",
      weighted_winner_before_threshold = winner,
      weighted_winner_mean_score = mean_scores[[winner]],
      runner_up_mean_score_label = runner_up,
      runner_up_mean_score = mean_scores[[runner_up]],
      top2_mean_score_margin = mean_scores[[winner]] - mean_scores[[runner_up]],
      mean_prediction_score_max = mean(prediction_score_max[index]),
      mean_normalized_entropy = mean(normalized_entropy[index]),
      mean_top2_cell_score_margin = mean(top2_cell_margin[index]),
      mean_mapping_score = mean(mapping_by_cell[Cells(query)[index]], na.rm = TRUE),
      median_mapping_score = median(mapping_by_cell[Cells(query)[index]], na.rm = TRUE),
      fraction_cells_above_legacy_threshold = mean(prediction_score_max[index] >= legacy_threshold),
      legacy_threshold = legacy_threshold,
      stringsAsFactors = FALSE
    )
    if (!any(is.finite(mapping_by_cell[Cells(query)[index]]))) {
      row$mean_mapping_score <- NA_real_
      row$median_mapping_score <- NA_real_
    }
    for (class_name in focus_classes) {
      row[[paste0("mean_score_", class_name)]] <- if (class_name %in% names(mean_scores)) {
        mean_scores[[class_name]]
      } else NA_real_
      row[[paste0("fraction_predicted_", class_name)]] <- mean(predicted_id[index] == class_name)
    }
    row
  })
  cluster_summaries[[as.character(k_weight)]] <- do.call(rbind, cluster_rows)

  target_index <- cluster_ids %in% sample_targets$seurat_cluster
  target_table <- data.frame(
    sample_id = sample_id, reference = reference_name, condition = condition,
    renormalize_10k = renormalize_10k, reduction = reduction,
    k_anchor = k_anchor, k_score = k_score, k_weight = k_weight,
    seurat_cluster = cluster_ids[target_index], cell_id = Cells(query)[target_index],
    predicted_id = predicted_id[target_index],
    thresholded_predicted_id = thresholded_id[target_index],
    prediction_score_max = prediction_score_max[target_index],
    top2_cell_score_margin = top2_cell_margin[target_index],
    normalized_entropy = normalized_entropy[target_index],
    mapping_score = unname(mapping_by_cell[Cells(query)[target_index]]),
    stringsAsFactors = FALSE
  )
  for (class_name in focus_classes) {
    target_table[[paste0("score_", class_name)]] <- if (class_name %in% colnames(score_matrix)) {
      score_matrix[target_index, class_name]
    } else NA_real_
  }
  target_table <- merge(
    target_table, sample_targets,
    by = c("sample_id", "seurat_cluster"), all.x = TRUE, sort = FALSE
  )
  target_cell_tables[[as.character(k_weight)]] <- target_table
}

cluster_summary <- do.call(rbind, cluster_summaries)
target_cells <- do.call(rbind, target_cell_tables)
mapping_summary <- do.call(rbind, lapply(unique(cluster_ids), function(cluster_id) {
  values <- mapping_by_cell[Cells(query)[cluster_ids == cluster_id]]
  data.frame(
    sample_id = sample_id, reference = reference_name, condition = condition,
    seurat_cluster = cluster_id, n_cells = length(values),
    mapping_available = any(is.finite(values)),
    mean_mapping_score = if (any(is.finite(values))) mean(values, na.rm = TRUE) else NA_real_,
    median_mapping_score = if (any(is.finite(values))) median(values, na.rm = TRUE) else NA_real_,
    q10_mapping_score = if (any(is.finite(values))) quantile(values, 0.10, na.rm = TRUE) else NA_real_,
    q90_mapping_score = if (any(is.finite(values))) quantile(values, 0.90, na.rm = TRUE) else NA_real_,
    mapping_error = mapping_error,
    stringsAsFactors = FALSE
  )
}))

provenance <- data.frame(
  field = c(
    "task_id", "sample_id", "reference", "condition", "query_path",
    "reference_path", "renormalize_10k", "normalization_scale_factor",
    "reduction", "random_seed", "query_cells", "full_reference_cells",
    "balanced_reference_cells", "shared_genes", "transfer_features",
    "pca_npcs", "dimensions", "k_anchor", "k_score", "k_weight_values", "threshold",
    "anchor_count", "mapping_score_available", "mapping_error",
    "reference_strata_fields", "production_effect"
  ),
  value = c(
    task_id, sample_id, reference_name, condition, query_path, reference_path,
    renormalize_10k, normalization_scale_factor, reduction,
    config$runtime$random_seed, ncol(query), ncol(reference),
    ncol(reference_balanced), length(shared_genes), length(transfer_features),
    pca_npcs, paste(transfer_dims, collapse = ","), k_anchor, k_score,
    paste(k_weight_values, collapse = "|"), legacy_threshold,
    nrow(anchor_matrix), any(is.finite(mapping_by_cell)), mapping_error,
    paste(strata_fields, collapse = "|"), "none; isolated tables only"
  ),
  stringsAsFactors = FALSE
)

dir.create(dirname(output_paths[["cluster_summary"]]), recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_summary, output_paths[["cluster_summary"]], row.names = FALSE, na = "")
write.csv(target_cells, output_paths[["target_cells"]], row.names = FALSE, na = "")
write.csv(anchor_composition, output_paths[["anchor_composition"]], row.names = FALSE, na = "")
write.csv(anchor_coverage, output_paths[["anchor_coverage"]], row.names = FALSE, na = "")
write.csv(mapping_summary, output_paths[["mapping_summary"]], row.names = FALSE, na = "")
write.csv(reference_balance, output_paths[["reference_balance"]], row.names = FALSE, na = "")
write.csv(reference_strata, output_paths[["reference_strata"]], row.names = FALSE, na = "")
write.csv(feature_coverage, output_paths[["feature_coverage"]], row.names = FALSE, na = "")
write.csv(provenance, output_paths[["provenance"]], row.names = FALSE, na = "")
message("Wrote isolated transfer diagnostics for task ", task_id, "/27.")
