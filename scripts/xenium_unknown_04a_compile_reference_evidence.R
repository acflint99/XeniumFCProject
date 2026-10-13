#!/usr/bin/env Rscript

# Compile one review row per sample-qualified Unknown population by joining the
# existing broad route, continuous reference-transfer score audits, and the
# whole-tissue top-20 marker workbooks. This diagnostic never modifies Seurat
# objects, labels, routes, thresholds, or QC decisions.

rm(list = ls())

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
production_paths <- get_production_output_paths(config)
route_paths <- unknown_reference_route_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[!args %in% valid_options]
if (length(unknown_options)) {
  stop(
    "Usage: Rscript scripts/xenium_unknown_04a_compile_reference_evidence.R ",
    "[--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "The 'readxl' package is required to inspect Unknown marker workbooks. ",
    "Restore it through the project renv environment; do not install it from ",
    "this production script."
  )
}

reference_keys <- unknown_reference_names()
reference_titles <- c(
  aldinger = "Aldinger",
  sepp = "Sepp",
  science = "Science"
)
if (!identical(names(reference_titles), reference_keys)) {
  stop("Internal reference-name mapping does not match the routing references.")
}

score_threshold <- as.numeric(config$label_transfer$prediction_score_threshold)
margin_threshold <- as.numeric(
  config$label_transfer$resolution5_all_samples$min_top2_mean_score_margin
)
if (!is.finite(score_threshold) || !is.finite(margin_threshold)) {
  stop("Configured reference score or margin threshold is not finite.")
}

routes_path <- route_paths$route_manifest
if (!file.exists(routes_path) && !dry_run) {
  stop("Unknown reference-route manifest not found: ", routes_path)
}

routes_for_paths <- if (file.exists(routes_path)) {
  read.csv(
    routes_path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
} else {
  data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
}
route_sample_ids <- unique(trimws(as.character(routes_for_paths$sample_id)))
if (anyNA(route_sample_ids) || any(!nzchar(route_sample_ids))) {
  stop("Reference-route manifest contains a blank or missing sample_id.")
}
if (any(!route_sample_ids %in% sample_ids)) {
  stop("Reference-route manifest contains sample IDs absent from config/samples.csv.")
}
route_sample_ids <- sample_ids[sample_ids %in% route_sample_ids]

reference_table_path <- function(sample_id, reference_key) {
  file.path(
    production_paths$label_transfer_dir,
    reference_key,
    "tables",
    paste0(
      sample_id, "_", reference_titles[[reference_key]],
      "_majority_vs_weighted.csv"
    )
  )
}

marker_table_dir <- file.path(
  production_paths$annotation_root,
  config$production$consensus_label_stage,
  "tables", "unknown_markers"
)
marker_workbook_path <- function(sample_id) {
  file.path(
    marker_table_dir,
    paste0(sample_id, "_Unknown_cluster_top20_markers.xlsx")
  )
}

reference_paths <- unlist(lapply(route_sample_ids, function(sample_id) {
  vapply(
    reference_keys,
    function(reference_key) reference_table_path(sample_id, reference_key),
    character(1)
  )
}), use.names = FALSE)
marker_paths <- vapply(route_sample_ids, marker_workbook_path, character(1))

output_dir <- file.path(route_paths$table_dir)
output_path <- file.path(
  output_dir,
  "Xenium_Unknown_Reference_Score_Marker_Evidence.csv"
)

input_paths <- c(routes_path, reference_paths, marker_paths)
if (dry_run) {
  compact_dry_run(
    "Compile Unknown reference-score and marker evidence",
    inputs = input_paths,
    outputs = output_path,
    checks = c(
      `readxl package available` = requireNamespace("readxl", quietly = TRUE),
      `score threshold finite` = is.finite(score_threshold),
      `margin threshold finite` = is.finite(margin_threshold)
    )
  )
  quit(save = "no", status = 0L)
}

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing evidence input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
if (file.exists(output_path) && !overwrite) {
  stop(
    "Refusing to overwrite existing evidence table:\n", output_path,
    "\nUse --overwrite only after reviewing the existing file."
  )
}

routes <- load_unknown_reference_routes(config)
routes$sample_id <- trimws(as.character(routes$sample_id))
routes$seurat_clusters <- trimws(as.character(routes$seurat_clusters))
routes$consensus_label <- trimws(as.character(routes$consensus_label))
route_keys <- paste(routes$sample_id, routes$consensus_label, sep = "\u001f")
if (anyDuplicated(route_keys)) {
  stop("Reference-route manifest contains duplicate sample/Unknown keys.")
}

reference_required <- c(
  "seurat_clusters", "n_cells", "cluster_weighted",
  "cluster_weighted_before_abstention", "weighted_winner_mean_score",
  "weighted_runner_up_label", "weighted_runner_up_mean_score",
  "weighted_top2_mean_score_margin", "weighted_winner_cell_fraction",
  "weighted_failed_score", "weighted_failed_margin", "weighted_abstained",
  "prediction_score_threshold", "min_top2_mean_score_margin"
)

one_logical <- function(value, field, context) {
  result <- suppressWarnings(as.logical(value))
  if (length(result) != 1L || is.na(result)) {
    stop(context, " has an invalid ", field, " value: ", value)
  }
  result
}

one_finite_number <- function(value, field, context) {
  result <- suppressWarnings(as.numeric(value))
  if (length(result) != 1L || !is.finite(result)) {
    stop(context, " has a non-finite ", field, " value: ", value)
  }
  result
}

format_abstention_reason <- function(failed_score, failed_margin, abstained) {
  if (!abstained) return("passed")
  reasons <- c(
    if (failed_score) "winner_score_below_threshold" else character(),
    if (failed_margin) "top2_margin_below_threshold" else character()
  )
  if (!length(reasons)) {
    stop("An abstained reference has neither a failed-score nor failed-margin flag.")
  }
  paste(reasons, collapse = ";")
}

population_key <- function(table) {
  paste(table$sample_id, table$consensus_label, sep = "\u001f")
}

reference_evidence <- vector("list", length(reference_keys))
names(reference_evidence) <- reference_keys

for (reference_key in reference_keys) {
  rows <- vector("list", nrow(routes))
  table_cache <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(routes))) {
    sample_id <- routes$sample_id[[i]]
    table_path <- reference_table_path(sample_id, reference_key)
    if (!exists(sample_id, envir = table_cache, inherits = FALSE)) {
      table <- read.csv(
        table_path, stringsAsFactors = FALSE, check.names = FALSE,
        na.strings = character()
      )
      missing_columns <- setdiff(reference_required, names(table))
      if (length(missing_columns)) {
        stop(
          "Reference table lacks required column(s): ", table_path, "\n- ",
          paste(missing_columns, collapse = "\n- ")
        )
      }
      table$seurat_clusters <- trimws(as.character(table$seurat_clusters))
      if (anyNA(table$seurat_clusters) || any(!nzchar(table$seurat_clusters)) ||
          anyDuplicated(table$seurat_clusters)) {
        stop("Reference table has invalid or duplicate cluster IDs: ", table_path)
      }
      assign(sample_id, table, envir = table_cache)
    }
    table <- get(sample_id, envir = table_cache, inherits = FALSE)
    cluster_id <- routes$seurat_clusters[[i]]
    index <- match(cluster_id, table$seurat_clusters)
    if (is.na(index)) {
      stop("Cluster ", cluster_id, " for ", sample_id, " is absent from ", table_path)
    }
    row <- table[index, , drop = FALSE]
    context <- paste(sample_id, routes$consensus_label[[i]], reference_key)

    observed_score_threshold <- one_finite_number(
      row$prediction_score_threshold[[1]], "prediction_score_threshold", context
    )
    observed_margin_threshold <- one_finite_number(
      row$min_top2_mean_score_margin[[1]], "min_top2_mean_score_margin", context
    )
    if (!isTRUE(all.equal(observed_score_threshold, score_threshold)) ||
        !isTRUE(all.equal(observed_margin_threshold, margin_threshold))) {
      stop(context, " does not use the active production transfer thresholds.")
    }

    failed_score <- one_logical(
      row$weighted_failed_score[[1]], "weighted_failed_score", context
    )
    failed_margin <- one_logical(
      row$weighted_failed_margin[[1]], "weighted_failed_margin", context
    )
    abstained <- one_logical(
      row$weighted_abstained[[1]], "weighted_abstained", context
    )
    if (!identical(abstained, failed_score || failed_margin)) {
      stop(context, " has internally inconsistent abstention flags.")
    }

    winner_score <- one_finite_number(
      row$weighted_winner_mean_score[[1]], "weighted_winner_mean_score", context
    )
    runner_up_score <- one_finite_number(
      row$weighted_runner_up_mean_score[[1]],
      "weighted_runner_up_mean_score", context
    )
    top2_margin <- one_finite_number(
      row$weighted_top2_mean_score_margin[[1]],
      "weighted_top2_mean_score_margin", context
    )
    if (!isTRUE(all.equal(
      winner_score - runner_up_score, top2_margin, tolerance = 1e-12
    ))) {
      stop(context, " has a top-two margin inconsistent with its class scores.")
    }

    expected_vote <- as.character(
      routes[[paste0(reference_key, "_reference_vote")]][[i]]
    )
    observed_vote <- trimws(as.character(row$cluster_weighted[[1]]))
    if (!identical(observed_vote, expected_vote)) {
      stop(
        context, " reference vote differs between its transfer table ('",
        observed_vote, "') and route manifest ('", expected_vote, "')."
      )
    }

    winner_fraction <- one_finite_number(
      row$weighted_winner_cell_fraction[[1]],
      "weighted_winner_cell_fraction", context
    )
    if (winner_fraction < 0 || winner_fraction > 1) {
      stop(context, " has weighted_winner_cell_fraction outside [0, 1].")
    }
    n_cells <- suppressWarnings(as.integer(row$n_cells[[1]]))
    if (length(n_cells) != 1L || is.na(n_cells) || n_cells <= 0L) {
      stop(context, " has an invalid n_cells value.")
    }

    rows[[i]] <- data.frame(
      sample_id = sample_id,
      consensus_label = routes$consensus_label[[i]],
      n_cells = n_cells,
      winner_before_abstention = trimws(as.character(
        row$cluster_weighted_before_abstention[[1]]
      )),
      reference_vote = observed_vote,
      winner_mean_score = winner_score,
      runner_up_label = trimws(as.character(row$weighted_runner_up_label[[1]])),
      runner_up_mean_score = runner_up_score,
      top2_mean_score_margin = top2_margin,
      winner_cell_fraction = winner_fraction,
      failed_score = failed_score,
      failed_margin = failed_margin,
      abstained = abstained,
      abstention_reason = format_abstention_reason(
        failed_score, failed_margin, abstained
      ),
      score_threshold = score_threshold,
      score_minus_threshold = winner_score - score_threshold,
      margin_threshold = margin_threshold,
      margin_minus_threshold = top2_margin - margin_threshold,
      reference_table_path = normalizePath(table_path, mustWork = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  evidence <- do.call(rbind, rows)
  prefix_columns <- setdiff(names(evidence), c("sample_id", "consensus_label"))
  names(evidence)[match(prefix_columns, names(evidence))] <- paste0(
    reference_key, "_", prefix_columns
  )
  reference_evidence[[reference_key]] <- evidence
}

evidence_table <- routes
for (reference_key in reference_keys) {
  reference_table <- reference_evidence[[reference_key]]
  index <- match(route_keys, population_key(reference_table))
  if (anyNA(index)) stop("Reference evidence is incomplete for ", reference_key, ".")
  add_columns <- setdiff(names(reference_table), c("sample_id", "consensus_label"))
  evidence_table[add_columns] <- reference_table[index, add_columns, drop = FALSE]
}

n_cell_columns <- paste0(reference_keys, "_n_cells")
n_cell_matrix <- as.matrix(evidence_table[n_cell_columns])
storage.mode(n_cell_matrix) <- "integer"
if (any(apply(n_cell_matrix, 1L, function(values) length(unique(values)) != 1L))) {
  stop("Reference tables disagree on cell counts for at least one Unknown population.")
}
evidence_table$n_cells <- n_cell_matrix[, 1L]
evidence_table[n_cell_columns] <- NULL

marker_rows <- vector("list", nrow(routes))
for (sample_id in unique(routes$sample_id)) {
  workbook_path <- marker_workbook_path(sample_id)
  workbook_sheets <- readxl::excel_sheets(workbook_path)
  sample_labels <- routes$consensus_label[routes$sample_id == sample_id]
  if (!setequal(workbook_sheets, sample_labels)) {
    stop(
      "Unknown marker workbook sheets do not match route populations for ",
      sample_id, ". Missing sheets: ",
      paste(setdiff(sample_labels, workbook_sheets), collapse = ", "),
      "; unexpected sheets: ",
      paste(setdiff(workbook_sheets, sample_labels), collapse = ", "), "."
    )
  }

  for (consensus_label in sample_labels) {
    marker_table <- as.data.frame(
      readxl::read_excel(
        workbook_path, sheet = consensus_label, .name_repair = "minimal"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    marker_required <- c(
      "rank", "gene", "unknown_label", "seurat_cluster", "p_val",
      "avg_log2FC", "pct.1", "pct.2", "p_val_adj"
    )
    missing_columns <- setdiff(marker_required, names(marker_table))
    if (length(missing_columns)) {
      stop(
        "Marker sheet lacks required column(s): ", workbook_path, " / ",
        consensus_label, "\n- ", paste(missing_columns, collapse = "\n- ")
      )
    }
    if (!nrow(marker_table)) {
      stop("Marker sheet is empty: ", workbook_path, " / ", consensus_label)
    }
    ranks <- suppressWarnings(as.integer(marker_table$rank))
    if (anyNA(ranks) || anyDuplicated(ranks) ||
        !identical(sort(ranks), seq_len(nrow(marker_table)))) {
      stop("Marker ranks are invalid: ", workbook_path, " / ", consensus_label)
    }
    marker_table <- marker_table[order(ranks), , drop = FALSE]
    marker_labels <- unique(trimws(as.character(marker_table$unknown_label)))
    if (!identical(marker_labels, consensus_label)) {
      stop("Marker sheet label does not match its worksheet: ", workbook_path)
    }
    route_index <- which(
      routes$sample_id == sample_id & routes$consensus_label == consensus_label
    )
    if (length(route_index) != 1L) {
      stop("Could not resolve exactly one route row for marker sheet: ", workbook_path)
    }
    marker_clusters <- unique(trimws(as.character(marker_table$seurat_cluster)))
    if (!identical(marker_clusters, routes$seurat_clusters[[route_index]])) {
      stop("Marker sheet cluster does not match its route row: ", workbook_path)
    }
    genes <- trimws(as.character(marker_table$gene))
    if (anyNA(genes) || any(!nzchar(genes))) {
      stop("Marker sheet contains blank gene names: ", workbook_path)
    }
    adjusted_p <- suppressWarnings(as.numeric(marker_table$p_val_adj))
    log2_fc <- suppressWarnings(as.numeric(marker_table$avg_log2FC))
    pct_1 <- suppressWarnings(as.numeric(marker_table$pct.1))
    pct_2 <- suppressWarnings(as.numeric(marker_table$pct.2))
    if (any(!is.finite(adjusted_p)) || any(!is.finite(log2_fc)) ||
        any(!is.finite(pct_1)) || any(!is.finite(pct_2))) {
      stop("Marker sheet contains non-finite statistics: ", workbook_path)
    }
    significant <- adjusted_p < 0.05
    marker_rows[[route_index]] <- data.frame(
      sample_id = sample_id,
      consensus_label = consensus_label,
      marker_top_n = nrow(marker_table),
      marker_adjusted_significant_n = sum(significant),
      marker_min_p_val_adj = min(adjusted_p),
      marker_max_avg_log2FC = max(log2_fc),
      marker_median_avg_log2FC = stats::median(log2_fc),
      marker_top_genes = paste(genes, collapse = ";"),
      marker_adjusted_significant_genes = paste(
        genes[significant], collapse = ";"
      ),
      marker_workbook_path = normalizePath(workbook_path, mustWork = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

if (any(vapply(marker_rows, is.null, logical(1)))) {
  stop("Marker evidence is incomplete for at least one route population.")
}
marker_evidence <- do.call(rbind, marker_rows)
marker_index <- match(route_keys, population_key(marker_evidence))
if (anyNA(marker_index)) stop("Could not match every route population to marker evidence.")
marker_columns <- setdiff(names(marker_evidence), c("sample_id", "consensus_label"))
evidence_table[marker_columns] <- marker_evidence[
  marker_index, marker_columns, drop = FALSE
]
evidence_table$route_manifest_path <- normalizePath(routes_path, mustWork = TRUE)

identity_columns <- c(
  "sample_id", "seurat_clusters", "consensus_label", "n_cells",
  unknown_reference_vote_columns(), paste0(reference_keys, "_broad_route"),
  "reference_route", "reference_route_support_n",
  "reference_route_nonabstain_n", "reference_route_status",
  "qc_decision", "qc_review_notes"
)
reference_columns <- unlist(lapply(reference_keys, function(reference_key) {
  grep(paste0("^", reference_key, "_"), names(evidence_table), value = TRUE)
}), use.names = FALSE)
marker_columns <- grep("^marker_", names(evidence_table), value = TRUE)
provenance_columns <- c(
  "consensus_table_path", "consensus_table_size", "consensus_table_mtime",
  "route_manifest_path"
)
output_columns <- unique(c(
  identity_columns, reference_columns, marker_columns, provenance_columns
))
if (any(!output_columns %in% names(evidence_table))) {
  stop("Internal evidence table lacks at least one requested output column.")
}
evidence_table <- evidence_table[output_columns]

if (nrow(evidence_table) != nrow(routes) ||
    anyDuplicated(population_key(evidence_table))) {
  stop("Evidence compilation changed population cardinality or created duplicates.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
partial_path <- paste0(output_path, ".", Sys.getpid(), ".part")
on.exit(if (file.exists(partial_path)) unlink(partial_path), add = TRUE)
write.csv(evidence_table, partial_path, row.names = FALSE, na = "")

backup_path <- NULL
if (file.exists(output_path)) {
  backup_path <- paste0(output_path, ".backup_", Sys.getpid())
  if (!file.rename(output_path, backup_path)) {
    stop("Could not preserve the existing evidence table as: ", backup_path)
  }
}
if (!file.rename(partial_path, output_path)) {
  if (!is.null(backup_path) && file.exists(backup_path)) {
    file.rename(backup_path, output_path)
  }
  stop("Could not move the completed evidence table into place: ", output_path)
}
if (!is.null(backup_path) && file.exists(backup_path)) unlink(backup_path)

message(
  "Wrote continuous reference-score and marker evidence for ",
  nrow(evidence_table), " Unknown populations: ", output_path
)
print(
  as.data.frame(
    table(
      evidence_table$reference_route,
      evidence_table$marker_adjusted_significant_n == 0L
    ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)
