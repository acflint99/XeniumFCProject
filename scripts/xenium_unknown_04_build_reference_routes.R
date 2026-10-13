#!/usr/bin/env Rscript

# Build the production VZ/RL routing table for reviewed, sample-qualified
# Unknown populations. This stage now consumes the explicitly reviewed
# broad-reference consensus, whose per-reference VZ/RL/Other scores are gated
# before two-of-three voting. It never modifies a Seurat object or detailed
# consensus label.

rm(list = ls())

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))
source(here("scripts", "R", "xenium_broad_reference.R"))

config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
production_paths <- get_production_output_paths(config)
route_paths <- unknown_reference_route_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) {
  stop(
    "Usage: Rscript scripts/xenium_unknown_04_build_reference_routes.R ",
    "[--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

qc_per_sample_dir <- file.path(
  production_paths$annotation_root, "04_unknown_qc_reviewed", "tables", "per_sample"
)
qc_input_paths <- file.path(
  qc_per_sample_dir, paste0(sample_ids, "_Unknown_QC_routing.csv")
)
consensus_table_paths <- file.path(
  production_paths$consensus_table_dir,
  paste0(sample_ids, "_comparison_merged.csv")
)
broad_consensus_table_paths <- file.path(
  here(config$project$outputs_dir), "xenium", "annotation",
  "05_broad_reference_consensus", "tables", "per_sample",
  paste0(sample_ids, "_Broad_Reference_Consensus.csv")
)
reviewed_qc_path <- resolve_config_path(config$manifests$reviewed_unknown_qc, config)
reviewed_threshold_path <- resolve_config_path(
  config$manifests$reviewed_broad_reference_thresholds, config
)
input_paths <- c(
  qc_input_paths, consensus_table_paths, broad_consensus_table_paths,
  reviewed_qc_path, reviewed_threshold_path
)
output_paths <- unlist(route_paths[c(
  "route_manifest", "input_manifest", "validation", "summary"
)], use.names = FALSE)

if (dry_run) {
  compact_dry_run(
    "Build reviewed Unknown reference routes",
    inputs = input_paths,
    outputs = output_paths,
    checks = c(
      reviewed_broad_thresholds = tryCatch({
        load_reviewed_broad_thresholds(config, require_approved = TRUE)
        TRUE
      }, error = function(e) FALSE),
      reference_support = identical(
        as.integer(config$regional_subsets$unknown_route_min_reference_support),
        2L
      )
    )
  )
  quit(save = "no", status = 0L)
}

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing reference-route input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite reference-route output(s):\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing files."
  )
}

input_manifest <- do.call(rbind, lapply(qc_input_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = character())
}))
inspection <- inspect_unknown_qc_routing_manifest(input_manifest, sample_ids)
if (!inspection$valid) {
  stop("Reviewed input manifest is invalid:\n- ",
       paste(inspection$errors, collapse = "\n- "))
}

recorded_review_path <- normalizePath(
  as.character(input_manifest$reviewed_qc_path), mustWork = TRUE
)
current_review_path <- normalizePath(reviewed_qc_path, mustWork = TRUE)
if (any(recorded_review_path != current_review_path)) {
  stop("Per-sample reviewed inputs were not created from the active reviewed QC manifest.")
}
review_info <- file.info(reviewed_qc_path)
if (any(!vapply(seq_len(nrow(input_manifest)), function(i) {
  isTRUE(all.equal(
    as.numeric(input_manifest$reviewed_qc_size[[i]]),
    as.numeric(review_info$size)
  )) && isTRUE(all.equal(
    as.numeric(input_manifest$reviewed_qc_mtime[[i]]),
    as.numeric(review_info$mtime)
  ))
}, logical(1)))) {
  stop("The reviewed QC manifest changed after per-sample application; rerun stage 03.")
}

review <- load_reviewed_unknown_qc(config)
validate_unknown_qc_review(review, sample_ids, require_complete = TRUE)
review$decision <- tolower(trimws(as.character(review$decision)))

all_route_tables <- vector("list", length(sample_ids))
for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  table_path <- consensus_table_paths[[i]]
  table <- read.csv(
    table_path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  broad_path <- broad_consensus_table_paths[[i]]
  broad <- read.csv(
    broad_path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = character()
  )
  broad_required <- c(
    "sample_id", "seurat_clusters", "consensus_label",
    paste0(unknown_reference_names(), "_broad_vote"),
    "broad_reference_consensus", "broad_reference_route",
    "broad_reference_support_n", "broad_reference_nonabstain_n",
    "broad_reference_status"
  )
  broad_missing <- setdiff(broad_required, names(broad))
  if (length(broad_missing)) {
    stop(broad_path, " lacks column(s): ", paste(broad_missing, collapse = ", "))
  }
  table$seurat_clusters <- trimws(as.character(table$seurat_clusters))
  broad$seurat_clusters <- trimws(as.character(broad$seurat_clusters))
  if (anyDuplicated(table$seurat_clusters) || anyDuplicated(broad$seurat_clusters) ||
      !setequal(table$seurat_clusters, broad$seurat_clusters)) {
    stop("Detailed and broad consensus cluster sets differ for ", sample_id, ".")
  }
  broad_index <- match(table$seurat_clusters, broad$seurat_clusters)
  if (!identical(
    trimws(as.character(table$consensus_label)),
    trimws(as.character(broad$consensus_label[broad_index]))
  )) {
    stop("Detailed consensus labels changed before broad routing for ", sample_id, ".")
  }
  for (reference_key in unknown_reference_names()) {
    broad_vote <- trimws(as.character(
      broad[[paste0(reference_key, "_broad_vote")]][broad_index]
    ))
    table[[paste0(reference_key, "_calibrated_broad_vote")]] <- broad_vote
    broad_route <- rep("abstain", length(broad_vote))
    nonabstain <- broad_vote != "Abstain"
    broad_route[nonabstain] <- broad_category_to_route(broad_vote[nonabstain])
    table[[paste0(reference_key, "_broad_route")]] <- broad_route
  }
  table$broad_reference_consensus <- broad$broad_reference_consensus[broad_index]
  table$reference_route <- broad$broad_reference_route[broad_index]
  table$reference_route_support_n <- as.integer(
    broad$broad_reference_support_n[broad_index]
  )
  table$reference_route_nonabstain_n <- as.integer(
    broad$broad_reference_nonabstain_n[broad_index]
  )
  table$reference_route_status <- broad$broad_reference_status[broad_index]
  table$sample_id <- sample_id
  table$consensus_table_path <- normalizePath(table_path, mustWork = TRUE)
  table_info <- file.info(table_path)
  table$consensus_table_size <- as.numeric(table_info$size)
  table$consensus_table_mtime <- as.numeric(table_info$mtime)
  table$broad_consensus_table_path <- normalizePath(broad_path, mustWork = TRUE)
  broad_info <- file.info(broad_path)
  table$broad_consensus_table_size <- as.numeric(broad_info$size)
  table$broad_consensus_table_mtime <- as.numeric(broad_info$mtime)
  all_route_tables[[i]] <- table
}
all_routes <- do.call(rbind, all_route_tables)

known <- all_routes[!is_unknown_consensus_label(all_routes$consensus_label), , drop = FALSE]
known$expected_route_from_consensus <- collapse_reference_vote(
  known$consensus_label, config
)
known$route_matches_consensus <- (
  known$reference_route == known$expected_route_from_consensus
)
# Known-label conflicts are retained as calibration evidence rather than being
# silently relabelled or treated as a pipeline error.

routes <- all_routes[
  is_unknown_consensus_label(all_routes$consensus_label), , drop = FALSE
]
route_keys <- unknown_qc_review_key(routes)
review_keys <- unknown_qc_review_key(review)
if (!setequal(route_keys, review_keys)) {
  stop(
    "Reference-route and reviewed-QC Unknown populations differ. Missing review ",
    "rows: ", length(setdiff(route_keys, review_keys)), "; unused review rows: ",
    length(setdiff(review_keys, route_keys)), "."
  )
}
review_index <- match(route_keys, review_keys)
if (anyNA(review_index)) stop("Internal error matching reference routes to QC review.")
routes$qc_decision <- review$decision[review_index]
routes$qc_review_notes <- as.character(review$review_notes[review_index])
routes <- routes[order(
  match(routes$sample_id, sample_ids),
  suppressWarnings(as.numeric(as.character(routes$seurat_clusters))),
  routes$seurat_clusters
), , drop = FALSE]
rownames(routes) <- NULL
validate_unknown_reference_route_manifest(routes, sample_ids)

summary_table <- as.data.frame(
  table(routes$qc_decision, routes$reference_route),
  stringsAsFactors = FALSE
)
names(summary_table) <- c("qc_decision", "reference_route", "n_populations")
summary_table <- summary_table[summary_table$n_populations > 0L, , drop = FALSE]

validation_columns <- c(
  "sample_id", "seurat_clusters", "consensus_label",
  unknown_reference_vote_columns(),
  paste0(unknown_reference_names(), "_calibrated_broad_vote"),
  paste0(unknown_reference_names(), "_broad_route"),
  "broad_reference_consensus",
  "reference_route", "reference_route_support_n",
  "reference_route_status", "expected_route_from_consensus",
  "route_matches_consensus"
)
route_columns <- c(
  "sample_id", "seurat_clusters", "consensus_label",
  unknown_reference_vote_columns(),
  paste0(unknown_reference_names(), "_calibrated_broad_vote"),
  paste0(unknown_reference_names(), "_broad_route"),
  "broad_reference_consensus",
  "reference_route", "reference_route_support_n",
  "reference_route_nonabstain_n", "reference_route_status",
  "qc_decision", "qc_review_notes", "consensus_table_path",
  "consensus_table_size", "consensus_table_mtime",
  "broad_consensus_table_path", "broad_consensus_table_size",
  "broad_consensus_table_mtime"
)

dir.create(route_paths$table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(routes[route_columns], route_paths$route_manifest, row.names = FALSE, na = "")
write.csv(input_manifest, route_paths$input_manifest, row.names = FALSE, na = "")
write.csv(known[validation_columns], route_paths$validation, row.names = FALSE, na = "")
write.csv(summary_table, route_paths$summary, row.names = FALSE, na = "")

message(
  "Wrote reference routes for ", nrow(routes),
  " reviewed sample-qualified Unknown populations."
)
print(summary_table, row.names = FALSE)
