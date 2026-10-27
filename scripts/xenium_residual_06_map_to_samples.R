#!/usr/bin/env Rscript

# Map explicitly reviewed residual identities back to every manifest-defined
# reviewed whole-tissue object. Samples without residual cells receive the same
# metadata columns filled with NA so downstream joins remain deterministic.

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
  library(future.apply)
})
here::i_am("scripts/xenium_residual_06_map_to_samples.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))
source(here("scripts", "R", "xenium_residual_qc_removal.R"))

config <- load_pipeline_config()
sample_ids <- as.character(load_sample_manifest(config)$sample_id)
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_residual_06_map_to_samples.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
qc_state <- load_residual_analysis_state(config)

output_root <- here(config$project$outputs_dir)
residual_root <- file.path(output_root, "xenium", "residual")
master_path <- file.path(
  residual_root, "05_subclusters", "rds", "Xenium_Residual_Subclusters.rds"
)
label_manifest_path <- file.path(
  residual_root, "05_subclusters", "tables",
  "Xenium_Residual_Subcluster_Labels.csv"
)
status_paths <- file.path(
  residual_root, "01_subsets", "tables",
  paste0(sample_ids, "_ResidualSubset_Status.csv")
)
output_dir <- file.path(residual_root, "06_mapped", "rds")
output_names <- paste0(sample_ids, "_Residual_Subclusters.rds")
output_paths <- file.path(output_dir, output_names)
mapping_manifest_path <- file.path(
  residual_root, "06_mapped", "tables", "Xenium_Residual_Mapping_Manifest.csv"
)

reviewed_input_state <- tryCatch(
  list(valid = TRUE, manifest = load_unknown_reviewed_input_manifest(config)),
  error = function(error) list(valid = FALSE, error = conditionMessage(error))
)
input_paths <- if (isTRUE(reviewed_input_state$valid)) {
  input_manifest <- reviewed_input_state$manifest
  if (anyDuplicated(input_manifest$sample_id) ||
      !setequal(as.character(input_manifest$sample_id), sample_ids)) {
    stop("Reviewed Unknown-QC input manifest does not map one-to-one to config/samples.csv.")
  }
  input_manifest <- input_manifest[match(sample_ids, input_manifest$sample_id), , drop = FALSE]
  as.character(input_manifest$mapping_input_path)
} else {
  rep(NA_character_, length(sample_ids))
}

all_input_paths <- c(master_path, label_manifest_path, status_paths, input_paths)
if (dry_run) {
  existing_inputs <- !is.na(all_input_paths) & file.exists(all_input_paths)
  compact_dry_run(
    "Residual label mapping",
    inputs = all_input_paths,
    outputs = c(output_paths, mapping_manifest_path),
    checks = c(
      reviewed_input_manifest = isTRUE(reviewed_input_state$valid),
      sample_manifest_valid = length(sample_ids) > 0L &&
        !anyNA(sample_ids) && !anyDuplicated(sample_ids) && all(nzchar(sample_ids))
    )
  )
  if (!isTRUE(reviewed_input_state$valid)) {
    cat("  Reviewed input manifest: ", reviewed_input_state$error, "\n")
  }
  quit(save = "no", status = if (all(existing_inputs) &&
    isTRUE(reviewed_input_state$valid) && length(sample_ids) > 0L &&
    !anyNA(sample_ids) && !anyDuplicated(sample_ids) && all(nzchar(sample_ids))) 0L else 1L)
}

missing_inputs <- all_input_paths[is.na(all_input_paths) | !file.exists(all_input_paths)]
if (length(missing_inputs)) {
  stop("Missing residual mapping input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
expected_outputs <- c(output_paths, mapping_manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing residual mapping outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

statuses <- do.call(rbind, lapply(status_paths, function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
if (nrow(statuses) != length(sample_ids) || anyDuplicated(statuses$sample_id) ||
    !setequal(as.character(statuses$sample_id), sample_ids)) {
  stop("Residual Stage-1 status files do not map one-to-one to config/samples.csv.")
}
statuses <- statuses[match(sample_ids, statuses$sample_id), , drop = FALSE]
required_status <- c("sample_id", "input_rds", "input_size", "input_mtime", "residual_cells")
missing_status <- setdiff(required_status, names(statuses))
if (length(missing_status)) {
  stop("Residual Stage-1 statuses lack: ", paste(missing_status, collapse = ", "))
}
if (!identical(normalizePath(input_paths, mustWork = TRUE),
               normalizePath(as.character(statuses$input_rds), mustWork = TRUE))) {
  stop("Residual Stage-1 status inputs differ from the current reviewed input manifest.")
}
input_info <- file.info(input_paths)
if (!isTRUE(all.equal(
      as.numeric(statuses$input_size), as.numeric(input_info$size)
    )) || !isTRUE(all.equal(
      as.numeric(statuses$input_mtime), as.numeric(input_info$mtime)
    ))) {
  stop("A reviewed whole-tissue input changed after residual Stage 1.")
}
expected_cells_by_sample <- suppressWarnings(as.integer(statuses$residual_cells))
if (anyNA(expected_cells_by_sample) || any(expected_cells_by_sample < 0L)) {
  stop("Residual Stage-1 statuses contain invalid residual cell counts.")
}
if (qc_state$active) {
  if (!identical(expected_cells_by_sample, as.integer(qc_state$manifest$input_cells))) {
    stop("Pre-QC Stage-1 counts differ from the QC-removal manifest.")
  }
  expected_cells_by_sample <- as.integer(qc_state$manifest$remaining_cells)
}

message("Loading reviewed residual subcluster object...")
master_obj <- readRDS(master_path)
validate_residual_qc_object(master_obj, qc_state)
label_manifest <- read.csv(
  label_manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_label_manifest <- c(
  "cluster", "final_label", "approved", "source_cluster_column",
  "review_manifest", "review_manifest_size", "review_manifest_mtime",
  "input_rds", "input_size", "input_mtime"
)
missing_label_manifest <- setdiff(required_label_manifest, names(label_manifest))
if (length(missing_label_manifest)) {
  stop(
    "Residual Stage-5 label manifest lacks: ",
    paste(missing_label_manifest, collapse = ", ")
  )
}
if (length(unique(label_manifest$input_rds)) != 1L ||
    length(unique(label_manifest$input_size)) != 1L ||
    length(unique(label_manifest$input_mtime)) != 1L ||
    !identical(normalizePath(label_manifest$input_rds[[1]], mustWork = TRUE),
               normalizePath(qc_state$input, mustWork = TRUE)) ||
    !residual_qc_file_matches(qc_state$input, label_manifest$input_size[[1]], label_manifest$input_mtime[[1]])) {
  stop("Stage-5 labels do not match the current residual QC analysis input.")
}
approved <- suppressWarnings(as.logical(label_manifest$approved))
if (!nrow(label_manifest) || anyNA(approved) || any(!approved) ||
    anyDuplicated(as.character(label_manifest$cluster))) {
  stop("Residual Stage-5 label manifest is incomplete, unapproved, or duplicated.")
}
cluster_columns <- unique(trimws(as.character(label_manifest$source_cluster_column)))
if (length(cluster_columns) != 1L || !nzchar(cluster_columns) ||
    !cluster_columns %in% colnames(master_obj[[]])) {
  stop("Residual Stage-5 label manifest has an invalid source_cluster_column.")
}
review_paths <- unique(as.character(label_manifest$review_manifest))
review_sizes <- unique(suppressWarnings(as.numeric(label_manifest$review_manifest_size)))
review_mtimes <- unique(suppressWarnings(as.numeric(label_manifest$review_manifest_mtime)))
if (length(review_paths) != 1L || length(review_sizes) != 1L ||
    length(review_mtimes) != 1L || !file.exists(review_paths) ||
    !is.finite(review_sizes) || !is.finite(review_mtimes)) {
  stop("Residual Stage-5 label manifest has invalid review-file provenance.")
}
review_info <- file.info(review_paths)
if (!isTRUE(all.equal(as.numeric(review_info$size), review_sizes)) ||
    !isTRUE(all.equal(as.numeric(review_info$mtime), review_mtimes))) {
  stop("The approved residual-label review file changed after Stage 5.")
}
required_master_metadata <- c(
  "consensus_label", "residual_reference_route", "residual_subset_reason",
  "residual_subcluster", "residual_annotation_confidence",
  "residual_annotation_status"
)
missing_master_metadata <- setdiff(required_master_metadata, colnames(master_obj[[]]))
if (length(missing_master_metadata)) {
  stop("Residual master object lacks metadata: ", paste(missing_master_metadata, collapse = ", "))
}
if (ncol(master_obj) == 0L || anyDuplicated(colnames(master_obj))) {
  stop("Residual master object is empty or contains duplicate cell IDs.")
}
if (any(!as.character(master_obj$residual_reference_route) %in%
        c("route_other", "route_unresolved"))) {
  stop("Residual master object contains a VZ/RL or invalid reference route.")
}
if (any(as.character(master_obj$residual_subset_reason) !=
        "retained_high_quality_other_or_unresolved")) {
  stop("Residual master object contains an invalid residual_subset_reason.")
}
if (any(!is_unknown_consensus_label(master_obj$consensus_label))) {
  stop("Residual master object contains cells that were not consensus Unknown.")
}
master_metadata <- master_obj[[]][, required_master_metadata, drop = FALSE]
rownames(master_metadata) <- colnames(master_obj)
if (anyNA(master_metadata$residual_subcluster) ||
    any(!nzchar(trimws(as.character(master_metadata$residual_subcluster))))) {
  stop("Residual master object contains a missing or blank reviewed label.")
}
master_clusters <- as.character(master_obj[[cluster_columns, drop = TRUE]])
label_lookup <- setNames(
  trimws(as.character(label_manifest$final_label)),
  trimws(as.character(label_manifest$cluster))
)
expected_master_labels <- unname(label_lookup[master_clusters])
if (anyNA(expected_master_labels) ||
    !identical(
      trimws(as.character(master_metadata$residual_subcluster)),
      expected_master_labels
    )) {
  stop("Residual master labels do not exactly match the approved Stage-5 manifest.")
}
expected_annotation_status <- ifelse(
  tolower(expected_master_labels) == "unknown", "retained_unknown", "reviewed_label"
)
if (!identical(
      as.character(master_metadata$residual_annotation_status),
      expected_annotation_status
    )) {
  stop("Residual master annotation status does not match its reviewed labels.")
}

sample_prefixes <- paste0(sample_ids, "_")
master_sample_index <- vapply(rownames(master_metadata), function(cell_id) {
  matches <- which(startsWith(cell_id, sample_prefixes))
  if (length(matches) == 1L) matches else NA_integer_
}, integer(1))
if (anyNA(master_sample_index)) {
  bad_cells <- rownames(master_metadata)[is.na(master_sample_index)]
  stop(
    "Residual master cell names must match exactly one configured sample prefix. Examples: ",
    paste(head(bad_cells, 10L), collapse = ", ")
  )
}
observed_cells_by_sample <- tabulate(master_sample_index, nbins = length(sample_ids))
if (!identical(as.integer(observed_cells_by_sample), expected_cells_by_sample)) {
  differences <- sample_ids[observed_cells_by_sample != expected_cells_by_sample]
  stop(
    "Residual master cell counts differ from Stage-1 statuses for: ",
    paste(differences, collapse = ", ")
  )
}
rm(master_obj, master_sample_index)
invisible(gc(verbose = FALSE))

workers <- suppressWarnings(as.integer(Sys.getenv("RESIDUAL_MAPPING_WORKERS", unset = "1")))
if (is.na(workers) || workers < 1L) {
  stop("RESIDUAL_MAPPING_WORKERS must be a positive integer when set.")
}
allocated_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(allocated_cpus) || allocated_cpus < 1L) stop("Invalid SLURM_CPUS_PER_TASK.")
workers <- min(workers, allocated_cpus)
options(future.globals.maxSize = 200 * 1024^3)
if (workers > 1L) plan(multisession, workers = workers) else plan(sequential)
on.exit(plan(sequential), add = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(mapping_manifest_path), recursive = TRUE, showWarnings = FALSE)
message("Mapping reviewed residual labels across all configured samples with ", workers, " worker(s)...")

mapping_rows <- future_lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  object <- readRDS(input_paths[[i]])
  required_input_metadata <- c("consensus_label", "PCW")
  missing_input_metadata <- setdiff(required_input_metadata, colnames(object[[]]))
  if (length(missing_input_metadata)) {
    stop(sample_id, " lacks metadata: ", paste(missing_input_metadata, collapse = ", "))
  }
  residual_columns <- intersect(required_master_metadata[-1L], colnames(object[[]]))
  if (length(residual_columns)) {
    stop(sample_id, " input already contains residual metadata: ",
         paste(residual_columns, collapse = ", "))
  }
  if (anyDuplicated(colnames(object))) stop(sample_id, " contains duplicate cell IDs.")
  pcw_values <- unique(trimws(as.character(object$PCW)))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) {
    stop(sample_id, " must contain exactly one nonblank PCW value.")
  }

  prefixed_cells <- paste0(sample_id, "_", colnames(object))
  qc_removed <- residual_qc_exclusions_for_sample(colnames(object), object$consensus_label, sample_id, qc_state)
  match_index <- match(prefixed_cells, rownames(master_metadata))
  mapped <- !is.na(match_index)
  if (any(qc_removed & mapped)) stop(sample_id, ": a QC-excluded cell regained a residual label.")
  if (sum(mapped) != expected_cells_by_sample[[i]]) {
    stop(
      sample_id, ": mapped ", sum(mapped), " residual cells; Stage 1 and the master contain ",
      expected_cells_by_sample[[i]], "."
    )
  }
  mapped_metadata <- master_metadata[match_index, required_master_metadata[-1L], drop = FALSE]
  rownames(mapped_metadata) <- colnames(object)
  expected_consensus <- as.character(master_metadata$consensus_label[match_index[mapped]])
  observed_consensus <- as.character(object$consensus_label[mapped])
  if (anyNA(observed_consensus) || !identical(observed_consensus, expected_consensus)) {
    stop(sample_id, " consensus labels changed between the reviewed input and residual master.")
  }
  object <- AddMetaData(object, metadata = mapped_metadata)
  object$residual_qc_status <- ifelse(qc_removed, "removed", ifelse(mapped, "retained", "not_residual"))
  object@misc$residual_qc_decision_md5 <- qc_state$decision_md5
  object@misc$residual_post_qc_manifest_md5 <- unname(tools::md5sum(qc_state$post_qc_paths$provenance))
  if (sum(!is.na(object$residual_subcluster)) != expected_cells_by_sample[[i]]) {
    stop(sample_id, ": residual metadata failed post-injection validation.")
  }
  if (any(!is.na(object$residual_subcluster) &
          !is_unknown_consensus_label(object$consensus_label))) {
    stop(sample_id, ": a residual label was assigned to a non-Unknown consensus cell.")
  }

  saveRDS(object, output_paths[[i]], compress = FALSE)
  data.frame(
    sample_id = sample_id,
    input_path = normalizePath(input_paths[[i]], mustWork = TRUE),
    output_path = output_paths[[i]],
    total_cells = ncol(object),
    mapped_residual_cells = sum(mapped),
    removed_residual_qc_cells = sum(qc_removed),
    mapped_reviewed_label_cells = sum(
      mapped & master_metadata$residual_annotation_status[match_index] == "reviewed_label",
      na.rm = TRUE
    ),
    mapped_retained_unknown_cells = sum(
      mapped & master_metadata$residual_annotation_status[match_index] == "retained_unknown",
      na.rm = TRUE
    ),
    PCW = pcw_values[[1]],
    stringsAsFactors = FALSE
  )
}, future.seed = TRUE)

plan(sequential)
mapping_manifest <- do.call(rbind, mapping_rows)
if (!identical(as.character(mapping_manifest$sample_id), sample_ids)) {
  stop("Residual mapping results are not in manifest order.")
}
if (sum(mapping_manifest$mapped_residual_cells) != nrow(master_metadata)) {
  stop("Mapped residual cell total does not equal the residual master object cell count.")
}
if (sum(mapping_manifest$removed_residual_qc_cells) != nrow(qc_state$removed)) {
  stop("Mapped QC-exclusion total differs from the approved removed-cell manifest.")
}
if (sum(mapping_manifest$mapped_reviewed_label_cells) +
    sum(mapping_manifest$mapped_retained_unknown_cells) != nrow(master_metadata)) {
  stop("Residual annotation-status totals do not equal the mapped residual cell total.")
}
write.csv(mapping_manifest, mapping_manifest_path, row.names = FALSE)
message("Mapped every reviewed residual identity across all configured samples.")
message("Saved residual mapping manifest: ", mapping_manifest_path)
