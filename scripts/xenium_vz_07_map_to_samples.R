#!/usr/bin/env Rscript

# Map refined VZ identities back to all manifest-defined whole-tissue objects.
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
  library(future.apply)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_vz_rl_qc_routing.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_vz_07_map_to_samples.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "vz", "04_qc", "rds")
output_dir <- file.path(output_root, "xenium", "vz", "07_mapped", "rds")
decision_manifest_path <- file.path(
  output_root, "xenium", "vz", "04_qc", "tables",
  "Xenium_VZ_QC_removal_manifest.csv"
)
master_path <- file.path(
  output_root, "xenium", "vz", "06_subclusters", "rds",
  "Xenium_VZ_subclusters.rds"
)
input_names <- paste0(sample_ids, "_VZ_QC.rds")
legacy_input_paths <- file.path(input_dir, input_names)
output_names <- paste0(sample_ids, "_VZ_QC_Subclusters.rds")
output_paths <- file.path(output_dir, output_names)
manifest_path <- file.path(output_dir, "Xenium_VZ_Mapping_manifest.csv")
decision_manifest <- if (file.exists(decision_manifest_path)) {
  read.csv(decision_manifest_path, stringsAsFactors = FALSE)
} else NULL
decision_inspection <- inspect_regional_qc_manifest(decision_manifest, sample_ids)
preview_input_paths <- if (isTRUE(decision_inspection$valid) &&
                           isTRUE(decision_inspection$has_routing)) {
  as.character(decision_manifest$mapping_input_path)
} else {
  legacy_input_paths
}
routing_result <- if (isTRUE(decision_inspection$valid)) {
  tryCatch(
    resolve_regional_qc_mapping_inputs(
      decision_manifest, sample_ids, legacy_input_paths
    ),
    error = function(error) error
  )
} else NULL
routing_ready <- is.list(routing_result) && !inherits(routing_result, "error")
legacy_routing <- isTRUE(decision_inspection$valid) &&
  isTRUE(decision_inspection$legacy)
observed_names <- if (legacy_routing && dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
unexpected_names <- if (legacy_routing) setdiff(observed_names, input_names) else character()

if (dry_run) {
  compact_dry_run(
    "VZ label mapping",
    inputs = c(master_path, decision_manifest_path, preview_input_paths),
    outputs = c(output_paths, manifest_path),
    checks = c(
      decision_manifest = isTRUE(decision_inspection$valid),
      routed_inputs_current = routing_ready,
      no_unexpected_legacy_inputs = !length(unexpected_names)
    )
  )
  if (length(unexpected_names)) {
    cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  }
  quit(save = "no", status = 0L)
}

if (!file.exists(master_path)) stop("Master VZ subcluster object not found: ", master_path)
if (!isTRUE(decision_inspection$valid)) {
  stop(
    "VZ Stage-4 decision manifest is missing or invalid: ",
    decision_manifest_path, "\n- ",
    paste(decision_inspection$errors, collapse = "\n- ")
  )
}
if (!routing_ready) {
  stop("Cannot resolve VZ mapping inputs: ", conditionMessage(routing_result))
}
if (length(unexpected_names)) {
  stop("Cannot map VZ labels; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
input_paths <- routing_result$input_paths
materialization_modes <- routing_result$materialization_mode
expected_outputs <- c(output_paths, manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing VZ mapping outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

message("Loading master VZ subcluster object...")
master_obj <- readRDS(master_path)
required_master_metadata <- c("consensus_label", "VZ_subset_reason")
missing_master_metadata <- setdiff(required_master_metadata, colnames(master_obj[[]]))
if (length(missing_master_metadata)) {
  stop("Master VZ object lacks metadata: ", paste(missing_master_metadata, collapse = ", "))
}
all_new_labels <- setNames(as.character(Idents(master_obj)), colnames(master_obj))
all_subset_reasons <- setNames(as.character(master_obj$VZ_subset_reason), colnames(master_obj))
if (anyNA(all_new_labels) || any(!nzchar(all_new_labels))) {
  stop("Master VZ object contains missing or blank active identities.")
}
if (anyDuplicated(names(all_new_labels))) stop("Master VZ object contains duplicate cell names.")
if (anyNA(all_subset_reasons) ||
    any(!all_subset_reasons %in% c("configured_broad_label", "reference_route_vz"))) {
  stop("Master VZ object contains missing or invalid VZ_subset_reason values.")
}
if (!all(
  (all_subset_reasons == "reference_route_vz") ==
    is_unknown_consensus_label(master_obj$consensus_label)
)) {
  stop("Master VZ subset provenance does not match consensus Unknown identities.")
}

sample_prefixes <- paste0(sample_ids, "_")
master_sample_index <- vapply(names(all_new_labels), function(cell_id) {
  matches <- which(startsWith(cell_id, sample_prefixes))
  if (length(matches) == 1L) matches else NA_integer_
}, integer(1))
if (anyNA(master_sample_index)) {
  bad_cells <- names(all_new_labels)[is.na(master_sample_index)]
  stop(
    "Master VZ cell names must match exactly one configured sample prefix. Examples: ",
    paste(head(bad_cells, 10L), collapse = ", ")
  )
}
master_cells_by_sample <- tabulate(master_sample_index, nbins = length(sample_ids))
if (any(master_cells_by_sample == 0L)) {
  stop(
    "Master VZ object contains no cells for configured sample(s): ",
    paste(sample_ids[master_cells_by_sample == 0L], collapse = ", ")
  )
}
rm(master_obj, master_sample_index)
gc()

workers <- suppressWarnings(as.integer(Sys.getenv("VZ_MAPPING_WORKERS", unset = "1")))
if (is.na(workers) || workers < 1L) stop("VZ_MAPPING_WORKERS must be a positive integer when set.")
allocated_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(allocated_cpus) || allocated_cpus < 1L) stop("Invalid SLURM_CPUS_PER_TASK.")
workers <- min(workers, allocated_cpus)
options(future.globals.maxSize = 200 * 1024^3)
if (workers > 1L) plan(multisession, workers = workers) else plan(sequential)
on.exit(plan(sequential), add = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Mapping refined VZ labels across all configured samples with ", workers, " worker(s)...")

mapping_rows <- future_lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  temp_obj <- readRDS(input_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(temp_obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  existing_regional_metadata <- intersect(
    c("VZ_subcluster", "VZ_subset_reason"), colnames(temp_obj[[]])
  )
  if (length(existing_regional_metadata)) {
    stop(sample_id, " input already contains metadata: ",
         paste(existing_regional_metadata, collapse = ", "))
  }
  pcw_values <- unique(as.character(temp_obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

  match_index <- match(paste0(sample_id, "_", colnames(temp_obj)), names(all_new_labels))
  relevant_labels <- all_new_labels[match_index]
  relevant_reasons <- all_subset_reasons[match_index]
  names(relevant_labels) <- colnames(temp_obj)
  names(relevant_reasons) <- colnames(temp_obj)
  mapped_cells <- sum(!is.na(relevant_labels))
  expected_cells <- master_cells_by_sample[[i]]
  if (mapped_cells != expected_cells) {
    stop(sample_id, ": mapped ", mapped_cells, " cells; master contains ", expected_cells, ".")
  }
  temp_obj <- AddMetaData(temp_obj, relevant_labels, col.name = "VZ_subcluster")
  temp_obj <- AddMetaData(temp_obj, relevant_reasons, col.name = "VZ_subset_reason")
  if (sum(!is.na(temp_obj$VZ_subcluster)) != expected_cells) {
    stop(sample_id, ": VZ_subcluster metadata failed post-injection validation.")
  }
  mapped_unknown_cells <- sum(temp_obj$VZ_subset_reason == "reference_route_vz", na.rm = TRUE)
  expected_unknown_cells <- sum(
    relevant_reasons == "reference_route_vz", na.rm = TRUE
  )
  if (mapped_unknown_cells != expected_unknown_cells) {
    stop(sample_id, ": VZ reference-routed Unknown provenance failed post-injection validation.")
  }
  saveRDS(temp_obj, output_paths[[i]], compress = FALSE)
  data.frame(
    sample_id = sample_id, input_path = input_paths[[i]], output_path = output_paths[[i]],
    input_materialization_mode = materialization_modes[[i]],
    total_cells = ncol(temp_obj), mapped_vz_cells = mapped_cells,
    mapped_vz_unknown_cells = mapped_unknown_cells,
    PCW = pcw_values[[1]], stringsAsFactors = FALSE
  )
}, future.seed = TRUE)

plan(sequential)
mapping_manifest <- do.call(rbind, mapping_rows)
if (!identical(mapping_manifest$sample_id, sample_ids)) stop("VZ mapping results are not in manifest order.")
if (sum(mapping_manifest$mapped_vz_cells) != length(all_new_labels)) {
  stop("Mapped VZ cell total does not equal the master VZ object cell count.")
}
if (sum(mapping_manifest$mapped_vz_unknown_cells) !=
    sum(all_subset_reasons == "reference_route_vz")) {
  stop("Mapped VZ reference-routed Unknown total does not equal the master VZ object count.")
}
write.csv(mapping_manifest, manifest_path, row.names = FALSE)
message("Mapped every master VZ identity across all configured samples.")
message("Saved VZ mapping manifest: ", manifest_path)
