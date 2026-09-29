#!/usr/bin/env Rscript

# Map refined RL identities back to all manifest-defined whole-tissue objects.
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(future)
  library(future.apply)
})
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
if (any(!args %in% valid_options)) {
  stop("Usage: Rscript scripts/xenium_rl_07_map_to_samples.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
input_dir <- file.path(output_root, "xenium", "rl", "04_qc", "rds")
output_dir <- file.path(output_root, "xenium", "rl", "07_mapped", "rds")
master_path <- file.path(
  output_root, "xenium", "rl", "06_subclusters", "rds",
  "Xenium_RL_subclusters.rds"
)
input_names <- paste0(sample_ids, "_RL_QC.rds")
input_paths <- file.path(input_dir, input_names)
output_names <- paste0(sample_ids, "_RL_QC_Subclusters.rds")
output_paths <- file.path(output_dir, output_names)
manifest_path <- file.path(output_dir, "Xenium_RL_Mapping_manifest.csv")
observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  compact_dry_run(
    "RL label mapping",
    inputs = c(master_path, input_paths),
    outputs = c(output_paths, manifest_path),
    checks = c(no_unexpected_inputs = !length(unexpected_names))
  )
  if (length(unexpected_names)) {
    cat("  Unexpected inputs: ", paste(unexpected_names, collapse = "; "), "\n")
  }
  quit(save = "no", status = 0L)
}

if (!file.exists(master_path)) stop("Master RL subcluster object not found: ", master_path)
if (length(missing_names)) {
  stop("Cannot map RL labels; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot map RL labels; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
expected_outputs <- c(output_paths, manifest_path)
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing RL mapping outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

message("Loading master RL subcluster object...")
master_obj <- readRDS(master_path)
required_master_metadata <- c("consensus_label", "RL_subset_reason")
missing_master_metadata <- setdiff(required_master_metadata, colnames(master_obj[[]]))
if (length(missing_master_metadata)) {
  stop("Master RL object lacks metadata: ", paste(missing_master_metadata, collapse = ", "))
}
all_new_labels <- setNames(as.character(Idents(master_obj)), colnames(master_obj))
all_subset_reasons <- setNames(as.character(master_obj$RL_subset_reason), colnames(master_obj))
if (anyNA(all_new_labels) || any(!nzchar(all_new_labels))) {
  stop("Master RL object contains missing or blank active identities.")
}
if (anyDuplicated(names(all_new_labels))) stop("Master RL object contains duplicate cell names.")
if (anyNA(all_subset_reasons) ||
    any(!all_subset_reasons %in% c("configured_broad_label", "shared_unknown"))) {
  stop("Master RL object contains missing or invalid RL_subset_reason values.")
}
if (!all(
  (all_subset_reasons == "shared_unknown") ==
    is_unknown_consensus_label(master_obj$consensus_label)
)) {
  stop("Master RL subset provenance does not match consensus Unknown identities.")
}

sample_prefixes <- paste0(sample_ids, "_")
master_sample_index <- vapply(names(all_new_labels), function(cell_id) {
  matches <- which(startsWith(cell_id, sample_prefixes))
  if (length(matches) == 1L) matches else NA_integer_
}, integer(1))
if (anyNA(master_sample_index)) {
  bad_cells <- names(all_new_labels)[is.na(master_sample_index)]
  stop(
    "Master RL cell names must match exactly one configured sample prefix. Examples: ",
    paste(head(bad_cells, 10L), collapse = ", ")
  )
}
master_cells_by_sample <- tabulate(master_sample_index, nbins = length(sample_ids))
if (any(master_cells_by_sample == 0L)) {
  stop(
    "Master RL object contains no cells for configured sample(s): ",
    paste(sample_ids[master_cells_by_sample == 0L], collapse = ", ")
  )
}
rm(master_obj, master_sample_index)
gc()

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
options(future.globals.maxSize = 200 * 1024^3)
plan(multisession, workers = workers)
on.exit(plan(sequential), add = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Mapping refined RL labels across all configured samples with ", workers, " worker(s)...")

mapping_rows <- future_lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  temp_obj <- readRDS(input_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(temp_obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  existing_regional_metadata <- intersect(
    c("RL_subcluster", "RL_subset_reason"), colnames(temp_obj[[]])
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
  temp_obj <- AddMetaData(temp_obj, relevant_labels, col.name = "RL_subcluster")
  temp_obj <- AddMetaData(temp_obj, relevant_reasons, col.name = "RL_subset_reason")
  if (sum(!is.na(temp_obj$RL_subcluster)) != expected_cells) {
    stop(sample_id, ": RL_subcluster metadata failed post-injection validation.")
  }
  mapped_unknown_cells <- sum(temp_obj$RL_subset_reason == "shared_unknown", na.rm = TRUE)
  expected_unknown_cells <- sum(
    relevant_reasons == "shared_unknown", na.rm = TRUE
  )
  if (mapped_unknown_cells != expected_unknown_cells) {
    stop(sample_id, ": RL shared-Unknown provenance failed post-injection validation.")
  }
  saveRDS(temp_obj, output_paths[[i]], compress = FALSE)
  data.frame(
    sample_id = sample_id, input_path = input_paths[[i]], output_path = output_paths[[i]],
    total_cells = ncol(temp_obj), mapped_rl_cells = mapped_cells,
    mapped_rl_unknown_cells = mapped_unknown_cells,
    PCW = pcw_values[[1]], stringsAsFactors = FALSE
  )
}, future.seed = TRUE)

plan(sequential)
mapping_manifest <- do.call(rbind, mapping_rows)
if (!identical(mapping_manifest$sample_id, sample_ids)) stop("RL mapping results are not in manifest order.")
if (sum(mapping_manifest$mapped_rl_cells) != length(all_new_labels)) {
  stop("Mapped RL cell total does not equal the master RL object cell count.")
}
if (sum(mapping_manifest$mapped_rl_unknown_cells) !=
    sum(all_subset_reasons == "shared_unknown")) {
  stop("Mapped RL shared-Unknown total does not equal the master RL object count.")
}
write.csv(mapping_manifest, manifest_path, row.names = FALSE)
message("Mapped every master RL identity across all configured samples.")
message("Saved RL mapping manifest: ", manifest_path)
