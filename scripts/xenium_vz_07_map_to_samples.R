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
master_path <- file.path(
  output_root, "xenium", "vz", "06_subclusters", "rds",
  "Xenium_VZ_subclusters_Res1.5.rds"
)
input_names <- paste0(sample_ids, "_Ald_VZ_QC.rds")
input_paths <- file.path(input_dir, input_names)
output_names <- paste0(sample_ids, "_Ald_VZ_QC_Subclusters.rds")
output_paths <- file.path(output_dir, output_names)
manifest_path <- file.path(output_dir, "Xenium_VZ_Mapping_manifest.csv")
observed_names <- if (dir.exists(input_dir)) {
  list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)
} else character()
missing_names <- input_names[!file.exists(input_paths)]
unexpected_names <- setdiff(observed_names, input_names)

if (dry_run) {
  cat("Expected VZ mapping inputs:", length(input_paths), "\n")
  cat("Existing VZ mapping inputs:", sum(file.exists(input_paths)), "\n")
  cat("Unexpected top-level RDS files:", length(unexpected_names), "\n")
  cat("Master VZ object:", master_path, "\n")
  cat("Master VZ object exists:", file.exists(master_path), "\n")
  cat("Mapping manifest:", manifest_path, "\n")
  cat("Mapping manifest exists:", file.exists(manifest_path), "\n")
  write.table(
    data.frame(
      sample_id = sample_ids, input = input_paths,
      input_exists = file.exists(input_paths), output = output_paths,
      output_exists = file.exists(output_paths)
    ),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  if (length(unexpected_names)) {
    cat("Unexpected files:\n- ", paste(unexpected_names, collapse = "\n- "), "\n")
  }
  quit(save = "no", status = 0L)
}

if (!file.exists(master_path)) stop("Master VZ subcluster object not found: ", master_path)
if (length(missing_names)) {
  stop("Cannot map VZ labels; missing files:\n- ", paste(missing_names, collapse = "\n- "))
}
if (length(unexpected_names)) {
  stop("Cannot map VZ labels; unexpected top-level RDS files:\n- ", paste(unexpected_names, collapse = "\n- "))
}
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
all_new_labels <- setNames(as.character(Idents(master_obj)), colnames(master_obj))
if (anyNA(all_new_labels) || any(!nzchar(all_new_labels))) {
  stop("Master VZ object contains missing or blank active identities.")
}
if (anyDuplicated(names(all_new_labels))) stop("Master VZ object contains duplicate cell names.")

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

workers <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))
if (is.na(workers) || workers < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer when set.")
options(future.globals.maxSize = 200 * 1024^3)
plan(multisession, workers = workers)
on.exit(plan(sequential), add = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Mapping refined VZ labels across all configured samples with ", workers, " worker(s)...")

mapping_rows <- future_lapply(seq_along(sample_ids), function(i) {
  sample_id <- sample_ids[[i]]
  temp_obj <- readRDS(input_paths[[i]])
  required_metadata <- c("consensus_label", "PCW")
  missing_metadata <- setdiff(required_metadata, colnames(temp_obj[[]]))
  if (length(missing_metadata)) stop(sample_id, " lacks metadata: ", paste(missing_metadata, collapse = ", "))
  if ("VZ_subcluster" %in% colnames(temp_obj[[]])) {
    stop(sample_id, " input already contains VZ_subcluster metadata.")
  }
  pcw_values <- unique(as.character(temp_obj$PCW))
  pcw_values <- pcw_values[!is.na(pcw_values) & nzchar(pcw_values)]
  if (length(pcw_values) != 1L) stop(sample_id, " must contain exactly one nonblank PCW value.")

  match_index <- match(paste0(sample_id, "_", colnames(temp_obj)), names(all_new_labels))
  relevant_labels <- all_new_labels[match_index]
  names(relevant_labels) <- colnames(temp_obj)
  mapped_cells <- sum(!is.na(relevant_labels))
  expected_cells <- master_cells_by_sample[[i]]
  if (mapped_cells != expected_cells) {
    stop(sample_id, ": mapped ", mapped_cells, " cells; master contains ", expected_cells, ".")
  }
  temp_obj <- AddMetaData(temp_obj, relevant_labels, col.name = "VZ_subcluster")
  if (sum(!is.na(temp_obj$VZ_subcluster)) != expected_cells) {
    stop(sample_id, ": VZ_subcluster metadata failed post-injection validation.")
  }
  saveRDS(temp_obj, output_paths[[i]], compress = FALSE)
  data.frame(
    sample_id = sample_id, input_path = input_paths[[i]], output_path = output_paths[[i]],
    total_cells = ncol(temp_obj), mapped_vz_cells = mapped_cells,
    PCW = pcw_values[[1]], stringsAsFactors = FALSE
  )
}, future.seed = TRUE)

plan(sequential)
mapping_manifest <- do.call(rbind, mapping_rows)
if (!identical(mapping_manifest$sample_id, sample_ids)) stop("VZ mapping results are not in manifest order.")
if (sum(mapping_manifest$mapped_vz_cells) != length(all_new_labels)) {
  stop("Mapped VZ cell total does not equal the master VZ object cell count.")
}
write.csv(mapping_manifest, manifest_path, row.names = FALSE)
message("Mapped every master VZ identity across all configured samples.")
message("Saved VZ mapping manifest: ", manifest_path)
