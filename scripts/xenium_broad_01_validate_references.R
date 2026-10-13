#!/usr/bin/env Rscript

# Validate that each production reference is restricted to the Xenium panel and
# that every retained reference identity has an explicit VZ/RL/Other mapping.
# This stage is read-only with respect to reference objects and production data.

rm(list = ls())

suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_broad_reference.R"))

config <- load_pipeline_config()
taxonomy <- load_reference_broad_taxonomy(config)
validate_taxonomy_against_config(taxonomy, config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) {
  stop("Usage: Rscript scripts/xenium_broad_01_validate_references.R [--dry-run|--overwrite]")
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

reference_paths <- c(
  Aldinger = resolve_config_path(config$inputs$references$aldinger, config),
  Sepp = resolve_config_path(config$inputs$references$sepp, config),
  Science = resolve_config_path(config$inputs$references$science, config)
)
panel_path <- resolve_config_path(config$inputs$xenium_panel_genes, config)
taxonomy_path <- attr(taxonomy, "path")

output_dir <- here(
  config$project$outputs_dir, "validation", "broad_reference_workflow",
  "01_reference_validation", "tables"
)
outputs <- c(
  summary = file.path(output_dir, "Xenium_Broad_Reference_Validation.csv"),
  labels = file.path(output_dir, "Xenium_Broad_Reference_Label_Taxonomy.csv"),
  features = file.path(output_dir, "Xenium_Broad_Reference_Feature_Audit.csv"),
  provenance = file.path(output_dir, "Xenium_Broad_Reference_Validation_Provenance.csv")
)

if (dry_run) {
  compact_dry_run(
    "Validate broad references, panel, and taxonomy",
    inputs = c(panel_path, taxonomy_path, reference_paths),
    outputs = outputs,
    checks = c(
      taxonomy_matches_config = tryCatch({
        validate_taxonomy_against_config(taxonomy, config)
        TRUE
      }, error = function(e) FALSE)
    )
  )
  quit(
    save = "no",
    status = if (all(file.exists(c(panel_path, taxonomy_path, reference_paths)))) 0L else 1L
  )
}

missing_inputs <- c(panel_path, taxonomy_path, reference_paths)[
  !file.exists(c(panel_path, taxonomy_path, reference_paths))
]
if (length(missing_inputs)) {
  stop("Missing reference-validation input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing_outputs <- outputs[file.exists(outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite broad-reference validation output(s):\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing the existing files."
  )
}

suppressPackageStartupMessages(library(Seurat))

panel_genes <- unique(trimws(as.character(readRDS(panel_path))))
if (!length(panel_genes) || anyNA(panel_genes) || any(!nzchar(panel_genes)) ||
    anyDuplicated(panel_genes)) {
  stop("Xenium panel gene vector must contain unique, nonblank gene names.")
}

summary_rows <- vector("list", length(reference_paths))
label_rows <- vector("list", length(reference_paths))
feature_rows <- vector("list", length(reference_paths))

for (i in seq_along(reference_paths)) {
  reference_name <- names(reference_paths)[[i]]
  reference_path <- reference_paths[[i]]
  message("Validating ", reference_name, " reference: ", reference_path)
  reference <- readRDS(reference_path)
  if (anyDuplicated(Cells(reference))) {
    stop(reference_name, " reference contains duplicate cell IDs.")
  }
  if (!identical(rownames(reference[[]]), Cells(reference))) {
    stop(reference_name, " reference metadata does not align exactly with cell IDs.")
  }
  transfer_assay <- if ("RNA" %in% Assays(reference)) "RNA" else Assays(reference)[[1]]
  reference_features <- rownames(reference[[transfer_assay]])
  if (!length(reference_features) || anyDuplicated(reference_features)) {
    stop(reference_name, " transfer assay has missing or duplicate features.")
  }
  outside_panel <- setdiff(reference_features, panel_genes)
  missing_panel <- setdiff(panel_genes, reference_features)
  if (length(outside_panel)) {
    stop(
      reference_name, " transfer assay contains ", length(outside_panel),
      " feature(s) outside the Xenium panel. First values: ",
      paste(head(sort(outside_panel), 10L), collapse = ", ")
    )
  }
  if (!"clusters_refined" %in% names(reference[[]])) {
    stop(reference_name, " reference lacks clusters_refined metadata.")
  }
  labels <- trimws(as.character(reference$clusters_refined))
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop(reference_name, " reference contains blank or missing clusters_refined values.")
  }
  label_counts <- as.data.frame(table(fine_label = labels), stringsAsFactors = FALSE)
  names(label_counts)[[2]] <- "n_cells"
  label_counts$fine_label <- as.character(label_counts$fine_label)
  taxonomy_index <- match(label_counts$fine_label, taxonomy$fine_label)
  if (anyNA(taxonomy_index)) {
    stop(
      reference_name, " identity label(s) are absent from the reviewed broad taxonomy: ",
      paste(sort(label_counts$fine_label[is.na(taxonomy_index)]), collapse = ", ")
    )
  }
  label_counts$reference <- reference_name
  label_counts$broad_label <- taxonomy$broad_label[taxonomy_index]
  label_counts$broad_route <- broad_category_to_route(label_counts$broad_label)
  label_counts <- label_counts[c(
    "reference", "fine_label", "broad_label", "broad_route", "n_cells"
  )]
  label_rows[[i]] <- label_counts

  feature_rows[[i]] <- data.frame(
    reference = reference_name,
    feature = panel_genes,
    present_in_reference = panel_genes %in% reference_features,
    stringsAsFactors = FALSE
  )
  summary_rows[[i]] <- data.frame(
    reference = reference_name,
    reference_path = normalizePath(reference_path, mustWork = TRUE),
    transfer_assay = transfer_assay,
    n_cells = ncol(reference),
    panel_genes = length(panel_genes),
    reference_panel_genes = length(reference_features),
    missing_panel_genes = length(missing_panel),
    outside_panel_genes = length(outside_panel),
    n_reference_labels = nrow(label_counts),
    all_labels_mapped = TRUE,
    stringsAsFactors = FALSE
  )
  rm(reference)
  invisible(gc(verbose = FALSE))
}

summary_table <- do.call(rbind, summary_rows)
label_table <- do.call(rbind, label_rows)
feature_table <- do.call(rbind, feature_rows)
provenance <- data.frame(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  panel_path = normalizePath(panel_path, mustWork = TRUE),
  panel_size = as.numeric(file.info(panel_path)$size),
  panel_mtime = as.numeric(file.info(panel_path)$mtime),
  taxonomy_path = taxonomy_path,
  taxonomy_size = as.numeric(file.info(taxonomy_path)$size),
  taxonomy_mtime = as.numeric(file.info(taxonomy_path)$mtime),
  broad_categories = paste(broad_reference_categories(), collapse = "|"),
  route_labels = paste(broad_category_to_route(broad_reference_categories()), collapse = "|"),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, outputs[["summary"]], row.names = FALSE)
write.csv(label_table, outputs[["labels"]], row.names = FALSE)
write.csv(feature_table, outputs[["features"]], row.names = FALSE)
write.csv(provenance, outputs[["provenance"]], row.names = FALSE)
message("Broad-reference panel and taxonomy validation passed for all three references.")
