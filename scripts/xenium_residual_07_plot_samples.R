#!/usr/bin/env Rscript

# Manifest-driven residual spatial/marker reports; no RDS outputs are written.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_07_plot_samples.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))
source(here("scripts", "R", "residual_spatial_reports.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
sample_ids <- as.character(samples$sample_id)
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--list")) {
  write.table(data.frame(task_id = seq_along(sample_ids), sample_id = sample_ids),
              row.names = FALSE, quote = FALSE, sep = "\t")
  quit(save = "no", status = 0L)
}
valid <- c("--combine", "--validate", "--dry-run", "--overwrite")
task_args <- args[!startsWith(args, "--")]
if (any(!args[startsWith(args, "--")] %in% valid)) stop("Unknown residual report option.")
combine <- "--combine" %in% args
validate <- "--validate" %in% args
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if ((combine && validate) || (dry_run && overwrite) || length(task_args) > 1L ||
    ((combine || validate) && length(task_args)) || (validate && overwrite)) {
  stop("Usage: Rscript scripts/xenium_residual_07_plot_samples.R ",
       "[--list | --validate | --combine [--dry-run|--overwrite] | [--dry-run|--overwrite] TASK_ID]")
}
root <- here(config$project$outputs_dir, "xenium", "residual")
reports <- file.path(root, "07_sample_reports")
input_paths <- file.path(root, "06_mapped", "rds", paste0(sample_ids, "_Residual_Subclusters.rds"))
mapping_path <- file.path(root, "06_mapped", "tables", "Xenium_Residual_Mapping_Manifest.csv")
review_path <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels, config)
review <- load_reviewed_residual_subcluster_labels(config)
# Load shared contracts before deriving residual panels. Do not use the broad
# global `markers` variable defined by the shared palette as a residual panel.
source(here("scripts", "color_palette.R"))
residual_report_markers <- residual_report_marker_sets(review, residual_shared_marker_panels())
marker_panel_md5 <- residual_marker_panel_signature(residual_report_markers)
review_md5 <- unname(tools::md5sum(review_path))
mapping <- if (file.exists(mapping_path)) read.csv(mapping_path, stringsAsFactors = FALSE) else NULL
required <- c("sample_id", "output_path", "total_cells", "mapped_residual_cells", "removed_residual_qc_cells", "PCW")
mapping_valid <- !is.null(mapping) && all(required %in% names(mapping)) &&
  identical(as.character(mapping$sample_id), sample_ids) &&
  identical(as.character(mapping$output_path), input_paths)
if (mapping_valid) {
  for (column in c("total_cells", "mapped_residual_cells", "removed_residual_qc_cells")) {
    values <- mapping[[column]]
    mapping_valid <- mapping_valid && is.numeric(values) && all(is.finite(values)) &&
      all(values >= 0 & values == floor(values))
  }
}
dependencies <- c(pdfinfo = nzchar(Sys.which("pdfinfo")))
if (combine || validate) {
  if (!mapping_valid) stop("Residual mapping manifest is missing, incomplete, or differs from config/samples.csv.")
  audited <- audit_residual_spatial_reports(reports, mapping, input_paths, review_md5, marker_panel_md5)
  if (validate) {
    cat("PASS | residual reports | samples=", nrow(audited), "/", length(sample_ids),
        " | residual cells=", sum(audited$mapped_residual_cells),
        " | no-residual samples=", sum(audited$mapped_residual_cells == 0), "\n", sep = "")
  } else {
    page_map <- build_stage8_pcw_page_map(config, samples)
    combine_stage8_sample_reports(page_map, reports, "Residual", dry_run, overwrite)
  }
  quit(save = "no", status = 0L)
}
# Validate taxonomy during the path-only preflight as well as the plotting run.
missing_taxonomy <- setdiff(names(residual_report_markers),
                            intersect(master_subcluster_order, names(subcluster_palette)))
dependencies <- c(dependencies,
  approved_residual_taxonomy = !length(missing_taxonomy),
  Seurat = requireNamespace("Seurat", quietly = TRUE),
  ggplot2 = requireNamespace("ggplot2", quietly = TRUE),
  patchwork = requireNamespace("patchwork", quietly = TRUE),
  Cairo = requireNamespace("Cairo", quietly = TRUE),
  raster_backend = requireNamespace("scattermore", quietly = TRUE) || requireNamespace("ggrastr", quietly = TRUE))
if (dry_run) {
  cat("Residual marker panels (from current shared production definitions):\n")
  write.table(data.frame(
    identity = names(residual_report_markers),
    source = unname(attr(residual_report_markers, "marker_sources")),
    genes = vapply(residual_report_markers, paste, collapse = ";", FUN.VALUE = character(1)),
    row.names = NULL), row.names = FALSE, quote = FALSE, sep = "\t")
}
if (dry_run && !length(task_args)) {
  all_outputs <- unlist(lapply(sample_ids, function(id) residual_report_paths(reports, id)), use.names = FALSE)
  compact_dry_run("Residual spatial/marker reports (all samples)",
    inputs = c(input_paths, mapping_path, review_path), outputs = all_outputs,
    checks = c(mapping_manifest_current = mapping_valid, dependencies))
  ready <- mapping_valid && all(dependencies) && all(file.exists(c(input_paths, mapping_path, review_path)))
  quit(save = "no", status = if (ready) 0L else 1L)
}
task <- if (length(task_args)) task_args else Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
task <- suppressWarnings(as.integer(task))
if (length(task) != 1L || is.na(task) || task < 1L || task > length(sample_ids)) stop("TASK_ID must match --list.")
sample_id <- sample_ids[[task]]
paths <- residual_report_paths(reports, sample_id)
if (dry_run) {
  compact_dry_run(paste0("Residual sample report ", task, "/", length(sample_ids), " [", sample_id, "]"),
    inputs = c(input_paths[[task]], mapping_path, review_path), outputs = paths,
    checks = c(mapping_manifest_current = mapping_valid, dependencies))
  ready <- mapping_valid && all(dependencies) && all(file.exists(c(input_paths[[task]], mapping_path, review_path)))
  quit(save = "no", status = if (ready) 0L else 1L)
}
if (!mapping_valid) stop("Residual mapping manifest is missing, incomplete, or stale.")
if (!file.exists(input_paths[[task]])) stop("Residual mapped input missing: ", input_paths[[task]])
if (!all(dependencies)) stop("Missing report dependencies: ", paste(names(dependencies)[!dependencies], collapse = ", "),
                            ". Use the established renv environment and Poppler; no packages are installed by this script.")
if (any(file.exists(paths)) && !overwrite) stop("Residual report outputs exist; review before --overwrite.")
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
if (length(missing_taxonomy)) stop("Approved residual identities lack shared taxonomy/colors: ", paste(missing_taxonomy, collapse = ", "))
generate_residual_spatial_report(sample_id, input_paths[[task]], paths, mapping[task, , drop = FALSE],
                                 unique(as.character(review$final_label)), residual_report_markers, review_md5)
message("Completed residual spatial and marker report: ", sample_id)
