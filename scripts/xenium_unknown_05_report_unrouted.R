#!/usr/bin/env Rscript

# Create a non-mutating inventory of retained Unknown populations that are not
# assigned to either production regional branch. This is a report, not a third
# integration or clustering workflow.

rm(list = ls())

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_unknown_qc_review.R"))
source(here("scripts", "R", "xenium_unknown_reference_routes.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)
route_paths <- unknown_reference_route_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown <- args[!args %in% valid_options]
if (length(unknown)) {
  stop(
    "Usage: Rscript scripts/xenium_unknown_05_report_unrouted.R ",
    "[--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

qc_evidence_path <- file.path(
  production_paths$annotation_root, "04_unknown_qc_review", "tables",
  "Xenium_Unknown_QC_Evidence.csv"
)
report_dir <- file.path(
  production_paths$annotation_root, "05_unrouted_unknown_report", "tables"
)
population_path <- file.path(report_dir, "Xenium_Unrouted_Unknown_Populations.csv")
sample_summary_path <- file.path(report_dir, "Xenium_Unrouted_Unknown_BySample.csv")
route_summary_path <- file.path(report_dir, "Xenium_Unrouted_Unknown_Summary.csv")
input_paths <- c(route_paths$route_manifest, qc_evidence_path)
output_paths <- c(population_path, sample_summary_path, route_summary_path)

if (dry_run) {
  compact_dry_run(
    "Report retained unrouted Unknown populations",
    inputs = input_paths,
    outputs = output_paths
  )
  quit(save = "no", status = 0L)
}
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing unrouted-report input(s):\n- ", paste(missing_inputs, collapse = "\n- "))
}
existing <- output_paths[file.exists(output_paths)]
if (length(existing) && !overwrite) {
  stop(
    "Refusing to overwrite unrouted-report output(s):\n- ",
    paste(existing, collapse = "\n- "),
    "\nUse --overwrite only after reviewing the existing files."
  )
}

routes <- load_unknown_reference_routes(config)
evidence <- read.csv(
  qc_evidence_path, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = character()
)
evidence_keys <- unknown_qc_review_key(evidence)
if (anyDuplicated(evidence_keys)) stop("Unknown-QC evidence contains duplicate keys.")
route_keys <- unknown_qc_review_key(routes)
evidence_index <- match(route_keys, evidence_keys)
if (anyNA(evidence_index)) {
  stop("Unknown-QC evidence is missing one or more reference-route populations.")
}
qc_columns <- setdiff(
  names(evidence), c("sample_id", "seurat_clusters", "consensus_label")
)
for (column in qc_columns) routes[[column]] <- evidence[[column]][evidence_index]

keep <- tolower(routes$qc_decision) == "keep" &
  routes$reference_route %in% c("route_other", "route_unresolved")
unrouted <- routes[keep, , drop = FALSE]
marker_dir <- file.path(
  production_paths$annotation_root, config$production$consensus_label_stage,
  "tables", "unknown_markers"
)
unrouted$marker_workbook_path <- file.path(
  marker_dir,
  paste0(unrouted$sample_id, "_Unknown_cluster_top20_markers.xlsx")
)
unrouted$marker_workbook_exists <- file.exists(unrouted$marker_workbook_path)

if (nrow(unrouted)) {
  sample_summary <- aggregate(
    unrouted$n_cells,
    by = list(
      sample_id = unrouted$sample_id,
      reference_route = unrouted$reference_route
    ),
    FUN = sum
  )
  names(sample_summary)[[3]] <- "n_cells"
  population_counts <- aggregate(
    unrouted$consensus_label,
    by = list(
      sample_id = unrouted$sample_id,
      reference_route = unrouted$reference_route
    ),
    FUN = length
  )
  names(population_counts)[[3]] <- "n_populations"
  sample_summary <- merge(
    population_counts, sample_summary,
    by = c("sample_id", "reference_route"), sort = FALSE
  )
  sample_summary <- sample_summary[order(
    match(sample_summary$sample_id, sample_ids), sample_summary$reference_route
  ), , drop = FALSE]

  route_summary <- aggregate(
    unrouted$n_cells,
    by = list(reference_route = unrouted$reference_route),
    FUN = sum
  )
  names(route_summary)[[2]] <- "n_cells"
  route_population_counts <- as.data.frame(
    table(unrouted$reference_route), stringsAsFactors = FALSE
  )
  names(route_population_counts) <- c("reference_route", "n_populations")
  route_summary <- merge(
    route_population_counts, route_summary,
    by = "reference_route", sort = FALSE
  )
} else {
  sample_summary <- data.frame(
    sample_id = character(), reference_route = character(),
    n_populations = integer(), n_cells = integer()
  )
  route_summary <- data.frame(
    reference_route = character(), n_populations = integer(),
    n_cells = integer()
  )
}

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(unrouted, population_path, row.names = FALSE, na = "")
write.csv(sample_summary, sample_summary_path, row.names = FALSE, na = "")
write.csv(route_summary, route_summary_path, row.names = FALSE, na = "")
message(
  "Reported ", nrow(unrouted),
  " retained populations outside the VZ and RL production branches."
)
