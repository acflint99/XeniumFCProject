#!/usr/bin/env Rscript

# Compare the independently reviewed VZ and RL subcluster annotations assigned
# to the same whole-tissue sample-qualified Unknown cells. This is a review
# stage only: it never changes an RDS object or chooses a final cell identity.

rm(list = ls())
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "consensus_labels.R"))
source(here("scripts", "R", "xenium_vz_rl_unknown_relabel.R"))

config <- load_pipeline_config()
sample_ids <- load_sample_manifest(config)$sample_id
production_paths <- get_production_output_paths(config)

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c("--dry-run", "--overwrite")
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}
if (any(!args %in% valid_options)) {
  stop(
    "Usage: Rscript scripts/xenium_vz_rl_diagnostic_compare_unknown_annotations.R ",
    "[--dry-run|--overwrite]"
  )
}
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

output_root <- here(config$project$outputs_dir)
base_input_dir <- production_paths$consensus_rds_dir
vz_input_dir <- file.path(output_root, "xenium", "vz", "07_mapped", "rds")
rl_input_dir <- file.path(output_root, "xenium", "rl", "07_mapped", "rds")
stage_dir <- file.path(output_root, "xenium", "vz_rl", "00_unknown_annotation_review")
table_dir <- file.path(stage_dir, "tables")
plot_dir <- file.path(stage_dir, "plots")

base_input_paths <- file.path(
  base_input_dir, paste0(sample_ids, "_Consensus_annotated.rds")
)
vz_input_paths <- file.path(
  vz_input_dir, paste0(sample_ids, "_VZ_QC_Subclusters.rds")
)
rl_input_paths <- file.path(
  rl_input_dir, paste0(sample_ids, "_RL_QC_Subclusters.rds")
)
input_paths <- c(base_input_paths, vz_input_paths, rl_input_paths)

cell_table_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_Cell_Annotations.csv"
)
population_summary_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_Population_Summary.csv"
)
pair_table_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_VZ_RL_Pairs.csv"
)
global_pair_table_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_VZ_RL_Pairs_Global.csv"
)
review_template_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_Relabel_Review_Template.csv"
)
comparison_manifest_path <- file.path(
  table_dir, "Xenium_VZRL_Unknown_Annotation_Comparison_manifest.csv"
)
comparison_plot_path <- file.path(
  plot_dir, "Xenium_VZRL_Unknown_VZ_RL_Comparison.pdf"
)
expected_outputs <- c(
  cell_table_path, population_summary_path, pair_table_path,
  global_pair_table_path, review_template_path, comparison_manifest_path,
  comparison_plot_path
)

expected_vz_names <- c(
  paste0(sample_ids, "_VZ_QC_Subclusters.rds"), "Xenium_VZ_Mapping_manifest.csv"
)
expected_rl_names <- c(
  paste0(sample_ids, "_RL_QC_Subclusters.rds"), "Xenium_RL_Mapping_manifest.csv"
)
observed_vz_names <- if (dir.exists(vz_input_dir)) {
  list.files(vz_input_dir, pattern = "\\.rds$|_manifest\\.csv$", full.names = FALSE)
} else character()
observed_rl_names <- if (dir.exists(rl_input_dir)) {
  list.files(rl_input_dir, pattern = "\\.rds$|_manifest\\.csv$", full.names = FALSE)
} else character()
unexpected_vz_names <- setdiff(observed_vz_names, expected_vz_names)
unexpected_rl_names <- setdiff(observed_rl_names, expected_rl_names)

if (dry_run) {
  compact_dry_run(
    "VZ/RL Unknown annotation comparison",
    inputs = input_paths,
    outputs = expected_outputs,
    checks = c(
      sample_manifest_nonempty = length(sample_ids) > 0L,
      no_duplicate_samples = !anyDuplicated(sample_ids),
      no_unexpected_vz_inputs = !length(unexpected_vz_names),
      no_unexpected_rl_inputs = !length(unexpected_rl_names)
    )
  )
  if (length(unexpected_vz_names)) {
    cat("  Unexpected VZ inputs: ", paste(unexpected_vz_names, collapse = "; "), "\n")
  }
  if (length(unexpected_rl_names)) {
    cat("  Unexpected RL inputs: ", paste(unexpected_rl_names, collapse = "; "), "\n")
  }
  quit(save = "no", status = 0L)
}

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Cannot compare VZ/RL annotations; missing inputs:\n- ", paste(missing_inputs, collapse = "\n- "))
}
if (length(unexpected_vz_names) || length(unexpected_rl_names)) {
  stop(
    "Cannot compare VZ/RL annotations; unexpected mapped-stage files:\n- ",
    paste(c(unexpected_vz_names, unexpected_rl_names), collapse = "\n- ")
  )
}
existing_outputs <- expected_outputs[file.exists(expected_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing Unknown-comparison outputs:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nUse --overwrite only after reviewing them."
  )
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_manifest <- rbind(
  data.frame(sample_id = sample_ids, input_type = "consensus", input_path = base_input_paths),
  data.frame(sample_id = sample_ids, input_type = "VZ", input_path = vz_input_paths),
  data.frame(sample_id = sample_ids, input_type = "RL", input_path = rl_input_paths)
)
input_manifest <- input_manifest[order(
  match(input_manifest$sample_id, sample_ids),
  match(input_manifest$input_type, c("consensus", "VZ", "RL"))
), , drop = FALSE]
input_manifest$input_path <- normalizePath(input_manifest$input_path, mustWork = TRUE)
input_info <- file.info(input_manifest$input_path)
input_manifest$input_size <- as.numeric(input_info$size)
input_manifest$input_mtime <- as.numeric(input_info$mtime)
rownames(input_manifest) <- NULL

cell_rows <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message("Comparing Unknown annotations ", i, " of ", length(sample_ids), ": ", sample_id)
  base_obj <- readRDS(base_input_paths[[i]])
  vz_obj <- readRDS(vz_input_paths[[i]])
  rl_obj <- readRDS(rl_input_paths[[i]])

  required_base <- "consensus_label"
  required_vz <- c("consensus_label", "VZ_subcluster", "VZ_subset_reason")
  required_rl <- c("consensus_label", "RL_subcluster", "RL_subset_reason")
  missing_metadata <- c(
    setdiff(required_base, colnames(base_obj[[]])),
    setdiff(required_vz, colnames(vz_obj[[]])),
    setdiff(required_rl, colnames(rl_obj[[]]))
  )
  if (length(missing_metadata)) {
    stop(sample_id, " lacks required metadata: ", paste(unique(missing_metadata), collapse = ", "))
  }
  if (anyDuplicated(colnames(base_obj)) || anyDuplicated(colnames(vz_obj)) ||
      anyDuplicated(colnames(rl_obj))) {
    stop(sample_id, " has duplicate cell IDs in a consensus or regional object.")
  }
  if (!all(colnames(vz_obj) %in% colnames(base_obj)) ||
      !all(colnames(rl_obj) %in% colnames(base_obj))) {
    stop(sample_id, " has regional cell IDs absent from the consensus object.")
  }

  base_consensus <- setNames(as.character(base_obj$consensus_label), colnames(base_obj))
  if (anyNA(base_consensus) || any(!nzchar(base_consensus))) {
    stop(sample_id, " has missing or blank consensus labels.")
  }
  if (!identical(
    unname(as.character(vz_obj$consensus_label)),
    unname(base_consensus[colnames(vz_obj)])
  ) || !identical(
    unname(as.character(rl_obj$consensus_label)),
    unname(base_consensus[colnames(rl_obj)])
  )) {
    stop(sample_id, " consensus labels changed between base and regional objects.")
  }

  unknown_cells <- names(base_consensus)[is_unknown_consensus_label(base_consensus)]
  if (!length(unknown_cells)) {
    rm(base_obj, vz_obj, rl_obj)
    gc()
    next
  }
  vz_labels <- setNames(as.character(vz_obj$VZ_subcluster), colnames(vz_obj))
  rl_labels <- setNames(as.character(rl_obj$RL_subcluster), colnames(rl_obj))
  vz_reasons <- setNames(as.character(vz_obj$VZ_subset_reason), colnames(vz_obj))
  rl_reasons <- setNames(as.character(rl_obj$RL_subset_reason), colnames(rl_obj))
  vz_values <- unname(vz_labels[unknown_cells])
  rl_values <- unname(rl_labels[unknown_cells])
  vz_reason_values <- unname(vz_reasons[unknown_cells])
  rl_reason_values <- unname(rl_reasons[unknown_cells])
  has_vz <- !is.na(vz_values) & nzchar(vz_values)
  has_rl <- !is.na(rl_values) & nzchar(rl_values)
  if (any(has_vz & (is.na(vz_reason_values) | vz_reason_values != "shared_unknown")) ||
      any(has_rl & (is.na(rl_reason_values) | rl_reason_values != "shared_unknown"))) {
    stop(sample_id, " has Unknown regional labels without shared_unknown provenance.")
  }

  availability <- ifelse(
    has_vz & has_rl, "both",
    ifelse(has_vz, "VZ_only", ifelse(has_rl, "RL_only", "neither"))
  )
  consensus_values <- unname(base_consensus[unknown_cells])
  cell_rows[[i]] <- data.frame(
    sample_id = sample_id,
    cell_id = unknown_cells,
    consensus_label = consensus_values,
    sample_qualified_unknown = sample_qualify_unknown_consensus_labels(
      consensus_values, sample_id
    ),
    VZ_subcluster = ifelse(has_vz, vz_values, ""),
    RL_subcluster = ifelse(has_rl, rl_values, ""),
    annotation_availability = availability,
    exact_subcluster_label_agreement = has_vz & has_rl & vz_values == rl_values,
    stringsAsFactors = FALSE
  )

  rm(base_obj, vz_obj, rl_obj)
  gc()
}

cell_rows <- cell_rows[!vapply(cell_rows, is.null, logical(1))]
if (!length(cell_rows)) {
  stop("No consensus Unknown cells were found in the manifest-defined samples.")
}
cell_annotations <- do.call(rbind, cell_rows)
if (anyDuplicated(paste(cell_annotations$sample_id, cell_annotations$cell_id, sep = "\u001f"))) {
  stop("Unknown cell comparison contains duplicate sample/cell IDs.")
}

dominant_label <- function(labels) {
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) return(c(label = "", cells = 0, fraction = NA_real_))
  counts <- sort(table(labels), decreasing = TRUE)
  winners <- sort(names(counts)[counts == max(counts)])
  winner <- winners[[1]]
  c(label = winner, cells = as.numeric(counts[[winner]]), fraction = counts[[winner]] / length(labels))
}

population_indices <- split(
  seq_len(nrow(cell_annotations)),
  paste(
    cell_annotations$sample_id,
    cell_annotations$consensus_label,
    sep = "\u001f"
  )
)
population_summary <- do.call(rbind, lapply(population_indices, function(index) {
  current <- cell_annotations[index, , drop = FALSE]
  vz_dominant <- dominant_label(current$VZ_subcluster)
  rl_dominant <- dominant_label(current$RL_subcluster)
  both <- current$annotation_availability == "both"
  data.frame(
    sample_id = current$sample_id[[1]],
    consensus_label = current$consensus_label[[1]],
    sample_qualified_unknown = current$sample_qualified_unknown[[1]],
    base_unknown_cells = nrow(current),
    VZ_annotated_cells = sum(current$annotation_availability %in% c("both", "VZ_only")),
    RL_annotated_cells = sum(current$annotation_availability %in% c("both", "RL_only")),
    shared_annotated_cells = sum(both),
    VZ_only_cells = sum(current$annotation_availability == "VZ_only"),
    RL_only_cells = sum(current$annotation_availability == "RL_only"),
    neither_annotated_cells = sum(current$annotation_availability == "neither"),
    fraction_shared_annotated = mean(both),
    exact_subcluster_label_agreement_cells = sum(current$exact_subcluster_label_agreement),
    exact_subcluster_label_agreement_fraction_of_shared = if (sum(both)) {
      sum(current$exact_subcluster_label_agreement) / sum(both)
    } else NA_real_,
    dominant_VZ_subcluster = vz_dominant[["label"]],
    dominant_VZ_cells = as.numeric(vz_dominant[["cells"]]),
    dominant_VZ_fraction = as.numeric(vz_dominant[["fraction"]]),
    dominant_RL_subcluster = rl_dominant[["label"]],
    dominant_RL_cells = as.numeric(rl_dominant[["cells"]]),
    dominant_RL_fraction = as.numeric(rl_dominant[["fraction"]]),
    stringsAsFactors = FALSE
  )
}))
population_summary <- population_summary[order(
  match(population_summary$sample_id, sample_ids),
  population_summary$consensus_label
), , drop = FALSE]
rownames(population_summary) <- NULL

shared_annotations <- cell_annotations[
  cell_annotations$annotation_availability == "both", , drop = FALSE
]
if (!nrow(shared_annotations)) {
  stop("No consensus Unknown cells have both VZ and RL annotations.")
}
pair_table <- aggregate(
  list(cells = rep(1L, nrow(shared_annotations))),
  by = shared_annotations[c(
    "sample_id", "consensus_label", "sample_qualified_unknown",
    "VZ_subcluster", "RL_subcluster"
  )],
  FUN = sum
)
pair_table <- merge(
  pair_table,
  population_summary[c(
    "sample_id", "consensus_label", "base_unknown_cells", "shared_annotated_cells"
  )],
  by = c("sample_id", "consensus_label"),
  all.x = TRUE,
  sort = FALSE
)
pair_table$fraction_of_base_unknown_population <-
  pair_table$cells / pair_table$base_unknown_cells
pair_table$fraction_of_shared_annotated_cells <-
  pair_table$cells / pair_table$shared_annotated_cells
pair_table <- pair_table[order(
  match(pair_table$sample_id, sample_ids), pair_table$consensus_label,
  -pair_table$cells, pair_table$VZ_subcluster, pair_table$RL_subcluster
), , drop = FALSE]
rownames(pair_table) <- NULL

global_pair_table <- aggregate(
  list(cells = pair_table$cells),
  by = pair_table[c("VZ_subcluster", "RL_subcluster")],
  FUN = sum
)
global_pair_table <- global_pair_table[order(
  -global_pair_table$cells,
  global_pair_table$VZ_subcluster,
  global_pair_table$RL_subcluster
), , drop = FALSE]
global_pair_table$fraction_of_all_shared_unknown_cells <-
  global_pair_table$cells / sum(global_pair_table$cells)
rownames(global_pair_table) <- NULL

review_template <- pair_table[c(
  "sample_id", "consensus_label", "VZ_subcluster", "RL_subcluster",
  "cells", "fraction_of_base_unknown_population",
  "fraction_of_shared_annotated_cells", "base_unknown_cells",
  "shared_annotated_cells"
)]
review_template$final_label <- ""
review_template$review_notes <- ""
review_template <- review_template[c(
  unknown_relabel_key_columns(), "final_label", "review_notes",
  "cells", "fraction_of_base_unknown_population",
  "fraction_of_shared_annotated_cells", "base_unknown_cells",
  "shared_annotated_cells"
)]
unknown_relabel_key(review_template)

write.csv(cell_annotations, cell_table_path, row.names = FALSE)
write.csv(population_summary, population_summary_path, row.names = FALSE)
write.csv(pair_table, pair_table_path, row.names = FALSE)
write.csv(global_pair_table, global_pair_table_path, row.names = FALSE)
write.csv(review_template, review_template_path, row.names = FALSE)
write.csv(input_manifest, comparison_manifest_path, row.names = FALSE)

plot_pair_heatmap <- function(data, title, subtitle = NULL) {
  ggplot(data, aes(x = .data$RL_subcluster, y = .data$VZ_subcluster)) +
    geom_tile(aes(fill = log10(.data$cells + 1)), color = "white", linewidth = 0.4) +
    geom_text(aes(label = .data$cells), size = 3) +
    scale_fill_viridis_c(name = "log10(cells + 1)") +
    labs(x = "RL subcluster", y = "VZ subcluster", title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

grDevices::cairo_pdf(comparison_plot_path, width = 12, height = 9)
pdf_open <- TRUE
on.exit({
  if (isTRUE(pdf_open)) grDevices::dev.off()
}, add = TRUE)
print(plot_pair_heatmap(
  global_pair_table,
  "VZ versus RL annotations for all shared consensus Unknown cells",
  paste0(
    format(sum(global_pair_table$cells), big.mark = ","),
    " cells with both regional annotations"
  )
))
for (qualified_label in population_summary$sample_qualified_unknown) {
  current <- pair_table[
    pair_table$sample_qualified_unknown == qualified_label, , drop = FALSE
  ]
  if (!nrow(current)) next
  summary_row <- population_summary[
    population_summary$sample_qualified_unknown == qualified_label, , drop = FALSE
  ]
  print(plot_pair_heatmap(
    current,
    qualified_label,
    paste0(
      summary_row$shared_annotated_cells, "/", summary_row$base_unknown_cells,
      " base Unknown cells have both annotations"
    )
  ))
}
grDevices::dev.off()
pdf_open <- FALSE

message("Compared ", nrow(cell_annotations), " consensus Unknown cells across ",
        nrow(population_summary), " sample-qualified populations.")
message("Shared VZ/RL annotations yielded ", nrow(pair_table),
        " reviewable annotation-pair decisions.")
message("Review template: ", review_template_path)
message(
  "After scientific review, copy its decision rows to ",
  resolve_config_path(config$manifests$reviewed_unknown_relabels, config),
  " and fill final_label plus review_notes."
)
