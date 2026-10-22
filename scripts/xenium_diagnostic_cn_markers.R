#!/usr/bin/env Rscript
# Descriptive CN identity evidence. Never changes labels or saves a Seurat object.
# Expected HPC runtime: 15-90 min merged, 10-60 min/regional branch; low confidence.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_diagnostic_cn_markers.R")
source(here("scripts", "R", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
branch_arg <- args[grepl("^--branch=", args)]
if (length(branch_arg) > 1L) stop("Specify --branch only once.")
branch <- if (length(branch_arg)) sub("^--branch=", "", branch_arg) else "nonspatial"
args <- args[!grepl("^--branch=", args)]
if (!branch %in% c("nonspatial", "vz", "rl", "residual") ||
    any(!args %in% c("--dry-run", "--overwrite")) ||
    all(c("--dry-run", "--overwrite") %in% args)) {
  stop("Usage: xenium_diagnostic_cn_markers.R [--branch=nonspatial|vz|rl|residual] [--dry-run|--overwrite]")
}
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
stage <- if (branch == "residual") "05_subclusters" else "06_subclusters"
rds_name <- switch(branch, vz = "Xenium_VZ_subclusters.rds",
                   rl = "Xenium_RL_subclusters.rds", residual = "Xenium_Residual_Subclusters.rds",
                   nonspatial = "Xenium_VZRL_merged.rds")
label_column <- switch(branch, vz = "VZ_subcluster", rl = "RL_subcluster",
                       residual = "residual_subcluster", nonspatial = "comb_subcluster")
cluster_column <- "Xenium_snn_res.0.8"
input <- here(config$project$outputs_dir, "xenium", branch, stage, "rds", rds_name)
root <- here(config$project$outputs_dir, "xenium", branch, "diagnostics", "canonical_cn_markers_gaba_ecn_mli")
if (branch == "nonspatial") {
  input <- here(config$project$outputs_dir, "xenium", "vz_rl", "nonspatial",
                "03_merged", "rds", "Xenium_VZRL_merged.rds")
  root <- here(config$project$outputs_dir, "xenium", "vz_rl", "nonspatial",
               "diagnostics", "canonical_cn_markers_gaba_ecn_mli")
  label_column <- "comb_subcluster"
  # Regional resolution-0.8 IDs are not propagated by the combined-label stage.
  # Never interpret a whole-tissue clustering column as a regional cluster ID.
  cluster_column <- "comb_subcluster"
}
# Existing broad GABA-family subclusters; PCs are a separate project lineage.
# This is a report selection only, never a change to cell identities or QC.
selected_groups <- c("eCN", "iCN", "MLIs", "Golgi Cells", "GABA Interneurons", "GABA Progenitors")
paths <- file.path(root, c("marker_availability.csv", "cluster_expression.csv",
                          "sample_cluster_expression.csv", "sample_balanced_expression.csv",
                          "coexpression.csv", "sample_cluster_counts.csv", "sample_coverage.csv",
                          "provenance.csv", "interpretation.txt", "CN_markers_DotPlot.pdf", "sessionInfo.txt",
                          "CN_identity_Focused_DotPlot.pdf", "CN_identity_Focused_Expression.csv"))
names(paths) <- c("availability", "pooled", "sample", "balanced", "pairs", "counts",
                  "coverage", "provenance", "notes", "plot", "session", "focused_plot", "focused_expression")

# Stronger CN candidates are separated from broadly shared programs. These are
# review panels, not classifiers; mouse evidence is not assumed human-specific.
panels <- list(
  "iCN candidates" = c("SOX14", "DMBX1"),
  "iCN supporting (shared)" = c("ZFHX3", "ZFHX4", "TOX", "SLC24A4", "CHRM2", "MEIS1", "TAC1"),
  "GABA transmission / lineage" = c("GAD1", "GAD2", "SLC32A1", "SLC6A1", "LHX1", "LHX5"),
  "eCN core" = c("SLC17A6", "MEIS2"),
  "eCN subtype / developmental" = c("TBR1", "LHX9", "LMX1A", "BARHL1", "NEUROD6"),
  "eCN supporting (shared)" = c("NTNG1", "ZNF804A", "NR4A2", "NRP1", "C1QL3", "NTNG2", "CREG2", "NEUROD2"),
  "MLI markers (stage / subtype dependent)" = c("PVALB", "SORCS3", "NXPH1", "KCNA2", "GRIK1", "KIT", "ERBB4", "GRM8", "CACNA1E"),
  "GABA subtype context (shared)" = c("PAX2", "TFAP2A", "TFAP2B", "SLC6A5", "NXPH2", "GRIK3"),
  "Pan-neuronal (not CN-specific)" = c("RBFOX3", "TUBB3", "NEFL", "NEFM")
)
availability <- data.frame(panel = rep(names(panels), lengths(panels)),
                           gene = unlist(panels, use.names = FALSE))
markers <- availability$gene
stopifnot(!anyDuplicated(markers))
pair_spec <- data.frame(
  gene1 = c("SOX14", "SOX14", "DMBX1", "SLC17A6", "SLC17A6", "SOX14", "ZFHX4", "SLC17A6", "SLC17A6"),
  gene2 = c("GAD1", "GAD2", "GAD1", "MEIS2", "GAD1", "SLC17A6", "GAD1", "NEUROD6", "NTNG1"),
  purpose = c("iCN support", "iCN support", "iCN support", "eCN support",
              "transmitter conflict review", "CN program conflict review", "iCN shared support",
              "eCN developmental support", "eCN shared support"))
packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2")
checks <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
checks <- c(checks, sample_manifest = nrow(samples) > 0L && !anyNA(samples$sample_id) &&
              all(nzchar(samples$sample_id)) && !anyDuplicated(samples$sample_id))
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("Canonical CN markers (paths only; no RDS loading)", input, paths, checks)
  cat("Branch: ", branch, "; manifest samples: ", nrow(samples),
      "; grouping: ", if (branch == "nonspatial") "existing label pooled across origins" else cluster_column,
      "\n", sep = "")
  if (branch == "nonspatial") cat("Selected labels (absent labels are reported): ", paste(selected_groups, collapse = "; "), "\n", sep = "")
  cat("Input: ", input, "\nOutput directory: ", root, "\n", sep = "")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!all(checks)) stop("Missing dependency or invalid manifest: ", paste(names(checks)[!checks], collapse = ", "))
if (!file.exists(input)) stop("Missing annotated input: ", input)
if (any(file.exists(paths)) && !"--overwrite" %in% args) stop("Outputs exist; review before using --overwrite.")
message(Sys.time(), ": loading ", input)
object <- readRDS(input)
metadata <- object[[]]
source_cells <- nrow(metadata)
source_sample_counts <- table(factor(metadata$orig.ident, levels = samples$sample_id))
required <- if (branch == "nonspatial") {
  c("orig.ident", "comb_subcluster", "regional_assignment")
} else c("orig.ident", label_column, cluster_column)
if (!nrow(metadata) || !all(required %in% names(metadata)) ||
    anyDuplicated(colnames(object)) || !identical(rownames(metadata), colnames(object))) {
  stop("Expected unique aligned cell IDs and metadata: ", paste(required, collapse = ", "))
}
for (column in required) {
  if (anyNA(metadata[[column]]) || any(!nzchar(trimws(as.character(metadata[[column]]))))) {
    stop("Missing/blank values in ", column)
  }
}
if (!all(metadata$orig.ident %in% samples$sample_id)) stop("Object sample IDs absent from active manifest.")
if (!any(metadata[[label_column]] %in% c("iCN", "eCN"))) stop("No iCN/eCN identities in this branch.")
if (branch == "vz" && !all(c("iCN", "eCN") %in% metadata[[label_column]])) stop("VZ requires both iCN and eCN labels.")
# A source cluster must map to exactly one existing label, with no metadata join.
if (branch == "nonspatial") {
  if (!setequal(unique(as.character(metadata$orig.ident)), samples$sample_id)) {
    stop("Merged object does not contain exactly the active manifest samples; rebuild/review the merge.")
  }
  if (!all(metadata$regional_assignment %in% c("VZ", "RL", "Residual", "unassigned_consensus"))) {
    stop("Invalid regional_assignment in merged object.")
  }
  branch_columns <- c(VZ = "VZ_subcluster", RL = "RL_subcluster", Residual = "residual_subcluster")
  if (!all(branch_columns %in% names(metadata))) stop("Merged object lacks source branch label columns.")
  present <- vapply(branch_columns, function(column) {
    !is.na(metadata[[column]]) & nzchar(trimws(as.character(metadata[[column]])))
  }, logical(nrow(metadata)))
  if (any(rowSums(present) > 1L)) stop("Overlapping source branch assignments.")
  for (origin in names(branch_columns)) {
    selected <- as.character(metadata$regional_assignment) == origin
    if (!identical(unname(present[, origin]), selected) ||
        any(as.character(metadata[[branch_columns[[origin]]]][selected]) !=
              as.character(metadata$comb_subcluster[selected]))) {
      stop("Combined label/source assignment mismatch for ", origin)
    }
  }
  if (!all(c("iCN", "eCN") %in% metadata$comb_subcluster) ||
      !any(metadata$regional_assignment == "Residual" & metadata$comb_subcluster == "eCN")) {
    stop("Expected iCN/eCN and Residual-derived eCN are missing; check merge provenance before interpreting.")
  }
  metadata$cluster <- as.character(metadata$comb_subcluster)
  # Restrict the in-memory report to the requested labels, pooling origins.
  metadata <- metadata[metadata$comb_subcluster %in% selected_groups, , drop = FALSE]
  missing_groups <- setdiff(selected_groups, unique(as.character(metadata$comb_subcluster)))
  if (length(missing_groups)) message("Requested labels absent from input: ", paste(missing_groups, collapse = ", "))
} else {
  metadata$cluster <- as.character(metadata[[cluster_column]])
}
mapping <- unique(data.frame(cluster = metadata$cluster,
                             label = as.character(metadata[[label_column]])))
if (anyDuplicated(mapping$cluster)) stop("Source cluster maps to multiple labels.")
mapping <- mapping[order(suppressWarnings(as.numeric(mapping$cluster)), mapping$cluster), ]
if (branch == "nonspatial") mapping <- mapping[order(match(mapping$label, selected_groups)), ]
metadata$label <- as.character(metadata[[label_column]])
metadata$sample_id <- as.character(metadata$orig.ident)
if (!"Xenium" %in% SeuratObject::Assays(object)) stop("Missing Xenium assay.")
availability$present <- markers %in% rownames(object[["Xenium"]])
features <- markers[availability$present]
if (!length(features)) stop("No requested markers in panel.")
if (any(!availability$present)) message("Unavailable genes: ", paste(markers[!availability$present], collapse = ", "))
SeuratObject::DefaultAssay(object) <- "Xenium"
# Report selection only; original RDS and all metadata remain unchanged.
object <- subset(object, cells = rownames(metadata), features = features)
matrices <- lapply(c("counts", "data"), function(layer) {
  layers <- SeuratObject::Layers(object[["Xenium"]], search = paste0("^", layer, "($|[.])"))
  if (!length(layers)) stop("Missing Xenium ", layer, " layer; no normalization is performed.")
  assay <- object[["Xenium"]]
  if (length(layers) > 1L) assay <- SeuratObject::JoinLayers(assay, layers = paste0("^", layer, "($|[.])"), new = layer)
  x <- SeuratObject::LayerData(assay, layer = layer)
  if (anyDuplicated(colnames(x)) || !setequal(colnames(x), rownames(metadata)) ||
      !all(features %in% rownames(x))) stop("Incomplete or duplicate layer cells/features: ", layer)
  x <- x[features, rownames(metadata), drop = FALSE]
  values <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (any(!is.finite(values)) || any(values < 0)) stop("Non-finite/negative ", layer, " values.")
  x
})
names(matrices) <- c("counts", "data")
rm(object)
gc(verbose = FALSE)
detected <- matrices$counts > 0
linear <- expm1(matrices$data)
if (any(!is.finite(if (inherits(linear, "sparseMatrix")) linear@x else as.vector(linear)))) stop("Non-finite linear expression.")

# Sparse indicator aggregation avoids cell-wise statistical pseudoreplication.
summarize_groups <- function(columns) {
  keys <- unique(metadata[, columns, drop = FALSE])
  # Match by integer factor codes, avoiding delimiters in real cell/sample IDs.
  make_key <- function(frame) do.call(paste, c(lapply(columns, function(n) match(frame[[n]], unique(metadata[[n]]))), sep = ":"))
  index <- match(make_key(metadata), make_key(keys))
  if (anyNA(index)) stop("Incomplete grouping map.")
  design <- Matrix::sparseMatrix(i = seq_len(nrow(metadata)), j = index, x = 1,
                                  dims = c(nrow(metadata), nrow(keys)))
  n <- as.numeric(Matrix::colSums(design))
  means <- as.matrix(linear %*% design)
  pct <- as.matrix(detected %*% design)
  result <- keys[rep(seq_len(nrow(keys)), each = length(features)), , drop = FALSE]
  result$gene <- rep(features, nrow(keys))
  result$cells <- rep(n, each = length(features))
  result$mean_linear <- as.vector(means) / result$cells
  result$log1p_mean <- log1p(result$mean_linear)
  result$pct_detected <- 100 * as.vector(pct) / result$cells
  if (any(!is.finite(as.matrix(result[, c("mean_linear", "log1p_mean", "pct_detected")])))) stop("Non-finite summaries.")
  result
}
pooled <- summarize_groups(c("cluster", "label"))
by_sample <- summarize_groups(c("sample_id", "cluster", "label"))
balanced <- aggregate(by_sample[, c("mean_linear", "pct_detected")],
                      by_sample[, c("cluster", "label", "gene")], mean)
replicates <- aggregate(sample_id ~ cluster + label + gene, by_sample, length)
if (!identical(balanced[, c("cluster", "label", "gene")], replicates[, c("cluster", "label", "gene")])) stop("Replicate summary alignment failed.")
balanced$samples_present <- replicates$sample_id
balanced$log1p_mean <- log1p(balanced$mean_linear)
counts <- unique(by_sample[, c("sample_id", "cluster", "label", "cells")])
if (sum(counts$cells) != nrow(metadata)) stop("Cell counts do not reconcile.")
pair_spec$available <- pair_spec$gene1 %in% features & pair_spec$gene2 %in% features
pair_rows <- lapply(seq_len(nrow(pair_spec)), function(i) {
  spec <- pair_spec[i, ]
  frame <- metadata[, c("sample_id", "cluster", "label")]
  if (spec$available) {
    a <- as.numeric(detected[spec$gene1, ]) > 0
    b <- as.numeric(detected[spec$gene2, ]) > 0
    frame$gene1_cells <- as.integer(a)
    frame$gene2_cells <- as.integer(b)
    frame$both_cells <- as.integer(a & b)
    frame$cells <- 1L
    out <- aggregate(frame[, c("gene1_cells", "gene2_cells", "both_cells", "cells")],
                     frame[, c("sample_id", "cluster", "label")], sum)
    out$pct_both <- 100 * out$both_cells / out$cells
  } else {
    out <- counts
    out$gene1_cells <- out$gene2_cells <- out$both_cells <- out$pct_both <- NA_real_
  }
  for (name in names(spec)) out[[name]] <- spec[[name]]
  out[, c("sample_id", "cluster", "label", "cells", "gene1_cells", "gene2_cells",
          "both_cells", "pct_both", "gene1", "gene2", "purpose", "available")]
})
pair_table <- do.call(rbind, pair_rows)
coverage <- data.frame(sample_id = samples$sample_id,
                       cells_in_branch = as.integer(table(factor(metadata$sample_id, levels = samples$sample_id))))
coverage$status <- ifelse(coverage$cells_in_branch > 0, "present", "absent from input branch; not treated as zero expression")
if (branch == "nonspatial") {
  names(coverage)[names(coverage) == "cells_in_branch"] <- "selected_cells"
  coverage$source_cells <- as.integer(source_sample_counts)
  coverage$status <- ifelse(coverage$selected_cells > 0, "selected cells present", "no requested cells; not treated as zero expression")
}
pooled$panel <- availability$panel[match(pooled$gene, availability$gene)]
pooled$display_group <- paste0(pooled$cluster, " | ", pooled$label)
display_levels <- paste0(mapping$cluster, " | ", mapping$label)
if (branch == "nonspatial") {
  pooled$display_group <- pooled$cluster
  display_levels <- mapping$cluster
}
pooled$display_group <- factor(pooled$display_group,
                               levels = rev(display_levels))
pooled$gene <- factor(pooled$gene, levels = features)
pooled$panel <- factor(pooled$panel, levels = names(panels))
plot <- ggplot2::ggplot(pooled, ggplot2::aes(gene, display_group)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_detected, color = log1p_mean)) +
  ggplot2::facet_grid(cols = ggplot2::vars(panel), scales = "free_x", space = "free_x") +
  ggplot2::scale_color_gradient(low = "grey90", high = "#B2182B", name = "log1p(mean\nnormalized expression)") +
  ggplot2::scale_size_area(max_size = 5, limits = c(0, 100), name = "% raw count > 0") +
  ggplot2::labs(title = paste(toupper(branch), "canonical CN marker evidence"),
                subtitle = if (branch == "nonspatial") "eCN and GABA subclusters pooled across all origins; unscaled expression" else "Regional subclusters; unscaled expression",
                x = NULL, y = if (branch == "nonspatial") "Combined subcluster" else "Source cluster | existing label") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                 strip.text.x = ggplot2::element_text(angle = 90), panel.grid = ggplot2::element_blank())
# A compact, sample-balanced view makes relative marker enrichment legible.
# Plot scaling changes presentation only; it never assigns identities.
focused_sets <- list(
  "eCN core" = c("SLC17A6", "MEIS2"),
  "eCN support" = c("NEUROD6", "NTNG1", "NR4A2", "TBR1", "LHX9"),
  "iCN core" = c("SOX14", "DMBX1"),
  "iCN support" = c("ZFHX4", "SLC24A4", "TOX"),
  "GABA" = c("GAD1", "GAD2", "SLC32A1"),
  "MLI" = panels[["MLI markers (stage / subtype dependent)"]]
)
focused_genes <- unlist(focused_sets, use.names = FALSE)
focused <- balanced[balanced$gene %in% focused_genes, , drop = FALSE]
if (!nrow(focused)) stop("No focused CN markers available; see requested marker panel.")
focused$relative_z <- ave(focused$log1p_mean, focused$gene, FUN = function(values) {
  deviation <- stats::sd(values)
  if (!is.finite(deviation) || deviation == 0) return(rep(0, length(values)))
  (values - mean(values)) / deviation
})
if (any(!is.finite(focused$relative_z))) stop("Non-finite focused expression scaling.")
focused$panel <- factor(rep(names(focused_sets), lengths(focused_sets))[
  match(focused$gene, focused_genes)], levels = names(focused_sets))
focused$gene <- factor(focused$gene, levels = focused_genes[focused_genes %in% features])
focused$display_group <- factor(if (branch == "nonspatial") focused$label else
                                 paste0(focused$cluster, " | ", focused$label), levels = rev(display_levels))
focused_plot <- ggplot2::ggplot(focused, ggplot2::aes(gene, display_group)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_detected, color = relative_z)) +
  ggplot2::facet_grid(cols = ggplot2::vars(panel), scales = "free_x", space = "free_x") +
  ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, name = "Relative expression\n(z-score per gene)") +
  ggplot2::scale_size_area(max_size = 7, limits = c(0, 100), name = "Mean sample\n% detected") +
  ggplot2::labs(title = "eCN, iCN and MLI identity evidence",
                subtitle = "Equal sample weight; color shows relative enrichment, size shows detection",
                x = NULL, y = NULL,
                caption = "Shared and developmental markers support review; relative enrichment alone does not establish identity.") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                 strip.background = ggplot2::element_rect(fill = "grey95"),
                 panel.grid = ggplot2::element_blank())
notes <- c(
  "Interpretation: marker evidence supports or challenges existing identities; it does not prove anatomy or projections.",
  if (branch == "nonspatial") {
    "Merged report groups comb_subcluster only: eCN, GABA Progenitors, GABA Interneurons, Golgi Cells, MLIs, iCN. Cells from VZ, RL and Residual are pooled within each label in plots, means, coexpression and sample-level tables. Other labels are excluded from this report only. Original numeric regional cluster IDs are not inferred."
  } else "Inspect source clusters separately (especially VZ iCN clusters 11 and 21), not only merged labels.",
  "iCN: inspect SOX14/DMBX1 together with GAD1/GAD2/SLC32A1; GABA markers alone cannot separate CN from cortical interneurons or PCs.",
  "SOX14 marks a subset of inhibitory projection CN in mouse; absence does not exclude every inhibitory CN subtype or developmental state.",
  "iCN supporting candidates ZFHX3/ZFHX4/TOX/SLC24A4 and project markers CHRM2/MEIS1/TAC1 are shared and require SOX14/DMBX1 plus transmitter context; they are not iCN-specific classifiers.",
  "eCN: inspect SLC17A6 with MEIS2. TBR1/LHX9/LMX1A/BARHL1/NEUROD6 are subtype/development-dependent, not universal requirements. NTNG1/ZNF804A/NR4A2/NRP1 and project C1QL3/NTNG2/CREG2/NEUROD2 are shared supporting evidence.",
  "The MLI panel replaces Purkinje contrast genes: PVALB, SORCS3, NXPH1, KCNA2, GRIK1, KIT, ERBB4, GRM8, CACNA1E. MLI markers vary with maturation/subtype and are not individually exclusive. SORCS3 and NXPH1 mark different MLI transcriptional types; do not require their coexpression or equate them directly to basket/stellate morphology.",
  "CN_identity_Focused_DotPlot includes eCN core/support, iCN core/support, GABA and MLI panels, using equal-sample mean expression/detection within each label. Color is a gene-wise z-score of log1p(sample-balanced mean linear expression) across the displayed groups; constant genes are zero. A red dot means relative enrichment, even for weak genes: consult dot size, samples_present and unscaled values in CN_identity_Focused_Expression.csv. The main CN_markers_DotPlot remains unscaled and cell-pooled. Neither plot tests identity or significance.",
  "NEFL/NEFM and other generic neuronal markers are insufficient to establish CN identity.",
  "Missing panel genes are unavailable evidence, not negative expression. A present gene with no detection is distinct from a missing gene.",
  "Coexpression is raw count >0 in the same segmented cell; it cannot distinguish biology from segmentation contamination/doublets.",
  "Sample-balanced means give each represented biological sample equal weight; no donor independence or adjusted statistical significance is assumed.",
  "Samples absent from a branch and sample-cluster combinations with no cells are not imputed as zero expression. Review sample coverage/counts.",
  "Normalized data are assumed existing log1p nonnegative normalized Xenium expression; report cells are selected by existing labels only. No normalization, QC filtering, relabeling or stochastic analysis is performed.",
  "No identity thresholds, pass/fail calls or cell-level significance tests are applied. Spatial nuclear localization and independent reference evidence strengthen an identity claim.",
  "Sources (marker specificity differs by species and developmental stage):",
  "Human developing GABA CN/SOX14: https://pubmed.ncbi.nlm.nih.gov/40249772/",
  "Mouse SOX14 inhibitory CN subset: https://pmc.ncbi.nlm.nih.gov/articles/PMC6706002/",
  "Mouse eCN MEIS2/SLC17A6: https://elifesciences.org/articles/50617/figures",
  "Human fetal eCN context: https://www.nature.com/articles/s41467-023-43568-6",
  "Expanded shared CN candidates (human developmental/organoid comparison): https://doi.org/10.1242/dmm.052290",
  "GABA subtype context (human development): https://pubmed.ncbi.nlm.nih.gov/40249772/",
  "MLI transcriptional types (mouse): https://www.nature.com/articles/s41586-021-03220-z",
  "MLI morphology and SORCS3/NXPH1/GRM8/CACNA1E (mouse): https://www.nature.com/articles/s41467-022-30977-2",
  "eCN developmental subtypes (mouse and comparative): https://www.nature.com/articles/s41586-023-06884-x",
  "Existing exploratory supporting/contrast panels: scripts/R/subcluster_review_markers.R and scripts/color_palette.R"
)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
for (name in c("availability", "pooled", "sample", "balanced", "pairs", "counts", "coverage")) {
  value <- switch(name, availability = availability, pooled = pooled, sample = by_sample,
                  balanced = balanced, pairs = pair_table, counts = counts, coverage = coverage)
  write.csv(value, paths[[name]], row.names = FALSE)
}
writeLines(notes, paths[["notes"]])
ggplot2::ggsave(paths[["plot"]], plot, device = "pdf", width = max(16, length(features) * 0.3),
                height = max(7, nrow(mapping) * 0.25 + 4), limitsize = FALSE)
ggplot2::ggsave(paths[["focused_plot"]], focused_plot, device = "pdf", width = max(12, length(unique(focused$gene)) * 0.55),
                height = max(4.5, nrow(mapping) * 0.3 + 2.5), limitsize = FALSE)
write.csv(focused, paths[["focused_expression"]], row.names = FALSE)
write.csv(data.frame(input_rds = normalizePath(input), input_size = file.info(input)$size,
                     input_mtime = as.character(file.info(input)$mtime), branch = branch,
                     cells = nrow(metadata), source_cells = source_cells, manifest_samples = nrow(samples),
                     requested_labels = if (branch == "nonspatial") paste(selected_groups, collapse = ";") else "all regional labels",
                     evaluated_labels = paste(unique(metadata$label), collapse = ";"),
                     grouping = if (branch == "nonspatial") "comb_subcluster pooled across origins" else cluster_column,
                     represented_samples = length(unique(metadata$sample_id)),
                     sample_manifest_md5 = unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))),
                     script_md5 = unname(tools::md5sum(here("scripts", "xenium_diagnostic_cn_markers.R"))),
                     cluster_column = cluster_column, label_column = label_column,
                     assay = "Xenium", layers = "counts;data", detection = "raw count > 0",
                     requested_markers = paste(markers, collapse = ";"),
                     available_markers = paste(features, collapse = ";"),
                     focused_markers = paste(focused_genes[focused_genes %in% features], collapse = ";"),
                     focused_scaling = "gene-wise z-score across displayed groups of log1p equal-sample mean expression; constant genes=0",
                     sample_weighting = "equal represented sample weight within each cluster; no inference",
                     created = as.character(Sys.time())), paths[["provenance"]], row.names = FALSE)
writeLines(capture.output(sessionInfo()), paths[["session"]])
if (any(is.na(file.info(paths)$size) | file.info(paths)$size <= 0)) stop("Incomplete diagnostic outputs.")
message("Completed CN marker evidence: ", root)
