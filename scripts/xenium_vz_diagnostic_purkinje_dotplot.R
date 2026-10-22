#!/usr/bin/env Rscript
# Read-only expression summary of existing Purkinje labels; no reclustering.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_vz_diagnostic_purkinje_dotplot.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "dotplot_helpers.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite")) ||
    all(c("--dry-run", "--overwrite") %in% args)) {
  stop("Use --dry-run or --overwrite, or no arguments.")
}
config <- load_pipeline_config()
samples <- load_sample_manifest(config)
# Fully legible genes only, left to right in the user-supplied cropped image.
markers <- c("VSTM2L", "COL5A1", "ITPR1", "CNTNAP4", "NDRG1", "KCNC3",
             "NEFL", "NXPH4", "PCDH10", "BCL11A", "NPTX2", "RFTN1",
             "EBF2", "EBF1", "RORB")
groups <- c("Maturing PCs", "Early-born PCs", "Late-born PCs", "Patterning PCs")
input <- here(config$project$outputs_dir, "xenium", "vz", "06_subclusters",
              "rds", "Xenium_VZ_subclusters.rds")
root <- here(config$project$outputs_dir, "xenium", "vz", "diagnostics", "purkinje_reference_markers")
prefix <- "Xenium_VZ_Purkinje_ReferenceMarker"
paths <- c(
  pdf = file.path(root, "plots", paste0(prefix, "_DotPlot.pdf")),
  tiff = file.path(root, "plots", paste0(prefix, "_DotPlot.tif")),
  availability = file.path(root, "tables", paste0(prefix, "_Availability.csv")),
  expression = file.path(root, "tables", paste0(prefix, "_Expression.csv")),
  counts = file.path(root, "tables", paste0(prefix, "_CellCounts.csv")),
  provenance = file.path(root, "tables", paste0(prefix, "_Provenance.csv"))
)
packages <- c("Seurat", "SeuratObject", "ggplot2", "scales", "Cairo")
checks <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
checks <- c(checks, cairo = capabilities("cairo"),
            sample_ids_valid = nrow(samples) > 0L && !anyNA(samples$sample_id) &&
              all(nzchar(samples$sample_id)) && !anyDuplicated(samples$sample_id))
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("Purkinje reference-marker DotPlot (paths only; no RDS loading)",
                        inputs = input, outputs = paths, checks = checks)
  cat("Groups: ", paste(groups, collapse = "; "), "\nMarkers: ",
      paste(markers, collapse = ", "), "\n", sep = "")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (!all(checks)) stop("Failed checks: ", paste(names(checks)[!checks], collapse = ", "))
if (!file.exists(input)) stop("Missing annotated VZ object: ", input)
if (any(file.exists(paths)) && !"--overwrite" %in% args) {
  stop("Diagnostic outputs already exist; review before using --overwrite.")
}
message(Sys.time(), ": loading annotated VZ object: ", input)
object <- readRDS(input)
metadata <- object[[]]
required <- c("VZ_subcluster", "orig.ident", "Xenium_snn_res.0.8")
if (!nrow(metadata) || !all(required %in% names(metadata)) ||
    anyDuplicated(colnames(object)) || !identical(rownames(metadata), colnames(object))) {
  stop("Invalid VZ object: expected unique, aligned cells and metadata: ", paste(required, collapse = ", "))
}
if (anyNA(metadata$VZ_subcluster) || any(!nzchar(as.character(metadata$VZ_subcluster))) ||
    anyNA(metadata$orig.ident) || !all(metadata$orig.ident %in% samples$sample_id)) {
  stop("Missing identities or sample IDs absent from the active sample manifest.")
}
if (!all(groups %in% metadata$VZ_subcluster)) {
  stop("Expected Purkinje groups missing: ", paste(setdiff(groups, metadata$VZ_subcluster), collapse = ", "))
}
if (!"Xenium" %in% names(object@assays)) stop("Missing Xenium assay.")
availability <- data.frame(order = seq_along(markers), gene = markers,
                           present = markers %in% rownames(object[["Xenium"]]))
features <- markers[availability$present]
if (!length(features)) stop("None of the reference markers exists in the Xenium assay.")
if (any(!availability$present)) {
  message("Absent panel genes (omitted, recorded in availability CSV): ",
          paste(markers[!availability$present], collapse = ", "))
}
pc_metadata <- metadata[metadata$VZ_subcluster %in% groups, required, drop = FALSE]
if (anyNA(pc_metadata$Xenium_snn_res.0.8) ||
    any(!nzchar(as.character(pc_metadata$Xenium_snn_res.0.8)))) stop("Missing PC source-cluster assignments.")
Seurat::DefaultAssay(object) <- "Xenium"
# Restrict the in-memory copy before joining normalized layers. Never save it.
object <- subset(object, cells = rownames(pc_metadata), features = features)
layers <- SeuratObject::Layers(object[["Xenium"]], search = "^data($|[.])")
if (!length(layers)) stop("Normalized Xenium data are missing; this diagnostic does not normalize.")
if (length(layers) > 1L) {
  object <- SeuratObject::JoinLayers(object, assay = "Xenium", layers = "^data($|[.])", new = "data")
}
normalized <- SeuratObject::LayerData(object, assay = "Xenium", layer = "data")
if (!all(features %in% rownames(normalized)) ||
    !setequal(colnames(normalized), rownames(pc_metadata))) stop("Normalized layer lacks requested genes/cells.")
rm(normalized)
object$VZ_subcluster <- factor(as.character(object$VZ_subcluster), levels = rev(groups))
Seurat::Idents(object) <- "VZ_subcluster"
grDevices::pdf(file = NULL)
plot <- workflow_dotplot(object, features = features, assay = "Xenium",
                         idents = rev(groups), cluster.idents = FALSE, scale = TRUE) +
  ggplot2::labs(title = "Purkinje subclusters", x = NULL, y = NULL,
                subtitle = "Reference markers available in Xenium; existing VZ labels")
if (any(!is.finite(as.matrix(plot$data[, c("avg.exp", "pct.exp", "avg.exp.scaled")])))) {
  stop("Non-finite DotPlot expression summaries; no outputs written.")
}
if (nrow(plot$data) != length(features) * length(groups)) stop("Incomplete gene-by-group expression summary.")
dimensions <- workflow_dotplot_dimensions(plot)
counts <- as.data.frame(table(pc_metadata$orig.ident, pc_metadata$VZ_subcluster,
                              pc_metadata$Xenium_snn_res.0.8), stringsAsFactors = FALSE)
names(counts) <- c("sample_id", "VZ_subcluster", "source_cluster", "cells")
counts <- counts[counts$cells > 0L, , drop = FALSE]
if (sum(counts$cells) != ncol(object)) stop("PC cell counts do not reconcile.")
for (directory in unique(dirname(paths))) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(paths[["pdf"]], plot, device = grDevices::cairo_pdf,
                width = dimensions[["width"]], height = dimensions[["height"]], bg = "white")
save_cairo_tiff(plot, paths[["tiff"]], dimensions[["width"]], dimensions[["height"]], background = "white")
grDevices::dev.off()
write.csv(availability, paths[["availability"]], row.names = FALSE)
write.csv(plot$data, paths[["expression"]], row.names = FALSE)
write.csv(counts, paths[["counts"]], row.names = FALSE)
write.csv(data.frame(input_rds = normalizePath(input), input_size = file.info(input)$size,
                     input_mtime = as.numeric(file.info(input)$mtime),
                     sample_manifest_md5 = unname(tools::md5sum(resolve_config_path(config$manifests$samples, config))),
                     source_cells = nrow(metadata), plotted_cells = ncol(object),
                     assay = "Xenium", layer = "data", group_by = "VZ_subcluster",
                     groups = paste(groups, collapse = ";"), markers = paste(markers, collapse = ";"),
                     marker_source = "15 fully visible genes in user reference image; cropped right edge excluded",
                     scaling = "Gene-wise z-score across four plotted PC groups, clipped to [-2.5,2.5]",
                     pdf_path = paths[["pdf"]], tiff_path = paths[["tiff"]]),
          paths[["provenance"]], row.names = FALSE)
if (any(is.na(file.info(paths)$size) | file.info(paths)$size <= 0)) stop("Incomplete diagnostic outputs.")
message("Completed: ", paths[["pdf"]])
