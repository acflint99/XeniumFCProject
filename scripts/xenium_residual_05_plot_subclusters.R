#!/usr/bin/env Rscript
# Plot the pooled, already labeled residual object without rewriting analysis RDS files.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_05_plot_subclusters.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "spatial_report_helpers.R"))
source(here("scripts", "R", "residual_spatial_reports.R"))
source(here("scripts", "color_palette.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--dry-run", "--overwrite")) ||
    all(c("--dry-run", "--overwrite") %in% args)) stop("Use --dry-run or --overwrite, or no arguments.")
root <- here(config$project$outputs_dir, "xenium", "residual", "05_subclusters")
input <- file.path(root, "rds", "Xenium_Residual_Subclusters.rds")
review_path <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels, config)
review <- load_reviewed_residual_subcluster_labels(config)
panels <- residual_report_marker_sets(review, residual_shared_marker_panels())
paths <- c(tiff = file.path(root, "plots", "Xenium_Residual_SubclusterMarker_DotPlot.tif"),
           pdf = file.path(root, "plots", "Xenium_Residual_SubclusterMarker_DotPlot.pdf"),
           markers = file.path(root, "tables", "Xenium_Residual_SubclusterMarker_Availability.csv"),
           provenance = file.path(root, "tables", "Xenium_Residual_DotPlot_Provenance.csv"))
checks <- c(Seurat = requireNamespace("Seurat", quietly = TRUE),
            ggplot2 = requireNamespace("ggplot2", quietly = TRUE),
            patchwork = requireNamespace("patchwork", quietly = TRUE),
            Cairo = requireNamespace("Cairo", quietly = TRUE),
            cairo = capabilities("cairo"),
            taxonomy = all(names(panels) %in% intersect(master_subcluster_order, names(subcluster_palette))))
if ("--dry-run" %in% args) {
  compact_dry_run("Pooled residual subcluster DotPlot", inputs = c(input, review_path), outputs = paths, checks = checks)
  quit(save = "no", status = if (all(checks) && all(file.exists(c(input, review_path)))) 0L else 1L)
}
if (!all(checks)) stop("Missing plotting dependencies or taxonomy: ", paste(names(checks)[!checks], collapse = ", "))
if (!file.exists(input)) stop("Missing labeled residual object: ", input)
if (any(file.exists(paths)) && !"--overwrite" %in% args) stop("DotPlot outputs exist; review before using --overwrite.")
object <- readRDS(input)
metadata <- object[[]]
required <- c("residual_subcluster", "Xenium_snn_res.0.8")
if (!nrow(metadata) || anyDuplicated(rownames(metadata)) || !all(required %in% names(metadata))) stop("Invalid labeled residual metadata.")
index <- match(as.character(metadata$Xenium_snn_res.0.8), as.character(review$cluster))
if (anyNA(index) || !setequal(as.character(metadata$Xenium_snn_res.0.8), as.character(review$cluster)) ||
    !identical(as.character(metadata$residual_subcluster), as.character(review$final_label[index]))) stop("Residual labels do not match the current approved review.")
if (!"Xenium" %in% names(object@assays)) stop("Labeled residual object lacks Xenium assay.")
Seurat::DefaultAssay(object) <- "Xenium"
layers <- SeuratObject::Layers(object[["Xenium"]], search = "^data($|\\.)")
if (!length(layers)) stop("Xenium normalized data missing; this plotter does not normalize data.")
if (length(layers) > 1L) object <- SeuratObject::JoinLayers(object, assay = "Xenium", layers = "^data($|[.])", new = "data")
present <- intersect(master_subcluster_order, unique(as.character(metadata$residual_subcluster)))
requested <- panels[present]
availability <- data.frame(residual_subcluster = rep(names(requested), lengths(requested)),
                           gene = unlist(requested, use.names = FALSE), stringsAsFactors = FALSE)
availability$marker_source <- unname(attr(panels, "marker_sources")[availability$residual_subcluster])
availability$present <- availability$gene %in% rownames(object[["Xenium"]])
features <- lapply(requested, intersect, y = rownames(object[["Xenium"]]))
if (any(!lengths(features))) stop("No available markers for: ", paste(names(features)[!lengths(features)], collapse = ", "))
object$residual_subcluster <- factor(metadata$residual_subcluster, levels = rev(present))
Seurat::Idents(object) <- "residual_subcluster"
# A null measurement device prevents grid conversions creating Rplots.pdf.
grDevices::pdf(file = NULL)
plot <- workflow_dotplot(object, features = features, assay = "Xenium", idents = rev(present),
                         cols = c("lightgrey", "red"), dot.scale = 6, cluster.idents = FALSE,
                         scale = length(present) > 1L) + Seurat::RotatedAxis() +
  ggplot2::ggtitle("Xenium Merged Residual Subclusters")
dimensions <- workflow_dotplot_dimensions(plot)
display <- patchwork::wrap_elements(full = workflow_dotplot_grob(plot, strip_colors = subcluster_palette[present]))
for (directory in unique(dirname(paths))) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
save_cairo_tiff(display, paths[["tiff"]], dimensions[["width"]], dimensions[["height"]], background = "white")
ggplot2::ggsave(paths[["pdf"]], display, device = grDevices::cairo_pdf,
                width = dimensions[["width"]], height = dimensions[["height"]], bg = "white")
grDevices::dev.off()
write.csv(availability, paths[["markers"]], row.names = FALSE)
if (any(is.na(file.info(paths[1:3])$size) | file.info(paths[1:3])$size <= 0)) stop("Incomplete DotPlot outputs.")
write.csv(data.frame(input_rds = normalizePath(input), input_size = file.info(input)$size,
  input_mtime = as.numeric(file.info(input)$mtime), cells = nrow(metadata), identities = length(present),
  review_md5 = unname(tools::md5sum(review_path)), marker_panel_md5 = residual_marker_panel_signature(panels),
  tiff_path = paths[["tiff"]], pdf_path = paths[["pdf"]], stringsAsFactors = FALSE),
  paths[["provenance"]], row.names = FALSE)
message("Completed pooled residual DotPlot: ", paths[["pdf"]])
