#!/usr/bin/env Rscript
# Repaint saved reference objects into their original plot paths. No RDS writes.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/refresh_reference_plots_original.R")
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "palette_refresh.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
flags <- c("--list", "--dry-run", "--overwrite")
if (any(!args %in% flags & !grepl("^[0-9]+$", args))) stop("Use [--list|--dry-run|--overwrite] [TASK_ID].")
if (sum(c("--list", "--dry-run", "--overwrite") %in% args) > 1L) stop("Choose one mode.")
references <- names(config$inputs$references)
if (!setequal(references, c("aldinger", "sepp", "science"))) stop("Review the reference export map before changing references.")
if ("--list" %in% args) {
  utils::write.table(data.frame(task_id=seq_along(references), reference=references), row.names=FALSE, quote=FALSE, sep="\t")
  quit(save="no", status=0L)
}
task_args <- args[grepl("^[0-9]+$", args)]
if (length(task_args)>1L) stop("Specify at most one TASK_ID.")
task <- if (length(task_args)) as.integer(task_args) else suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))
if (is.na(task) && !"--dry-run" %in% args) stop("TASK_ID is required; inspect --list.")
if (!is.na(task) && !task %in% seq_along(references)) stop("TASK_ID is outside --list.")
selected <- if (is.na(task)) references else references[[task]]

reference_specs <- function(reference) {
  panel <- resolve_config_path(config$inputs$references[[reference]], config)
  rds <- dirname(panel)
  plots <- file.path(dirname(rds), "plots")
  prefix <- switch(reference, aldinger="Aldinger", sepp="Sepp", science="Science")
  full <- file.path(rds, paste0(prefix, "_newClusters_newUMAPv1.rds"))
  original <- file.path(rds, paste0(prefix, "_newClusters.rds"))
  prepanel <- file.path(rds, paste0(prefix, if(reference=="science") "_newClusters_UMAPv1_5k.rds" else "_newClusters_newUMAPv1_5k.rds"))
  names <- switch(reference,
    aldinger=c("AldingerUMAP_newClusters.tiff", "AldingerUMAP_newClusters_newUMAPv1.tiff",
               "AldingerUMAP_newClusters_5k_newUMAP50v2.tiff", "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.tiff",
               "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"),
    sepp=c("SeppUMAP_FC_newClusters.tiff", "SeppUMAP_FC_newClusters_newUMAP50v1.tiff",
           "SeppUMAP_FC_newClusters_newUMAP50v2_5k.tiff", "SeppDotPlot_FC_newClusters_5k_newUMAP50v2_markers.tiff",
           "SeppDotPlot_FC_newClusters_5k_newUMAP50v2_markers.pdf"),
    science=c("ScienceUMAP_newClusters.tiff", "ScienceUMAP_newClusters_newUMAP50v1.tiff",
              "ScienceUMAP_newClusters_5k_newUMAP50v2.tiff", "ScienceDotPlot_newClusters_5k_newUMAP50v2_markers.tiff",
              "ScienceDotPlot_newclusters_5k_newUMAP50v2_markers.pdf"))
  entries <- list(
    list(input=original, reduction=if(reference=="sepp") "umap2d" else "umap", output=file.path(plots,names[[1]]), required=FALSE),
    list(input=full, reduction="umap", output=file.path(plots,names[[2]]), required=TRUE),
    list(input=prepanel, reduction="umap", output=file.path(plots,paste0(prefix,"UMAP_FC_newClusters_5k.tiff")), required=FALSE),
    list(input=panel, reduction="umap", output=file.path(plots,names[[3]]), required=TRUE))
  # Aldinger panel subset did not export a pre-recomputed UMAP.
  if (reference=="aldinger") entries <- entries[-3L]
  entries <- entries[vapply(entries, function(e) e$required || file.exists(e$output), logical(1))]
  aliases <- if(reference=="aldinger") file.path(plots,c("AldingerUMAP_newClusters_5k_newUMAP50v2.tif",
       "AldingerDotPlot_newclusters_5k_newUMAP50v2_markers.tif")) else character()
  aliases <- aliases[file.exists(aliases)]
  list(reference=reference, panel=panel, entries=entries, dot=file.path(plots,names[4:5]),
       aliases=aliases, manifest=file.path(plots,"Reference_Plot_Refresh_manifest.csv"))
}
specs <- lapply(selected,reference_specs)
paths <- function(spec) c(vapply(spec$entries, `[[`, "", "output"), spec$dot, spec$aliases, spec$manifest)
inputs <- unique(unlist(lapply(specs,function(spec) vapply(spec$entries, `[[`, "", "input"))))
outputs <- unlist(lapply(specs,paths),use.names=FALSE)
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("Reference plots: original paths and UMAP style", inputs, outputs)
  cat("  UMAP: 7 x 6 inches; original Seurat point-size/raster defaults; saved coordinates.\n")
  cat("  Existing plot replacement requires --overwrite; no RDS writes.\n")
  quit(save="no",status=if(ok)0L else 1L)
}
if (any(!file.exists(inputs))) stop("Missing saved reference inputs: ",paste(inputs[!file.exists(inputs)],collapse="; "))
if (any(file.exists(outputs)) && !"--overwrite" %in% args) stop("Existing reference plots are protected; use --overwrite for the requested replacement.")
dependencies <- c("Seurat","SeuratObject","ggplot2","patchwork","scales","Cairo")
missing <- dependencies[!vapply(dependencies,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing project dependencies: ",paste(missing,collapse=", "),". Resolve explicitly through renv; this script installs nothing.")
source(here("scripts","color_palette.R"))
set.seed(config$runtime$random_seed)
grDevices::pdf(file=NULL)
save_tiff <- function(plot,path,width,height) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  Cairo::CairoTIFF(path,width=width,height=height,units="in",res=600)
  tryCatch(print(plot),finally=grDevices::dev.off())
}
for (spec in specs) {
  provenance <- list()
  for (entry in spec$entries) {
    input_info <- file.info(entry$input)
    object <- readRDS(entry$input)
    cells <- colnames(object); metadata <- object[[]]
    if(!length(cells)||anyDuplicated(cells)||!identical(rownames(metadata),cells)||!"clusters_refined" %in% names(metadata)) stop("Invalid reference metadata: ",entry$input)
    values <- as.character(metadata$clusters_refined)
    if(anyNA(values)||any(!nzchar(values))) stop("Missing reference identities: ",entry$input)
    validate_palette(values)
    if(any(!values %in% celltype_order)) stop("Reference labels are outside the configured broad order: ",entry$input)
    if(!entry$reduction %in% Seurat::Reductions(object)) stop("Saved reduction missing: ",entry$reduction," in ",entry$input)
    levels <- intersect(celltype_order,unique(values))
    object$clusters_refined <- factor(values,levels=levels)
    Seurat::Idents(object) <- "clusters_refined"
    # Match the reference source scripts: no pt.size/raster/raster.dpi override.
    umap <- Seurat::DimPlot(object,reduction=entry$reduction,group.by="clusters_refined") +
      ggplot2::scale_color_manual(values=cluster_colors)
    panel_view <- identical(entry$input,spec$panel)
    primary_umap <- if(panel_view && spec$reference=="aldinger") umap+ggplot2::ggtitle("UMAP of Refined Clusters") else umap
    save_tiff(primary_umap,entry$output,7,6)
    produced <- entry$output
    if(panel_view) {
      umap_alias <- spec$aliases[grepl("UMAP",basename(spec$aliases),fixed=TRUE)]
      for(alias in umap_alias) save_tiff(umap,alias,7,6)
      produced <- c(produced,umap_alias)
      assay <- Seurat::DefaultAssay(object)
      layers <- refresh_normalized_layers(object[[assay]])
      if(!length(layers)) stop("Stored normalized reference data missing: ",entry$input)
      if(length(layers)>1L) object <- SeuratObject::JoinLayers(object,assay=assay,layers="^data($|[.])",new="data")
      Seurat::Idents(object) <- factor(values,levels=rev(levels))
      panels <- refresh_marker_panels("broad",levels)
      features <- lapply(panels,intersect,y=rownames(object[[assay]]))
      features <- features[lengths(features)>0L]
      if(!length(features)) stop("No configured broad markers in the reference assay.")
      dot <- workflow_dotplot(object,features=features,assay=assay,cluster.idents=FALSE) +
        ggplot2::ggtitle("High-Specificity Marker Expression by Cluster")
      dot <- standardize_broad_dotplot(dot,levels)
      size <- workflow_dotplot_dimensions(dot)
      display <- patchwork::wrap_elements(full=workflow_dotplot_grob(dot))
      save_tiff(display,spec$dot[[1]],size[["width"]],size[["height"]])
      ggplot2::ggsave(spec$dot[[2]],display,device=grDevices::cairo_pdf,width=size[["width"]],height=size[["height"]])
      dot_alias <- spec$aliases[grepl("DotPlot",basename(spec$aliases),fixed=TRUE)]
      for(alias in dot_alias) save_tiff(display,alias,size[["width"]],size[["height"]])
      produced <- c(produced,spec$dot,dot_alias)
    }
    if(!identical(file.info(entry$input)[c("size","mtime")],input_info[c("size","mtime")])) stop("Reference source changed while plotting.")
    if(!identical(colnames(object),cells)||!identical(as.character(object$clusters_refined),values)) stop("Reference cells or labels changed.")
    other_columns <- setdiff(names(metadata),"clusters_refined")
    if(!identical(object[[]][other_columns],metadata[other_columns])) stop("Reference metadata changed during plotting.")
    if(anyNA(file.info(produced)$size)||any(file.info(produced)$size<=0)) stop("Incomplete reference plot exports.")
    provenance[[length(provenance)+1L]] <- data.frame(input_rds=normalizePath(entry$input),
      input_size=input_info$size,input_mtime=as.numeric(input_info$mtime),
      reduction=entry$reduction,output_path=normalizePath(produced),
      palette_md5=unname(tools::md5sum(here("scripts","color_palette.R"))),
      umap_width=7,umap_height=6,point_size="Seurat default",raster="Seurat default",stringsAsFactors=FALSE)
    rm(object);gc()
  }
  utils::write.csv(do.call(rbind,provenance),spec$manifest,row.names=FALSE)
  message("Updated reference plots in original directory: ",dirname(spec$manifest))
}
grDevices::dev.off()
