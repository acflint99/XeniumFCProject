#!/usr/bin/env Rscript
# Redraw the six original combined/corrected spatial TIFFs from saved objects.
suppressPackageStartupMessages(library(here))
here::i_am("scripts/refresh_combined_sample_plots_original.R")
source(here("scripts","R","config.R"))
source(here("scripts","R","palette_refresh.R"))
source(here("scripts","R","xenium_original_plots.R"))
config <- load_pipeline_config()
spec <- original_combined_plot_spec(config)
args <- commandArgs(trailingOnly=TRUE)
if(any(!args %in% c("--list","--dry-run","--overwrite") & !grepl("^[0-9]+$",args))||sum(c("--list","--dry-run","--overwrite") %in% args)>1L) stop("Use [--list|--dry-run|--overwrite] [TASK_ID].")
if("--list" %in% args) {write.table(data.frame(task_id=seq_along(spec$ids),sample_id=spec$ids),row.names=FALSE,quote=FALSE,sep="\t");quit(save="no",status=0L)}
task_args <- args[grepl("^[0-9]+$",args)]
if(length(task_args)>1L) stop("One TASK_ID is required.")
task <- suppressWarnings(as.integer(if(length(task_args)) task_args else Sys.getenv("SLURM_ARRAY_TASK_ID","")))
if(is.na(task)||!task %in% seq_along(spec$ids)) stop("TASK_ID is outside --list.")
inputs <- c(spec$combined[[task]],spec$corrected[[task]])
outputs <- spec$outputs[[task]]
if("--dry-run" %in% args) {ok<-compact_dry_run(paste("Original combined plots",spec$ids[[task]]),inputs,outputs);quit(save="no",status=if(ok)0L else 1L)}
if(any(!file.exists(inputs))) stop("Missing saved combined/corrected object.")
if(any(file.exists(outputs))&&!"--overwrite" %in% args) stop("Existing plots require --overwrite.")
suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
if(!requireNamespace("Cairo",quietly=TRUE)) stop("Missing project Cairo dependency; use the established renv environment.")
source(here("scripts","color_palette.R"))
set.seed(config$runtime$random_seed)
for(index in seq_along(inputs)) {
  input_info <- file.info(inputs[[index]])
  object <- readRDS(inputs[[index]])
  cells <- colnames(object); metadata <- object[[]]
  columns <- if(index==1L)c("comb_subcluster","consensus_label") else "corrected_broad_label"
  if(!length(cells)||anyDuplicated(cells)||!identical(rownames(metadata),cells)||!all(columns %in% names(metadata))) stop("Invalid source metadata: ",inputs[[index]])
  for(column in columns) {
    labels <- as.character(metadata[[column]])
    if(anyNA(labels)||any(!nzchar(labels))) stop("Missing plotting identities.")
    preferred <- if(column=="comb_subcluster")unique(c(master_subcluster_order,celltype_order)) else celltype_order
    levels <- order_consensus_levels(labels,preferred)
    colors <- if(column=="comb_subcluster")subcluster_palette[levels] else resolve_cluster_colors(levels)
    if(column=="comb_subcluster") {unknown<-is_unknown_consensus_label(levels);colors[unknown]<-resolve_cluster_colors(levels[unknown])}
    if(anyNA(colors)) stop("Missing plotting colors.")
    object[[column]] <- factor(labels,levels=levels)
    indices <- switch(column,comb_subcluster=1:2,consensus_label=3:4,corrected_broad_label=5:6)
    sizes <- if(column=="consensus_label")c(.85,.6) else c(.6,.85)
    title <- paste(spec$ids[[task]],switch(column,comb_subcluster="Combined Subclusters",consensus_label="-Post QC (Consensus)",corrected_broad_label="Corrected Broad Labels"))
    for(j in seq_along(indices)) {
      plot <- ImageDimPlot(object,group.by=column,cols=colors,size=sizes[[j]])+ggtitle(title)
      path <- outputs[[indices[[j]]]]
      dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
      Cairo::CairoTIFF(path,width=if(column=="comb_subcluster")14 else 12,height=10,units="in",res=600,bg="black")
      tryCatch(print(plot),finally=grDevices::dev.off())
      message("Updated original plot: ",path)
    }
  }
  if(!identical(colnames(object),cells)||any(vapply(columns,function(column)!identical(as.character(object[[column]][,1]),as.character(metadata[[column]])),logical(1)))) stop("Source cell IDs or labels changed.")
  if(!identical(input_info[c("size","mtime")],file.info(inputs[[index]])[c("size","mtime")])) stop("Source RDS changed during redraw.")
  rm(object);gc()
}
if(anyNA(file.info(outputs)$size)||any(file.info(outputs)$size<=0)) stop("Partial original spatial plot set.")
message("PASS: six original TIFFs updated; no analysis RDS written.")
