#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_original_plot_plan.R")
source(here("scripts","R","config.R"))
source(here("scripts","R","palette_refresh.R"))
source(here("scripts","R","xenium_original_plots.R"))
args <- commandArgs(trailingOnly=TRUE)
stage_args <- args[startsWith(args,"--stage=")]
if(length(stage_args)!=1L||any(!args %in% c("--list","--dry-run") & !startsWith(args,"--stage=") & !grepl("^[0-9]+$",args))) stop("Use --stage=NAME [--list|--dry-run] [TASK_ID].")
if(all(c("--list","--dry-run") %in% args)) stop("Choose --list or --dry-run.")
stage <- sub("^--stage=","",stage_args)
plan <- original_xenium_plan(load_pipeline_config(),stage)
if("--list" %in% args) {
  write.table(plan$tasks,row.names=FALSE,quote=FALSE,sep="\t")
  quit(save="no",status=0L)
}
task <- args[grepl("^[0-9]+$",args)]
if(length(task)>1L) stop("At most one TASK_ID.")
indices <- if(length(task)) match(as.integer(task),plan$tasks$task_id) else seq_len(nrow(plan$tasks))
if(anyNA(indices)) stop("TASK_ID is outside the manifest mapping.")
input_paths <- unique(c(unlist(plan$inputs[indices],use.names=FALSE),here("scripts",plan$script),
  here("scripts",c("color_palette.R","R/dotplot_helpers.R","R/palette_refresh.R","R/xenium_original_plots.R"))))
ok <- compact_dry_run(paste("Original Xenium plots",stage),inputs=input_paths)
cat("  Tasks: ",length(indices),"; original stage output paths; plot replacement authorized by --submit; RDS files remain inputs only.\n",sep="")
quit(save="no",status=if(ok)0L else 1L)
