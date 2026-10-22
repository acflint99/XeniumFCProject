# Paths and task mappings only: no Seurat objects are loaded here.
original_manifest_ids <- function(config) {
  ids <- load_sample_manifest(config)$sample_id
  if(!length(ids)||anyNA(ids)||any(!nzchar(ids))||anyDuplicated(ids)) stop("Sample manifest requires unique, nonempty sample IDs.")
  ids
}

original_plot_directory <- function(root, sample_ids, suffixes) {
  expected <- as.vector(outer(sample_ids, suffixes, paste0))
  found <- if (dir.exists(root)) list.files(root, recursive=TRUE, full.names=TRUE) else character()
  found <- found[basename(found) %in% expected]
  directories <- unique(dirname(found))
  if(length(directories)>1L) stop("Plot family occurs in multiple original directories; resolve duplicate/mixed layouts: ",paste(directories,collapse="; "))
  if(length(directories)) directories[[1]] else root
}

original_combined_plot_spec <- function(config) {
  ids <- original_manifest_ids(config)
  root <- here::here(config$project$outputs_dir,"xenium","vz_rl")
  corrected <- refresh_corrected_broad_paths(here::here(config$project$outputs_dir),ids)
  combined_root <- file.path(root,"01_combined_labels","plots")
  fine_suffix <- c("_combSubcluster_GlobalSpatial_pt0.6.tif","_combSubcluster_GlobalSpatial_pt0.85.tif")
  broad_suffix <- c("_postQC_GlobalSpatial_pt0.85.tif","_postQC_GlobalSpatial_pt0.6.tif")
  corrected_suffix <- c("_CorrectedBroad_GlobalSpatial_pt0.6.tif","_CorrectedBroad_GlobalSpatial_pt0.85.tif")
  fine_dir <- original_plot_directory(combined_root,ids,fine_suffix)
  broad_dir <- original_plot_directory(combined_root,ids,broad_suffix)
  corrected_root <- file.path(dirname(dirname(corrected[[1]])),"plots")
  corrected_dir <- original_plot_directory(corrected_root,ids,corrected_suffix)
  list(ids=ids, combined=file.path(root,"01_combined_labels","rds",paste0(ids,"_VZ_RL_QC_Subclusters.rds")),
       corrected=unname(corrected),outputs=lapply(ids,function(id) c(
         file.path(fine_dir,paste0(id,fine_suffix)),file.path(broad_dir,paste0(id,broad_suffix)),
         file.path(corrected_dir,paste0(id,corrected_suffix)))))
}

original_xenium_plan <- function(config,stage) {
  ids <- original_manifest_ids(config)
  root <- here::here(config$project$outputs_dir)
  singleton <- function(script,input,args=character()) list(
    tasks=data.frame(task_id=1L,sample_id=stage),inputs=list(input),script=script,args=args,per_task=FALSE)
  sample_plan <- function(script,inputs,args=character()) list(
    tasks=data.frame(task_id=seq_along(ids),sample_id=ids),inputs=as.list(inputs),script=script,args=args,per_task=TRUE)
  result <- switch(stage,
    consensus=sample_plan("xenium_annotate_03_apply_consensus.R",
      file.path(get_production_output_paths(config)$consensus_rds_dir,paste0(ids,"_Consensus_annotated.rds")),
      c("--all-samples-res5","--weighted-2of3","--plots-only")),
    vz=sample_plan("xenium_vz_08_plot_samples.R",file.path(root,"xenium","vz","07_mapped","rds",paste0(ids,"_VZ_QC_Subclusters.rds"))),
    rl=sample_plan("xenium_rl_08_plot_samples.R",file.path(root,"xenium","rl","07_mapped","rds",paste0(ids,"_RL_QC_Subclusters.rds"))),
    residual=sample_plan("xenium_residual_07_plot_samples.R",file.path(root,"xenium","residual","06_mapped","rds",paste0(ids,"_Residual_Subclusters.rds"))),
    combined={
      spec <- original_combined_plot_spec(config)
      sample_plan("refresh_combined_sample_plots_original.R",lapply(seq_along(ids),function(i)c(spec$combined[[i]],spec$corrected[[i]])))
    },
    `vz-pooled`=singleton("xenium_vz_06_annotate_subclusters.R",file.path(root,"xenium","vz","06_subclusters","rds","Xenium_VZ_subclusters.rds"),c("--plots-only","--resolution=0.8")),
    `rl-pooled`=singleton("xenium_rl_06_annotate_subclusters.R",file.path(root,"xenium","rl","06_subclusters","rds","Xenium_RL_subclusters.rds"),c("--plots-only","--resolution=0.8")),
    `residual-pooled`=singleton("xenium_residual_05_plot_subclusters.R",file.path(root,"xenium","residual","05_subclusters","rds","Xenium_Residual_Subclusters.rds")),
    `counts-vz`=singleton("xenium_vz_09_plot_counts.R",file.path(root,"xenium","vz","06_subclusters","rds","Xenium_VZ_subclusters.rds")),
    `counts-rl`=singleton("xenium_rl_09_plot_counts.R",file.path(root,"xenium","rl","06_subclusters","rds","Xenium_RL_subclusters.rds")),
    nonspatial=singleton("xenium_vz_rl_04_05_nonspatial_process_and_plot.R",c(
      file.path(root,"xenium","vz_rl","nonspatial","04_processed","rds","Xenium_VZRL_processed.rds"),
      file.path(root,"xenium","vz_rl","nonspatial","04_processed","tables","Xenium_VZRL_processed_manifest.csv"),
      file.path(root,"xenium","vz_rl","nonspatial","03_merged","rds","Xenium_VZRL_merged.rds")),"--plots-only"),
    spatial=singleton("refresh_spatial_umaps_original.R",file.path(root,"xenium","vz_rl","spatial","04_integrated","rds","Xenium_VZRL_spatial_integrated.rds")),
    stop("Unknown native plot stage: ",stage))
  if(stage %in% c("residual","residual-pooled")) {
    extra <- resolve_config_path(config$manifests$reviewed_residual_subcluster_labels,config)
    if(stage=="residual") extra <- c(extra,file.path(root,"xenium","residual","06_mapped","tables","Xenium_Residual_Mapping_Manifest.csv"))
    result$inputs <- lapply(result$inputs,function(paths)c(paths,extra))
  }
  result
}
