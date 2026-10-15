#!/bin/bash
# Targeted new-sample preprocessing, resolution-5 annotation, and QC audit.
# --dry-run checks raw inputs/task maps only; --submit queues the dependency chain.
set -euo pipefail
if [[ $# -ne 1 || ( $1 != --dry-run && $1 != --submit ) ]]; then
  echo "Usage: bash scripts/submit_xenium_new_samples_20261014.sh [--dry-run|--submit]" >&2
  exit 2
fi
mode=$1
cd /home/acflint/R/Projects/XeniumFCProject
cpu_vendor=$(lscpu | awk -F: '/Vendor ID/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}')
[[ $cpu_vendor == GenuineIntel ]] || {
  echo "Run preflight on an Intel interactive node; the project library requires Intel." >&2
  exit 2
}
module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
cohort=scripts/config/new_samples_20261014.txt
[[ -f $cohort ]] || cohort=config/new_samples_20261014.txt
[[ -f $cohort ]] || { echo "Missing cohort file" >&2; exit 2; }
# Keep this helper protective even if overwrite flags were set in the shell.
export SINGLE_PREPROCESS_OVERWRITE=false SPLIT_PREPROCESS_OVERWRITE=false
export RES5_ALL_OVERWRITE=false RES5_ALL_ABT_OVERWRITE=false
export RES5_SELECTED_CONSENSUS_OVERWRITE=false RES5_ALL_W2OF3_OVERWRITE=false
export UNKNOWN_QC_AUDIT_OVERWRITE=false BROAD_SCORE_OVERWRITE=false
Rscript --max-ppsize=500000 scripts/validate_config.R --check-files
single_map=$(Rscript --max-ppsize=500000 scripts/xenium_preprocess_single_slides.R --list)
split_map=$(Rscript --max-ppsize=500000 scripts/xenium_preprocess_split_slides.R --list)
global_map=$(Rscript --max-ppsize=500000 scripts/xenium_preprocess_03b_resolution2_pilot.R --all-samples-res5 --list)
broad_map=$(Rscript --max-ppsize=500000 scripts/xenium_broad_02_score_clusters.R --list)
select_ids() {
  awk 'NR==FNR {if (NF) wanted[$1]=1; next} $2 in wanted {print $1}' "$cohort" -
}
single_ids=$(printf '%s\n' "$single_map" | select_ids)
split_ids=$(printf '%s\n' "$split_map" | select_ids)
global_ids=$(printf '%s\n' "$global_map" | select_ids)
broad_ids=$(printf '%s\n' "$broad_map" | select_ids)
[[ $(printf '%s\n' "$global_ids" | awk 'NF {n++} END {print n}') == 7 ]] || {
  echo "Expected exactly seven matching production tasks; inspect --list." >&2; exit 2;
}
[[ $(printf '%s\n' "$single_ids" | awk 'NF {n++} END {print n}') == 6 &&
   $(printf '%s\n' "$split_ids" | awk 'NF {n++} END {print n}') == 1 &&
   $(printf '%s\n' "$broad_ids" | awk 'NF {n++} END {print n}') == 21 ]] || {
  echo "Expected six single, one split, and 21 broad-score tasks." >&2; exit 2;
}
for task in $single_ids; do
  Rscript --max-ppsize=500000 scripts/xenium_preprocess_single_slides.R --dry-run "$task"
done
for task in $split_ids; do
  Rscript --max-ppsize=500000 scripts/xenium_preprocess_split_slides.R --dry-run "$task"
done
# Check crop cell lists and PCW joins without loading Xenium or RDS objects.
export NEW_SAMPLE_COHORT="$cohort"
Rscript --max-ppsize=500000 -e 'source("scripts/R/config.R"); cfg <- load_pipeline_config(); ids <- readLines(Sys.getenv("NEW_SAMPLE_COHORT")); manifest <- load_sample_manifest(cfg); meta <- readxl::read_excel(resolve_config_path(cfg$manifests$sample_metadata,cfg)); stopifnot(all(c("sample","PCW") %in% names(meta))); for(id in ids) {p <- resolve_sample_paths(get_sample_record(id,manifest,cfg),cfg); x <- read.csv(p$cell_stats_path,comment.char="#",check.names=FALSE); stopifnot("Cell ID" %in% names(x),nrow(x)>0,!anyNA(x[["Cell ID"]]),!anyDuplicated(x[["Cell ID"]]),all(nzchar(x[["Cell ID"]]))); area <- names(x)[grepl("area",names(x),ignore.case=TRUE)]; stopifnot(length(area)>0); key <- sub("_[0-9]+$","",id); hit <- which(as.character(meta$sample)==key); if(length(hit)!=1L || is.na(meta$PCW[hit]) || !nzchar(as.character(meta$PCW[hit]))) stop("Missing/duplicate PCW for ",key); cat(id," crop cells=",nrow(x)," PCW=",as.character(meta$PCW[hit]),"\n",sep=""); }'
printf 'Single tasks: %s\nSplit tasks: %s\nProduction tasks: %s\n' "$single_ids" "$split_ids" "$global_ids"
[[ $mode == --submit ]] || exit 0
mkdir -p logs
to_array() { printf '%s\n' "$1" | paste -sd, -; }
submit_job() {
  local result
  result=$(sbatch --parsable "$@")
  result=${result%%;*}
  [[ $result =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $result" >&2; exit 2; }
  printf '%s\n' "$result"
}
single_job=$(submit_job --array="$(to_array "$single_ids")%3" scripts/run_xenium_preprocess_single_slides.slurm)
echo "Single preprocessing: $single_job"
split_job=$(submit_job --array="$(to_array "$split_ids")%3" --dependency="afterok:$single_job" scripts/run_xenium_preprocess_split_slides.slurm)
echo "Split preprocessing: $split_job"
res5_job=$(submit_job --array="$(to_array "$global_ids")%4" --dependency="afterok:$split_job" scripts/run_xenium_preprocess_04_resolution5_clustered.slurm)
echo "Resolution 5: $res5_job"
transfers=()
for reference in Aldinger Sepp Science; do
  job=$(submit_job --array="$(to_array "$global_ids")%3" --dependency="afterok:$res5_job" --export="ALL,REFERENCE=$reference" scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
  transfers+=("$job")
  echo "$reference transfer: $job"
done
transfer_dependency=$(IFS=:; echo "${transfers[*]}")
# Export via the environment: sbatch --export comma syntax would truncate IDs.
export SELECTED_SAMPLE_IDS
SELECTED_SAMPLE_IDS=$(paste -sd, "$cohort")
build_job=$(submit_job --dependency="afterok:$transfer_dependency" scripts/run_xenium_annotate_02_build_consensus_resolution5_selected_samples_weighted2of3.slurm)
echo "Selected consensus tables: $build_job"
apply_job=$(submit_job --array="$(to_array "$global_ids")%6" --dependency="afterok:$build_job" scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
echo "Consensus objects/plots: $apply_job"
audit_job=$(submit_job --array="$(to_array "$global_ids")%3" --dependency="afterok:$apply_job" scripts/run_xenium_unknown_01_audit_qc.slurm)
echo "Unknown QC audit: $audit_job"
broad_job=$(submit_job --array="$(to_array "$broad_ids")%6" --dependency="afterok:$apply_job" scripts/run_xenium_broad_02_score_clusters.slurm)
echo "Broad reference scores: $broad_job"
echo "Stop for review before replacing shared QC/broad-routing tables or making regional subsets."
