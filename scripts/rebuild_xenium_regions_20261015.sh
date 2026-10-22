#!/bin/bash
# Rebuild the 41-sample VZ/RL branches through fresh QC review, never QC apply.
# --dry-run is read-only. --submit archives stages 02-09 before queuing jobs.
set -euo pipefail

mode=${1:-}
case "$mode" in
  --dry-run|--submit) [[ $# == 1 ]] || exit 2 ;;
  --stage) [[ $# == 4 ]] || exit 2 ;;
  *) echo "Usage: bash scripts/rebuild_xenium_regions_20261015.sh [--dry-run|--submit]" >&2; exit 2 ;;
esac
cd /home/acflint/R/Projects/XeniumFCProject
cpu_vendor=$(lscpu | awk -F: '/Vendor ID/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}')
[[ $cpu_vendor == GenuineIntel ]] || { echo "Use an Intel session for the project R library." >&2; exit 2; }
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

run_dry() {
  local result
  result=$(Rscript --max-ppsize=500000 "$@" 2>&1) || { printf '%s\n' "$result" >&2; return 2; }
  printf '%s\n' "$result"
  [[ $result == *"DRY-RUN PASS"* && $result != *"DRY-RUN FAIL"* ]]
}

# Pin execution modes: no inherited plots-only, old checkpoint, QC-removal,
# overwrite, or worker override may silently alter this fresh cohort rebuild.
export VZ_MERGE_OVERWRITE=false RL_MERGE_OVERWRITE=false
export VZ_INTEGRATE_OVERWRITE=false RL_INTEGRATE_OVERWRITE=false
export VZ_INTEGRATE_PLOTS_ONLY=false RL_INTEGRATE_PLOTS_ONLY=false
export VZ_INTEGRATE_INCLUDE_UNCORRECTED=false RL_INTEGRATE_INCLUDE_UNCORRECTED=false
export VZ_INTEGRATE_RESUME_CHECKPOINT=false RL_INTEGRATE_RESUME_CHECKPOINT=false
export VZ_INTEGRATE_CLUSTER_WORKERS=1 RL_INTEGRATE_CLUSTER_WORKERS=1
export VZ_QC_MODE=review RL_QC_MODE=review VZ_QC_OVERWRITE=false RL_QC_OVERWRITE=false
export VZ_HARMONY_AUDIT_OVERWRITE=false RL_HARMONY_AUDIT_OVERWRITE=false

stage_paths() {
  local branch=$1 stage=$2
  case "$branch" in vz|rl) ;; *) return 2 ;; esac
  case "$stage" in
    merge) driver="scripts/xenium_${branch}_02_merge_samples.R"; launcher="scripts/run_xenium_${branch}_02_merge_samples.slurm" ;;
    integrate) driver="scripts/xenium_${branch}_03_integrate.R"; launcher="scripts/run_xenium_${branch}_03_integrate.slurm" ;;
    qc) driver="scripts/xenium_${branch}_04_review_qc.R"; launcher="scripts/run_xenium_${branch}_04_review_qc.slurm" ;;
    audit) driver="scripts/xenium_${branch}_03a_audit_harmony.R"; launcher="scripts/run_xenium_${branch}_03a_audit_harmony.slurm" ;;
    *) return 2 ;;
  esac
}

if [[ $mode == --stage ]]; then
  [[ ${SLURM_JOB_ID:-} =~ ^[0-9]+$ ]] || { echo "Internal --stage mode must run in Slurm." >&2; exit 2; }
  stage_paths "$2" "$3"
  [[ -f $4 ]] || { echo "Missing submission record: $4" >&2; exit 2; }
  printf 'Run record: %s | branch: %s | stage: %s\n' "$4" "$2" "$3"
  if [[ $3 == qc ]]; then
    run_dry "$driver" --qc-only --dry-run
  else
    run_dry "$driver" --dry-run
  fi
  exec bash "$launcher"
fi

# Stop if any regional or combined-region job could still be using old files.
regional_names=$(awk '/^#SBATCH --job-name=/ {sub(/^#SBATCH --job-name=/, ""); print}' scripts/run_xenium_vz_*.slurm scripts/run_xenium_rl_*.slurm)
active_names=$(squeue -h -u acflint -o '%j')
while IFS= read -r name; do
  if [[ -n $name ]] && printf '%s\n' "$regional_names" | grep -Fxq "$name"; then
    echo "Regional job still active: $name. Wait before archiving or rebuilding." >&2
    exit 2
  fi
done <<< "$active_names"

preflight=$(Rscript --max-ppsize=500000 - <<'RSCRIPT'
source("scripts/R/config.R")
cfg <- load_pipeline_config()
samples <- load_sample_manifest(cfg)
ids <- samples$sample_id
if(length(ids)!=41L || anyDuplicated(ids)) stop("This dated rebuild expects 41 unique manifest samples.")
packages <- c("Seurat","harmony","dplyr","future","ggplot2","patchwork","scales","RANN")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing project dependencies (no packages installed): ",paste(missing,collapse=", "))
root <- normalizePath(here::here(cfg$project$outputs_dir),mustWork=TRUE)
for(branch in c("vz","rl")) {
  paths <- file.path(root,"xenium",branch,"01_subsets","rds",paste0(ids,"_",toupper(branch),"subset.rds"))
  if(any(!file.exists(paths)) || any(file.info(paths)$size<=0,na.rm=TRUE)) stop("Missing/empty ",branch," subset inputs.")
}
inputs <- file.path(get_production_output_paths(cfg)$consensus_rds_dir,paste0(ids,"_Consensus_annotated.rds"))
if(any(!file.exists(inputs)) || any(file.info(inputs)$size<=0,na.rm=TRUE)) stop("Missing/empty consensus inputs for later QC review.")
cat("PASS: 41 manifest samples, 82 nonempty regional subsets and 41 consensus inputs.\n")
cat("OUTPUT_ROOT\t",root,"\n",sep="")
RSCRIPT
)
printf '%s\n' "$preflight"
output_root=$(printf '%s\n' "$preflight" | awk -F '\t' '$1=="OUTPUT_ROOT" {print $2}')
[[ -d $output_root && $output_root != / ]] || { echo "Invalid configured output root." >&2; exit 2; }
run_record="$output_root/regional_rebuild_20261015.run"
[[ ! -e $run_record ]] || { echo "A rebuild was already started. Inspect $run_record and recorded jobs; do not submit twice." >&2; exit 2; }
Rscript --max-ppsize=500000 scripts/validate_config.R --check-files
run_dry scripts/xenium_vz_02_merge_samples.R --dry-run
run_dry scripts/xenium_rl_02_merge_samples.R --dry-run

for branch in vz rl; do
  for stage in merge integrate qc audit; do
    stage_paths "$branch" "$stage"
    bash -n "$launcher"
    grep -q 'Rscript --max-ppsize=500000' "$launcher"
    printf 'Launcher: %s\n' "$launcher"
  done
done
echo "Planned archive: existing 02_* through 09_* directories in both branches."
echo "Retained: all 01_subsets inputs, whole-tissue inputs, scientific settings and labels."
echo "Plan per branch: merge -> integrate -> QC review + Harmony audit. Stop before QC apply."
[[ $mode == --submit ]] || exit 0

# Renames stay on the configured output filesystem, preserving large RDS bytes
# and avoiding a second full copy. The persistent run record prevents retrying
# a partial submission and accidentally archiving new results.
lock="$output_root/.regional_rebuild_20261015.lock"
mkdir "$lock" || { echo "Another regional submission holds $lock" >&2; exit 2; }
trap 'rmdir "$lock"' EXIT
[[ ! -e $run_record ]] || exit 2
mkdir -p "$output_root/regional_rebuild_backups" logs
backup=$(mktemp -d "$output_root/regional_rebuild_backups/cohort41_20261015_XXXXXX")
printf '%s\n' "$backup" > "$run_record"
mkdir -p "$backup/provenance"
cp -p config/samples.csv config/config.yml scripts/color_palette.R "$backup/provenance/"
cp -p "$0" "$backup/provenance/"
printf 'old_path\tarchived_path\n' > "$backup/archive_manifest.tsv"
shopt -s nullglob
for branch in vz rl; do
  mkdir -p "$backup/$branch"
  for path in "$output_root/xenium/$branch"/0[2-9]_*; do
    [[ -d $path && ! -L $path ]] || { echo "Unexpected regional stage path: $path" >&2; exit 2; }
    destination="$backup/$branch/${path##*/}"
    mv "$path" "$destination"
    printf '%s\t%s\n' "$path" "$destination" >> "$backup/archive_manifest.tsv"
  done
done
echo "Previous regional results/checkpoints archived in: $backup"
echo "Combined regional outputs are now stale until both new branches are reviewed and mapped."
printf 'branch\tstage\tjob_id\n' > "$backup/job_ids.tsv"

submit_stage() {
  local branch=$1 stage=$2 dependency=${3:-} directive id wrapper
  local -a options=()
  stage_paths "$branch" "$stage"
  # Read resource/log/mail directives from the existing launcher verbatim.
  # Each current directive uses one whitespace-free --option=value token.
  while IFS= read -r directive; do
    directive=${directive#\#SBATCH }
    [[ $directive == --* && $directive != *[[:space:]]* ]] || { echo "Unsupported SBATCH directive: $directive" >&2; return 2; }
    options+=("$directive")
  done < <(awk '/^#SBATCH / {print}' "$launcher")
  [[ -z $dependency ]] || options+=("--dependency=afterok:$dependency")
  printf -v wrapper '%q ' bash scripts/rebuild_xenium_regions_20261015.sh --stage "$branch" "$stage" "$run_record"
  id=$(sbatch --parsable "${options[@]}" --export=ALL --wrap="$wrapper")
  id=${id%%;*}
  [[ $id =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $id" >&2; return 2; }
  printf '%s\t%s\t%s\n' "$branch" "$stage" "$id" >> "$backup/job_ids.tsv"
  echo "$branch $stage: $id" >&2
  printf '%s\n' "$id"
}
for branch in vz rl; do
  merge_id=$(submit_stage "$branch" merge)
  integrate_id=$(submit_stage "$branch" integrate "$merge_id")
  qc_id=$(submit_stage "$branch" qc "$integrate_id")
  audit_id=$(submit_stage "$branch" audit "$integrate_id")
done
echo "Submitted both branches. Job IDs: $backup/job_ids.tsv"
echo "Next: inspect new stage-4 QC reports and Harmony audits; no old removal IDs or annotations were applied."
