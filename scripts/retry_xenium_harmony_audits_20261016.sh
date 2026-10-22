#!/bin/bash
# Retry only the two diagnostics; never rerun regional merges or integrations.
set -euo pipefail
[[ $# == 1 && ( $1 == --dry-run || $1 == --submit ) ]] || {
  echo "Usage: bash scripts/retry_xenium_harmony_audits_20261016.sh [--dry-run|--submit]" >&2
  exit 2
}
mode=$1
cd /home/acflint/R/Projects/XeniumFCProject
vendor=$(lscpu | awk -F: '/Vendor ID/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}')
[[ $vendor == GenuineIntel ]] || { echo "Use an Intel terminal." >&2; exit 2; }
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

active=$(squeue -h -u acflint -o '%j')
if printf '%s\n' "$active" | grep -Eq '^Xen_(VZ|RL)_Harmony_audit$'; then
  echo "A Harmony audit is already active. Inspect it before retrying." >&2
  exit 2
fi
paths=$(Rscript --max-ppsize=500000 -e '
source("scripts/R/config.R")
cfg <- load_pipeline_config()
cat("OUTPUT_ROOT\t",normalizePath(here::here(cfg$project$outputs_dir),mustWork=TRUE),"\n",sep="")
')
root=$(printf '%s\n' "$paths" | awk -F '\t' '$1=="OUTPUT_ROOT" {print $2}')
[[ -d $root && $root != / ]] || exit 2
record="$root/regional_rebuild_20261015.run"
[[ -f $record ]] || { echo "Missing original regional rebuild run record: $record" >&2; exit 2; }
run_dir=$(cat "$record")
[[ -f $run_dir/job_ids.tsv ]] || { echo "Missing original job IDs." >&2; exit 2; }
ids=$(awk -F '\t' 'NR>1 {print $3}' "$run_dir/job_ids.tsv" | paste -sd, -)
echo "Original regional job outcomes (QC may need a separate investigation):"
sacct -X -j "$ids" --format=JobID,JobName%28,Elapsed,State,ExitCode
for branch in vz rl; do
  id=$(awk -F '\t' -v b="$branch" '$1==b && $2=="integrate" {print $3}' "$run_dir/job_ids.tsv")
  [[ $id =~ ^[0-9]+$ ]] || { echo "Missing unambiguous $branch integration job ID." >&2; exit 2; }
  state=$(sacct -X -n -P -j "$id" --format=JobID,State,ExitCode)
  if ! printf '%s\n' "$state" | awk -F '|' -v id="$id" '$1==id && $2=="COMPLETED" && $3=="0:0" {ok=1} END {exit !ok}'; then
    echo "$branch integration $id did not complete successfully. Diagnose that failure before an audit retry." >&2
    exit 2
  fi
  launcher="scripts/run_xenium_${branch}_03a_audit_harmony.slurm"
  bash -n "$launcher"
  grep -Fxq '#SBATCH --mem=256G' "$launcher" || { echo "Upload the patched 256-GB audit launcher first." >&2; exit 2; }
  result=$(Rscript --max-ppsize=500000 "scripts/xenium_${branch}_03a_audit_harmony.R" --dry-run 2>&1) || { printf '%s\n' "$result" >&2; exit 2; }
  printf '%s\n' "$result"
  [[ $result == *"DRY-RUN PASS"* && $result != *"DRY-RUN FAIL"* ]] || exit 2
done
[[ $mode == --submit ]] || exit 0

mkdir -p "$root/report_backups" logs
backup=$(mktemp -d "$root/report_backups/harmony_audit_retry_20261016_XXXXXX")
printf 'branch\tjob_id\n' > "$backup/job_ids.tsv"
for branch in vz rl; do
  dir="$root/xenium/$branch/03_integrated/tables/harmony_audit"
  if [[ -e $dir ]]; then
    [[ -d $dir && ! -L $dir ]] || { echo "Unexpected audit output path: $dir" >&2; exit 2; }
    mv "$dir" "$backup/${branch}_prior_audit"
  fi
done
echo "Previous/partial audit tables preserved in: $backup"
for branch in vz rl; do
  if [[ $branch == vz ]]; then export VZ_HARMONY_AUDIT_OVERWRITE=false; else export RL_HARMONY_AUDIT_OVERWRITE=false; fi
  id=$(sbatch --parsable --export=ALL "scripts/run_xenium_${branch}_03a_audit_harmony.slurm")
  id=${id%%;*}
  [[ $id =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $id. Inspect recorded jobs before retrying." >&2; exit 2; }
  printf '%s\t%s\n' "$branch" "$id" >> "$backup/job_ids.tsv"
  echo "$branch Harmony audit: $id"
done
echo "Submitted only the audits. Integrations, clusters, QC decisions and annotations were not rerun."
