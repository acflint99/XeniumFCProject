#!/usr/bin/env bash
set -euo pipefail

project_root=/home/acflint/R/Projects/XeniumFCProject
run_dir="$project_root/outputs/reference_label_refresh_runs/updated_reference_labels_20261022"
ledger="$run_dir/jobs.tsv"
cd "$project_root"

require_ledger_job() {
  local label=$1 expected_id=$2
  local line
  line=$(awk -F '\t' -v label="$label" '$1 == "annotation" && $2 == label {print}' "$ledger")
  if [[ $line != $'annotation\t'"$label"$'\t'"$expected_id"$'\t1' ]]; then
    printf 'Ledger differs from the reviewed Aldinger recovery state: %s\n' "$label" >&2
    exit 1
  fi
}

require_sacct() {
  local job_id=$1 state=$2 exit_code=$3
  if ! sacct -X -n -P -j "$job_id" --format=JobID,State,ExitCode |
      awk -F '|' -v job="$job_id" -v state="$state" -v code="$exit_code" \
        '$1 == job && $2 == state && $3 == code {found = 1} END {exit !found}'; then
    printf 'Unexpected accounting state for job %s. Stop and inspect.\n' "$job_id" >&2
    exit 1
  fi
}

require_ledger_job Aldinger_standardize 41015193
require_ledger_job Aldinger_panel 41015194
require_ledger_job Science_standardize 41049984
require_ledger_job Science_panel 41052256
require_ledger_job reference_validation 41015200

require_sacct 41015193 COMPLETED 0:0
require_sacct 41015194 COMPLETED 0:0
require_sacct 41049984 COMPLETED 0:0
require_sacct 41052256 COMPLETED 0:0
require_sacct 41015200 FAILED 1:0

active_conflicts=$(squeue -h -u acflint -o '%i|%j' |
  awk -F '|' '$2 ~ /^(Aldinger_|Xen_broad_ref_validate)/ {print}')
if [[ -n $active_conflicts ]]; then
  printf 'Conflicting active job(s) found; stop and inspect:\n%s\n' "$active_conflicts" >&2
  exit 1
fi

recovery_dir=$(mktemp -d "$run_dir/aldinger_vascular_recovery_XXXXXX")
cp -p "$ledger" "$recovery_dir/jobs.before.tsv"
printf '%s\n' \
  'Aldinger source label Endothelial was not converted in clusters_refined.' \
  'Replacement runs the corrected Aldinger relabeling and panel stages, then requeues validation under the same job ID.' \
  > "$recovery_dir/recovery_reason.txt"
printf 'Recovery record: %s\n' "$recovery_dir"

submit_held() {
  local label=$1 script=$2 dependency=${3:-}
  local args=(sbatch --parsable --hold --export=ALL)
  if [[ -n $dependency ]]; then
    args+=(--dependency="afterok:$dependency")
  fi
  local job
  job=$("${args[@]}" "scripts/$script")
  job=${job%%;*}
  if [[ ! $job =~ ^[0-9]+$ ]]; then
    printf 'Unexpected sbatch response for %s: %s\n' "$label" "$job" >&2
    exit 1
  fi
  printf '%s\t%s\n' "$label" "$job" >> "$recovery_dir/replacement_jobs.tsv"
  printf '%s: %s (held)\n' "$label" "$job" >&2
  printf '%s' "$job"
}

aldinger_standardize=$(submit_held Aldinger_standardize run_aldinger_02_standardize_cell_types.slurm)
aldinger_panel=$(submit_held Aldinger_panel run_aldinger_03_subset_gene_panel.slurm "$aldinger_standardize")

# Preserve downstream references to job 41015200 by rerunning that job under
# the same ID, held until the corrected Aldinger panel completes.
scontrol requeuehold 41015200
scontrol update JobId=41015200 Dependency="afterok:$aldinger_panel"

original_ledger=$(<"$ledger")
if [[ $original_ledger != $(<"$recovery_dir/jobs.before.tsv") ]]; then
  printf 'Ledger changed concurrently. New jobs and validation remain held; inspect %s\n' "$recovery_dir" >&2
  exit 1
fi

staged_ledger=$(mktemp "$run_dir/.jobs_aldinger_recovery.XXXXXX")
awk -F '\t' -v OFS='\t' -v standardize="$aldinger_standardize" -v panel="$aldinger_panel" '
  $1 == "annotation" && $2 == "Aldinger_standardize" {$3 = standardize}
  $1 == "annotation" && $2 == "Aldinger_panel" {$3 = panel}
  {print}
' "$ledger" > "$staged_ledger"
chmod --reference="$ledger" "$staged_ledger"
mv "$staged_ledger" "$ledger"

scontrol release 41015200
scontrol release "$aldinger_panel"
scontrol release "$aldinger_standardize"
printf 'Recovery submitted. Updated ledger: %s\n' "$ledger"
squeue -u acflint -o '%.24i %.36j %.10T %.12M %.100E %R'
