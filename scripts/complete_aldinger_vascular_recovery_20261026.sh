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
    printf 'Ledger differs from the reviewed completion state: %s\n' "$label" >&2
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

require_pending() {
  local job_id=$1 state
  state=$(squeue -h -j "$job_id" -o '%T')
  if [[ $state != PENDING ]]; then
    printf 'Expected held replacement job %s to be pending; found: %s\n' "$job_id" "$state" >&2
    exit 1
  fi
}

require_ledger_job Aldinger_standardize 41015193
require_ledger_job Aldinger_panel 41015194
require_ledger_job Science_standardize 41049984
require_ledger_job Science_panel 41052256
require_ledger_job reference_validation 41015200

require_sacct 41049984 COMPLETED 0:0
require_sacct 41052256 COMPLETED 0:0
require_sacct 41015200 FAILED 1:0
require_pending 41052564
require_pending 41052565

active_conflicts=$(squeue -h -u acflint -o '%i|%j' |
  awk -F '|' '$2 ~ /^(Aldinger_|Xen_broad_ref_validate)/ && $1 !~ /^4105256[45]$/ {print}')
if [[ -n $active_conflicts ]]; then
  printf 'Conflicting active job(s) found; stop and inspect:\n%s\n' "$active_conflicts" >&2
  exit 1
fi

recovery_dir=$(mktemp -d "$run_dir/aldinger_vascular_completion_XXXXXX")
cp -p "$ledger" "$recovery_dir/jobs.before.tsv"
printf '%s\n' \
  'Validation job 41015200 had aged out of the live Slurm controller and could not be requeued.' \
  'A new validation job will replace it; only the three transfer arrays are repointed to the new job.' \
  > "$recovery_dir/recovery_reason.txt"
printf 'Recovery record: %s\n' "$recovery_dir"

validation=$(sbatch --parsable --hold \
  --export=ALL,BROAD_REFERENCE_VALIDATION_OVERWRITE=true \
  --dependency=afterok:41052565 \
  scripts/run_xenium_broad_01_validate_references.slurm)
validation=${validation%%;*}
if [[ ! $validation =~ ^[0-9]+$ ]]; then
  printf 'Unexpected validation submission response: %s\n' "$validation" >&2
  exit 1
fi
printf 'reference_validation\t%s\n' "$validation" > "$recovery_dir/replacement_jobs.tsv"
printf 'reference_validation: %s (held)\n' "$validation"

for transfer in 41015201 41015202 41015203; do
  if [[ $(squeue -h -j "$transfer" -o '%T') != PENDING ]]; then
    printf 'Transfer job %s is no longer pending; new validation remains held.\n' "$transfer" >&2
    exit 1
  fi
  scontrol update JobId="$transfer" Dependency="afterok:$validation"
done

if [[ $(<"$ledger") != $(<"$recovery_dir/jobs.before.tsv") ]]; then
  printf 'Ledger changed concurrently. All replacements remain held; inspect %s\n' "$recovery_dir" >&2
  exit 1
fi

staged_ledger=$(mktemp "$run_dir/.jobs_aldinger_completion.XXXXXX")
awk -F '\t' -v OFS='\t' -v standardize=41052564 -v panel=41052565 -v validation="$validation" '
  $1 == "annotation" && $2 == "Aldinger_standardize" {$3 = standardize}
  $1 == "annotation" && $2 == "Aldinger_panel" {$3 = panel}
  $1 == "annotation" && $2 == "reference_validation" {$3 = validation}
  {print}
' "$ledger" > "$staged_ledger"
chmod --reference="$ledger" "$staged_ledger"
mv "$staged_ledger" "$ledger"

scontrol release "$validation"
scontrol release 41052565
scontrol release 41052564
printf 'Recovery submitted. Updated ledger: %s\n' "$ledger"
squeue -u acflint -o '%.24i %.36j %.10T %.12M %.100E %R'
