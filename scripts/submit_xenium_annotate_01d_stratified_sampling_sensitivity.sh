#!/bin/bash
# Submit preflight -> nine-task stratified array -> three-mode combine. Every
# downstream stage is released only when the preceding stage succeeds.

set -euo pipefail

if [[ "$#" -gt 1 || ( "$#" -eq 1 && "$1" != "--dry-run" ) ]]; then
  echo "Usage: bash scripts/submit_xenium_annotate_01d_stratified_sampling_sensitivity.sh [--dry-run]" >&2
  exit 2
fi

cd /home/acflint/R/Projects/XeniumFCProject

Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R --list
Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R \
  --inspect-strata --dry-run

if [[ "${1:-}" == "--dry-run" ]]; then
  echo "DRY-RUN PASS: reference paths and the nine-task mapping are valid; no jobs submitted."
  echo "Each array task will run its path-only dry-run after the preflight creates the strata manifest."
  exit 0
fi

preflight_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01d_stratified_sampling_preflight.slurm)
preflight_job_id="${preflight_submit%%;*}"
array_submit=$(sbatch --parsable --dependency="afterok:${preflight_job_id}" \
  scripts/run_xenium_annotate_01d_stratified_sampling_sensitivity.slurm)
array_job_id="${array_submit%%;*}"
combine_submit=$(sbatch --parsable --dependency="afterok:${array_job_id}" \
  scripts/run_xenium_annotate_01d_stratified_sampling_sensitivity_combine.slurm)
combine_job_id="${combine_submit%%;*}"

echo "Submitted age/stage preflight job: ${preflight_job_id}"
echo "Submitted dependent stratified array job: ${array_job_id}"
echo "Submitted dependent three-mode combine job: ${combine_job_id}"
