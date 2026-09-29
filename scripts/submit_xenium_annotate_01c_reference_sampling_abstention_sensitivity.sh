#!/bin/bash
# Validate and submit the isolated 18-task sensitivity array, followed by the
# combine job only after every array task succeeds.

set -euo pipefail

if [[ "$#" -gt 1 || ( "$#" -eq 1 && "$1" != "--dry-run" ) ]]; then
  echo "Usage: bash scripts/submit_xenium_annotate_01c_reference_sampling_abstention_sensitivity.sh [--dry-run]" >&2
  exit 2
fi

cd /home/acflint/R/Projects/XeniumFCProject

Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R --list
for task_id in $(seq 1 18); do
  Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R \
    --dry-run "${task_id}"
done

if [[ "${1:-}" == "--dry-run" ]]; then
  echo "DRY-RUN PASS: all 18 task mappings and inputs are valid; no jobs submitted."
  exit 0
fi

array_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01c_reference_sampling_abstention_sensitivity.slurm)
array_job_id="${array_submit%%;*}"
combine_submit=$(sbatch --parsable --dependency="afterok:${array_job_id}" \
  scripts/run_xenium_annotate_01c_reference_sampling_abstention_sensitivity_combine.slurm)
combine_job_id="${combine_submit%%;*}"

echo "Submitted sensitivity array job: ${array_job_id}"
echo "Submitted dependent combine job: ${combine_job_id}"
