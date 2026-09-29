#!/bin/bash
# Validate and submit the isolated 27-task cap-convergence array, followed by the
# combine job only after every array task succeeds.

set -euo pipefail

if [[ "$#" -gt 1 || ( "$#" -eq 1 && "$1" != "--dry-run" ) ]]; then
  echo "Usage: bash scripts/submit_xenium_annotate_01e_reference_cap_convergence.sh [--dry-run]" >&2
  exit 2
fi

cd /home/acflint/R/Projects/XeniumFCProject

existing_combined_dir="outputs/validation/resolution5_method_development/reference_sampling_abstention_sensitivity/combined_tables"
required_existing_tables=(
  "${existing_combined_dir}/sampling_abstention_cluster_metrics.csv"
  "${existing_combined_dir}/sampling_abstention_reference_votes.csv"
  "${existing_combined_dir}/sampling_abstention_reference_class_balance.csv"
)
for required_table in "${required_existing_tables[@]}"; do
  if [[ ! -f "${required_table}" ]]; then
    echo "Missing completed cap1000/full input: ${required_table}" >&2
    exit 1
  fi
done

Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R --list
for task_id in $(seq 1 27); do
  Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R \
    --dry-run "${task_id}"
done

if [[ "${1:-}" == "--dry-run" ]]; then
  echo "DRY-RUN PASS: all 27 task mappings and inputs are valid; no jobs submitted."
  exit 0
fi

array_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01e_reference_cap_convergence.slurm)
array_job_id="${array_submit%%;*}"
combine_submit=$(sbatch --parsable --dependency="afterok:${array_job_id}" \
  scripts/run_xenium_annotate_01e_reference_cap_convergence_combine.slurm)
combine_job_id="${combine_submit%%;*}"

echo "Submitted cap-convergence array job: ${array_job_id}"
echo "Submitted dependent combine job: ${combine_job_id}"
