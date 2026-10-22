#!/bin/bash
# Preflight all mapped samples, submit their report array, then merge afterok.
# Expected preflight: under one minute, file/CSV metadata only; no RDS loads.
set -euo pipefail
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != --dry-run ) ]]; then
  echo "Usage: bash scripts/submit_xenium_residual_reports.sh [--dry-run]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "${script_dir}/.."
task_map="$(Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R --list)"
task_count="$(printf '%s\n' "${task_map}" | awk -F '\t' '
  $1 ~ /^[0-9]+$/ {count++; if ($1 != count || NF != 2 || $2 == "" || seen[$2]++) bad=1}
  END {if (!count || bad) exit 2; print count}
')"
printf '%s\n' "${task_map}"
Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R --dry-run
if [[ "${1:-}" == --dry-run ]]; then
  echo "Would submit ${task_count} sample tasks at concurrency 3, followed by one dependent combine task."
  exit 0
fi
sample_job="$(sbatch --parsable "--array=1-${task_count}%3" \
  --export=ALL,RESIDUAL_PLOTS_MODE=samples scripts/run_xenium_residual_07_plot_samples.slurm)"
sample_job="${sample_job%%;*}"
if [[ ! "${sample_job}" =~ ^[0-9]+$ ]]; then
  echo "Could not identify the residual sample array job: ${sample_job}" >&2
  exit 2
fi
echo "Submitted residual sample array: ${sample_job}"
combined_job="$(sbatch --parsable --array=1 "--dependency=afterok:${sample_job}" \
  --export=ALL,RESIDUAL_PLOTS_MODE=combine scripts/run_xenium_residual_07_plot_samples.slurm)"
echo "Submitted dependent residual combined report: ${combined_job}"
