#!/bin/bash
# Derive sample-array bounds from the selected R driver's manifest task map.
# Expected preflight runtime: under one minute; no RDS objects are loaded.
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash scripts/submit_xenium_downstream_array.sh {vz-plots|rl-plots|residual-plots|residual-subset} [--dry-run]" >&2
  exit 2
fi
stage="$1"
if [[ $# -eq 2 && "$2" != "--dry-run" ]]; then
  echo "The only optional argument is --dry-run." >&2
  exit 2
fi
case "${stage}" in
  vz-plots)
    driver=xenium_vz_08_plot_samples
    if [[ "${VZ_PLOTS_MODE:-samples}" != samples ]]; then
      echo "This helper submits sample mode only; submit report combine mode separately with --array=1." >&2
      exit 2
    fi
    ;;
  rl-plots)
    driver=xenium_rl_08_plot_samples
    if [[ "${RL_PLOTS_MODE:-samples}" != samples ]]; then
      echo "This helper submits sample mode only; submit report combine mode separately with --array=1." >&2
      exit 2
    fi
    ;;
  residual-subset) driver=xenium_residual_01_subset_samples ;;
  residual-plots)
    driver=xenium_residual_07_plot_samples
    if [[ "${RESIDUAL_PLOTS_MODE:-samples}" != samples ]]; then
      echo "This helper submits sample mode only; submit residual combine mode separately with --array=1." >&2
      exit 2
    fi
    ;;
  *) echo "Unknown stage: ${stage}" >&2; exit 2 ;;
esac
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "${script_dir}/.."
task_map="$(Rscript --max-ppsize=500000 "scripts/${driver}.R" --list)"
# Startup notices are allowed; numbered TSV rows must be contiguous and nonempty.
task_count="$(printf '%s\n' "${task_map}" | awk -F '\t' '
  $1 ~ /^[0-9]+$/ {
    count++
    if ($1 != count || NF != 2 || $2 == "" || seen[$2]++) bad=1
  }
  END {if (!count || bad) exit 2; print count}
')"
printf '%s\n' "${task_map}"
submission=(sbatch "--array=1-${task_count}%3" "scripts/run_${driver}.slurm")
printf 'Submission: '
printf '%q ' "${submission[@]}"
printf '\n'
if [[ "${2:-}" != --dry-run ]]; then
  "${submission[@]}"
fi
