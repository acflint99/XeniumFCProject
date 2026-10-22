#!/bin/bash
# Run after loading the project's R module, from any working directory.
# Preflight: <1 minute, paths only; analysis: 1-5 min/task (low confidence).
set -euo pipefail
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != --dry-run ) ]]; then
  echo "Usage: bash scripts/submit_xenium_sex_marker_expression.sh [--dry-run]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "${script_dir}/.."
driver=scripts/xenium_sex_marker_expression.R
launcher=scripts/run_xenium_sex_marker_expression.slurm
concurrency="${SEX_MARKER_CONCURRENCY:-3}"
if [[ ! "${concurrency}" =~ ^[1-9][0-9]*$ ]]; then
  echo "SEX_MARKER_CONCURRENCY must be a positive integer." >&2
  exit 2
fi
export SEX_MARKER_OUTPUT_DIR="${SEX_MARKER_OUTPUT_DIR:-outputs/xenium/sex_marker_expression}"
task_map="$(Rscript --max-ppsize=500000 "${driver}" --list)"
task_count="$(printf '%s\n' "${task_map}" | awk -F '\t' '
  $1 ~ /^[0-9]+$/ {
    count++
    if ($1 != count || NF != 2 || $2 == "" || seen[$2]++) bad=1
  }
  END {if (!count || bad) exit 2; print count}
')"
export SEX_MARKER_MAP_HASH="$(printf '%s\n' "${task_map}" | sha256sum | awk '{print $1}')"
printf '%s\n' "${task_map}"
Rscript --max-ppsize=500000 "${driver}" --dry-run --output-dir "${SEX_MARKER_OUTPUT_DIR}"
array=(sbatch --parsable "--array=1-${task_count}%${concurrency}" "${launcher}")
printf 'Sample array: '
printf '%q ' "${array[@]}"
printf '\nCombine: afterok dependency, 4 GB, 15-minute limit, same launcher.\n'
printf 'Email: acflint@uab.edu; array FAIL; combine END,FAIL,INVALID_DEPEND.\n'
if [[ "${1:-}" == --dry-run ]]; then exit 0; fi
mkdir -p logs
export SEX_MARKER_MODE=samples
job_id="$("${array[@]}")"
job_id="${job_id%%;*}"
if [[ ! "${job_id}" =~ ^[0-9]+$ ]]; then
  echo "Could not parse array job ID: ${job_id}" >&2
  exit 2
fi
echo "Sample array submitted: ${job_id}"
export SEX_MARKER_MODE=combine
sbatch --dependency="afterok:${job_id}" --mem=4G --time=00:15:00 \
  --mail-type=END,FAIL,INVALID_DEPEND "${launcher}"
