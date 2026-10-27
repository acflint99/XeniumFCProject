#!/bin/bash
# Expected 1-5 min/sample; 15-70 min for 41 tasks at concurrency 3.
# Low confidence: RDS I/O estimate, no matching sacct measurement yet.
# 64 GB is provisional headroom; collect MaxRSS after the first run.
#SBATCH --job-name=Xen_sex_markers
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_sex_markers_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_sex_markers_%A_%a.err
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject

driver=scripts/xenium_sex_marker_expression.R
task_map="$(Rscript --max-ppsize=500000 "${driver}" --list)"
map_hash="$(printf '%s\n' "${task_map}" | sha256sum | awk '{print $1}')"
if [[ -z "${SEX_MARKER_MAP_HASH:-}" || "${map_hash}" != "${SEX_MARKER_MAP_HASH}" ]]; then
  echo "Task mapping changed or is unverified. Use submit_xenium_sex_marker_expression.sh." >&2
  exit 2
fi
options=(--output-dir "${SEX_MARKER_OUTPUT_DIR:-outputs/xenium/sex_marker_expression}")
if [[ "${SEX_MARKER_MODE:-samples}" == combine ]]; then
  options+=(--combine)
elif [[ "${SEX_MARKER_MODE:-samples}" == samples ]]; then
  options+=(--task "${SLURM_ARRAY_TASK_ID:?Array task ID required}")
else
  echo "Invalid SEX_MARKER_MODE." >&2
  exit 2
fi
Rscript --max-ppsize=500000 "${driver}" "${options[@]}"
