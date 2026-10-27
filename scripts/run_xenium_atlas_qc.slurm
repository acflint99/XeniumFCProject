#!/bin/bash
# Expected 2-10 min/sample; 28-140 min for 41 tasks at concurrency 3.
# Low confidence: RDS I/O estimate, no matching sacct measurement yet.
# 64 GB is provisional headroom; collect MaxRSS after the first run.
#SBATCH --job-name=Xen_atlas_qc
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_atlas_qc_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_atlas_qc_%A_%a.err
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

driver=scripts/xenium_atlas_qc.R
task_map="$(Rscript --max-ppsize=500000 "${driver}" --list)"
map_hash="$(printf '%s\n' "${task_map}" | sha256sum | awk '{print $1}')"
if [[ -z "${ATLAS_QC_MAP_HASH:-}" || "${map_hash}" != "${ATLAS_QC_MAP_HASH}" ]]; then
  echo "Task mapping changed or is unverified. Use submit_xenium_atlas_qc.sh." >&2
  exit 2
fi
options=(--output-dir "${ATLAS_QC_OUTPUT_DIR:-outputs/xenium/atlas_qc}")
if [[ "${ATLAS_QC_MODE:-samples}" == combine ]]; then
  options+=(--compile)
  if [[ "${ATLAS_QC_ALLOW_INCOMPLETE:-false}" == true ]]; then options+=(--allow-incomplete); fi
elif [[ "${ATLAS_QC_MODE:-samples}" == samples ]]; then
  options+=(--task "${SLURM_ARRAY_TASK_ID:?Array task ID required}")
else
  echo "Invalid ATLAS_QC_MODE." >&2
  exit 2
fi
Rscript --max-ppsize=500000 "${driver}" "${options[@]}"

if [[ "${ATLAS_QC_MODE:-samples}" == combine ]]; then
  Rscript --max-ppsize=500000 scripts/xenium_atlas_plot_qc.R --output-dir "${ATLAS_QC_OUTPUT_DIR:-outputs/xenium/atlas_qc}"
fi
