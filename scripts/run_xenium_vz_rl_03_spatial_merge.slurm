#!/bin/bash
#SBATCH --job-name=Xen_SpatialMerge
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_SpatialMerge_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_SpatialMerge_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=600G
#SBATCH --partition=medium
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=()
spatial_batch_size="${SPATIAL_MERGE_BATCH_SIZE:-6}"
if [[ ! "${spatial_batch_size}" =~ ^([2-9]|10)$ ]]; then
  echo "SPATIAL_MERGE_BATCH_SIZE must be an integer from 2 through 10." >&2
  exit 2
fi
driver_options+=("--batch-size=${spatial_batch_size}")
if [[ "${SPATIAL_MERGE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
spatial_preflight="$(Rscript --max-ppsize=500000 \
  scripts/xenium_vz_rl_03_spatial_merge.R --dry-run "--batch-size=${spatial_batch_size}" 2>&1)" || {
  printf '%s\n' "${spatial_preflight}" >&2
  exit 2
}
printf '%s\n' "${spatial_preflight}"
if [[ "${spatial_preflight}" != *"DRY-RUN PASS"* ]]; then
  echo "Corrected-broad spatial merge input preflight failed; no analysis outputs were written." >&2
  exit 2
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_03_spatial_merge.R "${driver_options[@]}"
