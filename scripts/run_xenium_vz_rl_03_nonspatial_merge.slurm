#!/bin/bash
#SBATCH --job-name=Xen_CleanMerge
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CleanMerge_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CleanMerge_%j.err
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
if [[ "${CLEAN_MERGE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
merge_preflight="$(Rscript --max-ppsize=500000 \
  scripts/xenium_vz_rl_03_nonspatial_merge.R --dry-run 2>&1)" || {
  printf '%s\n' "${merge_preflight}" >&2
  exit 2
}
printf '%s\n' "${merge_preflight}"
if [[ "${merge_preflight}" != *"DRY-RUN PASS"* ]]; then
  echo "Corrected-broad merge input preflight failed; no analysis outputs were written." >&2
  exit 2
fi
Rscript --max-ppsize=500000 'scripts/xenium_vz_rl_03_nonspatial_merge.R' "${driver_options[@]}"
