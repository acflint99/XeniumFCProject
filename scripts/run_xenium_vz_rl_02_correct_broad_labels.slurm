#!/bin/bash
#SBATCH --job-name=Xen_CorrectBroad
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CorrectBroad_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CorrectBroad_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

preflight="$(Rscript --max-ppsize=500000 \
  scripts/xenium_vz_rl_02_correct_broad_labels.R --dry-run 2>&1)" || {
  printf '%s\n' "${preflight}" >&2
  exit 2
}
printf '%s\n' "${preflight}"
if [[ "${preflight}" != *"DRY-RUN PASS"* ]]; then
  echo "Corrected-broad input preflight failed; no objects were written." >&2
  exit 2
fi
driver_options=()
if [[ "${CORRECTED_BROAD_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 \
  scripts/xenium_vz_rl_02_correct_broad_labels.R "${driver_options[@]}" || exit 1
Rscript --max-ppsize=500000 scripts/validate_xenium_vz_rl_02_correct_broad.R
