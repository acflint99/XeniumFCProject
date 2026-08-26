#!/bin/bash
#SBATCH --job-name=Xen_VZ_Merge
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Merge_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Merge_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=medium
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
if [[ "${VZ_MERGE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_vz_02_merge_samples.R "${driver_options[@]}"
