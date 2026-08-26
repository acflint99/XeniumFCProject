#!/bin/bash
#SBATCH --job-name=Xen_VZ_Subclusters
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Subclusters_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Subclusters_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
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
if [[ "${VZ_SUBCLUSTER_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_vz_06_annotate_subclusters.R "${driver_options[@]}"
