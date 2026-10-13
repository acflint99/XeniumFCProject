#!/bin/bash
#SBATCH --job-name=Xen_broad_apply
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript scripts/xenium_broad_04_apply_reviewed_thresholds.R --dry-run || exit 2

driver_options=()
if [[ "${BROAD_CONSENSUS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_broad_04_apply_reviewed_thresholds.R "${driver_options[@]}"
