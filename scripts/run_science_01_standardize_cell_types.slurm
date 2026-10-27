#!/bin/bash
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --job-name=Science_standardize
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Science_standardize_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Science_standardize_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
# 128G is required through scaling for the current reviewed Science reference.
#SBATCH --mem=128G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript --max-ppsize=500000 scripts/science_01_standardize_cell_types.R --dry-run || exit 2

driver_options=()
if [[ "${SCIENCE_REFERENCE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/science_01_standardize_cell_types.R "${driver_options[@]}"
