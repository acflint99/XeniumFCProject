#!/bin/bash
#SBATCH --job-name=Aldinger_standardize
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Aldinger_standardize_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Aldinger_standardize_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

# Expected 5-60 min, low-confidence static estimate; no matching measured job.
# Use the established Sepp standardization resource envelope; collect sacct.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
Rscript --max-ppsize=500000 scripts/aldinger_02_standardize_cell_types.R
