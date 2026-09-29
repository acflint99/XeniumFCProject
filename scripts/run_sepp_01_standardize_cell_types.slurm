#!/bin/bash
#SBATCH --job-name=Sepp_standardize
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Sepp_standardize_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Sepp_standardize_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript scripts/sepp_01_standardize_cell_types.R
