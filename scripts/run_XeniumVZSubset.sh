#!/bin/bash
#SBATCH --job-name=Xenium_VZSubsets
#SBATCH --array=1 #1-15 for nearly all samples
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G             
#SBATCH --partition=short

# 1. Load Modules
module purge
module load rc-base
module load R/4.4.0-foss-2022b

# 2. Move into project directory (CRITICAL for here() and renv)
cd ~/R/Projects/XeniumFCProject

# 3. Disable renv sandbox (Smart move—prevents slow symlinking on HPC)
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

# 4. Run the R script
# We skip 'renv::restore' inside the array to prevent 15 jobs 
# from fighting over the library lock simultaneously.
Rscript scripts/Xenium_VZ_Subset.R $SLURM_ARRAY_TASK_ID
