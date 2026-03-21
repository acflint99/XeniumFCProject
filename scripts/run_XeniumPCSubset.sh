#!/bin/bash
#SBATCH --job-name=Xenium_PCSubsets
#SBATCH --array=1-16 #1-x# for all samples (based on how many listed in Xenium_PC_Subset.R)
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/PCsubset_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/PCsubset_%a.err
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
Rscript scripts/Xenium_PC_Subset.R $SLURM_ARRAY_TASK_ID
