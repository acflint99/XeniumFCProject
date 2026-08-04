#!/bin/bash
#SBATCH --job-name=Xen_res1.5_ABT
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/xen_ABT_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/xen_ABT_%a.err
#SBATCH --array=1-9
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                 # Stay at 128G if you can; 64G is risky for res 1.5
#SBATCH --partition=medium

# 1. Load Modules
module purge
module load rc-base
module load GCC/13.2.0    # Add this to provide the updated libstdc++.so.6
module load R/4.4.0-foss-2022b

# 2. Move into project directory
cd ~/R/Projects/XeniumFCProject

# 3. Disable renv sandbox 
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

# --------------------------
# Run your R script using the project library
# --------------------------
Rscript scripts/AnchorBasedTransfer_RPCA_res1.5.R $SLURM_ARRAY_TASK_ID