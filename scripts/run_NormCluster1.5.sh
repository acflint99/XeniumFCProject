#!/bin/bash
#SBATCH --job-name=Xen_res1.5
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/xen_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/xen_%a.err
#SBATCH --array=1-16
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                 # Stay at 128G if you can; 64G is risky for res 1.5
#SBATCH --partition=medium

# 1. Load Modules
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# 2. Move into project directory
cd ~/R/Projects/XeniumFCProject

# 3. Disable renv sandbox 
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

# 4. Run the R script
# Make sure the filename below matches your .R file exactly!
Rscript scripts/XeniumNormCluster_res1.5.R $SLURM_ARRAY_TASK_ID
