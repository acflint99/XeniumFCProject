#!/bin/bash
#SBATCH --job-name=VZSamplePlots_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZSamplePlots_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZSamplePlots_job_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=short

# Load R module first
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# 1. Project Directory
cd ~/R/Projects/XeniumFCProject

# 2. Tell R where your library is
export R_LIBS_USER=~/R/Projects/XeniumFCProject/renv/library/linux-rhel-7.9/R-4.4/x86_64-pc-linux-gnu

# 3. Use 'renv::refresh' to fix version mismatches automatically
Rscript -e 'renv::refresh()'

# 4. Run the main analysis
Rscript scripts/XeniumVZSamplePlots_Loop_res1.5.R