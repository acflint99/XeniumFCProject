#!/bin/bash
#SBATCH --job-name=RL_Analysis_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_Analysis_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_Analysis_job_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=short

# Load R module first
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# Move into the project directory (CRITICAL)
cd ~/R/Projects/XeniumFCProject

# Restore the renv environment (installs missing packages if needed)
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org/"); library(renv); renv::restore(prompt = FALSE)'

# Run your R script using the project library
Rscript scripts/Xenium_RL_Merge_postQC_Processing_res1.5.R
