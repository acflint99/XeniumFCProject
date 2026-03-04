#!/bin/bash
#SBATCH --job-name=Sepp_R_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Sepp_R_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Sepp_R_job_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=short

# Load R module first
module purge
module load rc-base
module load R/4.4.0-foss-2022b

# Move into the project directory (CRITICAL)
cd ~/R/Projects/XeniumFCProject

# Restore the renv environment (installs missing packages if needed)
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org/"); library(renv); renv::restore(prompt = FALSE)'

# Run your R script using the project library
Rscript scripts/SeppGenePanelSubset.R
