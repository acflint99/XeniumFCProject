#!/bin/bash
#SBATCH --job-name=Science_R_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Science_R_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Science_R_job_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=short

# Load R module first
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# Move into the project directory (CRITICAL)
cd ~/R/Projects/XeniumFCProject

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

# Run the analysis using the existing project library
Rscript scripts/science_02_subset_gene_panel.R
