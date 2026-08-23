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

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

# Run the analysis using the existing project library
Rscript scripts/Xenium_RL_Merge_postQC_Processing_res1.5.R
