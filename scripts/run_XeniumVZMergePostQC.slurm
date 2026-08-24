#!/bin/bash
#SBATCH --job-name=VZ_MergePostQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_MergePostQC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_MergePostQC_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd ~/R/Projects/XeniumFCProject

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

Rscript scripts/Xenium_VZ_Merge_postQC_Processing_res1.5.R
