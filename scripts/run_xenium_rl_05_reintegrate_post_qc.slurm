#!/bin/bash
#SBATCH --job-name=RL_MergePostQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd ~/R/Projects/XeniumFCProject

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

Rscript scripts/xenium_rl_05_reintegrate_post_qc.R
