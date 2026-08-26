#!/bin/bash
#SBATCH --job-name=RL_MergePostQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=medium
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

driver_options=()
if [[ "${RL_POST_QC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_rl_05_reintegrate_post_qc.R "${driver_options[@]}"
