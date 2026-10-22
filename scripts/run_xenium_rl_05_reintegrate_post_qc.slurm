#!/bin/bash
#SBATCH --job-name=RL_MergePostQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_MergePostQC_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=600G
#SBATCH --partition=medium
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# Default to sequential clustering after cohort expansion and parallel OOMs.
# An explicit worker override is bounded by the allocated CPU count.
export RL_POST_QC_CLUSTER_WORKERS="${RL_POST_QC_CLUSTER_WORKERS:-1}"

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

driver_options=()
if [[ "${RL_POST_QC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_rl_05_reintegrate_post_qc.R "${driver_options[@]}"
