#!/bin/bash
#SBATCH --job-name=VZ_MergePostQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_MergePostQC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_MergePostQC_%j.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=4
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
export VZ_POST_QC_CLUSTER_WORKERS="${VZ_POST_QC_CLUSTER_WORKERS:-${SLURM_CPUS_PER_TASK}}"

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

driver_options=()
if [[ "${VZ_POST_QC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
if [[ "${VZ_POST_QC_RESUME_CHECKPOINT:-false}" == "true" ]]; then
  driver_options+=(--resume-checkpoint)
fi

Rscript scripts/xenium_vz_05_reintegrate_post_qc.R "${driver_options[@]}"
