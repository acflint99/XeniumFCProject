#!/bin/bash
#SBATCH --job-name=Xen_broad_scores
#SBATCH --array=1-123%6
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --export=ALL
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=FAIL,INVALID_DEPEND
module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript --max-ppsize=500000 scripts/xenium_broad_02_score_clusters.R \
  --dry-run "${SLURM_ARRAY_TASK_ID}" || exit 2

driver_options=()
if [[ "${BROAD_SCORE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_broad_02_score_clusters.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
