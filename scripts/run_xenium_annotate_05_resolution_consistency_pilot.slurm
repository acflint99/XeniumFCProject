#!/bin/bash
# Compare weighted-2-of-3 consensus labels across completed resolution-3,
# resolution-4, and resolution-5 pilot objects for the three manifest samples.
# This audit does not modify any Seurat object.

#SBATCH --job-name=Xen_res_consistency
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res_consistency_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res_consistency_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
if [[ "${RESOLUTION_CONSISTENCY_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
