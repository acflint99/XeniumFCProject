#!/bin/bash
#SBATCH --job-name=Xen_unknown_QC_apply
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_QC_apply_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_QC_apply_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

driver_options=()
if [[ "${UNKNOWN_QC_APPLY_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_unknown_03_apply_qc_review.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
