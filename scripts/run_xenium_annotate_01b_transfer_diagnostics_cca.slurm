#!/bin/bash
# Nine isolated matched-10k CCA diagnostics: three resolution-5 samples and
# three references. Each task evaluates k.weight 25 and 50 from one
# k.anchor=10 AnchorSet. No production output is replaced.

#SBATCH --job-name=Xen_res5_transfer_diag_cca
#SBATCH --array=19-27%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --export=ALL
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=FAIL

module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript scripts/xenium_annotate_01b_transfer_diagnostics.R \
  --dry-run "${SLURM_ARRAY_TASK_ID}" || exit 2

driver_options=()
if [[ "${TRANSFER_DIAGNOSTICS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 \
  scripts/xenium_annotate_01b_transfer_diagnostics.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
