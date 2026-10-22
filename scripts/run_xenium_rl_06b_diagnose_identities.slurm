#!/bin/bash
#SBATCH --job-name=Xen_RL_identity_review
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
driver_options=(scripts/xenium_rl_06b_diagnose_identities.R)
if [[ "${RL_IDENTITY_DIAGNOSTICS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_rl_06b_diagnose_identities.R --dry-run
Rscript --max-ppsize=500000 "${driver_options[@]}"
