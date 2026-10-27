#!/bin/bash
#SBATCH --job-name=Xen_residual_postQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
set -eo pipefail
# Expected 15-120 min after exclusions; 1-10 min for no-removal reuse (~10k cells; low confidence).
# Reuse the existing integration envelope: requested 12 h / 128G / 4 CPUs, not measured post-QC requirements.
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
driver_options=()
inspection_options=()
case "${RESIDUAL_QC_ACTION:?Set remove or none after review}" in
  none) driver_options+=(--no-removals); inspection_options+=(--no-removals) ;;
  remove) ;;
  *) echo "RESIDUAL_QC_ACTION must be remove or none" >&2; exit 2 ;;
esac
if [[ "${RESIDUAL_POST_QC_OVERWRITE:-false}" == "true" ]]; then driver_options+=(--overwrite); fi
Rscript --max-ppsize=500000 scripts/xenium_residual_05_reintegrate_post_qc.R --dry-run "${inspection_options[@]}" || exit 2
Rscript --max-ppsize=500000 scripts/xenium_residual_05_reintegrate_post_qc.R "${driver_options[@]}"
