#!/bin/bash
#SBATCH --job-name=Xen_RL_marker_evidence
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
driver_options=()
if [[ "${RL_MARKER_EVIDENCE_OVERWRITE:-false}" == "true" ]]; then driver_options+=(--overwrite); fi
Rscript --max-ppsize=500000 scripts/xenium_rl_06a_build_replicate_marker_evidence.R --dry-run || exit 2
Rscript --max-ppsize=500000 scripts/xenium_rl_06a_build_replicate_marker_evidence.R "${driver_options[@]}"
