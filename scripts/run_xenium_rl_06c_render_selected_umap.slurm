#!/bin/bash
#SBATCH --job-name=RL_UMAP_Selected
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_UMAP_Selected_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_UMAP_Selected_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=express
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
driver_options=()
mode_options=(--render-selected)
case "${RL_UMAP_SCOPE:-production}" in
  production) ;;
  upstream) mode_options+=(--upstream-only) ;;
  *) echo "RL_UMAP_SCOPE must be 'production' or 'upstream'." >&2; exit 2 ;;
esac
if [[ "${RL_UMAP_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R "${mode_options[@]}" --dry-run
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R "${mode_options[@]}" "${driver_options[@]}"
