#!/bin/bash
#SBATCH --job-name=Xen_MetaNeighbor
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_MetaNeighbor_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_MetaNeighbor_%A_%a.err
#SBATCH --array=1-3%1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Expected 0.5-4 h/task, 1.5-12 h for all three at concurrency 1; low confidence.
# Queue wait excluded. Requested limit is 12 h PER TASK, not expected runtime.
# No prior MetaNeighbor sacct evidence. Original merged read/aggregation reached
# ~35.3 GiB in partial reference-similarity job 41010260; nearby VZ work reached
# 152.09 GiB (history job 40966278). Dense rank allocations add substantial memory.
# Keep a provisional 256-GiB first-run envelope and measure with time/sacct.
# CPU short partition limit verified at https://docs.rc.uab.edu/cheaha/hardware/
# on 2026-10-22. Fixed array maps references, never samples; driver validates
# the complete query sample set against the active authoritative manifest.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
case "${SLURM_ARRAY_TASK_ID:?Submit this launcher as its configured array}" in
  1) reference=aldinger ;;
  2) reference=sepp ;;
  3) reference=science ;;
  *) echo "Invalid reference array task." >&2; exit 2 ;;
esac
out_root="${METANEIGHBOR_OUT:-outputs/xenium/metaneighbor_original}"
panel="${REFERENCE_PANEL_JSON:-data/FCXeniumProject/FB124_X_G/gene_panel.json}"
options=("--reference=${reference}" "--out=${out_root}" "--panel-json=${panel}")
echo "MetaNeighborUS: job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; reference=${reference}; start=$(date -Is)"
Rscript --max-ppsize=500000 scripts/xenium_compare_metaneighbor.R --list
Rscript --max-ppsize=500000 scripts/xenium_compare_metaneighbor.R "${options[@]}" --dry-run
# Tiny synthetic engine test (~10-60 seconds) before reading original objects.
Rscript --max-ppsize=500000 scripts/test_xenium_metaneighbor_engine.R
/usr/bin/time -v Rscript --max-ppsize=500000 scripts/xenium_compare_metaneighbor.R "${options[@]}"
echo "MetaNeighborUS completed: reference=${reference}; end=$(date -Is)"
