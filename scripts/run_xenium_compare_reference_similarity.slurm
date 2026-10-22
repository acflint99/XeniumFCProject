#!/bin/bash
#SBATCH --job-name=Xen_RefSimilarity
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RefSimilarity_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RefSimilarity_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Expected runtime: 1-6 hours or longer, low confidence; queue wait excluded.
# Requested limit: 12 hours. No array; both query levels and all references run once.
# First-run memory envelope, not a measured requirement for this analysis:
# merged-object expression job 41005351 peaked at 22.86 GiB; related VZ
# job 40966278 peaked at 152.09 GiB. Allow headroom for original references
# and validation allocations; refine from this job's time/sacct measurements.
# short's 12-hour limit checked at https://docs.rc.uab.edu/cheaha/hardware/
# on 2026-10-21. Preserve existing R module and project renv activation.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Optional environment overrides select inputs/output paths, never overwrite.
panel_json="${REFERENCE_PANEL_JSON:-data/FCXeniumProject/FB124_X_G/gene_panel.json}"
driver_options=(--references=aldinger,sepp,science --bootstrap=100 "--panel-json=${panel_json}")
if [[ -n "${REFERENCE_SIMILARITY_OUT:-}" ]]; then
  driver_options+=("--out=${REFERENCE_SIMILARITY_OUT}")
fi

echo "Reference similarity: job=${SLURM_JOB_ID:-manual}; start=$(date -Is)"
echo "Panel: ${panel_json}; output: ${REFERENCE_SIMILARITY_OUT:-outputs/xenium/reference_similarity}"
# Lightweight preflight: paths, dependencies and small original panel metadata only.
Rscript --max-ppsize=500000 scripts/xenium_compare_reference_similarity.R "${driver_options[@]}" --dry-run
# One process loads the originals sequentially; time reports peak resident memory
# to stderr even if the R process is terminated. No package installation/restore.
/usr/bin/time -v Rscript --max-ppsize=500000 scripts/xenium_compare_reference_similarity.R "${driver_options[@]}"
echo "Reference similarity completed: end=$(date -Is)"
