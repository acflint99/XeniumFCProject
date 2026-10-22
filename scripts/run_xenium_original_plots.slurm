#!/bin/bash
# Resource/array overrides are supplied by submit_xenium_original_plots.sh.
#SBATCH --job-name=Xen_OriginalPlots
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_OriginalPlots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_OriginalPlots_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
stage="${XENIUM_ORIGINAL_STAGE:-consensus}"
if [[ "${XENIUM_ORIGINAL_OVERWRITE:-false}" != true ]]; then
  echo "This redraw replaces original plot files; set XENIUM_ORIGINAL_OVERWRITE=true." >&2
  exit 2
fi
set --
task="${SLURM_ARRAY_TASK_ID:-1}"
case "$stage" in
  consensus) driver=xenium_annotate_03_apply_consensus.R; set -- --all-samples-res5 --weighted-2of3 --plots-only "$task" ;;
  vz) driver=xenium_vz_08_plot_samples.R; set -- "$task" ;;
  rl) driver=xenium_rl_08_plot_samples.R; set -- "$task" ;;
  residual) driver=xenium_residual_07_plot_samples.R; set -- "$task" ;;
  combined) driver=refresh_combined_sample_plots_original.R; set -- "$task" ;;
  vz-pooled) driver=xenium_vz_06_annotate_subclusters.R; set -- --plots-only --resolution=0.8 ;;
  rl-pooled) driver=xenium_rl_06_annotate_subclusters.R; set -- --plots-only --resolution=0.8 ;;
  residual-pooled) driver=xenium_residual_05_plot_subclusters.R ;;
  counts-vz) driver=xenium_vz_09_plot_counts.R ;;
  counts-rl) driver=xenium_rl_09_plot_counts.R ;;
  nonspatial) driver=xenium_vz_rl_04_05_nonspatial_process_and_plot.R; set -- --plots-only ;;
  spatial) driver=refresh_spatial_umaps_original.R ;;
  vz-report) driver=xenium_vz_08_plot_samples.R; set -- --combine ;;
  rl-report) driver=xenium_rl_08_plot_samples.R; set -- --combine ;;
  residual-report) driver=xenium_residual_07_plot_samples.R; set -- --combine ;;
  combined-report) driver=xenium_vz_rl_reports.R; set -- --labels=both --point-size=0.6 ;;
  combined-report085) driver=xenium_vz_rl_reports.R; set -- --labels=both --point-size=0.85 ;;
  consensus-pages|consensus-report)
    driver=xenium_annotate_03d_plot_report.R
    Rscript --max-ppsize=500000 "scripts/${driver}" --all-samples-res5 --weighted-2of3 --dry-run
    if [[ "$stage" == consensus-pages ]]; then
      Rscript --max-ppsize=500000 "scripts/${driver}" --all-samples-res5 --weighted-2of3 --render-page --overwrite
    else
      Rscript --max-ppsize=500000 "scripts/${driver}" --all-samples-res5 --weighted-2of3 --merge-pages --overwrite
    fi
    exit 0 ;;
  *) echo "Unknown native plot stage: ${stage}" >&2; exit 2 ;;
esac
Rscript --max-ppsize=500000 "scripts/${driver}" "$@" --dry-run
# Consensus plot-only mode itself authorizes replacing plots and rejects
# --overwrite. Other drivers retain their explicit overwrite safeguard.
if [[ "$stage" != consensus ]]; then
  set -- "$@" --overwrite
fi
Rscript --max-ppsize=500000 "scripts/${driver}" "$@"
