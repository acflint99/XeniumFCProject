#!/bin/bash
#SBATCH --job-name=Xen_RL_Plots
#SBATCH --array=1-34%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Plots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Plots_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

plot_mode="${RL_PLOTS_MODE:-samples}"
driver_options=()
if [[ "${RL_PLOTS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

case "${plot_mode}" in
  samples)
    Rscript --max-ppsize=500000 scripts/xenium_rl_08_plot_samples.R \
      "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
    ;;
  combine)
    if [[ "${SLURM_ARRAY_TASK_COUNT:-1}" -ne 1 ]]; then
      echo "RL combine mode must be submitted with --array=1." >&2
      exit 2
    fi
    Rscript --max-ppsize=500000 scripts/xenium_rl_08_plot_samples.R \
      --combine "${driver_options[@]}"
    ;;
  *)
    echo "RL_PLOTS_MODE must be 'samples' or 'combine'; received '${plot_mode}'." >&2
    exit 2
    ;;
esac
