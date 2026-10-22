#!/usr/bin/env bash
# Install the staged VZ-only DotPlot palette after the current VZ/RL workflow,
# then refresh pooled and per-sample plots without rerunning labels or merges.
set -euo pipefail

if [[ $# -ne 1 || ( "$1" != "--submit" && "$1" != "--install" ) ]]; then
  echo "Usage: bash scripts/submit_vz_marker_plot_refresh_after_jobs.sh --submit" >&2
  exit 2
fi

cd /home/acflint/R/Projects/XeniumFCProject
staged=scripts/color_palette_vz_markers_pending.R
backup=scripts/color_palette_before_vz_marker_prune_20261020.R
active_name_pattern='^(Xen_VZ_|Xen_Combine$|Xen_CorrectBroad$|Xen_CleanMerge$|Xen_Comb|Xen_Spatial|Xen_VZRL_)'

if [[ "$1" == "--install" ]]; then
  [[ -n "${SLURM_JOB_ID:-}" ]] || { echo "--install must run in Slurm." >&2; exit 2; }
  module purge
  module load rc-base
  module load R/4.4.1-foss-2023b

  # Also catch new running jobs submitted after the dependency snapshot.
  running="$(squeue -h -u acflint -t R,CG -o '%F %j' |
    awk -v self="${SLURM_JOB_ID}" -v pattern="${active_name_pattern}" \
      '$1 != self && $2 ~ pattern {print $1, $2}')"
  if [[ -n "${running}" ]]; then
    printf 'Another VZ/RL job is running; palette not installed:\n%s\n' "${running}" >&2
    exit 2
  fi
  [[ -f "${staged}" ]] || { echo "Missing staged palette: ${staged}" >&2; exit 2; }
  [[ ! -e "${backup}" ]] || { echo "Backup already exists: ${backup}" >&2; exit 2; }

  # Both preflights inspect paths/metadata; neither loads a large RDS object.
  Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
    --plots-only --resolution=0.8 --dry-run
  Rscript --max-ppsize=500000 scripts/xenium_vz_08_plot_samples.R --dry-run 1
  Rscript --max-ppsize=500000 -e '
    source("scripts/color_palette_vz_markers_pending.R")
    stopifnot(identical(vz_other_markers[["Granule"]],
                        c("NEUROD1", "CBLN1", "RELN")))
    cat("VZ marker panel PASS\n")
  '

  cp -p scripts/color_palette.R "${backup}"
  cp -p "${staged}" scripts/color_palette.R
  cmp "${staged}" scripts/color_palette.R
  echo "Installed VZ plot palette; backup: ${backup}"
  exit 0
fi

[[ -f "${staged}" ]] || { echo "Missing staged palette: ${staged}" >&2; exit 2; }
[[ ! -e "${backup}" ]] || { echo "Backup already exists: ${backup}" >&2; exit 2; }
if squeue -h -u acflint -n Xen_VZ_MarkerInstall -o '%i' | grep -q .; then
  echo "A VZ marker installation job is already queued; refusing duplicate submission." >&2
  exit 2
fi

# Derive array bounds from the authoritative sample manifest through the R driver.
task_map="$(Rscript --max-ppsize=500000 scripts/xenium_vz_08_plot_samples.R --list)"
task_count="$(printf '%s\n' "${task_map}" | awk -F '\t' '
  $1 ~ /^[0-9]+$/ {
    count++
    if ($1 != count || NF != 2 || $2 == "" || seen[$2]++) bad = 1
  }
  END {if (!count || bad) exit 2; print count}
')"

# %F is the array base ID, or the ordinary job ID for a non-array job.
mapfile -t active_ids < <(squeue -h -u acflint -o '%F %j' |
  awk -v pattern="${active_name_pattern}" \
    '$1 ~ /^[0-9]+$/ && $2 ~ pattern {print $1}' | sort -u)
install_dependency=()
if ((${#active_ids[@]})); then
  dependency_ids="$(printf '%s:' "${active_ids[@]}")"
  install_dependency=("--dependency=afterok:${dependency_ids%:}")
fi

submit_job() {
  local raw job_id
  raw="$(sbatch --parsable "$@")"
  job_id="${raw%%;*}"
  [[ "${job_id}" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: ${raw}" >&2; return 1; }
  printf '%s' "${job_id}"
}

printf 'Waiting for successful completion of VZ/RL job IDs: %s\n' \
  "${active_ids[*]:-none}"
printf 'VZ sample plot array: %s manifest-mapped tasks at concurrency 3.\n' \
  "${task_count}"

install_job="$(submit_job "${install_dependency[@]}" \
  --job-name=Xen_VZ_MarkerInstall \
  --output=logs/Xen_VZ_MarkerInstall_%j.out \
  --error=logs/Xen_VZ_MarkerInstall_%j.err \
  --partition=express --time=00:20:00 --cpus-per-task=1 --mem=8G \
  --wrap='bash scripts/submit_vz_marker_plot_refresh_after_jobs.sh --install')"
pooled_job="$(submit_job --dependency="afterok:${install_job}" \
  --export=ALL,VZ_SUBCLUSTER_MODE=refresh-plots,VZ_SUBCLUSTER_RESOLUTION=0.8,VZ_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm)"
sample_job="$(submit_job --dependency="afterok:${pooled_job}" \
  --array="1-${task_count}%3" \
  --export=ALL,VZ_PLOTS_MODE=samples,VZ_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_08_plot_samples.slurm)"
combine_job="$(submit_job --dependency="afterok:${sample_job}" \
  --array=1 --export=ALL,VZ_PLOTS_MODE=combine,VZ_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_08_plot_samples.slurm)"

printf 'Submitted: install=%s pooled=%s samples=%s combine=%s\n' \
  "${install_job}" "${pooled_job}" "${sample_job}" "${combine_job}"
