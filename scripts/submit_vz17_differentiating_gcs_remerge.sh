#!/usr/bin/env bash
# Submit the approved VZ cluster-17 relabel and every label-dependent merge.
# Run after older VZ/RL jobs have finished and the numbered output layout is in place.
set -euo pipefail

if [[ $# -ne 1 || ( "$1" != "--dry-run" && "$1" != "--submit" ) ]]; then
  echo "Usage: bash scripts/submit_vz17_differentiating_gcs_remerge.sh --dry-run|--submit" >&2
  exit 2
fi
mode="$1"
cd /home/acflint/R/Projects/XeniumFCProject

for path in \
  outputs/xenium/vz_rl/01_combined_labels \
  outputs/xenium/vz_rl/02_corrected_broad \
  outputs/xenium/vz_rl/nonspatial/03_merged \
  outputs/xenium/vz_rl/spatial/03_merged; do
  if [[ ! -d "${path}" ]]; then
    echo "Missing numbered output stage: ${path}; finish the layout migration first." >&2
    exit 2
  fi
done

active_names='^(Xen_VZ_Subclusters|Xen_VZ_Map|Xen_VZ_Plots|Xen_VZ_Counts|Xen_Combine|Xen_CorrectBroad|Xen_CleanMerge|Xen_CombProcess|Xen_CombPlots|Xen_SpatialMerge|Xen_SpatialProcess|Xen_VZRL_Reports)$'
queue_names="$(squeue -h -u acflint -o '%j')"
active_jobs="$(printf '%s\n' "${queue_names}" | grep -E "${active_names}" || true)"
if [[ -n "${active_jobs}" ]]; then
  printf 'Wait for these older VZ/RL jobs before relabeling:\n%s\n' "${active_jobs}" >&2
  exit 2
fi

if ! grep -Fq '"17" = "Differentiating GCs"' scripts/xenium_vz_06_annotate_subclusters.R ||
   ! grep -Fq 'cluster 17 = Differentiating GCs' scripts/validate_xenium_regional_mapping.R; then
  echo "The VZ label map and mapping validator must both be updated before submission." >&2
  exit 2
fi
Rscript --max-ppsize=500000 -e '
  source("scripts/color_palette.R")
  source("scripts/R/corrected_broad_labels.R")
  stopifnot(
    "Differentiating GCs" %in% vz_subcluster_order,
    !"Mature UBCs" %in% vz_subcluster_order,
    corrected_broad_from_subcluster[["Differentiating GCs"]] == "Granule"
  )
'

# Lightweight checks only: the label preflight reads paths and small CSVs;
# the regional validator reads manifests and filesystem metadata, not RDS data.
label_preflight="$(Rscript --max-ppsize=500000 \
  scripts/xenium_vz_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run 2>&1)" || {
  printf '%s\n' "${label_preflight}" >&2
  exit 2
}
printf '%s\n' "${label_preflight}"
if [[ "${label_preflight}" != *"DRY-RUN PASS"* ]]; then
  echo "VZ label preflight did not pass." >&2
  exit 2
fi
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R RL RESIDUAL

task_map="$(Rscript --max-ppsize=500000 scripts/xenium_vz_08_plot_samples.R --list)"
task_count="$(printf '%s\n' "${task_map}" | awk -F '\t' '
  $1 ~ /^[0-9]+$/ {
    count++
    if ($1 != count || NF != 2 || $2 == "" || seen[$2]++) bad = 1
  }
  END {if (!count || bad) exit 2; print count}
')"
printf 'VZ sample plot array: %s manifest-mapped tasks, concurrency 3.\n' "${task_count}"
printf 'Will replace existing VZ annotation/map/plots, combined labels/reports, corrected broad labels, both merges, and both processed outputs.\n'
if [[ "${mode}" == "--dry-run" ]]; then
  echo "DRY-RUN PASS: no jobs submitted. Run again with --submit when ready."
  exit 0
fi

submit_job() {
  local raw job_id
  raw="$(sbatch --parsable "$@")"
  job_id="${raw%%;*}"
  if [[ ! "${job_id}" =~ ^[0-9]+$ ]]; then
    echo "Unexpected sbatch response: ${raw}" >&2
    return 1
  fi
  printf '%s' "${job_id}"
}

vz_apply="$(submit_job \
  --export=ALL,VZ_SUBCLUSTER_MODE=apply,VZ_SUBCLUSTER_RESOLUTION=0.8,VZ_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm)"
vz_map="$(submit_job --dependency="afterok:${vz_apply}" \
  --export=ALL,VZ_MAPPING_WORKERS=1,VZ_MAPPING_OVERWRITE=true \
  scripts/run_xenium_vz_07_map_to_samples.slurm)"
vz_plots="$(submit_job --dependency="afterok:${vz_map}" \
  --array="1-${task_count}%3" --export=ALL,VZ_PLOTS_MODE=samples,VZ_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_08_plot_samples.slurm)"
vz_plot_combine="$(submit_job --dependency="afterok:${vz_plots}" \
  --array=1 --export=ALL,VZ_PLOTS_MODE=combine,VZ_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_08_plot_samples.slurm)"
vz_counts="$(submit_job --dependency="afterok:${vz_apply}" \
  --export=ALL,VZ_COUNT_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_09_plot_counts.slurm)"

combine="$(submit_job --dependency="afterok:${vz_map}" \
  --export=ALL,COMBINE_OVERWRITE=true \
  scripts/run_xenium_vz_rl_01_combine_labels.slurm)"
corrected="$(submit_job --dependency="afterok:${combine}" \
  --export=ALL,CORRECTED_BROAD_OVERWRITE=true \
  scripts/run_xenium_vz_rl_02_correct_broad_labels.slurm)"
# Both-label reports read the corrected-broad TIFFs as well as the combined
# TIFFs. Waiting for corrected also waits transitively for combine.
combined_reports="$(submit_job --dependency="afterok:${corrected}" \
  --export=ALL,VZRL_REPORT_LABELS=both,VZRL_REPORT_POINT_SIZE=0.6,VZRL_REPORT_OVERWRITE=true \
  scripts/run_xenium_vz_rl_reports.slurm)"
nonspatial_merge="$(submit_job --dependency="afterok:${corrected}" \
  --export=ALL,CLEAN_MERGE_OVERWRITE=true \
  scripts/run_xenium_vz_rl_03_nonspatial_merge.slurm)"
spatial_merge="$(submit_job --dependency="afterok:${corrected}" \
  --export=ALL,SPATIAL_MERGE_BATCH_SIZE=6,SPATIAL_MERGE_OVERWRITE=true \
  scripts/run_xenium_vz_rl_03_spatial_merge.slurm)"
nonspatial_process="$(submit_job --dependency="afterok:${nonspatial_merge}" \
  --export=ALL,COMBINED_PROCESS_OVERWRITE=true \
  scripts/run_xenium_vz_rl_04_nonspatial_process.slurm)"
nonspatial_plots="$(submit_job --dependency="afterok:${nonspatial_process}" \
  --export=ALL,COMBINED_PLOTS_OVERWRITE=true \
  scripts/run_xenium_vz_rl_05_nonspatial_plots.slurm)"
spatial_process="$(submit_job --dependency="afterok:${spatial_merge}" \
  --export=ALL,SPATIAL_PROCESS_OVERWRITE=true \
  scripts/run_xenium_vz_rl_04_spatial_integrate.slurm)"

printf 'Submitted: VZ apply=%s map=%s plots=%s plot-combine=%s counts=%s\n' \
  "${vz_apply}" "${vz_map}" "${vz_plots}" "${vz_plot_combine}" "${vz_counts}"
printf 'Submitted: combined=%s reports=%s corrected-broad=%s\n' \
  "${combine}" "${combined_reports}" "${corrected}"
printf 'Submitted: nonspatial merge=%s process=%s plots=%s; spatial merge=%s integrate=%s\n' \
  "${nonspatial_merge}" "${nonspatial_process}" "${nonspatial_plots}" \
  "${spatial_merge}" "${spatial_process}"
