#!/bin/bash
# Rebuild the Sepp Xenium-panel reference after a reviewed update, then rerun
# only the resolution-5 all-sample weighted 2-of-3 workflow through the final
# consensus reports. Existing downstream outputs are replaced explicitly;
# clustering and unchanged reference transfers are intentionally not submitted.

set -Eeuo pipefail

dry_run=false
downstream_only=false
current_stage="argument validation"
jobs_submitted=false

handle_failure() {
  local exit_code="$1"
  local line_number="${2:-unknown}"
  local result_label="PREFLIGHT FAILED"
  trap - ERR

  if [[ "${dry_run}" == "true" ]]; then
    result_label="DRY-RUN FAILED"
  elif [[ "${jobs_submitted}" == "true" ]]; then
    result_label="SUBMISSION FAILED"
  fi

  printf '\n========================================\n' >&2
  printf '%s\n' "${result_label}" >&2
  printf 'Failed stage: %s\n' "${current_stage}" >&2
  printf 'Exit code: %s (script line %s)\n' "${exit_code}" "${line_number}" >&2
  if [[ "${jobs_submitted}" == "true" ]]; then
    printf 'Some Slurm jobs were already submitted. Inspect them with squeue -u acflint.\n' >&2
  else
    printf 'Slurm jobs submitted: 0\n' >&2
  fi
  printf 'Review the error immediately above this summary before retrying.\n' >&2
  printf '========================================\n' >&2
  exit "${exit_code}"
}

trap 'handle_failure $? "$LINENO"' ERR

for option in "$@"; do
  case "${option}" in
    --dry-run)
      dry_run=true
      ;;
    --downstream-only)
      downstream_only=true
      ;;
    *)
      echo "Usage: bash scripts/submit_sepp_reference_resolution5_consensus.sh [--dry-run] [--downstream-only]" >&2
      handle_failure 2 "${LINENO}"
      ;;
  esac
done

current_stage="project directory"
project_root="/home/acflint/R/Projects/XeniumFCProject"
if ! cd "${project_root}"; then
  echo "Could not enter project directory: ${project_root}" >&2
  handle_failure 1 "${LINENO}"
fi

# Avoid printing the same renv synchronization notice once for every one of the
# 34 path checks. This does not restore or change the environment.
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

current_stage="34-sample manifest validation"
manifest_path="config/samples.csv"
if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  handle_failure 2 "${LINENO}"
fi

manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 manifest rows; found ${manifest_count}. Review the manifest and Slurm array bounds." >&2
  handle_failure 2 "${LINENO}"
fi

# Lightweight structural and path-only checks. These do not load large RDS
# objects and do not create or replace analysis outputs. sepp_02 does not have
# a dry-run mode, so its two required inputs are checked directly.
current_stage="configuration validation"
Rscript scripts/validate_config.R
if [[ "${downstream_only}" == "false" ]]; then
  current_stage="Sepp panel-subset input validation"
  for required_input in \
    outputs/references/sepp/rds/Sepp_newClusters_newUMAPv1.rds \
    inputs/xenium_5k_genes.rds; do
    if [[ ! -f "${required_input}" ]]; then
      echo "Required Sepp panel-subset input not found: ${required_input}" >&2
      handle_failure 2 "${LINENO}"
    fi
  done
fi

current_stage="Sepp transfer task-map validation"
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
  --all-samples-res5 --list >/dev/null
printf 'DRY-RUN PASS | Sepp transfer task map | samples=%s | concurrency=6\n' "${manifest_count}"

for task_id in $(seq 1 "${manifest_count}"); do
  current_stage="Sepp transfer dry-run task ${task_id}/${manifest_count}"
  Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
    --all-samples-res5 --dry-run Sepp "${task_id}"
done

# The remaining checks are deferred because their inputs will be replaced by
# dependency jobs. Each production launcher performs its own path-only dry-run.
printf '\nDeferred preflight checks for consensus build, consensus application, and report generation until their dependency jobs run.\n'

if [[ "${dry_run}" == "true" ]]; then
  printf '\n========================================\n'
  printf 'DRY-RUN SUCCESS\n'
  printf 'Passed: configuration, 34-sample manifest, task map, and %s/%s Sepp transfer input checks.\n' \
    "${manifest_count}" "${manifest_count}"
  printf 'Informational: the outputs-existing field only counts files already present; this is not a failure.\n'
  printf 'Deferred: consensus build, consensus application, and report checks run inside dependency jobs after new inputs exist.\n'
  printf 'Slurm jobs submitted: 0\n'
  printf 'Environment: repeated renv synchronization notices were suppressed; run renv::status() before production.\n'
  printf '========================================\n'
  exit 0
fi

if [[ "${downstream_only}" == "false" ]]; then
  current_stage="submit Sepp panel-subset job"
  printf 'Submitting sepp_02; it writes the reviewed stable Sepp RDS and plot paths.\n'
  reference_panel_submit=$(sbatch --parsable \
    scripts/run_sepp_02_subset_gene_panel.slurm)
  reference_panel_job=${reference_panel_submit%%;*}
  jobs_submitted=true

  current_stage="submit Sepp transfer array"
  sepp_transfer_submit=$(sbatch --parsable \
    --dependency="afterok:${reference_panel_job}" \
    --array="1-${manifest_count}%6" \
    --job-name="Xen_res5_all_Sepp" \
    --export=ALL,REFERENCE=Sepp,RES5_ALL_ABT_OVERWRITE=true \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
else
  current_stage="submit Sepp transfer array"
  sepp_transfer_submit=$(sbatch --parsable \
    --array="1-${manifest_count}%6" \
    --job-name="Xen_res5_all_Sepp" \
    --export=ALL,REFERENCE=Sepp,RES5_ALL_ABT_OVERWRITE=true \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
fi
sepp_transfer_job=${sepp_transfer_submit%%;*}
jobs_submitted=true

current_stage="submit consensus-table build job"
consensus_build_submit=$(sbatch --parsable \
  --dependency="afterok:${sepp_transfer_job}" \
  --export=ALL,RES5_ALL_CONSENSUS_OVERWRITE=true \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_build_job=${consensus_build_submit%%;*}

current_stage="submit consensus-application array"
consensus_apply_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_build_job}" \
  --export=ALL,RES5_ALL_W2OF3_OVERWRITE=true \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_apply_job=${consensus_apply_submit%%;*}

current_stage="submit consensus-report page array"
report_pages_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_apply_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
report_pages_job=${report_pages_submit%%;*}

current_stage="submit consensus-report merge job"
report_merge_submit=$(sbatch --parsable \
  --dependency="afterok:${report_pages_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
report_merge_job=${report_merge_submit%%;*}

current_stage="submission summary"
printf '\nSubmitted Sepp reference rerun through the consensus report:\n'
if [[ "${downstream_only}" == "false" ]]; then
  printf 'Sepp panel subset:      %s\n' "${reference_panel_job}"
else
  printf 'Sepp panel subset:      skipped (using existing output)\n'
fi
printf 'Sepp label transfer:    %s (34 tasks, concurrency 6)\n' "${sepp_transfer_job}"
printf 'Consensus table build: %s\n' "${consensus_build_job}"
printf 'Consensus label apply: %s\n' "${consensus_apply_job}"
printf 'Consensus report pages:%s (all-cell and known-only)\n' "${report_pages_job}"
printf 'Consensus report merge:%s (two PDFs)\n\n' "${report_merge_job}"

squeue_job_ids="${sepp_transfer_job},${consensus_build_job},${consensus_apply_job},${report_pages_job},${report_merge_job}"
if [[ "${downstream_only}" == "false" ]]; then
  squeue_job_ids="${reference_panel_job},${squeue_job_ids}"
fi

current_stage="Slurm queue display"
squeue -j "${squeue_job_ids}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"

printf '\n========================================\n'
printf 'SUBMISSION SUCCESS\n'
printf 'All dependency jobs were submitted. The workflow is now managed by Slurm.\n'
printf '========================================\n'
