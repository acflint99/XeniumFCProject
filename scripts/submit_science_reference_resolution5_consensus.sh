#!/bin/bash
# Rebuild the Science reference after a reviewed update, then rerun only the
# resolution-5 all-sample weighted 2-of-3 workflow through the final consensus
# report. Existing outputs are replaced explicitly; clustering and regional
# analyses are intentionally not submitted.

set -euo pipefail

dry_run=false
downstream_only=false
for option in "$@"; do
  case "${option}" in
    --dry-run)
      dry_run=true
      ;;
    --downstream-only)
      downstream_only=true
      ;;
    *)
      echo "Usage: bash scripts/submit_science_reference_resolution5_consensus.sh [--dry-run] [--downstream-only]" >&2
      exit 2
      ;;
  esac
done

project_root="/home/acflint/R/Projects/XeniumFCProject"
cd "${project_root}" || exit 1

manifest_path="config/samples.csv"
if [[ ! -f "${manifest_path}" ]]; then
  echo "Could not find ${manifest_path} from $(pwd)." >&2
  exit 2
fi

manifest_count="$(awk 'END {print NR - 1}' "${manifest_path}")"
if [[ "${manifest_count}" != "34" ]]; then
  echo "Expected 34 manifest rows; found ${manifest_count}. Review the manifest and Slurm array bounds." >&2
  exit 2
fi

# Lightweight structural and path-only checks. These do not load large RDS
# objects and do not create or replace analysis outputs.
Rscript scripts/validate_config.R
if [[ "${downstream_only}" == "false" ]]; then
  Rscript scripts/science_01_standardize_cell_types.R --dry-run
  Rscript scripts/science_02_subset_gene_panel.R --dry-run
fi
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
  --all-samples-res5 --list >/dev/null
printf 'DRY-RUN | Science transfer task map | samples=%s\n' "${manifest_count}"

for task_id in $(seq 1 "${manifest_count}"); do
  Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
    --all-samples-res5 --dry-run Science "${task_id}"
done

# The remaining stages depend on outputs created by the Science transfer above.
# Checking them before submission would incorrectly report those future outputs
# as missing and stop the dependency chain. Their Slurm jobs run only after the
# preceding stage succeeds, so those input checks are intentionally deferred.
printf '\nDeferred preflight checks for consensus build, consensus application, and report generation until their dependency jobs run.\n'

if [[ "${dry_run}" == "true" ]]; then
  echo "Dry-run complete; no Slurm jobs were submitted."
  exit 0
fi

if [[ "${downstream_only}" == "false" ]]; then
  reference_standardize_submit=$(sbatch --parsable \
    --export=ALL,SCIENCE_REFERENCE_OVERWRITE=true \
    scripts/run_science_01_standardize_cell_types.slurm)
  reference_standardize_job=${reference_standardize_submit%%;*}

  reference_panel_submit=$(sbatch --parsable \
    --dependency="afterok:${reference_standardize_job}" \
    --export=ALL,SCIENCE_REFERENCE_OVERWRITE=true \
    scripts/run_science_02_subset_gene_panel.slurm)
  reference_panel_job=${reference_panel_submit%%;*}

  science_transfer_submit=$(sbatch --parsable \
    --dependency="afterok:${reference_panel_job}" \
    --array="1-${manifest_count}%6" \
    --job-name="Xen_res5_all_Science" \
    --export=ALL,REFERENCE=Science,RES5_ALL_ABT_OVERWRITE=true \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
else
  science_transfer_submit=$(sbatch --parsable \
    --array="1-${manifest_count}%6" \
    --job-name="Xen_res5_all_Science" \
    --export=ALL,REFERENCE=Science,RES5_ALL_ABT_OVERWRITE=true \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
fi
science_transfer_job=${science_transfer_submit%%;*}

consensus_build_submit=$(sbatch --parsable \
  --dependency="afterok:${science_transfer_job}" \
  --export=ALL,RES5_ALL_CONSENSUS_OVERWRITE=true \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_build_job=${consensus_build_submit%%;*}

consensus_apply_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_build_job}" \
  --export=ALL,RES5_ALL_W2OF3_OVERWRITE=true \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
consensus_apply_job=${consensus_apply_submit%%;*}

report_pages_submit=$(sbatch --parsable \
  --dependency="afterok:${consensus_apply_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
report_pages_job=${report_pages_submit%%;*}

report_merge_submit=$(sbatch --parsable \
  --dependency="afterok:${report_pages_job}" \
  --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=true \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
report_merge_job=${report_merge_submit%%;*}

printf '\nSubmitted Science reference rerun through the consensus report:\n'
if [[ "${downstream_only}" == "false" ]]; then
  printf 'Science standardization: %s\n' "${reference_standardize_job}"
  printf 'Science panel subset:    %s\n' "${reference_panel_job}"
else
  printf 'Science standardization: skipped (using existing outputs)\n'
  printf 'Science panel subset:    skipped (using existing outputs)\n'
fi
printf 'Science label transfer:  %s\n' "${science_transfer_job}"
printf 'Consensus table build:   %s\n' "${consensus_build_job}"
printf 'Consensus label apply:   %s\n' "${consensus_apply_job}"
printf 'Consensus report pages:  %s\n' "${report_pages_job}"
printf 'Consensus report merge:  %s\n\n' "${report_merge_job}"

squeue_job_ids="${science_transfer_job},${consensus_build_job},${consensus_apply_job},${report_pages_job},${report_merge_job}"
if [[ "${downstream_only}" == "false" ]]; then
  squeue_job_ids="${reference_standardize_job},${reference_panel_job},${squeue_job_ids}"
fi

squeue -j "${squeue_job_ids}" \
  -o "%.18i %.32j %.10T %.10M %.10l %R"
