#!/bin/bash
# Rebuild both PCW-ordered consensus PDFs and the Aldinger TIFF-based PDF.
# Existing PDFs/cached pages are archived before refresh; RDS files are read only.
set -euo pipefail
if [[ $# -ne 1 || ( $1 != --dry-run && $1 != --submit ) ]]; then
  echo "Usage: bash scripts/refresh_xenium_consensus_reports_20261015.sh [--dry-run|--submit]" >&2
  exit 2
fi
mode=$1
cd /home/acflint/R/Projects/XeniumFCProject
cpu_vendor=$(lscpu | awk -F: '/Vendor ID/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}')
[[ $cpu_vendor == GenuineIntel ]] || { echo "Use an Intel session." >&2; exit 2; }
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

# Job-name check prevents moving caches used by an active report build.
active_report_jobs=$(squeue -h -u acflint -o '%j' | awk '$0 ~ /^Xen_res5_w2of3_(pages|merge)$/ || $0 == "Xen_res5_Ald_pdf"')
[[ -z $active_report_jobs ]] || {
  echo "Existing consensus/Aldinger report jobs are active; let them finish before refreshing." >&2
  exit 2
}
page_map=$(Rscript --max-ppsize=500000 scripts/xenium_annotate_03d_plot_report.R --all-samples-res5 --weighted-2of3 --list)
printf '%s\n' "$page_map"
page_ids=$(printf '%s\n' "$page_map" | awk '$1 ~ /^[0-9]+$/ {print $1}' | paste -sd, -)
[[ -n $page_ids ]] || { echo "Report task map is empty." >&2; exit 2; }
for driver in consensus aldinger; do
  if [[ $driver == consensus ]]; then
    result=$(Rscript --max-ppsize=500000 scripts/xenium_annotate_03d_plot_report.R --all-samples-res5 --weighted-2of3 --dry-run)
  else
    result=$(Rscript --max-ppsize=500000 scripts/xenium_annotate_01b_build_aldinger_resolution5_plot_report.R --dry-run)
  fi
  printf '%s\n' "$result"
  [[ $result == *"DRY-RUN PASS"* && $result != *"DRY-RUN FAIL"* ]] || exit 2
done
Rscript --max-ppsize=500000 scripts/audit_xenium_merged_outputs.R --dry-run
[[ $mode == --submit ]] || exit 0

mkdir -p outputs/report_backups logs
backup_dir=$(mktemp -d outputs/report_backups/consensus_41_samples_20261015_XXXXXX)
consensus_plots=outputs/xenium/annotation/03_consensus_labels_weighted_2of3/plots
aldinger_plots=outputs/xenium/annotation/01_label_transfer/aldinger/plots
for file in \
  "$consensus_plots/Xenium_Consensus_Res5_Weighted2of3_All_Samples_Report.pdf" \
  "$consensus_plots/Xenium_Consensus_Res5_Weighted2of3_All_Samples_KnownOnly_Report.pdf" \
  "$aldinger_plots/Xenium_Aldinger_Res5_All_Samples_Plot_Report.pdf"; do
  if [[ -f $file ]]; then mv "$file" "$backup_dir/"; fi
done
# New PCW entries change page numbers for older samples too. Starting fresh
# avoids keeping old page_XX_sample.png files that fail the exact page-set check.
for directory in report_pages_res5_weighted2of3 report_pages_res5_weighted2of3_known_only; do
  if [[ -d $consensus_plots/$directory ]]; then
    mv "$consensus_plots/$directory" "$backup_dir/"
  fi
done
printf '%s\n' "$page_map" > "$backup_dir/new_page_map.tsv"
echo "Archived previous reports/pages in $backup_dir"

submit_job() {
  local id
  id=$(sbatch --parsable "$@")
  id=${id%%;*}
  [[ $id =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $id" >&2; exit 2; }
  printf '%s\n' "$id"
}
pages=$(submit_job --array="${page_ids}%6" --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=false scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
echo "Both consensus variants, page rendering: $pages"
merge=$(submit_job --dependency="afterok:$pages" --export=ALL,RES5_W2OF3_REPORT_OVERWRITE=false scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
echo "Both 41-page consensus PDFs, merge: $merge"
aldinger=$(submit_job --export=ALL,ALDINGER_RES5_REPORT_OVERWRITE=false scripts/run_xenium_annotate_01b_build_aldinger_resolution5_plot_report.slurm)
echo "41-page Aldinger report: $aldinger"
printf 'pages=%s\nmerge=%s\naldinger=%s\n' "$pages" "$merge" "$aldinger" > "$backup_dir/job_ids.txt"
echo "Submitted. QC/broad/routing/validation reports are separate; see MERGED_OUTPUTS_20261015.md."
