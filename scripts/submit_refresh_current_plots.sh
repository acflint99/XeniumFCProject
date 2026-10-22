#!/bin/bash
# Defaults to a path-only preflight. --submit submits plotting jobs, not analysis.
set -euo pipefail
cd /home/acflint/R/Projects/XeniumFCProject
mode=dry-run
tag=palette_20261019
selected_groups=references,samples,pooled,combined
for option in "$@"; do
  case "$option" in
    --dry-run) mode=dry-run ;;
    --submit) mode=submit ;;
    --tag=*) tag="${option#--tag=}" ;;
    --groups=*) selected_groups="${option#--groups=}" ;;
    *) echo "Use [--dry-run|--submit] [--tag=NEW_SNAPSHOT_NAME] [--groups=references,samples,pooled,combined]." >&2; exit 2 ;;
  esac
done
if [[ ! "$tag" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]]; then echo "Invalid snapshot tag." >&2; exit 2; fi
if [[ ! "$selected_groups" =~ ^(references|samples|pooled|combined)(,(references|samples|pooled|combined))*$ ]]; then
  echo "Invalid --groups list." >&2; exit 2
fi
IFS=',' read -r -a groups <<< "$selected_groups"
seen_groups=,
for group in "${groups[@]}"; do
  if [[ "$seen_groups" == *",${group},"* ]]; then echo "Repeated group: ${group}" >&2; exit 2; fi
  seen_groups+="${group},"
done

# Check all required inputs before submitting any group. No RDS is loaded.
preflight_failed=false
for group in "${groups[@]}"; do
  if ! Rscript --max-ppsize=500000 scripts/refresh_current_plots.R "--group=${group}" "--tag=${tag}" --dry-run; then
    preflight_failed=true
  fi
done
if [[ "$preflight_failed" == true ]]; then
  echo "Preflight failed; no jobs were submitted. Review the failures above." >&2
  exit 2
fi
for group in "${groups[@]}"; do
  task_map="$(Rscript --max-ppsize=500000 scripts/refresh_current_plots.R "--group=${group}" "--tag=${tag}" --list)"
  # renv may print a startup status line before the TSV header. Count only
  # numeric task IDs, not the header or any startup messages, and require a
  # contiguous 1..N mapping before using it as a Slurm array bound.
  if ! task_count="$(awk -F '\t' '
    $1 ~ /^[1-9][0-9]*$/ && NF >= 2 {
      if (!seen[$1]++) count++
      if ($1 > maximum) maximum = $1
    }
    END {
      if (count < 1 || maximum != count) exit 2
      print count
    }' <<< "$task_map")"; then
    echo "Invalid/noncontiguous --list task mapping for ${group}." >&2
    exit 2
  fi
  if [[ ! "$task_count" =~ ^[1-9][0-9]*$ ]]; then echo "Invalid --list task count for ${group}." >&2; exit 2; fi
  case "$group" in
    references) concurrency=1; memory=64G; limit=01:00:00; partition=express ;;
    samples) concurrency=3; memory=64G; limit=02:00:00; partition=express ;;
    pooled) concurrency=1; memory=256G; limit=12:00:00; partition=short ;;
    combined) concurrency=1; memory=600G; limit=24:00:00; partition=medium ;;
  esac
  echo "${group}: ${task_count} tasks, concurrency ${concurrency}, ${memory}, limit ${limit}, ${partition}"
  if [[ "$mode" == submit ]]; then
    sbatch --parsable --array="1-${task_count}%${concurrency}" --mem="$memory" \
      --time="$limit" --partition="$partition" \
      --export="ALL,PLOT_REFRESH_GROUP=${group},PLOT_REFRESH_TAG=${tag}" \
      scripts/run_refresh_current_plots.slurm
  fi
done
