#!/bin/bash
# Remove the redundant resolution5_all_samples wrapper now that resolution 5
# is the only active whole-tissue annotation workflow.

set -euo pipefail

mode=${1:---dry-run}
if [[ "${mode}" != "--dry-run" && "${mode}" != "--move" ]]; then
  echo "Usage: bash scripts/flatten_resolution5_annotation.sh [--dry-run|--move]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
annotation_root="${project_root}/outputs/xenium/annotation"
source_root="${annotation_root}/resolution5_all_samples"
manifest_path="${project_root}/outputs/validation/resolution5_annotation_flatten_2026-09-29.tsv"
sample_manifest="${project_root}/config/samples.csv"

stages=(
  "01_label_transfer"
  "02_consensus_weighted_2of3"
  "03_consensus_labels_weighted_2of3"
)

[[ -f "${sample_manifest}" ]] || {
  echo "Sample manifest not found: ${sample_manifest}" >&2
  exit 1
}
expected_samples=$(awk -F, '
  NR > 1 {
    gsub(/\r/, "", $1)
    gsub(/^"|"$/, "", $1)
    if (length($1)) count++
  }
  END { print count + 0 }
' "${sample_manifest}")
[[ "${expected_samples}" -eq 34 ]] || {
  echo "Expected 34 biological samples; found ${expected_samples}." >&2
  exit 1
}

if [[ ! -d "${source_root}" ]]; then
  all_flattened=true
  for stage in "${stages[@]}"; do
    [[ -d "${annotation_root}/${stage}" ]] || all_flattened=false
  done
  if [[ "${all_flattened}" == true ]]; then
    echo "Already flattened: ${annotation_root}"
    exit 0
  fi
  echo "Source annotation wrapper not found: ${source_root}" >&2
  exit 1
fi

unexpected=()
while IFS= read -r entry; do
  name=$(basename "${entry}")
  allowed=false
  for stage in "${stages[@]}"; do
    [[ "${name}" == "${stage}" ]] && allowed=true
  done
  [[ "${allowed}" == true ]] || unexpected+=("${entry}")
done < <(find "${source_root}" -mindepth 1 -maxdepth 1 -print | LC_ALL=C sort)

if [[ "${#unexpected[@]}" -gt 0 ]]; then
  echo "STOP: move method-development outputs first; unexpected wrapper contents:" >&2
  printf '  %s\n' "${unexpected[@]}" >&2
  exit 1
fi

for stage in "${stages[@]}"; do
  source_stage="${source_root}/${stage}"
  destination_stage="${annotation_root}/${stage}"
  [[ -d "${source_stage}" ]] || {
    echo "Missing production stage: ${source_stage}" >&2
    exit 1
  }
  [[ ! -e "${destination_stage}" ]] || {
    echo "STOP: destination already exists: ${destination_stage}" >&2
    exit 1
  }
done

count_files() {
  local directory=$1
  local pattern=$2
  find "${directory}" -maxdepth 1 -type f -name "${pattern}" -print | wc -l
}

for reference in Aldinger Sepp Science; do
  reference_key=$(printf '%s' "${reference}" | tr '[:upper:]' '[:lower:]')
  count=$(count_files \
    "${source_root}/01_label_transfer/${reference_key}/rds" \
    "*_${reference}_annotated.rds")
  [[ "${count}" -eq "${expected_samples}" ]] || {
    echo "Incomplete ${reference} transfer RDS set: ${count}/${expected_samples}" >&2
    exit 1
  }
done

consensus_tables=$(count_files \
  "${source_root}/02_consensus_weighted_2of3/tables" \
  '*_comparison_merged.csv')
consensus_rds=$(count_files \
  "${source_root}/03_consensus_labels_weighted_2of3/rds" \
  '*_Consensus_annotated.rds')
[[ "${consensus_tables}" -eq "${expected_samples}" ]] || {
  echo "Incomplete consensus table set: ${consensus_tables}/${expected_samples}" >&2
  exit 1
}
[[ "${consensus_rds}" -eq "${expected_samples}" ]] || {
  echo "Incomplete consensus RDS set: ${consensus_rds}/${expected_samples}" >&2
  exit 1
}

printf 'Mode: %s\n' "${mode}"
for stage in "${stages[@]}"; do
  du -sh "${source_root}/${stage}"
  printf '  %s -> %s\n' \
    "outputs/xenium/annotation/resolution5_all_samples/${stage}" \
    "outputs/xenium/annotation/${stage}"
done

if [[ "${mode}" == "--dry-run" ]]; then
  echo "Dry run only. Re-run with --move after reviewing these exact paths."
  exit 0
fi

mkdir -p "$(dirname "${manifest_path}")"
printf 'moved_at\tsource\tdestination\n' > "${manifest_path}"
for stage in "${stages[@]}"; do
  mv -- "${source_root}/${stage}" "${annotation_root}/${stage}"
  printf '%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" \
    "outputs/xenium/annotation/resolution5_all_samples/${stage}" \
    "outputs/xenium/annotation/${stage}" >> "${manifest_path}"
done
rmdir -- "${source_root}"

echo "Flattening complete. Manifest: ${manifest_path}"
