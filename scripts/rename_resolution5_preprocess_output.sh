#!/bin/bash
# Rename the completed resolution-5 preprocessing directory without rewriting data.

set -euo pipefail

mode="${1:---dry-run}"
if [[ "${mode}" != "--dry-run" && "${mode}" != "--move" ]]; then
  echo "Usage: bash scripts/rename_resolution5_preprocess_output.sh [--dry-run|--move]" >&2
  exit 2
fi

project_root="/home/acflint/R/Projects/XeniumFCProject"
manifest_path="${project_root}/config/samples.csv"
preprocess_root="${project_root}/outputs/xenium/preprocess"
source_path="${preprocess_root}/03i_resolution5_all_samples"
target_path="${preprocess_root}/04_resolution5_clustered"

[[ -f "${manifest_path}" ]] || {
  echo "STOP: sample manifest is missing: ${manifest_path}" >&2
  exit 1
}

expected_samples=$(awk 'NR > 1 && length($1) {count++} END {print count + 0}' \
  FS=, "${manifest_path}")
[[ "${expected_samples}" -eq 34 ]] || {
  echo "STOP: expected 34 biological samples; found ${expected_samples}." >&2
  exit 1
}

validate_output_set() {
  local stage_path=$1
  local missing=0
  local sample_id

  [[ -d "${stage_path}/rds" ]] || {
    echo "STOP: missing RDS directory: ${stage_path}/rds" >&2
    return 1
  }

  while IFS=, read -r sample_id _; do
    sample_id=${sample_id//$'\r'/}
    sample_id=${sample_id#\"}
    sample_id=${sample_id%\"}
    [[ -n "${sample_id}" ]] || continue
    if [[ ! -s "${stage_path}/rds/${sample_id}_whole_tissue_Res5.0.rds" ]]; then
      echo "MISSING: ${stage_path}/rds/${sample_id}_whole_tissue_Res5.0.rds" >&2
      missing=$((missing + 1))
    fi
  done < <(tail -n +2 "${manifest_path}")

  [[ "${missing}" -eq 0 ]] || {
    echo "STOP: ${missing} required resolution-5 RDS file(s) are missing." >&2
    return 1
  }

  local observed
  observed=$(find "${stage_path}/rds" -maxdepth 1 -type f \
    -name '*_whole_tissue_Res5.0.rds' -print | wc -l)
  [[ "${observed}" -eq "${expected_samples}" ]] || {
    echo "STOP: expected ${expected_samples} resolution-5 RDS files; found ${observed}." >&2
    return 1
  }

  printf 'Validated resolution-5 clustering RDS: %s/%s\n' \
    "${observed}" "${expected_samples}"
}

if [[ -d "${target_path}" && ! -e "${source_path}" ]]; then
  validate_output_set "${target_path}"
  echo "Already renamed: ${target_path}"
  exit 0
fi

[[ -d "${source_path}" ]] || {
  echo "STOP: source directory is missing: ${source_path}" >&2
  exit 1
}
[[ ! -e "${target_path}" ]] || {
  echo "STOP: target already exists; refusing to merge or overwrite: ${target_path}" >&2
  exit 1
}

validate_output_set "${source_path}"

echo "Source: ${source_path}"
echo "Target: ${target_path}"
if [[ "${mode}" == "--dry-run" ]]; then
  echo "Dry run only. Re-run with --move to rename this exact directory."
  exit 0
fi

mv -- "${source_path}" "${target_path}"
validate_output_set "${target_path}"
echo "Renamed resolution-5 preprocessing directory. No analysis files were rewritten."
