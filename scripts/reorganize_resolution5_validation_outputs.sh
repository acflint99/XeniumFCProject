#!/bin/bash
# Move completed resolution-selection and method-development evidence out of
# production preprocessing/annotation trees into outputs/validation/.

set -euo pipefail

mode=${1:---dry-run}
if [[ "${mode}" != "--dry-run" && "${mode}" != "--move" ]]; then
  echo "Usage: bash scripts/reorganize_resolution5_validation_outputs.sh [--dry-run|--move]" >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "${script_dir}/.." && pwd)
output_root="${project_root}/outputs"
validation_root="${output_root}/validation"
manifest_path="${validation_root}/resolution5_output_reorganization_2026-09-29.tsv"

[[ -d "${output_root}" ]] || {
  echo "Output root not found: ${output_root}" >&2
  exit 1
}

sources=(
  "outputs/xenium/preprocess/03f_resolution_consistency"
  "outputs/xenium/annotation/resolution5_all_samples/04_transfer_threshold_calibration_audit"
  "outputs/xenium/annotation/resolution5_all_samples/k_anchor_sensitivity"
  "outputs/xenium/annotation/resolution5_all_samples/reference_sampling_abstention_sensitivity"
  "outputs/xenium/annotation/resolution5_all_samples/transfer_diagnostics_k10"
)
destinations=(
  "outputs/validation/resolution_selection"
  "outputs/validation/resolution5_method_development/04_transfer_threshold_calibration_audit"
  "outputs/validation/resolution5_method_development/k_anchor_sensitivity"
  "outputs/validation/resolution5_method_development/reference_sampling_abstention_sensitivity"
  "outputs/validation/resolution5_method_development/transfer_diagnostics_k10"
)

for index in "${!sources[@]}"; do
  source_path="${project_root}/${sources[$index]}"
  destination_path="${project_root}/${destinations[$index]}"
  if [[ -e "${source_path}" && -e "${destination_path}" ]]; then
    echo "STOP: source and destination both exist:" >&2
    echo "  ${source_path}" >&2
    echo "  ${destination_path}" >&2
    exit 1
  fi
done

printf 'Mode: %s\n' "${mode}"
for index in "${!sources[@]}"; do
  source_path="${project_root}/${sources[$index]}"
  destination_path="${project_root}/${destinations[$index]}"
  if [[ -d "${source_path}" ]]; then
    size=$(du -sh "${source_path}" | cut -f1)
    printf '%-8s  %s -> %s\n' "${size}" "${sources[$index]}" "${destinations[$index]}"
  elif [[ -d "${destination_path}" ]]; then
    printf 'MOVED     %s\n' "${destinations[$index]}"
  else
    printf 'MISSING   %s\n' "${sources[$index]}"
  fi
done

if [[ "${mode}" == "--dry-run" ]]; then
  echo "Dry run only. Re-run with --move after reviewing these exact paths."
  exit 0
fi

mkdir -p "${validation_root}"
printf 'moved_at\tsource\tdestination\n' > "${manifest_path}"

for index in "${!sources[@]}"; do
  source_path="${project_root}/${sources[$index]}"
  destination_path="${project_root}/${destinations[$index]}"
  [[ -d "${source_path}" ]] || continue
  mkdir -p "$(dirname "${destination_path}")"
  mv -- "${source_path}" "${destination_path}"
  printf '%s\t%s\t%s\n' \
    "$(date --iso-8601=seconds)" \
    "${sources[$index]}" \
    "${destinations[$index]}" >> "${manifest_path}"
done

echo "Reorganization complete. Manifest: ${manifest_path}"
