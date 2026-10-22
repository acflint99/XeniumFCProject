#!/usr/bin/env python3
"""Move completed VZ/RL outputs to the numbered, branch-specific layout.

Default is a read-only preflight. Run --apply only after all VZ/RL Slurm jobs
have left the queue. RDS contents are never loaded or rewritten.
"""

import argparse
import csv
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


STAGES = (
    ("01_corrected_broad", "02_corrected_broad",
     "tables/Xenium_VZRL_corrected_broad_manifest.csv"),
    ("02_merged", "nonspatial/03_merged",
     "tables/Xenium_VZRL_merge_manifest.csv"),
    ("03_processed", "nonspatial/04_processed",
     "tables/Xenium_VZRL_processed_manifest.csv"),
    ("spatial/01_merged", "spatial/03_merged",
     "tables/Xenium_VZRL_spatial_merged_manifest.csv"),
    ("spatial/02_integrated", "spatial/04_integrated",
     "tables/Xenium_VZRL_spatial_integrated_manifest.csv"),
)
ACTIVE_JOB_NAMES = {
    "Xen_VZ_Subclusters", "Xen_VZ_Map", "Xen_VZ_Plots", "Xen_VZ_Counts",
    "Xen_Combine", "Xen_CorrectBroad", "Xen_CleanMerge",
    "Xen_CombProcess", "Xen_CombPlots", "Xen_SpatialMerge",
    "Xen_SpatialProcess", "Xen_VZRL_Reports", "Xen_Unknown_VZRL",
}
RENAMED_SCRIPTS = {
    "xenium_vz_rl_00_compare_unknown_annotations.R": "xenium_vz_rl_diagnostic_compare_unknown_annotations.R",
    "run_xenium_vz_rl_00_compare_unknown_annotations.slurm": "run_xenium_vz_rl_diagnostic_compare_unknown_annotations.slurm",
    "xenium_vz_rl_01b_plot_spatial.R": "xenium_vz_rl_optional_consensus_spatial.R",
    "xenium_vz_rl_01c_correct_broad_labels.R": "xenium_vz_rl_02_correct_broad_labels.R",
    "run_xenium_vz_rl_01c_correct_broad_labels.slurm": "run_xenium_vz_rl_02_correct_broad_labels.slurm",
    "validate_xenium_vz_rl_01c_correct_broad.R": "validate_xenium_vz_rl_02_correct_broad.R",
    "xenium_vz_rl_02_merge_samples.R": "xenium_vz_rl_03_nonspatial_merge.R",
    "run_xenium_vz_rl_02_merge_samples.slurm": "run_xenium_vz_rl_03_nonspatial_merge.slurm",
    "xenium_vz_rl_03_process_and_plot.R": "xenium_vz_rl_04_05_nonspatial_process_and_plot.R",
    "run_xenium_vz_rl_03_process.slurm": "run_xenium_vz_rl_04_nonspatial_process.slurm",
    "run_xenium_vz_rl_03_plots_only.slurm": "run_xenium_vz_rl_05_nonspatial_plots.slurm",
    "xenium_vz_rl_03b_plot_cluster_counts.R": "xenium_vz_rl_optional_cluster_counts.R",
    "xenium_vz_rl_spatial_01_merge.R": "xenium_vz_rl_03_spatial_merge.R",
    "run_xenium_vz_rl_spatial_01_merge.slurm": "run_xenium_vz_rl_03_spatial_merge.slurm",
    "xenium_vz_rl_spatial_02_integrate.R": "xenium_vz_rl_04_spatial_integrate.R",
    "run_xenium_vz_rl_spatial_02_integrate.slurm": "run_xenium_vz_rl_04_spatial_integrate.slurm",
    "xenium_vz_rl_plot_rl_vz_marker_dotplot.R": "xenium_vz_rl_diagnostic_rl_vz_marker_dotplot.R",
    "run_xenium_vz_rl_rl_vz_marker_dotplots.slurm": "run_xenium_vz_rl_diagnostic_rl_vz_marker_dotplots.slurm",
}
EXPECTED_MANIFEST_PREFIXES = {
    "01_corrected_broad": {"output_path": "01_corrected_broad"},
    "02_merged": {"input_path": "01_corrected_broad", "clean_path": "02_merged"},
    "03_processed": {"input_path": "02_merged"},
    "spatial/01_merged": {"input_path": "01_corrected_broad"},
    "spatial/02_integrated": {"input_path": "spatial/01_merged"},
}


def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if not fields:
        raise ValueError(f"CSV has no header: {path}")
    return fields, rows


def check_jobs():
    if not shutil.which("squeue"):
        raise ValueError("squeue is unavailable; apply this migration on Cheaha")
    result = subprocess.run(
        ["squeue", "-h", "-u", "acflint", "-o", "%j"],
        check=True, text=True, capture_output=True,
    )
    active = sorted(set(result.stdout.splitlines()) & ACTIVE_JOB_NAMES)
    if active:
        raise ValueError("VZ/RL jobs are queued or running: " + ", ".join(active))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="perform the checked migration")
    parser.add_argument("--project-root", type=Path,
                        default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    root = args.project_root.resolve()
    out = root / "outputs/xenium/vz_rl"
    manifest = root / "config/samples.csv"
    if not out.is_dir() or not manifest.is_file():
        raise ValueError(f"Expected output tree and sample manifest under {root}")
    _, sample_rows = read_csv(manifest)
    sample_ids = [row["sample_id"] for row in sample_rows]
    if not sample_ids or len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample manifest is empty or has duplicate sample IDs")
    scripts = root / "scripts"
    missing_scripts = [new for new in RENAMED_SCRIPTS.values() if not (scripts / new).is_file()]
    if not (scripts / "R/vz_rl_batch_merge.R").is_file():
        missing_scripts.append("R/vz_rl_batch_merge.R")
    if args.apply and missing_scripts:
        raise ValueError("Sync the renamed scripts before migrating outputs; missing: " +
                         ", ".join(missing_scripts[:3]))
    old_scripts = [(scripts / old, scripts / new)
                   for old, new in RENAMED_SCRIPTS.items() if (scripts / old).exists()
                   and (scripts / new).is_file()]
    old_script_archive = scripts / "legacy_vz_rl_names_20261019"
    for old, _ in old_scripts:
        if (old_script_archive / old.name).exists():
            raise ValueError(f"Old-script archive collision: {old_script_archive / old.name}")

    moves = []
    manifests = []
    prefixes = [(str(out / old), str(out / new)) for old, new, _ in STAGES]
    prefixes.append((str(out / "03_processed/plots"),
                     str(out / "nonspatial/05_plots")))
    prefixes.sort(key=lambda item: len(item[0]), reverse=True)
    def rewrite_path(value):
        if not value:
            return value
        normalized = str(Path(value).resolve())
        for old_prefix, new_prefix in prefixes:
            if normalized == old_prefix or normalized.startswith(old_prefix + "/"):
                return new_prefix + normalized[len(old_prefix):]
        return value

    for old, new, manifest_rel in STAGES:
        source, target = out / old, out / new
        if source.exists() and target.exists():
            raise ValueError(f"Both old and new stage directories exist: {source}, {target}")
        if source.exists():
            if not source.is_dir():
                raise ValueError(f"Stage is not a directory: {source}")
            path = source / manifest_rel
            if not path.is_file():
                raise ValueError(f"Stage is incomplete; missing manifest: {path}")
            fields, rows = read_csv(path)
            if "sample_id" in fields and [r["sample_id"] for r in rows] != sample_ids:
                raise ValueError(f"Stage manifest sample order differs from config/samples.csv: {path}")
            if "sample_id" not in fields and len(rows) != 1:
                raise ValueError(f"Expected one summary row: {path}")
            required = EXPECTED_MANIFEST_PREFIXES[old]
            for row in rows:
                for key, stage_prefix in required.items():
                    value = row.get(key, "")
                    expected = str(out / stage_prefix)
                    normalized = str(Path(value).resolve()) if value else ""
                    if not (normalized == expected or normalized.startswith(expected + "/")):
                        raise ValueError(f"Unexpected {key} in {path}: {value}")
                for key in ("input_path", "output_path", "clean_path"):
                    if key in row:
                        row[key] = rewrite_path(row[key])
            if old == "01_corrected_broad":
                required_files = [source / "rds" / f"{sample}_VZ_RL_QC_CorrectedBroad.rds"
                                  for sample in sample_ids]
                required_files += [source / "plots" / f"{sample}_CorrectedBroad_GlobalSpatial_pt{point}.tif"
                                   for sample in sample_ids for point in ("0.6", "0.85")]
            elif old == "02_merged":
                required_files = [source / "rds" / f"clean_{sample}.rds" for sample in sample_ids]
                required_files.append(source / "rds/Xenium_VZRL_merged.rds")
            elif old == "03_processed":
                required_files = [source / "rds/Xenium_VZRL_processed.rds"]
            elif old == "spatial/01_merged":
                required_files = [source / "rds/Xenium_VZRL_spatial_merged.rds"]
            else:
                required_files = [source / "rds/Xenium_VZRL_spatial_integrated.rds"]
            missing = [str(file) for file in required_files if not file.is_file() or file.stat().st_size == 0]
            if missing:
                raise ValueError(f"Stage is incomplete; missing or empty outputs: {missing[:3]}")
            moves.append((source, target))
            manifests.append((target / manifest_rel, fields, rows))
        elif target.exists():
            raise ValueError(f"Target already exists; migration may be partial: {target}")

    old_plots = out / "03_processed/plots"
    new_plots = out / "nonspatial/05_plots"
    if old_plots.exists():
        if new_plots.exists():
            raise ValueError(f"Both old and new nonspatial plot directories exist: {old_plots}, {new_plots}")
        moves.append((out / "nonspatial/04_processed/plots", new_plots))

    combined = out / "01_combined_labels/plots/combined"
    flat = combined.parent
    plots = []
    if combined.exists():
        plots = sorted(combined.iterdir())
        if len(plots) != 4 * len(sample_ids) or any(p.suffix.lower() != ".tif" or not p.is_file() for p in plots):
            raise ValueError(f"Expected exactly four combined TIFFs per manifest sample in {combined}")
        for path in plots:
            if (flat / path.name).exists():
                raise ValueError(f"Combined plot collision: {flat / path.name}")

    if args.apply or shutil.which("squeue"):
        check_jobs()
    if args.apply:
        for path, _, _ in manifests:
            backup = path.with_name(path.stem + ".pre_vzrl_layout_20261019.csv")
            if backup.exists():
                raise ValueError(f"Manifest backup already exists: {backup}")
        for source, target in moves:
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
            print(f"MOVED {source} -> {target}")
        for path in plots:
            path.rename(flat / path.name)
        if plots:
            combined.rmdir()
            print(f"MOVED {len(plots)} combined plots into {flat}")
        for path, fields, rows in manifests:
            backup = path.with_name(path.stem + ".pre_vzrl_layout_20261019.csv")
            if backup.exists():
                raise ValueError(f"Manifest backup already exists: {backup}")
            shutil.copy2(path, backup)
            with tempfile.NamedTemporaryFile("w", newline="", dir=path.parent,
                                             prefix=path.name + ".", delete=False) as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
                temp = Path(handle.name)
            temp.replace(path)
            print(f"UPDATED {path}; original saved at {backup}")
        if old_scripts:
            old_script_archive.mkdir(exist_ok=True)
            for old, _ in old_scripts:
                old.rename(old_script_archive / old.name)
            print(f"ARCHIVED {len(old_scripts)} old script names in {old_script_archive}")
    else:
        for source, target in moves:
            print(f"WOULD MOVE {source} -> {target}")
        if plots:
            print(f"WOULD MOVE {len(plots)} combined plots into {flat}")
        for path, _, _ in manifests:
            print(f"WOULD CHECK/UPDATE manifest paths in {path}")
        if old_scripts:
            print(f"WOULD ARCHIVE {len(old_scripts)} old script names in {old_script_archive}")
        print("DRY-RUN PASS; no files changed")


if __name__ == "__main__":
    try:
        main()
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(f"MIGRATION STOPPED: {error}", file=sys.stderr)
        sys.exit(1)
