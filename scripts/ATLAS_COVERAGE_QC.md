# Atlas coverage, QC, and figure panels

This reporting layer reads existing production objects and raw QC tables. It does
not filter cells, change annotations/parameters, alter task order, or write RDS
objects. Every output is protected: use a **new report directory** to regenerate.
No packages are installed. Use the established project renv and R/4.4.1 on Cheaha.

Sync the two atlas drivers, `R/atlas_qc.R`, the two `test_xenium_atlas_qc*.R`
files, the atlas Slurm launcher/submission helper, this guide, and the three
`config/atlas_*.csv` companions into the corresponding locations under HPC
`scripts/`. Existing `R/config.R` is reused. Do not replace HPC project paths or
production manifests with local copies merely to use this reporting layer.

## Metadata and scientific scope

The authoritative `config/samples.csv` remains unchanged (41 entries). Edit these
companion tables under `scripts/config/`:

- `atlas_sample_metadata.csv`: explicit cerebellum/donor grouping, PCW, section,
  physical slide/run, panel, software and segmentation provenance. Shared-prefix
  sections are grouped according to the user's confirmation. The FB124/GZFB124
  and FB328/GZFB328 merges are **provisional user-authorized assumptions**, not
  confirmed biological replication. This gives 12 working cerebellum groups;
  ten groups have confirmed grouping and two remain provisional. Do not equate
  the 41 entries or 32 raw-input directories with donor counts.
- PCW is initialized from the available 2026-10-14 metadata inspection snapshot;
  extraction/compilation compares it with final-object PCW and fails on mismatch.
  Replace or verify these ages against current specimen records before publication.
- `section_order_lateral_to_medial` is ordinal: section 1 lateral, increasing
  numbers toward medial. It is not physical distance. Unnumbered sections have
  no assumed position. The supplied `GZFB12_5.tif` supports a sagittal/parasagittal
  interpretation for that series only, recorded as inferred in
  `section_plane_status`. Confirm the plane rather than propagating it to all samples.
- `atlas_anatomical_coverage.csv`: one row per sample/compartment. Status is
  `present`, `absent`, `uncertain`, or `not_assessed`. All analyzed samples are
  user-confirmed cerebellar crops; other compartments remain unassessed. Retain
  a row for an absent compartment. Do not delete missing sections from the manifest.
  Measured `area_mm2` requires `area_method` and `source`; use reviewed calibrated
  tissue masks/polygons, not summed cell areas or bounding boxes. Areas can overlap
  across compartments and must not be summed indiscriminately.
- `atlas_segmentation_review.csv`: one row per reviewed field, documenting its
  compartment, field-selection basis, image/evidence path, reviewer/date,
  `pending`/`reviewed` status, and findings. Include fields across age, density,
  VZ/RL, large neurons, tissue edges and artifacts. One recorded field is a
  completeness indicator, not certification of adequate segmentation validation.

Raw sections include brainstem and surrounding tissue. All analysis summaries
start with the **cerebellar crop**. Excluded surrounding tissue is not a QC failure.
Raw instrument metrics are acquisition/input-directory metrics and may include
extracerebellar tissue. They are never relabeled as cerebellum-specific metrics.

## Lightweight coverage report (no RDS load)

Run from `/home/acflint/R/Projects/XeniumFCProject`, after loading the established
R module/environment. Compilation and plotting are each expected to take 1–5
minutes, low confidence from small-table processing and figure export. Preflights
take under one minute and inspect paths only. A draft explicitly preserves gaps:

```bash
Rscript --max-ppsize=500000 scripts/xenium_atlas_qc.R \
  --compile --allow-incomplete --dry-run --output-dir outputs/xenium/atlas_qc_review1
Rscript --max-ppsize=500000 scripts/xenium_atlas_qc.R \
  --compile --allow-incomplete --output-dir outputs/xenium/atlas_qc_review1
Rscript --max-ppsize=500000 scripts/xenium_atlas_plot_qc.R \
  --dry-run --output-dir outputs/xenium/atlas_qc_review1
Rscript --max-ppsize=500000 scripts/xenium_atlas_plot_qc.R \
  --output-dir outputs/xenium/atlas_qc_review1
```

Draft compilation permits unavailable extraction bundles and metadata gaps; it
does not permit duplicate/foreign sample keys, conflicting donor ages, corrupted
count reconciliation, or changed extraction sources. Strict compilation without
`--allow-incomplete` refuses missing required evidence. A complete report is not
scientific approval. A partially extracted donor is excluded from donor-level
population summaries, even in a draft.

## Extract actual QC on the HPC

Each task sequentially reads cropped, initial-QC, and final combined-label RDS
files, retaining only metadata in memory between reads. The cropped object's
sparse counts/control assays are used to reconstruct QC metrics because the crop
predates initial QC. It reconstructs the
**stored** initial QC filters and checks exact retained-cell membership, then checks
that final cell IDs are a subset. It reads the cerebellar crop list and, when
available, raw `cells.csv.gz`/`cells.csv` and `metrics_summary.csv` in configured raw
directories. Crop lists from the same acquisition must be disjoint; overlapping
membership fails before writing. Missing raw metrics remain explicitly unavailable. Raw CSV metrics
support either one wide summary row or a two-column key/value layout; unsupported
schemas fail rather than guessing units. `experiment.xenium` is inventoried as
provenance; its software fields are not automatically interpreted.

Expected extraction runtime: **2–10 minutes/task; 28–140 minutes for 41 tasks at
concurrency three**, excluding queue wait, plus 1–5 minutes each for compile and
plots. Low confidence: metadata extraction with large saved-object/raw-table I/O,
not a measured benchmark. Requested limits are **64 GB, one CPU, one hour/task**;
combine uses **8 GB, one CPU, 30 minutes**. These reuse existing per-sample launcher
conventions and are not expected runtimes. Time one representative task and collect
peak memory before scaling if input sizes are substantially larger than existing
per-sample objects.

```bash
export ATLAS_QC_OUTPUT_DIR=outputs/xenium/atlas_qc_hpc_review1
export ATLAS_QC_ALLOW_INCOMPLETE=true
bash scripts/submit_xenium_atlas_qc.sh --dry-run
bash scripts/submit_xenium_atlas_qc.sh
```

The helper derives and hashes the task map from `--list`, preflights required
RDS paths without loading them, submits a bounded array, then an `afterok`
compile/plot job. Metadata/templates must be synced beside the scripts. To run
one HPC pilot instead, expected 2–10 minutes after a <1-minute path-only check:

```bash
Rscript --max-ppsize=500000 scripts/xenium_atlas_qc.R \
  --task 1 --dry-run --output-dir outputs/xenium/atlas_qc_pilot1
Rscript --max-ppsize=500000 scripts/xenium_atlas_qc.R \
  --task 1 --output-dir outputs/xenium/atlas_qc_pilot1
```

Run these only on an allocated HPC compute node, not a login node. The current
final input is `outputs/xenium/vz_rl/01_combined_labels/rds/` with unchanged native
filenames. QC counts start after crop and use a linear cropped → initial QC →
final sequence; VZ/RL/residual branches are not falsely presented as sequential
filters. Individual failure counts overlap; exclusive combinations reconcile
to the cropped total. Removed cells receive no inferred final identity.

## Tables, panels, and limitations

Per-sample bundles and compiled CSVs include retention, QC quantiles, overlapping
and exclusive failures, final population counts, imported-control availability,
stored thresholds, raw morphology/crop provenance audits, and raw instrument
metrics. Missing control assays mean unavailable evidence, not measured zero.
Constant area 20 is a flag only; direct raw-area comparison uses 0.011 square
micrometer tolerance to allow two-decimal rounding. No findings trigger removals.
Provenance records paths/sizes/times and MD5 for files ≤10 MiB; large RDS files
are not reread just for hashing. Compilation rejects changed/missing recorded
sources. Existing production validators/reports are inventoried, **not rerun or
treated as passing**. Existing QC text reports are inventoried; detailed retention
comes from the new reconciled extraction rather than fragile text parsing.

Plot exports include PDF, 600-dpi Cairo TIFF, and small PNG previews:

- evidence completeness and anatomical coverage;
- technical input-directory coverage and lateral-to-medial section positions;
- confirmed donor-age coverage and separately labeled provisional working coverage;
- age/run confounding when run metadata are supplied;
- cell retention, per-section QC medians, QC failure reasons, control availability;
- raw morphology diagnostics and density when reviewed cerebellar areas exist;
- final-population donor detection, excluding unconfirmed groups/incomplete donors;
- descriptive leave-one-donor-out composition when at least three eligible donors exist.

`plots/plot_manifest.csv` records generated and blocked panels. Missing biology is
not replaced by synthetic results. Donor detection is ≥1 observed cell and is not
proof of identity, sufficient coverage, or statistical power. Composition pools
sampled sections within donor; unequal anatomical sampling remains a limitation.
Leave-one-out output is descriptive whole-cohort composition sensitivity, **not**
a developmental hypothesis test. The pipeline does not perform donor-level DEG,
formal spatial enrichment, or automatic image-based segmentation validation.

For local lightweight review, `--pdf-only` uses base vector PDF plus PNG, omitting
Cairo TIFFs. Production defaults retain Cairo PDF/TIFF. Local RDS extraction is
not authorized by this guide. Missing dependencies should be addressed through
the explicit project renv workflow, not installed by scripts/jobs.

## Validation and checkpoint

Parsing should take <10 seconds/script. Both synthetic tests should take <1 minute
each and use tiny fabricated inputs; the integration test exercises extraction,
raw-area/crop matching, compilation and plot export without real data. Configuration
structural validation takes <1 minute; `--check-files` on HPC is filesystem-dependent.
Shell syntax checks take <5 seconds/file.

From scripts/:

```bash
Rscript --max-ppsize=500000 test_xenium_atlas_qc.R
Rscript --max-ppsize=500000 test_xenium_atlas_qc_integration.R
bash -n run_xenium_atlas_qc.slurm
bash -n submit_xenium_atlas_qc.sh
```

Before a Git checkpoint, review
status/diff and the new configuration templates, validate, and commit only when
ready. This local folder is not a Git checkout; no commit/push was made. No package
environment changes were made and no renv snapshot is needed solely for these edits.
After HPC jobs complete, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`,
`AllocCPUS` and `ExitCode` using sacct to refine estimates.
