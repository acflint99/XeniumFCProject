# Apply approved QC decisions and finish the seven RL/VZ splits

## Quick start: two uploads, two commands

Upload `reviewed_unknown_qc_APPROVED.csv` to the HPC project root and
`finish_xenium_new_samples_20261015.sh` to its `scripts/` directory.
From an Intel interactive terminal in the project root:

```bash
bash scripts/finish_xenium_new_samples_20261015.sh --dry-run
# Only after all preflight checks pass:
bash scripts/finish_xenium_new_samples_20261015.sh --submit
```

Preflight takes about 1–3 minutes and checks paths, CSVs, approved thresholds,
and output absence without loading large RDS objects. Submission backs up the
shared reports/tables, installs the approved review, and queues the stages
below with success dependencies. It prints all job IDs and returns immediately
after submission; the jobs then run independently of the terminal. Do not also
run the manual steps below for the same attempt.

The final check confirms that all 14 new subset files exist and are nonempty;
scientific identity and cell-ID checks are performed within the subset jobs.
The small QC routing-combine and final file-check jobs each request 1 CPU,
4 GB, and 15 minutes on `express` with the Intel constraint; expected runtime
is under 1 minute each from static inspection. Other stages reuse their
existing launchers and requested resources.

Shell syntax, both embedded R expressions, task mappings, and compliance with
the R startup-option rule were checked locally. HPC preflight and analysis
jobs have not been run by Codex. If a submission fails partway through,
inspect the printed IDs before retrying; new RDS overwrite protection stays on.

Approved on 2026-10-15: retain all 19 new Unknown populations (15,222 cells).
The approved manifest has 91 decisions, preserving the previous 72 exactly.
Unknown QC retention does not assign a regional lineage. Existing reviewed
broad-reference thresholds and two-of-three routing remain authoritative.

Upload `outputs/qc_review_new_samples_20261015/reviewed_unknown_qc_APPROVED.csv`
to `/home/acflint/R/Projects/XeniumFCProject/` on Cheaha. The following commands
are for the Intel interactive terminal already used for this analysis.
No package installation, environment restore, commit, or push is included.

## 1. Install the manifest and preserve previous shared outputs

Expected runtime: under 2 minutes, depending on filesystem speed. This copies
small tables and the combined PDF, never large analysis RDS files.

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

backup_root="outputs/review_gate_backups/approved_new_samples_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$backup_root"
cp -p config/reviewed_unknown_qc.csv "$backup_root/reviewed_unknown_qc_before.csv"
cp -p config/reviewed_broad_reference_thresholds.csv "$backup_root/reviewed_broad_reference_thresholds_before.csv"

for path in \
  outputs/xenium/annotation/04_unknown_qc_review/tables \
  outputs/xenium/annotation/04_unknown_qc_review/plots/Xenium_Unknown_QC_All_Samples.pdf \
  outputs/xenium/annotation/04_unknown_qc_reviewed/tables \
  outputs/validation/broad_reference_workflow/03_calibration/tables \
  outputs/xenium/annotation/05_broad_reference_consensus/tables \
  outputs/xenium/annotation/04_unknown_reference_routing/tables; do
  if [[ -e "$path" ]]; then
    mkdir -p "$backup_root/$(dirname "$path")"
    cp -a "$path" "$backup_root/$path"
  fi
done
```

The updated helper treats the active HPC review as authoritative. It imports
only the 19 approved decisions for the new cohort, preserves all existing HPC
rows, notes, and extra columns, and rejects a conflicting existing decision.
Its preflight is read-only; it writes the merged manifest only after backup
in `--submit` mode. Repeated installation does not rewrite an unchanged review.
Do not replace the HPC manifest with the full uploaded local copy.

Use the Quick start helper above to perform installation and queue the rest.
The remaining sections document the jobs it runs; do not submit them a second
time if the helper has already submitted them.

## 2. Refresh combined QC and broad calibration reports

These commands deliberately replace the shared reports/tables backed up above.
They do not replace individual QC PDFs or consensus RDS files. Keep the reviewed
broad thresholds unchanged. Run each command only after the preceding one succeeds.

```bash
Rscript --max-ppsize=500000 scripts/xenium_unknown_02_compile_qc_review.R --dry-run
sbatch --wait --export=ALL,UNKNOWN_QC_COMPILE_OVERWRITE=true scripts/run_xenium_unknown_02_compile_qc_review.slurm

Rscript --max-ppsize=500000 scripts/xenium_broad_03_compile_calibration.R --dry-run
sbatch --wait --export=ALL,BROAD_CALIBRATION_OVERWRITE=true scripts/run_xenium_broad_03_compile_calibration.slurm
```

Dry-runs: under 1 minute each, paths/small tables only. QC compilation is
estimated at 1–5 minutes (15-minute requested limit); broad calibration at
2–10 minutes (30-minute requested limit). Estimates are low confidence from
static inspection; the completed 41-second selected-consensus job is a different
stage and is not a measurement of these compilation jobs. Queue wait is excluded.

The generated QC review template is blank by design: **do not copy it over the
approved `config/reviewed_unknown_qc.csv`**.

## 3. Refresh QC-application provenance for all 41 samples

The reference-route builder checks that every per-sample table records the
current review CSV's size and modification time. Appending 19 decisions changes
that file, so rerun all 41 QC-application tasks. Because all decisions are `keep`,
these tasks read the consensus objects, retain every cell, and replace only
their small QC-application tables. They fail if a stale filtered RDS exists;
do not delete or override that error without investigation.

```bash
Rscript --max-ppsize=500000 scripts/xenium_unknown_03_apply_qc_review.R --dry-run 1
Rscript --max-ppsize=500000 scripts/xenium_unknown_03_apply_qc_review.R --dry-run 35
sbatch --wait --array=1-41%3 --export=ALL,UNKNOWN_QC_APPLY_OVERWRITE=true scripts/run_xenium_unknown_03_apply_qc_review.slurm

Rscript --max-ppsize=500000 scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --dry-run
Rscript --max-ppsize=500000 scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --overwrite
```

Dry-runs and the combine operation read paths/small CSVs, not large RDS files;
each should take under 1 minute. Actual QC application reads an RDS per task:
estimate 0.5–2 minutes per task, roughly 7–30 minutes at concurrency three,
excluding queue wait. This is low confidence based on the 26–43-second new-sample
QC audits and static code inspection; request remains 30 minutes per task.

## 4. Apply existing broad thresholds and build routes

```bash
Rscript --max-ppsize=500000 scripts/xenium_broad_04_apply_reviewed_thresholds.R --dry-run
sbatch --wait --export=ALL,BROAD_CONSENSUS_OVERWRITE=true scripts/run_xenium_broad_04_apply_reviewed_thresholds.slurm

Rscript --max-ppsize=500000 scripts/xenium_unknown_04_build_reference_routes.R --dry-run
sbatch --wait --export=ALL,UNKNOWN_ROUTE_OVERWRITE=true scripts/run_xenium_unknown_04_build_reference_routes.slurm
```

Dry-runs are path/small-table checks, expected under 1 minute each. Each actual
stage is estimated at 1–5 minutes (low confidence, static inspection). Requested
limits are 30 minutes for broad application and 15 minutes for routes.

## 5. Create only the seven new VZ/RL subsets

```bash
for task in 35 36 37 38 39 40 41; do
  Rscript --max-ppsize=500000 scripts/xenium_vz_01_subset_samples.R --dry-run "$task"
  Rscript --max-ppsize=500000 scripts/xenium_rl_01_subset_samples.R --dry-run "$task"
done

sbatch --array=35-41%3 --export=ALL,VZ_SUBSET_OVERWRITE=false scripts/run_xenium_vz_01_subset.slurm
sbatch --array=35-41%3 --export=ALL,RL_SUBSET_OVERWRITE=false scripts/run_xenium_rl_01_subset.slurm
```

These dry-runs read route manifests and inspect paths/provenance, not large RDS
objects: approximately 1–3 minutes total for 14 separate R startups. Verify each
reports PASS and no existing output before submitting. Do not use overwrite
flags for these new subset RDS files.
Actual subsetting is estimated at 0.5–3 minutes per task, approximately 2–10
minutes per branch at concurrency three (30-minute requested limits). This is
low confidence, based on the completed new-sample 26–43-second audits,
1:26–2:49 consensus jobs, and static inspection; subsetting itself is not timed.

Both branches may run concurrently. Inspect all 14 task states/exit codes and
their validation logs, including unique cell IDs, nonempty retained populations,
cell counts, and route agreement. Expected output files are:

```text
outputs/xenium/vz/01_subsets/rds/<sample>_VZsubset.rds
outputs/xenium/rl/01_subsets/rds/<sample>_RLsubset.rds
```

The samples are GZFB328 sections 1–6 and GZFB124 section 2. This endpoint does
not rerun merged regional integration or change reviewed subcluster identities.
Collect JobID, Elapsed, State, MaxRSS, ReqMem, AllocCPUS, and ExitCode from `sacct`
for the new jobs. Overall active-time estimate for this continuation: roughly
15–60 minutes, low confidence, excluding queue waits between stages.

## Provenance and rollback

The backup directory stores the prior review, unchanged broad thresholds,
and replaced shared tables/report. Restore those together if rolling back;
do not mix old route provenance with the new review file. Original RDS inputs
and the previous 34 VZ/RL subsets are not overwritten by these commands.
For a Git checkpoint, review the approved manifest and this guide, validate,
and commit only when ready. The backup directory and generated report outputs
do not belong in the source commit. No package environment changed.
