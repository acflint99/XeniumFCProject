# Seven new samples: crop through RL/VZ split

The manifest contains 41 analyzed samples. Original rows 1–34 are unchanged.
New samples are `GZFB328_1_X_G_1` through `_6`, followed by
`GZFB124_1_X_G_2`. Only section 2 is included from
`GZFB124_1_X_G_3__2__1`; no section-1 or section-3 outputs are requested.

Confirmed crop files:

- Six GZFB328 directories: `cerebellum_cells_stats.csv` in each directory.
- GZFB124 shared directory: `GZFB124_2_cerebellum_cells_stats.csv`.

The pipeline uses the exported `Cell ID` list as the crop. Check the export
visually in Xenium Explorer to confirm that it contains only the desired tissue.
No QC, clustering, transfer, consensus, or regional identity settings were changed.

## Copy the update to Cheaha

Upload `new_samples_20261014_update.tar.gz` through the Cheaha file browser to
`/home/acflint/R/Projects/XeniumFCProject/`. It contains only configuration,
scripts, and this guide. In a Cheaha terminal, preserve the previous versions
before extracting:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
backup_dir="scripts_backup_new_samples_$(date +%Y%m%d_%H%M%S)"
mkdir "$backup_dir"
tar -tzf new_samples_20261014_update.tar.gz | while IFS= read -r path; do
  if [[ -f "$path" ]]; then
    mkdir -p "$backup_dir/$(dirname "$path")"
    cp -p "$path" "$backup_dir/$path"
  fi
done
tar -xzf new_samples_20261014_update.tar.gz
```

The HPC must already have the current underlying R scripts and project `renv`
library. This update does not install packages or update `renv.lock`.

## Start preprocessing and annotation

Use an existing Intel interactive R session/terminal, or request an allocation:

```bash
srun --partition=interactive --constraint=intel --cpus-per-task=1 --mem=8G --time=00:30:00 --pty bash
cd /home/acflint/R/Projects/XeniumFCProject
bash scripts/submit_xenium_new_samples_20261014.sh --dry-run
```

Preflight checks the full manifest/input files, selects tasks from each driver's
`--list`, reads the seven crop CSVs, and verifies unique PCW metadata matches.
It does not load Xenium or RDS objects. Expected runtime: 1–3 minutes, low
confidence from static inspection and filesystem dependence; allocation limit
30 minutes. Missing packages or files must be resolved explicitly before submission.
The metadata lookup must match `GZFB328_1_X_G` and `GZFB124_1_X_G` in
`metadata/samples_meta.xlsx`; enter verified PCW values if either is missing.

When preflight succeeds:

```bash
bash scripts/submit_xenium_new_samples_20261014.sh --submit | tee "logs/new_samples_submission_$(date +%Y%m%d_%H%M%S).txt"
```

The helper submits only the seven new samples, preserves overwrite protection,
prints job IDs immediately, and queues these dependencies:

1. Single-slide preprocessing tasks 20–25 at concurrency 3.
2. Split-slide task 16 after single-slide preprocessing finishes, keeping total
   heavy preprocessing concurrency at three or fewer.
3. Production resolution-5 clustering tasks 35–41 at concurrency 4.
4. Three reference-transfer arrays for tasks 35–41, concurrency 3 per reference.
5. Consensus tables for the seven selected samples only, then consensus
   objects/plots for tasks 35–41 at concurrency 6.
6. Unknown QC audits for tasks 35–41 at concurrency 3, plus broad reference
   scores for tasks 103–123 at concurrency 6.

The broad-score stage requires the existing broad reference validation table.
If that table is missing, run the protected reference-validation stage described
in README first (5–20 minutes expected; 45-minute requested limit).
If any submission fails partway through, use the printed job IDs to inspect
what was queued before retrying. Do not blindly resubmit completed stages.

| Stage | Expected per task | Cohort wall time at the configured concurrency | Requested task limit |
|---|---|---|---|
| Crop/QC/resolution-1.5 preprocessing | 1–6 h, low confidence | 3–18 h, including the sequential split task | 12 h |
| Resolution-5 reclustering | 4–10 min | 8–20 min at 4 | 1 h |
| Reference transfer | 7–47 min historically | roughly 25–150 min with three arrays, each at 3 | 4 h |
| Selected consensus tables | 1–5 min | 1–5 min | 30 min |
| Consensus objects/plots | 10–60 min, low confidence | 20–120 min at 6 | 2 h |
| Unknown QC audit | 1–3 min, low confidence | 3–9 min at 3 | 30 min |
| Broad scores | 2–15 min, low confidence | 8–60 min for 21 tasks at 6 | 2 h |

Allow roughly 4–24 hours of active processing to reach review, excluding queue
wait. Transfer estimates use `slurm_resource_history.psv`; other ranges use
README and static code inspection. QC audit also uses the earlier 42–75-second
per-sample run, allowing additional time for raw morphology reads. New cell counts and actual `sacct` readings
will improve the estimates. Existing CPU/memory/time/partition requests were retained.

## Review gate and completion of the split

Wait until audit and broad-score jobs finish successfully. The existing global
QC summaries, calibration tables, broad consensus, and route manifests contain
the original cohort. They must be extended to 41 samples before regional
subsetting can validate provenance. The commands below are the remaining
sequence, **to be run only after explicitly approving replacement of those
shared tables/reports and saving their previous versions**. They do not authorize
replacing existing Seurat objects or plots from other stages.

From the same Intel session (path-only dry-runs generally under 1 minute each):

```bash
Rscript scripts/xenium_unknown_02_compile_qc_review.R --dry-run
Rscript scripts/xenium_broad_03_compile_calibration.R --dry-run
```

After authorizing the shared report/table replacements, compile the 41-sample
QC evidence and broad calibration using the matching launchers' protected
overwrite flags. Review the new Unknown populations, append explicit `keep`
or `remove` decisions and notes to `config/reviewed_unknown_qc.csv`, and retain
all previous decisions. Do not automatically mark new populations as reviewed.
Keep `config/reviewed_broad_reference_thresholds.csv` unchanged unless a
scientific threshold change is explicitly approved. Compilation is approximately
2–10 minutes for broad scores and 2–15 minutes for QC reports, low confidence.

The exact protected commands for this reviewed step are:

```bash
sbatch --wait --export=ALL,UNKNOWN_QC_COMPILE_OVERWRITE=true scripts/run_xenium_unknown_02_compile_qc_review.slurm
sbatch --wait --export=ALL,BROAD_CALIBRATION_OVERWRITE=true scripts/run_xenium_broad_03_compile_calibration.slurm
# Review and append the new Unknown decisions here before proceeding.
Rscript scripts/validate_config.R --check-files
Rscript scripts/xenium_unknown_03_apply_qc_review.R --dry-run 35
sbatch --wait --array=35-41%3 scripts/run_xenium_unknown_03_apply_qc_review.slurm
Rscript scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --dry-run
Rscript scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --overwrite
sbatch --wait --export=ALL,BROAD_CONSENSUS_OVERWRITE=true scripts/run_xenium_broad_04_apply_reviewed_thresholds.slurm
Rscript scripts/xenium_unknown_04_build_reference_routes.R --dry-run
sbatch --wait --export=ALL,UNKNOWN_ROUTE_OVERWRITE=true scripts/run_xenium_unknown_04_build_reference_routes.slurm
```

Unknown-QC application should take 1–10 minutes per new task, approximately
3–30 minutes at concurrency 3. Combining its small routing tables is under
1 minute; broad-threshold application and route building are roughly 2–10
minutes each. These are low-confidence static estimates. Wait for successful
completion at each step; do not continue after a failure.

Once the 41-sample route manifests validate:

```bash
for task in 35 36 37 38 39 40 41; do
  Rscript scripts/xenium_vz_01_subset_samples.R --dry-run "$task"
  Rscript scripts/xenium_rl_01_subset_samples.R --dry-run "$task"
done
sbatch --array=35-41%3 scripts/run_xenium_vz_01_subset.slurm
sbatch --array=35-41%3 scripts/run_xenium_rl_01_subset.slurm
```

Regional dry-runs read small manifests and check paths/provenance, not large
RDS objects. Subsetting should take 1–10 minutes per task, approximately 3–30
minutes per branch at concurrency 3; each launcher requests 30 minutes.
Estimates are low confidence until actual cell counts and `sacct` are available.

The requested endpoint is seven new files in each directory:

- `outputs/xenium/vz/01_subsets/rds/<sample>_VZsubset.rds`
- `outputs/xenium/rl/01_subsets/rds/<sample>_RLsubset.rds`

Inspect each subset job's validation messages for cell counts, exact cell IDs,
review provenance, and route agreement. Empty regional identities are errors.
Do not rerun cohort merge/integration/QC/subcluster mapping for this request.
Existing merged and integrated regional objects still represent the original
cohort; including the seven samples would require regeneration and fresh cluster
review, since old manual subcluster mappings cannot be assumed to remain valid.
Rollback uses the saved manifests/shared tables and previous scripts; the original
34 per-sample objects are retained.

Collect resource evidence after completion (replace IDs with printed job IDs):

```bash
sacct -j JOB_IDS --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Local validation passed structural configuration, all 11 shell syntax checks,
and four R driver task-map listings. Local path-only preprocessing/reclustering
dry-runs stopped on missing raw data/RDS inputs, as expected for this scripts-only
checkout. HPC input checks and jobs have not been run. This local directory is not a Git
checkout; on the HPC review `git status` and `git diff`, validate the completed
logical change, and commit only when ready. No package environment changed,
so this update does not call for `renv::snapshot()`.
