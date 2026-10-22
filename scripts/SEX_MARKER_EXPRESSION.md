# Sample sex-marker expression

`xenium_sex_marker_expression.R` summarizes **SRY, DDX3Y, RPS4Y1, KDM5D,
UTY, EIF1AY, USP9Y, and XIST** in each
manifest sample's whole-tissue cerebellum QC object:
`outputs/xenium/preprocess/02_qc/rds/<sample>_CB_QC.rds`.
It uses raw Xenium counts and all cells retained by that stage. It does not
repeat QC or clustering, change metadata, or write RDS files. The current
manifest has 41 samples; task bounds are derived from the manifest.

## Interpretation

The six core Y markers are DDX3Y, RPS4Y1, KDM5D, UTY, EIF1AY, and USP9Y.
SRY is retained as supporting evidence. XIST is associated with X inactivation and is
commonly useful as a female-associated marker. SRY need not be expressed in
cerebellum, so absent SRY expression does not establish that a sample is female.
DDX3Y and XIST are more useful together than SRY alone. See the
[XIST/DDX3Y expression study](https://pubmed.ncbi.nlm.nih.gov/32792067/) and
[multi-tissue sex-chromosome expression study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7374973/).
These studies do not validate the thresholds below for this Xenium dataset.
The added Y markers also have support from a
[human fetal brain study](https://www.nature.com/articles/s42003-025-08070-3).
Only features available in the saved Xenium assay can be summarized; adding
a marker name cannot recover genes absent from the targeted panel. The
`marker_coverage.csv` output reports availability across samples.

Default **exploratory** criteria:

- A supported marker has at least 10 cells with raw count >0, comprising at
  least 1% of the sample's cells. Both conditions must hold.
- At least 100 cells are required for a provisional call.
- A missing feature is reported as `measured=FALSE`, with expression metrics
  `NA`, rather than zero.
- **male-like:** one or more core Y markers supported without XIST support.
  XIST need not be measured. `evidence_pattern` distinguishes `single_Y_marker`
  from corroboration by `multiple_Y_markers`. Single-marker evidence needs review.
- **female-like:** XIST supported, at least one core Y marker measured, and
  no Y markers (including SRY) supported.
- **unresolved:** insufficient cells, conflicting XIST/Y support, SRY-only
  support, insufficient coverage, or neither pattern supported. In particular,
  low Y expression with missing XIST is labelled `low_Y_XIST_unavailable`:
  compatible with female origin, but not enough for a female-like call.

This expanded version removes the original requirement for XIST coverage
before calling a Y-positive sample male-like. It does not lower the existing
marker-support thresholds. Original three-marker tables cannot be combined
with expanded tables: the gene list and `marker_set_version` are checked.
Regenerate only the sex-marker reports in a fresh directory; no downstream
Seurat objects become stale. Retain the original report directory for comparison
or rollback.

These are provisional sample expression patterns, not confirmed sex or
karyotype. Low Y counts can reflect background; mixed signals can reflect
sample mixture or other biological/technical causes. Review the continuous
metrics, panel coverage, and available donor records before assigning sex.
No calls are written back to sample metadata. Different criteria change these
new reports only; existing pipeline outputs remain valid.

## Parallel HPC workflow

Copy all three new code files into the HPC project's `scripts/` directory.
From the project root, load the established R module and inspect the preflight:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export SEX_MARKER_OUTPUT_DIR=outputs/xenium/sex_marker_expression_expanded_v2
bash scripts/submit_xenium_sex_marker_expression.sh --dry-run
```

Preflight takes under 1 minute, checks paths and protected output names, prints
the task map, and loads no RDS objects. After reviewing a passing preflight:

```bash
bash scripts/submit_xenium_sex_marker_expression.sh
```

The helper submits one sample per array task, with concurrency **3** by default,
and a dependent combine job that runs only if every sample task succeeds.
Notifications go to `acflint@uab.edu`: the array emails on failure, and the
final combine job emails on completion, failure, or an invalid dependency.
Per-task completion emails are not enabled. These settings apply to new
submissions after copying the updated launcher and submission helper to HPC;
editing local files does not change already-submitted jobs.
It checks the task-map checksum in every job to detect changes after submission.
Do not edit the manifest or these scripts while the array is running.
To use a different concurrency, for example six:

```bash
SEX_MARKER_CONCURRENCY=6 bash scripts/submit_xenium_sex_marker_expression.sh --dry-run
SEX_MARKER_CONCURRENCY=6 bash scripts/submit_xenium_sex_marker_expression.sh
```

For a new run, choose a fresh directory by setting `SEX_MARKER_OUTPUT_DIR` to a
project-relative or absolute path. Existing outputs are protected; the helper
will refuse to replace them. If an array fails, inspect its logs and retry
only failed sample tasks using the R command below with their task IDs.
After all tasks succeed, combine manually. A partially written CSV also stays
protected; use a fresh output directory instead of deleting/replacing it.

## Outputs

The example above writes to `outputs/xenium/sex_marker_expression_expanded_v2/`.
Without an explicit override, the original default remains
`outputs/xenium/sex_marker_expression/`; existing files stay protected.

- `samples/<sample>_sex_markers.csv`: eight marker rows per sample, measured
  status, positive-cell counts/percentages, total marker transcripts, mean raw
  counts per cell, pseudobulk transcripts per 10,000 Xenium transcripts,
  support criteria, provisional call/reason, single/multiple Y evidence pattern,
  supported/missing marker lists, input path/size/mtime, and versions.
- `sample_sex_evidence.csv`: one row per manifest sample, with calls and key metrics.
- `gene_expression.csv`: the complete combined marker table.
- `marker_coverage.csv`: one row per marker, with measured/missing sample
  counts and the IDs of samples missing that marker.
- `sex_marker_expression.pdf`: dot plot with positive-cell percentage as size
  and log1p(mean raw counts) as color. Missing markers are explicitly labelled;
  measured zero expression has no dot. Color is descriptive and is not
  normalized for cell composition or differences in panel design.
- `provenance.txt`: interpretation notes, manifest path, output paths, and session.

The combine step requires every manifest sample exactly once with the same
support criteria and manifest task mapping. Missing or inconsistent tables fail.

## Individual tasks and custom criteria

These commands run analysis and belong on HPC compute nodes. Expected runtime
is 1–5 minutes per sample; combining takes 1–3 minutes.

```bash
Rscript --max-ppsize=500000 scripts/xenium_sex_marker_expression.R --list
Rscript --max-ppsize=500000 scripts/xenium_sex_marker_expression.R --dry-run --task 1 --output-dir outputs/xenium/sex_marker_expression_expanded_v2
Rscript --max-ppsize=500000 scripts/xenium_sex_marker_expression.R --task 1 --output-dir outputs/xenium/sex_marker_expression_expanded_v2
Rscript --max-ppsize=500000 scripts/xenium_sex_marker_expression.R --combine --output-dir outputs/xenium/sex_marker_expression_expanded_v2
```

`--min-positive-cells`, `--min-positive-pct`, and `--min-cells` adjust the
exploratory criteria for direct R invocations. Use the same values and a fresh
`--output-dir` for every sample and the combine step. The submission helper
uses the documented defaults.

## Resources and validation

Low-confidence expected runtime from static inspection: **1–5 minutes/task**;
for 41 samples, **14–70 minutes at concurrency 3**, or **7–35 minutes at
concurrency 6**, plus **1–3 minutes** for combine, excluding queue wait.
Each sample requests **1 CPU, 64 GB, and a 1-hour limit**. Combine requests
**1 CPU, 4 GB, and a 15-minute limit**. Requested limits are not runtime estimates.
The express partition supports limits up to 2 hours according to
[current Cheaha documentation](https://docs.rc.uab.edu/cheaha/hardware/).

There is no matching sex-marker `sacct` history. Existing single-sample
resolution-5 launchers request 32 GB, while reference-transfer history peaks
around 44–49 GB and includes additional reference work. The 64-GB request is
provisional headroom, not a measured requirement. Time and peak memory for the
first completed sample would improve these estimates.

Completed local validation: R syntax parse, both shell syntax checks, project
launcher stack-limit check, and 41 tiny synthetic Seurat objects covering
single/multiple Y support, missing XIST, low Y without XIST, female-like,
conflicting/absent/zero/low-support/few-cell patterns, cohort combination,
coverage reporting, incompatible marker sets, and overwrite refusal.
The repeatable synthetic check is `test_xenium_sex_marker_expression.R`;
expected runtime is under 1 minute and it uses only temporary synthetic inputs.
Local R is 4.5.1; HPC production
remains R 4.4.1 with the existing renv environment. No packages were installed.
The real-input dry-run found 0/41 QC RDS files in this local checkout; repeat
the preflight on HPC. No real-data analysis or job submission was performed.

After jobs complete, collect resource evidence (normally under a minute):

```bash
sacct -j ARRAY_JOB_ID --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Before a Git checkpoint in a Git-enabled checkout, review `git status` and
`git diff`, run the validation above (under 1 minute without real analysis),
and commit this logical unit when ready. The current local directory is not
a Git repository. No package environment change or renv snapshot is needed.
