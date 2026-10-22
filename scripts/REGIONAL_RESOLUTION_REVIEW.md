# VZ and RL resolution review (Stage 5c)

Run this diagnostic after Stage 5 and before selecting a partition for Stage 6.
For existing completed analyses, start directly from their saved Stage-5 objects.
There is no need to rerun integration, QC, annotation, or mapping to add this review.
The current resolution-0.8 selection and all production RDS files remain unchanged.

## Method and scope

`xenium_regional_05c_resolution_stability.R` provides both branches through
`--region=VZ|RL`. It reads the resolutions recorded in the Stage-5 provenance:
VZ full path 0.3/0.5/0.8/1.0, RL full path 0.3/0.5/0.8; zero-removal reuse can
contain fewer resolutions. Missing recorded partitions are errors, not silently
skipped. No new resolution is introduced by the diagnostic.

1. Draw a clustree using full-object saved assignments at all recorded resolutions.
2. In each repeat, sample 80% of cells without replacement **within every sample**
   (floor the count, minimum one). Biological samples retain their relative cell
   contributions; this does not equalize sample sizes. The active manifest currently
   has 41 samples; the code reads it rather than assuming 34 or 41.
3. Hold the production Harmony representation fixed, rebuild the SNN graph on the
   retained embeddings, and run Louvain at each recorded resolution. Dimensions,
   graph k, pruning, distance/neighbor settings and clustering settings come from
   the object's saved Seurat command records. Older zero-removal objects can use
   k=20, while full reintegration currently uses k=30. Ambiguous or incomplete command
   provenance fails early; settings are never silently substituted.
4. Restrict original cluster memberships to the retained cells and match each
   cluster to the recovered cluster with maximum Jaccard overlap:
   `intersection / union`. Matching is many-to-one; merges and splits both reduce
   recovery. Record absent clusters as NA with explicit absence counts.
5. Report per-cluster mean/median, standard deviation, empirical 10th/90th
   percentiles, and all repeat scores. Percentiles describe repeat variability;
   they are not confidence intervals. Compare cluster-wise results as well as
   equal-cluster and cell-weighted resolution summaries. No cutoff or automatic
   resolution winner is imposed.

Each repeat uses configured seed + repeat number for sampling, graph construction
and clustering. The same subset and clustering seed are used at all resolutions
within a repeat. This measures cell-sampling and stochastic clustering robustness
conditional on the current integration. It does not validate integration itself
or establish reproducibility across biological specimens. Five repeats are a
runtime pilot; twenty are exploratory, and convergence needs review before final
reporting. A 20-repeat run reproduces the first five seeds from the pilot when
inputs, packages and other parameters are unchanged, but recomputes those repeats.

This is an adaptation of cluster-wise stability assessment, **not chooseR**:

- Zappia & Oshlack (2018), clustering trees: https://doi.org/10.1093/gigascience/giy083
- Hennig (2007), cluster-wise stability: https://doi.org/10.1016/j.csda.2006.11.025
- Related subsampling framework, Patterson-Cross et al. (2021):
  https://doi.org/10.1186/s12859-021-03957-4

## Dependencies and preflight

Use the existing R/4.4.1-foss-2023b module and project renv library. From the HPC
project root, the checks below should take under one minute per region. `--list`
and `--dry-run` inspect configuration, the sample manifest, paths and output
protection only; they do **not** load an RDS. They do not establish internal object
validity. The real job checks object provenance, IDs, sample coverage, saved graph
settings and finite embeddings before graph construction.

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R --region=VZ --check-packages
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R --region=VZ --list
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R --region=RL --list
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R --region=VZ --dry-run
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R --region=RL --dry-run
```

`--check-packages` also renders a tiny clustree PDF from 12 artificial cells,
then removes that temporary file; estimate 5–30 seconds for rendering plus package
loading (usually under one minute total). Every real job performs this same check
before loading a production RDS. The driver supplies ggraph's edge-colour guide
as an explicit object to avoid the `Unknown guide: edge_colourbar` lookup error
seen in failed VZ job 41007543. No package upgrade is required by this code fix.
The same rendering helper is used for the real clustree PDF.

Required packages: here, yaml, Seurat, SeuratObject, Matrix, future, ggplot2,
ggraph and clustree. Clustree was missing from the local development library; HPC availability
must be checked. If missing there, perform an explicit interactive project-environment
update **before** submission, not inside a production job:

```r
renv::status()
renv::install("clustree")
renv::status()
# After reviewing the dependency changes and approving the lockfile update:
renv::snapshot()
```

Package setup is estimated at 2–20 minutes, potentially longer with compilation
or network delays (low confidence). Review dependency changes before accepting
updates; no installation, package upgrade or lockfile change is performed by this
diagnostic. The snapshot is appropriate only if the package environment changed.

## Submit the five-repeat review for both branches

These are two independent pooled jobs, not sample arrays. Each launcher uses one
CPU and sequential clustering to avoid the earlier parallel-memory failures.
The launcher performs a path-only dry-run before starting a fresh diagnostic.

```bash
sbatch scripts/run_xenium_vz_05c_resolution_stability.slurm
sbatch scripts/run_xenium_rl_05c_resolution_stability.slurm
```

Low-confidence expected active runtimes: **VZ 4–20 hours**, **RL 2–12 hours**
for five 80% repeats across the existing resolutions. With both jobs running
concurrently, approximately **4–20 hours** total; sequentially **6–32 hours**.
Queue wait is excluded. These replace the earlier preliminary conversation estimates.
Requested limits are **36 hours VZ / 24 hours RL**, both on `medium` with 600 GB,
matching the corresponding Stage-5 launchers. The current Cheaha documentation
lists a 50-hour medium-partition limit:
https://docs.rc.uab.edu/cheaha/hardware/ .

Basis: recorded full RL reintegration took 3h32m in the earlier cohort; the current
VZ sequential Stage-5 run took 6h14m for roughly 1.5 million cells. Current RL has
roughly 827,000 cells. Those jobs include additional work and do not time this new
subsampling workload. Each diagnostic loads the large RDS once, releases the object
after extracting embeddings/labels, then rebuilds one graph and clusters sequentially
per repeat. The existing 600-GB allocation is retained conservatively, not claimed
as measured diagnostic demand. Per-repeat elapsed seconds are saved; use the pilot
to refine time and memory before expanding. A timeout leaves staged partial evidence,
not a published complete report, and rerunning starts from the beginning.

After completion, collect actual resource evidence (normally under one minute):

```bash
sacct -j <VZ_JOB_ID>,<RL_JOB_ID> --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Twenty repeats require roughly four times the repeat computation, plus a single
input load. Based on the unmeasured pilot range, this could exceed the default
limits or even the 50-hour medium limit. Do not submit twenty repeats with guessed
resources. Once measured, select justified partition/time overrides and use fresh
output names via `--export=ALL,STABILITY_REPEATS=20,STABILITY_RUN_ID=review20`.
Use `STABILITY_FRACTION` only for an explicitly documented sensitivity comparison.

## Outputs and review

Each branch writes beneath
`outputs/xenium/<vz-or-rl>/05_post_qc/diagnostics/resolution_stability/pilot5/`:

- `clustree.pdf`: full-data cluster relationships across resolutions.
- `stability.pdf`: cluster-wise recovery distributions by resolution.
- `cluster_summary.csv`, `recovery_by_repeat.csv`, `resolution_summary.csv`:
  all recovery metrics, including rare/absent clusters.
- `sample_composition.csv`: full cluster-by-sample counts, including zeros.
- `sampling_by_repeat.csv`: original/retained counts, seeds and repeat elapsed time.
- `provenance.csv`, `samples.csv`, `sessionInfo.txt`, `METHODS.txt`: input paths,
  size/time and provenance checksums, parameters, package versions and limitations.
- `outputs.csv`: expected published files, sizes and checksums. Its presence marks
  the complete set; failed runs retain only `.staging_*` directories.

Existing run directories are protected. Prefer a new `STABILITY_RUN_ID` when
rerunning, including after a Stage-5 input refresh. Explicitly reviewed diagnostic
replacement uses `STABILITY_OVERWRITE=true`; the previous diagnostic directory
is archived alongside the new one. This does not permit overwriting production
objects. A run-id lock prevents simultaneous publication to the same destination;
after a killed job, verify that it has stopped before manually removing its empty
`.lock_<run-id>` directory.

For failures handled by R (including the plotting failure in job 41007543), the
run-id lock is released automatically. Failed `.staging_*` directories are left
as evidence and do not block resubmission. Once the rendering preflight passes,
rerun both launchers with the same `pilot5` ID if neither published directory
exists; no `STABILITY_OVERWRITE` is needed. An RL failure notification alone does
not establish the same error; inspect its stderr if it fails after this fix.

Review the clustree and individual recovery distributions with the existing Stage-6
marker evidence, developmental/sample composition, QC and spatial maps. Stable
technical artifacts are possible; a lower global score can also reflect real
developmental continua or rare cells. No automatic merger or relabeling is applied.

If a different production resolution is approved, revisit its Stage-6 label map
before applying labels, regenerate Stage-7 mapped objects, and refresh dependent
Stage-8/9 plots, combined VZ/RL/residual objects, and downstream spatial/count/export
outputs. Do not reuse the resolution-0.8 integer-to-label map for another partition.
No downstream outputs become stale merely because this read-only review is run.

## Lightweight development validation

Syntax parsing should take under 10 seconds per R file; the synthetic test should
take about 5–30 seconds and uses 100 artificial cells. Shell checks take under
5 seconds per launcher, and the project-wide launcher audit usually under a minute.
These are local checks, not scientific validation or a production data run.

```bash
Rscript --max-ppsize=500000 -e 'for (f in c("scripts/R/regional_cluster_stability.R", "scripts/xenium_regional_05c_resolution_stability.R", "scripts/test_regional_cluster_stability.R")) parse(f)'
Rscript --max-ppsize=500000 scripts/test_regional_cluster_stability.R
bash -n scripts/run_xenium_vz_05c_resolution_stability.slurm
bash -n scripts/run_xenium_rl_05c_resolution_stability.slurm
cd scripts
python3 update_r_launchers_stack_limit.py --check
```

Before an authorized Git checkpoint in the actual repository, review `git status`
and `git diff`, run the checks, and commit this logical unit when ready. The local
scripts directory used for development is not a Git checkout; no commit or push
was created.
