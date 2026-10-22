# Plot-only refresh with the revised palette and VZ marker panels

This workflow creates a **new production figure snapshot** under
`outputs/plot_refresh/<tag>/`. It does not replace the old plots, edit metadata,
write analysis RDS files, renormalize counts, rerun PCA/UMAP/integration, or run
marker discovery. Existing snapshots are protected: choose a new tag for a rerun.
An interrupted, partial snapshot is an error, not a successful result.

The shared palette now contains the latest requested colors. The existing
`Microglia` palette keys remain as compatibility aliases for saved objects.
The snapshot shows them as `Immune` in temporary plotting labels only. This
broadens the displayed category; it does not change the saved marker evidence
or establish additional immune-cell identities. BG Progenitors retains BG's
violet until a separate replacement color is chosen. Retired Neuroblasts and
Cycling VZPs keys remain compatible with historical objects; they are not
introduced into current samples.

## Coverage

| Group | Tasks from `--list` | Inputs and refreshed figures |
|---|---|---|
| references | 6 | Full standardized and Xenium-panel Aldinger, Sepp, Science objects; all saved UMAP reductions, broad DotPlots, counts/proportions, PDF reports |
| samples | Current manifest count (41 at preparation) | Six saved objects per sample: consensus, VZ, RL, residual, combined subclusters, corrected broad; saved UMAPs, DotPlots, global/faceted spatial plots, counts/proportions, PDF reports |
| pooled | 3 | Pooled VZ, RL, residual subcluster objects; saved UMAPs where available, DotPlots, counts/proportions, PDF reports |
| combined | 2 | Saved nonspatial processed and spatial integrated objects; consensus and combined-label views, UMAPs, DotPlots, counts/proportions, PDF reports |

The snapshot uses consistent new filenames inside its separate output tree.
Corrected-broad sample inputs are detected in either `02_corrected_broad/rds/`
or the earlier `01_corrected_broad/rds/`. Files in both layouts are an ambiguity
error; partial sets are never stitched together. The chosen directory and
sample count are printed without reading RDS contents.
It does **not** re-create every historical filename or diagnostic: original
published-label comparisons, gene-specific plots, QC violins, marker-discovery
heatmaps, cross-study comparisons, and batch-only diagnostics remain separate.
It does not re-create historical whole-tissue resolutions. The samples group
uses `get_production_output_paths()` and the current manifest/configuration.
Merged spatial objects receive saved-UMAP views; whole-slide spatial maps are
drawn only from individual sample objects so different fields are not overlaid.

Each view writes an exact output manifest, source path/size/mtime, palette and
code hashes, marker availability, duplicate-marker notes, counts, and a PDF
report. `--validate` requires all expected views and files to exist, be nonempty,
and match the input and code provenance. A nonempty output directory cannot be
reused silently, even after a partial job. Use a fresh tag to retry that task.

## Broad DotPlot markers and reference formatting

Production broad plots for references, consensus samples, and merged/corrected
broad objects draw from `markers` in `color_palette.R`. Existing code intersects
that list with genes present in each assay, so missing genes can produce
different columns. The snapshot records these differences explicitly. VZ, RL,
and residual subcluster panels are intentionally different from broad panels.
The original Aldinger figure-cluster and Sepp cell-type diagnostics use the
same broad marker source but different identity axes. Sepp's original-cell-type
diagnostic calculates its summaries directly from counts, so it is not identical
to plotting the standardized, saved reference object.

The broad list repeats MKI67 under RL and Cycling Cells. The current shared
DotPlot helper rejects duplicated features. The snapshot keeps one displayed
copy in the first group (RL) and writes a duplicate-marker note; the scientific
marker definitions and expression values are unchanged. Combined VZ/RL panels
likewise show shared genes once in the first panel that requests them.

Reference snapshot DotPlots use dynamic width/height, colored group headers,
italic gene labels, top legends, the shared -2.5 to 2.5 scaled-expression range,
and a fixed 0-100% expression-size scale. Gene scaling remains within each
plotted dataset: equal colors across datasets are not absolute expression
comparisons. No expression normalization or label reassignment is performed.

VZ snapshot DotPlots use the new exact Purkinje, GABA and glial display marker
orders and move GABA Interneurons after MLIs. Missing required VZ display genes
are errors. The whole-tissue cell IDs, labels and all source metadata columns
are checked for preservation. Mapped regional NA labels are intentional
nonbranch background; a sample with no mapped regional cells receives an
explicit placeholder rather than a silently missing output.

## Run on Cheaha

Copy the code bundle into `/home/acflint/R/Projects/XeniumFCProject/scripts/`
and extract it **from that scripts directory**. The archive contains code only,
including the shared plotting helpers. Use the existing R/4.4.1 project renv
environment. Do not rerun reference standardization/panel-building scripts for
a figure refresh: those scripts can recompute embeddings and save RDS files.

```bash
cd /home/acflint/R/Projects/XeniumFCProject/scripts
tar -xzf plot_palette_refresh_20261019_v2.tar.gz
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# Path-only preflight for every group; no RDS objects are loaded.
bash scripts/submit_refresh_current_plots.sh --dry-run --tag=palette_20261019

# Submit only after the preflight passes and the task mapping is reviewed.
bash scripts/submit_refresh_current_plots.sh --submit --tag=palette_20261019
```

Inspect any full task map before submitting:

```bash
Rscript --max-ppsize=500000 scripts/refresh_current_plots.R --group=samples --list
Rscript --max-ppsize=500000 scripts/refresh_current_plots.R --group=references --list
```

To submit only references, independently of unavailable downstream objects:

```bash
bash scripts/submit_refresh_current_plots.sh --dry-run \
  --groups=references --tag=palette_20261019
bash scripts/submit_refresh_current_plots.sh --submit \
  --groups=references --tag=palette_20261019
```

The reference array currently contains six reference/scope rows. The submission
helper derives every array bound from the R driver's mapping and reports all
selected groups' preflight failures before stopping without submitting jobs.
Use `--groups=references,samples` or another explicit subset to run ready groups
independently. Each submitted task repeats its own path-only preflight.

After all jobs finish:

```bash
for group in references samples pooled combined; do
  Rscript --max-ppsize=500000 scripts/refresh_current_plots.R \
    "--group=${group}" --tag=palette_20261019 --validate || break
done

# Replace JOBIDS with the IDs returned by sbatch.
sacct -j JOBIDS --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Validation reads completion CSVs and file metadata, not large RDS objects.
Inspect representative reference, VZ and large-sample PDFs before treating the
new snapshot as the final figure set. Keep the original figures as a rollback.

## Estimated wall-clock runtime, excluding queue wait

| Group | Expected per task | Array concurrency | Approximate group wall time | Requested limit / memory per task |
|---|---|---|---|---|
| references | 5-20 min | 1 | 30-120 min for 6 tasks | 1 h / 64 GB |
| samples | 15-90 min for six sequential object views | 3 | 3.5-21 h for 41 samples | 2 h / 64 GB |
| pooled | VZ 10-60 min; RL 5-20 min; residual 1-5 min | 1 | 16-85 min | 12 h / 256 GB |
| combined | 15-120 min | 1 | 30-240 min for 2 tasks | 24 h / 600 GB |

These are **low-confidence estimates** for a new plotting-only workflow, based
on static inspection, 1.45-million-cell pooled VZ scale, the recorded 9.33-minute
RL plot refresh, and the 0.95-4-minute completed residual sample reports. Each
sample task here does substantially more work than one residual report. There
is no measured reference/combined snapshot run yet; measure the largest sample
and full reference to refine time/memory. The resource envelopes come from the
existing reference, per-sample report, pooled, and combined launchers; they are
not new memory measurements. The sample tasks load one RDS at a time. The
four groups can run independently; if resources are available concurrently,
the sample group is likely to determine total elapsed runtime (3.5-21 h).

Local syntax checks are expected under 10 seconds per R file and 5 seconds per
shell file. Path-only preflight and post-run manifest validation should take
under 1 minute locally, but HPC metadata latency can extend that. The small
synthetic plotting test is expected under 1 minute and establishes formatting
and control-flow behavior, not scientific correctness or large-object memory.
The local test passed for broad, VZ, RL, residual and combined panels using
base-PDF surrogate devices because the local R installation cannot load its
Cairo/X11 libraries. Production TIFF/Cairo rendering must be checked on Cheaha;
no dependencies were installed locally. All four path-only dry-runs reached
their expected mappings and failed because project analysis RDS inputs are
not available in this local checkout. No production analysis was run.

After reviewing the completed code and validation, make an appropriate Git
checkpoint in the tracked repository if ready. This local scripts directory
does not currently expose Git metadata. No package versions or renv lockfile
were changed, so this update alone does not call for `renv::snapshot()`.
