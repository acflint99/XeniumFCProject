# Replace Xenium plots in their original locations

The new submission helper redraws the original Xenium plot filenames under
`outputs/xenium/`, using native stage plotting code and saved objects. It does
not create a `plot_refresh/<tag>` tree. Existing original plot files are
explicitly replaced by `--submit`, as requested. Existing snapshot files are
left alone. No annotation, filtering, integration, normalization, or RDS-writing
stage is submitted.

UMAP point-size defaults, raster settings, dimensions, labels, and saved
reductions come from each native stage. The generic snapshot UMAP exporter is
not used. Updated colors and marker panels come from `color_palette.R`.
Consensus DotPlots now use the shared dynamic sizing and colored group headers,
including the per-sample PDF pages. The native consensus UMAP code is unchanged.

## Install and run on Cheaha

Transfer `xenium_original_plots_20261020.tar.gz` to the HPC scripts directory.
The archive contains code and documentation only, including the native plot
scripts used by the launcher; it does not contain configuration, data, RDS
objects, logs, or generated scientific plots.

```bash
cd /home/acflint/R/Projects/XeniumFCProject/scripts
tar -xzf xenium_original_plots_20261020.tar.gz
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b

bash scripts/submit_xenium_original_plots.sh --dry-run
# After all selected input checks pass:
bash scripts/submit_xenium_original_plots.sh --submit
```

Preflight checks paths and task mappings without loading RDS objects. Expect
1–3 minutes for the default ten groups, based on ten R startups plus mapping
reads. Every selected group is checked before any submission. Per-job native
preflights add their stage-specific checks before rendering. Path checks do
not establish that saved labels are current.

The default groups are consensus, VZ samples, RL samples, residual samples,
combined/corrected-broad samples, pooled VZ, pooled RL, pooled residual, and
VZ/RL count plots. Sample arrays use the manifest's current 41 samples at
concurrency 3. Bounds are computed, not taken from the launcher's fallback
SBATCH array. Combined reports wait for successful sample jobs; consensus
report pages and merged PDFs are also chained with `afterok`. Combined spatial
reports are regenerated at both 0.6 and 0.85 point sizes.

To run a subset, use the same group selection for preflight and submission:

```bash
bash scripts/submit_xenium_original_plots.sh --dry-run --groups=consensus,rl,residual,combined,rl-pooled,residual-pooled
bash scripts/submit_xenium_original_plots.sh --submit --groups=consensus,rl,residual,combined,rl-pooled,residual-pooled
```

This subset has six groups: its path preflight is expected to take 30–120
seconds. Do not duplicate already-running plotting jobs. Wait for any active
annotation/mapping jobs that supply the selected objects. Stale VZ objects
containing the former `Mature UBCs` identity are not relabeled by this redraw;
complete the approved upstream mapping first. A failed dependency should be
diagnosed, not released forcibly.

Processed combined objects can be redrawn separately once they exist:

```bash
bash scripts/submit_xenium_original_plots.sh --dry-run --groups=nonspatial,spatial
bash scripts/submit_xenium_original_plots.sh --submit --groups=nonspatial,spatial
```

That preflight should take 10–40 seconds. The nonspatial driver uses its existing
`--plots-only` path, with provenance checks. The spatial driver reads the saved
integrated object and stored UMAP; it does not rerun integration. Neither group
is enabled by default because these saved inputs were absent in the earlier
HPC preflight. Missing inputs are errors, not instructions to rerun processing.
Reference replacement remains covered by
`RESTORE_ORIGINAL_PLOT_FORMAT_20261020.md` and its separate reference launcher.

## Locations and scope

Native stage directories and filenames are preserved. For the six combined
sample spatial TIFFs, the helper detects the existing combined plot-family
directory and either corrected-broad RDS layout (`01_corrected_broad` or
`02_corrected_broad`), then redraws into the corresponding original plot paths.
Mixed/duplicate layouts are rejected instead of choosing silently. If no
original files exist for a family, its native default directory is used.
The `.6`/`.85` spatial point sizes and 12×10 or 14×10 canvases are preserved.
The saved spatial-integration UMAPs retain their native 14×14 canvases.

Pooled native plot stages may also refresh display marker tables and heatmaps.
They do not replace analysis RDS objects in the selected plot-only modes. QC,
reference transfer, cell identities, and scientific parameters are unchanged.
New colors affect visual interpretation, not the underlying annotations.

## Runtime and resources

Expected runtimes exclude queue wait. These are provisional ranges: prior
successful jobs inform RL/residual estimates; other ranges rely on static
inspection and the number of plots/RDS reads. Cell counts and filesystem speed
can substantially change runtime. Concurrency is per group, so simultaneous
groups can contend for resources.

| Group | Expected per task | Tasks / concurrency | Approximate group wall time | Requested limit / memory |
|---|---:|---:|---:|---:|
| Consensus sample plots | 5–20 min | 41 / 3 | 70–280 min | 2 h / 64 GB |
| VZ or RL sample plots, each | 5–30 min | 41 / 3 | 70–420 min | 2 h / 64 GB |
| Residual sample plots | 1–10 min | 41 / 3 | 14–140 min | 2 h / 64 GB |
| Combined/corrected spatial TIFFs | 3–15 min | 41 / 3 | 42–210 min | 2 h / 64 GB |
| Pooled VZ | 30–90 min | 1 / 1 | 30–90 min | 12 h / 256 GB |
| Pooled RL | 5–20 min | 1 / 1 | 5–20 min | 12 h / 256 GB |
| Pooled residual | 1–5 min | 1 / 1 | 1–5 min | 1 h / 32 GB |
| VZ or RL counts, each | 5–20 min | 1 / 1 | 5–20 min | 2 h / 128 GB |
| Consensus report pages | 5–20 min | 41 / 3 | 70–280 min, after samples | 2 h / 64 GB |
| Each merged report job | 5–30 min | 1 / 1 | 5–30 min, after inputs | 2 h / 64 GB |
| Optional nonspatial plot-only | 30–180 min | 1 / 1 | 30–180 min | 24 h / 600 GB |
| Optional saved spatial UMAPs | 10–60 min | 1 / 1 | 10–60 min | 24 h / 600 GB |

Low-confidence planning range for the default workflow including dependent
reports: roughly 2–10 hours once jobs can run, excluding queue wait. A prior
residual sample task took about 47 minutes; this unresolved outlier means its
array can exceed the typical range above. Pooled RL previously took about
9 minutes; the old snapshot pooled RL/residual tasks took about 2.6/0.8 minutes
but rendered fewer outputs. Limits are resource ceilings, not expected times.
The larger pooled/combined memory requests retain established stage envelopes.
Collect measurements after completion to improve estimates:

```bash
sacct -j JOBID --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

## Validation and follow-up

Local validation parsed changed R files, checked shell syntax and all R startup
flags, and exercised task mapping, preflight failure, duplicate-map rejection,
report dependencies, and every worker dispatch with mocked submission tools.
Synthetic Seurat objects produced six combined spatial and four saved-UMAP
exports with the expected filenames/dimensions; source RDS hashes stayed
unchanged. Local Cairo/X11 is unavailable, so these tests substituted base PDF
devices for TIFF output. Actual TIFF rendering remains an HPC check.

Local real-input dry-runs correctly fail because the HPC RDS files are absent
on this computer. No production job was submitted locally. R parsing takes
under 10 seconds per file; shell checks take under 5 seconds per file and are
syntax checks, not scientific validation.

After HPC completion, inspect the original plot files and report page counts,
check failed tasks/logs, and collect the accounting fields above. Keep the RDS
files as the rollback source for any further figure redraw. If this checkout
is tracked on the HPC, review `git status` and `git diff`, then commit this
plotting update when ready. No package-environment change was made; this does
not require `renv::snapshot()`.
