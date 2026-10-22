# Nonspatial plot refresh — 2026-10-21

The plotting driver now generates 24 figures, a provenance manifest, and a heatmap marker summary from the
existing processed object. Sync these files together to the HPC project:

- `scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R`
- `scripts/R/nonspatial_merge_plots.R` (new)
- `scripts/color_palette.R`
- `scripts/R/dotplot_helpers.R` (existing helper, included for compatibility)

The existing `scripts/run_xenium_vz_rl_05_nonspatial_plots.slurm` is unchanged.
Do not rerun processing or integration for this refresh. No packages or package
versions changed. Production continues to use its existing R 4.4.1/renv setup.

## Figures

Every family has a version with Unknowns and a NoUnknown version:

| Family | Grouping | Files per pair |
| --- | --- | --- |
| Batch comparison, pre/post Harmony | Sample; one shared bottom key | 2 TIFF |
| Original broad UMAP | `consensus_label` | 2 TIFF |
| Corrected broad UMAP | `corrected_broad_label` | 2 TIFF |
| Subcluster UMAP | `comb_subcluster` | 2 TIFF |
| Broad marker dotplot | Corrected broad | 2 TIFF + 2 PDF |
| Purkinje marker dotplot | Four established PC subclusters | 2 TIFF + 2 PDF |
| Top-five-marker heatmap | Corrected broad | 2 TIFF + 2 PDF |
| Signaling-gene violins | Selected corrected broad groups | 2 TIFF + 2 PDF |

The previous consensus-label dotplot used the original assignments. It is no
longer generated; use the corrected-broad dotplot. Original assignments remain
available in the original broad UMAP.

Broad display order is RL, UBC, Granule, eCN, VZ, Purkinje, GABA, Glia, OPC,
Meninges, Endothelial, Pericytes, Immune, Cycling Cells, Unresolved, Unknown(s).
Subcluster keys use this broad-family sequence and the existing within-family
order. This changes display ordering only; colors and biological labels are
preserved. Dotplot section fills match the registered broad colors.
The heatmap's redundant Identity legend is removed, while its Z-score key
remains. UMAPs retain ordered cluster keys and repel labels.

Unknown removal follows the displayed UMAP label column. Other cohort panels
use corrected broad Unknown status. Unresolved remains visible. Purkinje and
selected signaling panels already contain no Unknown identities, so both
versions are intentionally identical. The manifest records these exceptions.

The heatmap and signaling violins now use corrected broad groups rather than
original consensus groups. Heatmap marker-test settings are unchanged. Marker
tests run once on the complete object; both heatmaps share known marker genes,
full-cohort scaling, and sampled known cells. This refresh affects plotted
summaries only; it does not save cell filtering, labels, normalization or
embeddings back into the RDS. Earlier figures remain a visual rollback.

## Unresolved-marker correction

Job 41006727 failed because Unresolved had no qualifying positive markers and
the first version incorrectly required markers for every displayed group.
Unresolved is a heterogeneous retained category: zero qualifying markers is now
an explicitly supported plotting outcome for this category only. Its cells
remain in both applicable heatmaps and are shown using genes selected for other
groups; no gene block or marker evidence is invented for Unresolved. A visible
subtitle, the plot manifest notes, stderr, and
`Xenium_Merged_Heatmap_marker_summary.csv` report the absence. Other groups with
missing markers still fail, and an entirely empty marker result still fails.
Marker thresholds, identities, input RDS and all cell-selection rules are
unchanged. The auxiliary summary is protected against overwriting like figures.
Use the new v2 output folder below to preserve the failed run's partial exports.
The `renv` out-of-sync startup notice did not cause this explicit plotting error;
no package or lockfile changes are part of this fix.

## Run on Cheaha

From `/home/acflint/R/Projects/XeniumFCProject`, first inspect the input paths
and new output destination. Estimated runtime: under one minute; path checks
only, no RDS loading.

```bash
COMBINED_PLOT_OUTPUT_DIR=outputs/xenium/vz_rl/nonspatial/05_plots_review_20261021_v2 \
  Rscript --max-ppsize=500000 scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R \
  --plots-only --dry-run
```

If using an archived merge/processed run, preserve its existing
`COMBINED_MERGE_INPUT_DIR` and `COMBINED_NONSPATIAL_OUTPUT_DIR` exports as well.
After the dry-run reports PASS, submit the plotting launcher to a fresh review
folder. Existing outputs remain protected; choose another fresh folder if this
one already contains results.

```bash
sbatch --export=ALL,COMBINED_PLOTS_OVERWRITE=false,COMBINED_PLOT_OUTPUT_DIR=outputs/xenium/vz_rl/nonspatial/05_plots_review_20261021_v2 \
  scripts/run_xenium_vz_rl_05_nonspatial_plots.slurm
```

Estimated runtime: **30–90 minutes**, excluding queue wait, with low confidence.
Basis: job 41006727 reached the end of marker testing in 24 minutes 53 seconds
and then stopped on the zero-marker Unresolved check. Remaining heatmap and
violin rendering time is not yet measured. Its peak recorded batch RSS was
182.34 GiB; this is partial-run evidence, not the final successful job peak.
The single job has no array. Its existing request remains 24 hours,
600 GB, one CPU, medium partition; no resources were changed. Collect measured
runtime and memory after completion:

```bash
sacct -j JOB_ID --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Check that the job completed, all 24 files are nonempty, and the plot manifest
contains 24 unique filenames. Check the separate heatmap marker summary for
any explicitly retained zero-marker Unresolved group. Review UMAP label placement on the real cohort,
the 41-sample batch key, broad order, section colors, and absence of Unknowns
in each appropriate NoUnknown panel. The default protections prevent accidental
replacement; promoting reviewed figures to an existing destination requires
explicit overwrite review.

## Checks completed locally

- The plotting helper, driver, and lightweight regression script parse successfully
  (under ten seconds each). `test_nonspatial_heatmap_markers.R` verifies the
  zero-marker Unresolved exception and that missing named cell types still fail;
  runtime is under ten seconds, using in-memory marker tables only.
- Synthetic rendering with absent Unresolved marker results exercises all 24
  outputs, retained Unresolved cells, paired cell filtering, strict
  broad ordering, colored dotplot strips, and output completeness. Native TIFF
  and PDF devices substitute only for the unavailable local Cairo devices;
  production Cairo export still needs validation on the HPC. Synthetic checks
  take approximately 10–30 seconds here and do not use biological data.
- The path-only plotting dry-run ran and reported missing HPC input files;
  it could not pass locally because the processed object is not available.
- The existing plotting launcher passes shell syntax validation. All 159
  current shell files/267 Rscript invocations pass the stack-limit checker.
- Local R is 4.5.1, whereas the unchanged HPC launcher uses R 4.4.1. Local
  checks are lightweight compatibility checks, not production validation.

No analysis RDS or supplied plot was overwritten. A failed local Cairo check
unexpectedly replaced the incidental `scripts/Rplots.pdf`; no backup of that
file was available. Further rendering was confined to a temporary directory.

This local scripts folder is not a Git checkout. On the HPC checkout, review
`git status` and `git diff`, validate, then commit this completed plotting unit
if ready. No commit, push, package installation or lockfile update was made.
