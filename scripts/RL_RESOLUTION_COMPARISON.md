# RL resolutions 0.5, 0.6, 0.7 and 0.8

`xenium_rl_05d_compare_resolutions.R` compares full-data partitions after Stage 5
without changing production objects, labels, QC, integration, neighbors or UMAP.
It preserves the saved 0.5 and 0.8 assignments and computes only 0.6 and 0.7 on
the same `Xenium_snn` graph. Saved Louvain settings, including its original
clustering seed, are read from the production object's Seurat command record.
This seed can differ from the configured marker-test seed; both are recorded.
No clustering labels are transferred between resolutions by integer cluster ID.

This is a full-data comparison, not a new subsampling-stability run. First inspect
the marker and transition evidence; assess stability of any shortlisted intermediate
partitions separately. No resolution is automatically selected or annotated.

## Transfer and run

Copy the new driver, launcher and `R/regional_resolution_comparison.R` to their
matching locations under HPC `scripts/`. Keep the already patched
`R/regional_cluster_stability.R` (explicit ggraph guide),
`R/subcluster_review_markers.R`, `R/dotplot_helpers.R` and `color_palette.R` current.
The existing focused marker panels and palette marker additions are reused.

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b

if Rscript --max-ppsize=500000 scripts/xenium_rl_05d_compare_resolutions.R --check-packages &&
   Rscript --max-ppsize=500000 scripts/xenium_rl_05d_compare_resolutions.R --list &&
   Rscript --max-ppsize=500000 scripts/xenium_rl_05d_compare_resolutions.R --dry-run
then
  sbatch scripts/run_xenium_rl_05d_compare_resolutions.slurm
fi
```

Preflight should take under one minute per check. `--check-packages` renders the
same small synthetic clustree PDF used by the stability review and removes the
temporary file. `--list` and `--dry-run` inspect manifest/configuration/paths only;
they do not load an RDS. Actual job startup validates the provenance size/time,
cell/sample IDs and order, all active manifest samples (currently 41), graph values,
UMAP coordinates and saved clustering parameters. Required target clusters
1/3/5/6/8/17 must be present in the 0.8 anchor; otherwise the review stops rather
than silently analyzing different populations.

Expected HPC runtime is **2–8 hours** for this single pooled job, excluding queue
wait (low confidence). Basis: RL pilot repeats required 40–46 minutes for rebuilding
an 80% graph and three partitions; this comparison reuses the full saved graph,
computes two partitions, then runs four marker/DotPlot reviews. Earlier RL Stage-6
combined review completed in 12m12s with 80.32 GiB MaxRSS, but evidence reuse was
not confirmed. The launcher retains the established RL sensitivity allocation of
**12 hours, one CPU, 256 GB, short partition**, with sequential execution.
No package installation or renv restore runs inside the job. If dependencies are
missing, review `renv::status()` and perform an explicit project-environment setup;
snapshot only after reviewing and approving actual environment changes.

After completion, collect resource evidence (normally under one minute):

```bash
sacct -j <JOB_ID> --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

## Outputs

Default directory:

```text
/home/acflint/R/Projects/XeniumFCProject/outputs/xenium/rl/05_post_qc/diagnostics/resolution_comparison_05_06_07_08/review1/
```

- `clustree.pdf`: full-data membership relationships across the four resolutions.
- `res<resolution>_UMAP.pdf`: identical saved coordinates with each partition.
- `res<resolution>_cluster6_UMAP.pdf`: only original cluster-6 cells highlighted by
  their destination partition; caption says **identity unresolved**.
- `res<resolution>_top50_markers.csv`: descriptive positive-marker rankings.
- `res<resolution>_focused_markers.pdf` and `.tif`: established focused panels,
  including RL VZ, RL SVZ/UBC, glia and competing lineages.
- `res<resolution>_marker_availability.csv`: genes absent from the targeted panel.
- `res<resolution>_marker_status.csv`: marker counts and explicit not-testable
  status for clusters with fewer than three cells. These cells are never removed.
- `production0.8_transitions.csv`: exact cell-count transitions from every current
  0.8 cluster into 0.5, 0.6 and 0.7, with source-retention and destination-composition
  fractions. These quantify mergers without inferring identities.
- `focus_population_transitions.csv`: the same comparisons restricted to current
  clusters 1/3/5 (UBC review), 6 (unresolved), and 8/17 (RL SVZ review).
- `all_pair_transitions.csv`: transitions for all six pairs of resolutions.
- `assignments.csv.gz`: unique cell IDs, sample IDs and all four partitions.
- `cluster_summary.csv`, `sample_composition.csv`, `input_label_composition.csv`:
  counts, sample coverage and pooled input-label composition, retaining zeros.
- `provenance.csv`, `sessionInfo.txt`, `METHODS.txt`, `outputs.csv`: provenance,
  parameters, package versions and complete-file sizes/checksums.

Marker settings match the current RL Stage-5b diagnostic: Wilcoxon,
`only.pos=TRUE`, `min.pct=0.25`, `logfc.threshold=0.25`, and up to 1000 cells per
identity, using the configured marker seed. These are descriptive cell-level
rankings; P-values do not establish biological-replicate significance. DotPlot
colours are scaled within each plotted dataset, so inspect gene combinations and
expression fractions instead of interpreting between-plot colour differences as
absolute expression changes. No RL VZ or BG label is newly assigned.

Existing output directories are protected. A fresh run uses, for example,
`sbatch --export=ALL,RL_RESOLUTION_COMPARE_RUN_ID=review2 scripts/run_xenium_rl_05d_compare_resolutions.slurm`.
Explicit reviewed replacement uses `RL_RESOLUTION_COMPARE_OVERWRITE=true`; the old
diagnostic directory is archived first. Partial files remain in `.staging_*` and
the complete directory is published only after all expected outputs exist.
No production RDS is written. If a new production resolution is later approved,
review its Stage-6 label map before annotation and regenerate affected Stage-7
mappings and dependent combined, spatial, count and export outputs.

## Lightweight development checks

Parse every changed R file (under 10 seconds each). The synthetic test takes about
5–30 seconds and checks transition fractions, duplicate IDs, graph preservation and
deterministic clustering on 100 artificial cells; no project analysis object is read.
Shell syntax takes under 5 seconds per file; the launcher audit under one minute.

```bash
Rscript --max-ppsize=500000 -e 'for (f in c("scripts/R/regional_resolution_comparison.R", "scripts/xenium_rl_05d_compare_resolutions.R", "scripts/test_regional_resolution_comparison.R")) parse(f)'
Rscript --max-ppsize=500000 scripts/test_regional_resolution_comparison.R
bash -n scripts/run_xenium_rl_05d_compare_resolutions.slurm
cd scripts
python3 update_r_launchers_stack_limit.py --check
```

Before a Git checkpoint in the actual repository, review `git status` and `git diff`
and the checks, then commit the logical unit when ready. This local scripts folder
is not a Git checkout; no commit, push, package installation or lockfile update is
performed by this addition.
