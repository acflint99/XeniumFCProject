# Residual post-QC processing aligned with VZ/RL

The residual branch now follows the same stage order as VZ/RL:

**Subset → merge → integration → QC review/application → post-QC processing → clustree/stability review → marker/label review → label application → mapping.**

After approved exclusions, the new `xenium_residual_05_reintegrate_post_qc.R` repeats normalization, variable-gene selection, scaling, PCA, Harmony, UMAP, neighbor construction and clustering on retained cells. An explicitly reviewed no-removal decision uses `--no-removals` to reuse the original reductions, graph and partitions. The removal stage itself still produces a filtered snapshot with original coordinates; that snapshot is now an intermediate input, not the final analysis source.

## Scientific settings and provenance

This adds the previously absent residual refitting step. It preserves the established residual settings: 2,000 variable features, configured 30 dimensions and 30 neighbors, Harmony by `orig.ident`, resolutions 0.5/0.8/1 and current reviewed-label resolution 0.8. The configured stochastic seed is set before processing; the existing residual `FindClusters` defaults, including its command-recorded seed, are retained. VZ/RL parameters, QC decisions, consensus/routing rules, and reference labels are unchanged. No configuration or approval manifest is edited automatically.

All original metadata values and cell/sample IDs are checked for exact preservation. Post-QC cluster columns use `residual_postQC_Xenium_snn_res.<resolution>`; original `Xenium_snn_res.*` values remain available for provenance. The graph and reductions describe the new fit when removals occurred. Do not use the old cluster columns to interpret the new graph; label review, clustree and pooled plotting explicitly use the new columns.

New protected outputs:

- `outputs/xenium/residual/05_post_qc/rds/Xenium_Residual_postQC.rds`
- `outputs/xenium/residual/05_post_qc/tables/Xenium_Residual_postQC_manifest.csv`
- `outputs/xenium/residual/05_post_qc/tables/Xenium_Residual_postQC_sample_counts.csv`
- Three post-QC UMAP TIFFs under `05_post_qc/plots/`.
- Fresh marker evidence and an unapproved label template under `05_post_qc/marker_review/tables/`.

Post-QC provenance binds the source file, QC decision, sample manifest, output file, seed and coverage table. Sample counts cover every manifest row, including explicitly empty residual samples. Partial or stale inputs fail before downstream analysis. Output replacement requires the existing protected `--overwrite` pattern.

## HPC execution

Use the existing project R 4.4.1 / renv environment. No package is installed or upgraded by these scripts. From the project root, for an already applied approved removal:

```bash
Rscript --max-ppsize=500000 scripts/xenium_residual_05_reintegrate_post_qc.R --dry-run
sbatch --export=ALL,RESIDUAL_QC_ACTION=remove scripts/run_xenium_residual_05_reintegrate_post_qc.slurm
```

The dry-run inspects paths only and should take under one minute. Expected processing is **15–120 minutes** for a residual population near the prior 10,000 cells (low confidence; static inspection, no measured post-QC integration job). This is a single pooled job, not an array. The launcher reuses the existing residual integration envelope: **12-hour requested limit, 128 GB, four CPUs, short**, with one R worker. These are conservative inherited requests, not measured requirements.

For an explicitly reviewed no-removal decision, the removal-decision file and prior filtered outputs must first be reconciled/archived through the existing reviewed refresh workflow. `--no-removals` refuses to run when an approved removal decision is still active; it never ignores that decision. Then use the same dry-run with `--no-removals` and submit with `RESIDUAL_QC_ACTION=none`. Expected reuse runtime is **1–10 minutes**, low confidence from small-object I/O and plotting. It retains the same requested resource envelope.

The existing refresh helper now automatically schedules residual post-QC processing between QC application and clustree/stability. Its VZ/RL and reference chains are unchanged. Do not blindly resubmit an already recorded phase; follow the existing ledger recovery instructions or run just this new residual stage after reviewing completed inputs.

After the new stage completes, use the existing residual clustree and marker-review launchers. Clustree/stability is estimated at **5–60 minutes**, five 80% repeats on ~10,000 cells; its requested limit remains four hours. Fresh marker review is estimated at **5–60 minutes**, low confidence from static inspection, with its existing four-hour limit. Those estimates exclude queue wait. Path-only preflights each take under one minute and do not load RDS files.

## Fresh label approval is mandatory

Refitting can change cluster boundaries and integer cluster IDs even though the retained cell set stays fixed. The generated label template now carries `source_input_md5` and `source_cluster_column`. Populate and approve that fresh template in the active `config/reviewed_residual_subcluster_labels.csv`; do not copy old labels by matching cluster number. Label application rejects a review bound to any other input, even if its integer IDs match. Check markers, sample support, QC and spatial evidence before approval.

After approval, regenerate residual labels and all sample mappings. Existing stable label/mapping output names are retained, so replacement must be reviewed and explicitly enabled. Label application is estimated at **1–10 minutes**; mapping all manifest-defined samples is one pooled job estimated at **10–40 minutes**, based partly on the prior 22.1-minute residual mapping job. Regenerate pooled/per-sample residual reports as needed. The combined-label stage checks that mapped residual objects carry the current post-QC manifest fingerprint.

Every previous residual label object, mapping, combined/corrected-broad object, both merges, and dependent processed/spatial/plot/count outputs becomes stale for this new workflow. VZ/RL analysis does not need repeating solely because residual processing changed. Preserve old outputs and their source snapshot for rollback; reverting scripts alone does not make old labels valid on new clusters. Existing source/filtered objects are preserved, and no analysis output was changed during implementation.

After completed HPC jobs, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`, `AllocCPUS`, and `ExitCode` with `sacct` to refine the estimates.

## Validation and checkpoint

Local validation uses the already installed R 4.5.1 for syntax and tiny synthetic fixtures; it does not establish compatibility or scientific validity of the production R 4.4.1 run. All changed R files are parsed (normally under 10 seconds each). The synthetic representation check is expected to take 20–90 seconds; the QC-to-label/mapping fixture check is expected under one minute. Neither reads production RDS objects. Launcher syntax checks take under five seconds/file. The stack-limit checker verifies all current launcher calls.

Production path-only dry-runs cannot pass locally because the HPC analysis objects are absent. Re-run them in the established HPC environment before submission. Clustree/ggraph rendering also needs that environment; those packages were not available in the local library and were not installed.

Review `git status` and `git diff` in the tracked checkout after validation, then commit this logical unit when ready. This local scripts directory is not a Git checkout. No commit, push, package installation, or lockfile update was performed; no renv snapshot is needed solely for these code changes.

Implementation checks completed: all 12 changed/new R sources parsed; both synthetic test suites passed; both shell files passed `bash -n`; the launcher checker passed across 172 shell files and 299 Rscript calls; the authoritative manifest contains 41 unique sample IDs. Path-only production preflights were run and correctly failed for missing HPC inputs; the residual diagnostic `--list` resolved the new post-QC source. Local clustree/ggraph availability checks confirmed both packages are absent, so actual clustree rendering remains an HPC check. The updated SVG/PNG was rendered and visually inspected.
