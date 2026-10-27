# Review corrections — 22 October 2026

The corrections are saved locally. Regression tests were rerun successfully on 26 October 2026.

- The submission dependency test now requires residual QC (when removing cells) -> residual post-QC processing -> stability -> marker review, including the explicit no-removal branch.
- RL marker review generates a blank label-approval template. Label application requires an explicitly approved copy at `scripts/reviews/RL_subcluster_labels_Res0.8.csv`, bound to the source file metadata, review manifest and ordered cell-to-cluster fingerprint. Every cluster needs a label and notes matching the configured map. Regenerating evidence alone cannot approve historical labels for a new fit. Existing plot-only refresh bypasses label application.
- Workflow documentation now correctly describes residual refitting after approved removals and representation reuse for an explicitly reviewed no-removal decision.

Existing labels, cell/sample IDs, QC thresholds, clustering settings and output RDS names were preserved. Fresh VZ/RL/residual biological label decisions remain manual; no new mapping or approval was invented.

## Validation

All three changed/new R files parsed successfully. Synthetic RL approval tests passed, including the actual driver's path-only preflights using invalid RDS placeholders to establish that inspection does not load objects. Changed partitions with identical cell totals/cluster IDs, changed cell IDs, stale source/review fingerprints, missing/duplicate decisions, unapproved rows and different configured labels were rejected. Existing plot-only inspection passed without an approved label table. The updated submission dependency test passed for both residual QC choices and downstream merge/validation dependencies.

Structural configuration validation passed for 41 samples and 7 multi-sample slides. The global stack safeguard passed for 172 shell files / 299 Rscript calls during the correction. No production launcher resource settings were changed.

Production RL review/apply dry-runs ran from the project root and failed on unavailable local HPC inputs, as expected. Apply also reported the missing explicit label approval. They inspected paths and small tables only; no production RDS was loaded. HPC preflight and full-data scientific validation remain necessary.

Runtime estimates: parsing under 10 seconds per file; synthetic approval tests under 10 seconds; dependency tests under one second; structural config under one minute; production path-only preflight under one minute plus environment startup. Fingerprinting the already loaded ~827k-cell assignments is estimated to add under 10 seconds to label application (low confidence, static inspection); it adds no RDS read or refitting. Existing RL review/apply Slurm limits remain 12 hours. Full RL review is roughly 0.5–3 hours and label application 10–60 minutes, low confidence from related job history; queue waits are excluded.

No production outputs were replaced, jobs submitted, packages installed, lockfile updated, or Git commit/push performed. These are local changes; the Cheaha copies must be synchronized before use. On Cheaha, review `git status` and `git diff`, repeat the relevant lightweight checks, then commit this correction as a separate logical unit when ready. No renv snapshot is required for these code-only changes.
