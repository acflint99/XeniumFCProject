# Vascular reference and downstream label update

User-approved convention: **Vascular = Endothelial + Pericytes** at the broad-label level. The former standalone Endothelial identity is removed from current downstream configuration, taxonomy, corrected broad mappings, plotting orders/colors, count groups, and marker-panel headings. Pericytes remains an approved fine residual subcluster identity; its corrected broad label is now Vascular.

Aldinger already maps the original `13-Endothelial` and `16-Pericytes` source keys to Vascular. Sepp already maps its combined mural/endothelial source population to Vascular. Science standardization now explicitly harmonizes older source aliases into `clusters_refined=Vascular`, preserving source metadata and cells. Its supported legacy numeric mapping now also emits Vascular. The Science reviewed-label validator continues to reject obsolete review labels instead of accepting them as new biological identities.

Remaining literal Endothelial strings in active code are intentional ingestion keys, rejection checks, and synthetic test fixtures. Original reference metadata and scientific review notes are retained as provenance. Historical folders and the out-of-scope Giotto script are untouched. Original-reference validation still reports native published annotations; these are source evidence, not assigned Xenium labels.

The Vascular broad DotPlot marker group combines the existing gene sets: CLDN5, PECAM1, PDGFRB, ABCC9 and KCNJ8. No new marker genes were introduced. Former Endothelial contrast/plot headings are renamed Vascular without changing their genes. The existing red Vascular color is retained, while Pericytes keeps its fine-subcluster color.

## Scientific effect and outputs requiring regeneration

The corrected broad mapping changes Pericytes from a standalone broad category to Vascular. Fine subcluster identities, cell/sample IDs, QC filters, clustering/integration settings, reference sampling, transfer thresholds, consensus votes and routing categories are unchanged. Vascular stays nonregional/Other, so this name change does not itself move cells into VZ/RL.

Existing saved objects are not edited by changing these scripts. Rebuild standardized reference/panel outputs when their assigned labels change. Rebuild affected label transfers and consensus if a reference used for annotation changes; then revisit downstream routing/QC/subcluster decisions against the fresh source rather than reusing numeric IDs automatically. The existing refresh helper provides those manual gates. Already canonical Vascular reference objects do not need recomputing solely for this rename.

Regenerate corrected broad objects/counts, both merges and their dependent processing, figures and reference-consistency reports to adopt the new Pericytes-to-Vascular broad grouping. Update RL plots/counts if their source identities change. Existing combined objects containing the obsolete standalone Endothelial fine identity are rejected by the corrected mapping rather than silently creating that broad group; regenerate their upstream reviewed labels/mappings first. Retain old outputs as a rollback snapshot and review protected overwrite switches before replacement. No RDS, plot or analysis output was replaced during implementation.

## Validation performed

- Parsed all 12 changed/new R files; expected under 10 seconds each, completed successfully.
- `test_vascular_labels.R` passed: source-label preservation, harmonization, no duplicate palette/marker names, valid broad/fine identity distinction, expected genes/colors, and rejection of obsolete downstream labels. Runtime under 10 seconds; no RDS loaded.
- `validate_config.R` passed: 41 biological samples, structural checks only, under one minute.
- Broad taxonomy/config consistency check passed, under 10 seconds.
- Science standardization and corrected-broad `--dry-run` checks ran. They inspect paths only and do not load RDS files; both correctly reported missing HPC inputs locally. Science's existing driver prints DRY-RUN FAIL even though its legacy process exit code is zero; check the printed status as the submission helper does.
- No package changes, full analysis runs, Slurm submissions, commits, pushes or renv lockfile changes.

## Runtime expectations before any HPC rebuild

These estimates exclude queue wait and human review. They are expected runtimes, not requested Slurm limits. Use the established HPC R 4.4.1 / renv environment. Local parsing and synthetic checks used existing R 4.5.1 and do not replace HPC validation.

| Stage | Expected runtime | Basis/confidence |
|---|---|---|
| Any relevant path-only preflight | Under 1 minute per driver | Static inspection; no large RDS load |
| Structural config validation | Under 1 minute | Small YAML/CSV checks |
| `validate_config.R --check-files` on HPC | Approximately 1–5 minutes | Filesystem-dependent, low confidence |
| Reference standardization/panel preparation if needed | 5–60 minutes per stage | Existing refresh guide, static inspection, low confidence |
| Label transfer if reference labels changed | 8–50 minutes/task; 41 tasks/reference at concurrency 3: approximately 2–12 hours/reference | Prior arrays and refresh guide; three reference arrays may overlap |
| Corrected broad labels | 10–30 minutes for the pooled 41-sample job | Prior 12.93-minute job, low–moderate confidence |
| Nonspatial / spatial merge | 0.5–4 hours / 1–12 hours per pooled job | Existing guide, low confidence |

Other downstream stage estimates and requested limits are documented in `UPDATED_REFERENCES_TO_MERGES.md` and `RESIDUAL_POST_QC_PROCESSING.md`. Collect JobID, Elapsed, State, MaxRSS, ReqMem, AllocCPUS and ExitCode with sacct after jobs complete to refine estimates.

Review git status and git diff in the tracked checkout, rerun HPC preflights, and checkpoint this logical unit with a commit when ready. This local scripts directory is not a Git checkout. No renv snapshot is needed solely for label/code changes; no package environment was changed.
