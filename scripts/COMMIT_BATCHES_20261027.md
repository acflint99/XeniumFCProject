# Commit the changes shown in the 27 October Cheaha status

The local scripts folder is not a Git checkout. The helper groups explicit files
from the supplied status, informed by the current local sources and workflow
notes. Exact Cheaha diffs against HEAD have **not** been inspected. No commits
have been created. Review each batch's actual diff on Cheaha before committing.
Whole-file grouping cannot isolate unrelated hunks in shared files; if a batch
contains unrelated changes, use `git add -p` and commit that batch manually.

Helper validation passed locally: `bash -n`; 62 unique existing source paths;
preview without staging; untracked-file diff display; refusal of an occupied
index; seven commits in a disposable fixture repository; rerun skipping; and
preservation of an unrelated untracked file. Fixture checks took about eight
seconds. These checks validate the batching mechanics, not the analysis code or
the Cheaha diffs.

Copy `commit_batches_20261027.sh` into `scripts/` in the Cheaha checkout using
your usual synchronization method. Run from the `XeniumFCProject` root:

```bash
bash scripts/commit_batches_20261027.sh
bash scripts/commit_batches_20261027.sh --diff --batch 1
# After inspecting the diff and completing the relevant checks:
bash scripts/commit_batches_20261027.sh --commit --batch 1
```

Repeat the review/commit pair for batches 2 through 7 in order. A full reviewed
set can also be committed with `bash scripts/commit_batches_20261027.sh --commit`.
Estimated Git runtime is seconds to two minutes for all batches, excluding hooks
and human review, based on about sixty explicit paths. Preview does not stage.
Commit mode refuses a nonempty index, skips already committed files, checks
whitespace, and stops on errors. A failed staged whitespace check or commit hook
can leave the current batch staged; inspect `git diff --cached` before proceeding.
Stop editing these files during execution. Nothing is pushed or submitted.

| Batch | Logical unit | Review focus |
|---|---|---|
| 1 | Vascular reference/broad labels, palette, markers and counts | Scientific identity changes, reviewed Science CSV, shared config changes |
| 2 | Residual post-QC processing and downstream provenance | Cell exclusions, refit/reuse choice, ID preservation, source-bound labels and mapping |
| 3 | Regional/residual resolution diagnostics | Input coverage, source selection, diagnostic-only behavior |
| 4 | RL approval and selected UMAP handling | Approval binding, partition fingerprints, labels versus plot-only behavior |
| 5 | Reference similarity and MetaNeighbor | Original labels, grouping, output protection and reference task mapping |
| 6 | Refresh dependencies and workflow documentation | QC → post-QC → stability → review → mapping, post-merge comparisons |
| 7 | Reference preparation resources and recovery utilities | Requested resources, fixed run IDs/job IDs and hard-coded HPC paths |

These are related checkpoints, not a guarantee that every intermediate commit
works independently: shared config/marker/documentation files can contain changes
used by later groups. Batch 1 includes the explicitly referenced Science review
CSV, rather than recursively adding `scripts/reviews/`. Inspect its decisions
against the actual source before committing. No new biological approval is made
by this helper.

## Checks before committing on Cheaha

Use the established R/4.4.1-foss-2023b and project renv environment. Syntax checks
are lightweight and do not establish scientific correctness. Expected times:
R parsing under 10 seconds/file; structural config check under one minute;
shell parsing under 5 seconds/file; global stack check a few seconds. These
estimates are based on static inspection and the project's documented checks.

For each selected R file, use this pattern (substitute its explicit path):

```bash
Rscript --max-ppsize=500000 -e 'parse(file="scripts/science_01_standardize_cell_types.R"); invisible(NULL)'
Rscript --max-ppsize=500000 scripts/validate_config.R
python3 scripts/update_r_launchers_stack_limit.py --check
bash -n scripts/submit_updated_references_to_merges.sh
```

Parse every R file being committed, including newly added tests. Run `bash -n`
on every selected `.slurm`/`.sh` file and inspect all SBATCH directives and R
commands. The example is a pattern, not a claim these checks passed on Cheaha.
Use `validate_config.R --check-files` there when inputs are available (roughly
1–5 minutes, low confidence, filesystem-dependent). Verify the actual manifest;
local documentation reports 41 samples despite the historical 34-sample note.

Relevant existing synthetic tests are `test_vascular_labels.R` (under 10 seconds),
`test_residual_post_qc_processing.R` (20–90 seconds),
`test_xenium_residual_qc_removal.R` (under one minute), and
`test_rl_label_approval.R` (under 10 seconds). Run from the project root with
`Rscript --max-ppsize=500000 scripts/<test-name>`. They use synthetic fixtures,
not production RDS objects. `python3 scripts/test_workflow_dependencies.py`
takes under a second and records simulated submissions without Slurm calls.
Estimates follow the existing validation notes; dependencies must already exist.

Before any later production submission, use the documented supported path-only
`--dry-run`/`--list` checks for the affected stages (usually under one minute per
driver; inspect printed failures as well as exit codes). The configuration,
residual-processing and workflow guides document the relevant options. Do not
rerun full analyses just to make these commits.

## Deliberately left unstaged

- `config/backups/`, `plot_backups/`, `git_commit_batches_20261022/`,
  `scripts/legacy_vz_rl_names_20261019/` and all `scripts/._*` files.
- `scripts/color_palette_before_vz_marker_prune_20261020.R` and
  `scripts/color_palette_vz_markers_pending.R`: backups/drafts.
- `reviewed_unknown_qc_APPROVED.csv`: inspect against the configured authoritative
  decision file before deciding whether and where it belongs; never infer approval
  solely from its filename.
- Root-level `test_residual_post_qc_processing.R` and
  `test_xenium_residual_qc_removal.R`: compare with `scripts/` copies first.
- `scripts/subset_purk+gran.R` and `scripts/run_subset_purk+gran.slurm`: absent
  locally and therefore unreviewed. Review together as a separate eighth commit;
  check identities, outputs, overwrite protections, R parsing and shell syntax.
- Any other files inside `scripts/reviews/` and the new batching helper itself.

The dated recovery utilities in batch 7 contain fixed run/job IDs and can modify
job dependencies or ledgers when executed. Committing them records provenance;
the batching helper never executes them. Their absolute paths are intentional
run-specific details and should be reviewed as such.

Scientific consequences documented in the existing notes: Vascular broad labels
can stale corrected broad objects, merges and figures; changed reference labels
can also stale transfer/consensus and dependent reviews. Residual refitting
invalidates old marker evidence and numeric-cluster approvals, then dependent
mapping/merges/comparisons. Keep old outputs for rollback and review fresh
source-bound decisions before regeneration. Committing code does not regenerate
or replace any output. No package environment change or renv snapshot is needed
solely for these code changes.
