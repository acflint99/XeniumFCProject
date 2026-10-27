# Script and diagram alignment — 22 October 2026

The active workflow now includes clustree/subsampling stability review for VZ, RL
and residual before marker/label approval, and Spearman plus MetaNeighborUS
validation after the nonspatial merge. The refresh submitter wires these jobs
into its existing dependencies. No jobs were submitted by this change.

## Scientific behavior

QC thresholds, production resolutions, integration parameters, references,
consensus voting, identities, cell/sample IDs, metadata and stable RDS names are
unchanged. All three branches share the same review sequence; this does not
make their numerical settings identical. VZ/RL retain their existing post-QC
processing. Residual refits normalization, PCA, Harmony, UMAP, neighbors and
clusters after approved removals; an explicitly reviewed no-removal decision
reuses its representations. The filtered snapshot still retains the original
coordinates and is an intermediate input. See `RESIDUAL_POST_QC_PROCESSING.md`.
Refitting invalidates the prior residual marker evidence and numeric-cluster
label approvals. Clustree/stability is diagnostic: it
never selects a resolution, applies labels, writes RDS files or changes filtering.
Resolution 0.8 remains the reviewed annotation partition. If review leads to a
new production partition, obtain explicit approval, regenerate affected marker
and label maps, then map/combine/correct/merge and rerun dependent validation and
figures. Do not map integer cluster IDs across different partitions.

Residual checks use existing QC source binding and explicit per-sample counts.
Known zero-cell samples are included in coverage; missing nonempty samples,
unexpected samples, invalid counts, duplicate IDs, non-finite embeddings, stale
QC decisions or incomplete source provenance fail. Legacy residual integration
provenance fingerprints the original merge rather than the integrated output;
approved-removal provenance additionally binds integrated and filtered RDS files.
Source files/sidecars are checked for replacement while the diagnostic runs.

Spearman and MetaNeighborUS reuse the existing original-reference comparison
scripts. Each compares corrected broad and combined fine labels to original
Aldinger, Sepp and Science labels. The references also contributed to annotation;
these assess transcriptional consistency, not independently measured accuracy.
No new validation threshold or automatic relabeling is introduced.

## Scripts and entry points

- `xenium_regional_05c_resolution_stability.R`: common engine, now accepts Residual.
- `xenium_residual_05c_resolution_stability.R`: dedicated residual entry point.
- `R/regional_resolution_inputs.R`: common source/coverage checks.
- `xenium_residual_05_reintegrate_post_qc.R`: refit/reuse before residual stability.
- `R/rl_label_approval.R`: explicit source-bound RL label decisions, separate
  from automatically generated marker evidence.
- `run_xenium_residual_05c_resolution_stability.slurm`: residual HPC launcher.
- `submit_updated_references_to_merges.sh`: now queues all three stability jobs
  before their marker reviews; queues Spearman and MetaNeighbor after the
  nonspatial merge. Spatial merge remains a separate parallel branch.

The new residual entry point supports `--list`, path-only `--dry-run`,
`--check-packages`, `--repeats`, `--fraction`, `--run-id`, and protected `--overwrite`,
with the same defaults as VZ/RL: five repeats, 80%, configured seed.
No packages are installed by scripts or jobs.

## Runtime and resources

Expected active runtimes exclude queue waits and human review. They are estimates,
not requested time limits. All are low confidence for a refreshed cohort.

| Stage | Expected runtime | Requested limit / resources | Basis |
|---|---|---|---|
| Path-only list/dry-run | Under 1 min per driver | Existing allocation | Paths/configuration/small CSVs; no RDS reads |
| VZ clustree/stability | 4–20 h, single pooled job | Existing 36 h / 600 GB / 1 CPU / medium | Related Stage-5 history and static inspection |
| RL clustree/stability | 2–12 h, single pooled job | Existing 24 h / 600 GB / 1 CPU / medium | Related Stage-5 history and static inspection |
| Residual clustree/stability | 5–60 min, single pooled job | New 4 h / 96 GB / 1 CPU / short | Prior ~10k cells; reuse residual marker-review resource envelope |
| Spearman | 1–6 h, single job | Existing 12 h / 256 GB / 1 CPU / short | Existing comparison launcher |
| MetaNeighborUS | 0.5–4 h per reference; 1.5–12 h for three tasks at concurrency one | Existing 12 h per task / 256 GB / 1 CPU / short | Existing comparison launcher; no measured full-run history |

With the three stability jobs overlapping, allow about 4–20 hours at the prior
sizes. Spearman and MetaNeighbor can overlap after the nonspatial merge: allow
about 1.5–12 hours for both validations. Node caps/memory availability can limit
actual overlap. No sample-array bounds or sample concurrency were changed.
MetaNeighbor's separate reference map was verified as tasks 1/2/3 =
Aldinger/Sepp/Science; it is not a sample array. The authoritative manifest
currently contains 41 samples.

## Safe execution and completion review

Use the existing four-phase guide in `UPDATED_REFERENCES_TO_MERGES.md`. New
refresh runs use diagnostic run ID `<REFERENCE_REFRESH_RUN_ID>_postqc`, with
fresh-output protection enabled. The helper does not retroactively add jobs to
an older phase ledger. Existing completed analyses can use the standalone
launchers without rerunning annotation/QC/merges solely to add these diagnostics.
Production HPC R remains R/4.4.1-foss-2023b with existing renv activation.

Stability reports appear under each branch's
`05_post_qc/diagnostics/resolution_stability/<run_id>/`. Require `outputs.csv`
and all its listed nonempty files, including clustree/stability PDFs, provenance,
source inputs and coverage. Review marker and stability results together at the
existing label-approval gate. Fresh source objects still require fresh QC and
label decisions; no old removal ID/map is approved automatically.
RL marker review writes a blank resolution-specific approval template. Complete
its labels, notes and approved flags only after review, then save it as
`scripts/reviews/RL_subcluster_labels_Res0.8.csv` for production resolution 0.8.
Its provenance fields must remain unchanged; application checks the configured
map against these decisions and fingerprints the exact cell-to-cluster
assignments. Plot-only refresh does not require a new label decision.

Post-merge validation writes to fresh directories:
`outputs/reference_label_refresh_runs/<run_id>/cluster_validation/spearman/` and
`.../metaneighbor/<aldinger|sepp|science>/`. Require Spearman's `COMPLETED.txt`
and all three MetaNeighbor `COMPLETED.txt` files; inspect correlation/AUROC,
gene-coverage and best-match tables alongside scientific evidence. Their existing
drivers refuse to reuse an output directory. No existing RDS or figure output
was overwritten by this implementation.

After HPC completion, collect JobID, Elapsed, State, MaxRSS, ReqMem, AllocCPUS and
ExitCode from sacct (normally under one minute) and reconcile the job ledger.
Use measured cell counts and per-repeat timings to refine future requests.

## Local validation and limits

R parsing is expected under ten seconds per changed R file and does not establish
scientific correctness. Synthetic removal/mapping and coverage fixtures take
under one minute; dependency-graph tests take under one second and never submit
jobs. Shell syntax checking takes under five seconds per launcher. Structural
configuration validation takes under one minute. Local R is the existing 4.5.1
installation; production 4.4.1 remains unchanged and needs HPC validation.

Local production dry-runs report missing HPC RDS/reference/raw-panel paths. They
inspect paths or small tables and do not load large RDS objects. Clustree/ggraph
and MetaNeighbor are absent locally; no installation or lockfile update was
performed. On HPC, run the existing `--check-packages`/`--dry-run` checks before
submission (usually under one minute each). If dependencies are missing there,
use an explicit project renv maintenance session, review `renv::status()`, install
only the missing packages, and snapshot only after approving environment changes.
Allow 2–20 minutes for package setup, low confidence and potentially longer with
compilation/network latency. Do not install or restore inside production jobs.

This local scripts directory is not a Git checkout. In the tracked checkout,
review `git status` and `git diff`, repeat relevant validation and checkpoint this
logical unit when ready. No commit, push or renv.lock change was performed.

### Results from this implementation

The original checks below describe the earlier alignment revision. The
22 October correction adds residual post-QC dependency assertions and RL
source-bound approval tests. Current correction results are recorded in
`REVIEW_CORRECTIONS_20261022.md`; the earlier shell-call totals below are
historical and must not be read as a new validation run.

- All four changed/new R scripts parsed successfully.
- Synthetic residual removal/marker/label/mapping checks passed, including the new
  input adapter for no-removal/removal states and explicit empty-sample coverage.
- Shared stability tests passed; synthetic PDF rendering was skipped because
  clustree/ggraph are absent locally.
- Submission dependency tests passed for both residual QC choices, all three
  stability-to-marker dependencies and both merge-to-validation dependencies.
- Shell syntax checks passed for the submitter and all relevant launchers.
- Global stack-limit check passed: 171 shell files / 297 Rscript invocations.
- Structural configuration validation passed: 41 samples, seven multi-sample
  slide directories. No config/manifest edits were needed.
- Supported path-only dry-runs ran for VZ/RL/residual stability and both
  reference validators. They correctly failed on unavailable HPC inputs;
  MetaNeighbor additionally failed its local package/API check.
- Residual package check confirmed missing clustree/ggraph. No dependency
  installation or full-data execution was attempted.
- The updated diagram was rendered and visually inspected.
