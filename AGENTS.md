# Project Instructions

## Project scope

- Treat `README.md` as the primary guide to workflow order, expected inputs, and expected outputs, while verifying its claims against the current script before making a change.
- Work on current scripts by default. Treat `OLD/` as read-only historical reference unless the user explicitly requests work there.
- Giotto and SpaTrack are outside the current work scope unless the user explicitly brings them back into scope.
- Prefer focused changes that preserve the established pipeline structure and stable output names.

## Scientific safeguards

- Treat `config/samples.csv` as the authoritative sample manifest. It currently contains 34 biological samples; verify all array bounds and task mappings against the manifest rather than relying only on a hard-coded count.
- Treat `config/config.yml` as the source of truth for stages that have already been migrated to configuration. Do not force unmigrated stages into it without checking the intended scientific workflow.
- Do not change QC thresholds, clustering resolution, dimensionality, integration settings, reference paths or labels, label-transfer thresholds, consensus voting, or VZ/RL identity definitions without explicit user approval.
- Preserve cell IDs, sample IDs, and existing metadata columns. Validate joins for one-to-one or otherwise explicitly intended cardinality before writing outputs.
- Explain the likely scientific effect of any proposed parameter change, including which downstream objects would need to be regenerated.

## Data and output safety

- Use `--list` and/or `--dry-run` before expensive or multi-sample stages when the script supports them.
- Never overwrite, delete, or replace existing RDS files, plots, manifests, or analysis outputs without explicit user authorization. Prefer the script's protected `--overwrite` mechanism after review.
- Do not add raw Xenium data, large RDS objects, generated outputs, or Slurm logs to Git unless the user explicitly requests it and confirms the repository policy.
- Do not modify files under `OLD/` as part of cleanup or refactoring unless explicitly requested.

## R and reproducibility conventions

- Use `here::here()` or the project configuration utilities for new project paths instead of adding absolute local or HPC paths.
- Use the existing project `renv` environment and the R version established by the active Slurm launchers. Do not upgrade R or packages implicitly.
- Do not install packages from analysis scripts or Slurm jobs. Report missing dependencies and propose an explicit `renv` workflow.
- Set and document random seeds for stochastic analysis steps when reproducible results are expected.
- Fail early with informative messages when required files, samples, metadata columns, identities, or expected output counts are missing.
- Preserve provenance for important results by recording input paths, parameters, sample/task mappings, and output paths in logs, manifests, or summaries where practical.

## Runtime estimates

- Whenever recommending, creating, editing, reviewing, or asking the user to run R code, an R script, or a Slurm job, always provide an estimated wall-clock runtime.
- Give a realistic range and state the basis for the estimate, such as prior `sacct` history, sample or cell count, array size and concurrency, or static code inspection.
- For Slurm arrays, report both the expected runtime per task and the approximate total wall-clock runtime at the configured concurrency. Exclude queue wait unless it is explicitly estimated.
- Distinguish the expected runtime from the requested Slurm time limit.
- If there is insufficient evidence, label the estimate as low-confidence and state what measurement would improve it.

## Validation requirements

- After editing R code, parse every changed R script before handoff. A syntax parse is expected to take under 10 seconds per script and does not establish scientific correctness.
- Run the most relevant supported `--dry-run` after syntax checks. State whether it only inspects paths or loads large RDS objects, and provide an estimated runtime before recommending or running it.
- After changing configuration or manifests, run `Rscript validate_config.R`; use `--check-files` on the HPC when referenced project inputs are available. The structural check is expected to take under 1 minute; full file checking is filesystem-dependent and must be estimated separately.
- Validate changed Slurm launchers with `bash -n`, expected to take under 5 seconds per file, and inspect every `#SBATCH` directive and invoked R command.
- Use stage-specific validation scripts or output manifests after completed jobs. Do not run a full expensive analysis merely as a test without explicit authorization.
- Treat local checks as lightweight validation. Run memory-intensive Seurat/Xenium stages on the HPC unless the user explicitly requests and authorizes a local run.

## Slurm conventions

- Verify array bounds, task-to-sample mapping, and concurrency limits against the active manifest or the R driver's `--list` output before submission.
- Review requested memory and time against `slurm_resource_history.psv` and any relevant `sacct` evidence. Exclude Giotto and SpaTrack jobs from this review while they remain out of scope.
- Keep the requested time limit safely above the expected runtime while avoiding unsupported guesses. Explain the evidence and confidence behind resource recommendations.
- Preserve the established project working directory, module setup, `renv` activation strategy, and job/array identifiers in log filenames unless there is a documented reason to change them.
- Never add package installation or `renv::restore()` to the body of a production Slurm job.
- After important jobs complete, suggest collecting `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`, `AllocCPUS`, and `ExitCode` with `sacct` so future requests can be refined.

## Code Review Rules

- Flag any change that could silently alter sample-to-array-task mapping, cell filtering or removal, cell/sample identifiers, metadata joins, label-transfer thresholds, consensus-label calculation, VZ/RL identity selection, output filenames, or overwrite behavior.
- For a scientifically meaningful behavior change, identify the affected stage, the downstream outputs that become stale, and a safe validation or rollback path.
- Flag hard-coded sample lists, dated input filenames, absolute paths, and array bounds when a manifest or configuration value should be used instead.
- Treat unexpected missing samples, duplicate cell IDs, unmatched metadata, empty identities, non-finite values, and partial output sets as errors unless the workflow explicitly documents them as acceptable.

## Git and renv follow-up

- After major code changes, always suggest an appropriate Git checkpoint: review `git status` and `git diff`, run relevant validation, and then commit the completed logical unit if the user is ready.
- After major changes involving R packages, package versions, or the R environment, always suggest checking `renv::status()` and updating `renv.lock` with `renv::snapshot()` when appropriate.
- Do not recommend an `renv::snapshot()` solely because analysis code changed when the package environment did not change.
- Clearly distinguish suggested follow-up commands from actions already performed. Never create a Git commit, push changes, or update `renv.lock` without the user's authorization.
