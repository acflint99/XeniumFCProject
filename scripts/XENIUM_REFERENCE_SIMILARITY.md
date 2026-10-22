# Comparing corrected Xenium populations with original reference clusters

`xenium_compare_reference_similarity.R` compares both `corrected_broad_label`
and `comb_subcluster` with every retained original reference identity. It writes
correlation heatmaps, all-pairs scores, best/runner-up matches, and conditional
sample-bootstrap stability. It never changes Xenium annotations or source files.

## Inputs and the approved Science cohort

The default query is the corrected stage-03 nonspatial merge:
`outputs/xenium/vz_rl/nonspatial/03_merged/rds/Xenium_VZRL_merged.rds`.
The complete sample set is validated against the active `config/samples.csv`
(currently 41 samples), and corrected broad labels against the existing approved
fine-to-broad map. Raw counts come from the Xenium assay, grouped by orig.ident.

Original reference inputs are resolved under `config$inputs$published_data_root`:

| Reference | Original file | Counts | Original identity |
|---|---|---|---|
| Aldinger | `AldingerFC/Aldinger_seurat.rds` | RNA counts | figure_clusters |
| Sepp | `SeppFC/Sepp_hum_sce_final.rds` | umi, read directly from SCE | cell_type |
| Science | `ScienceBraunFC/Science_HumanFetalBrainPool.h5` | /shoji/Expression | /shoji/CellClass (default) |

The user approved broader fetal coverage for Aldinger and Sepp on 2026-10-21:
Aldinger age 9–21 PCW inclusive (69,174 cells in the inspected source), and Sepp
original fetal Stage values 7,8,9,11,17,20 wpc (100,375 cells before annotation
checks). Sepp postnatal/adult cells are excluded; its five reviewed Tissue
values are all cerebellar dissections. Unknown age/tissue values cause errors.
Published annotation values are retained, and source-to-selected cell audits are
written. Science uses the cohort explicitly
approved by the user on 2026-10-21: `Region == "Cerebellum"` and **9 ≤ Age ≤ 15**.
As approved 2026-10-22, original Class/CellClass names define the Science
comparison groups. Original numeric cluster IDs remain in the cell and lookup audits. No reference is reclustered, reintegrated, converted to
Seurat, upgraded, or relabeled. No gene-panel-subset reference object is saved.

The completed metadata audit predicts the following Science cohort:

| Age | Cells | Donor | Libraries |
|---|---:|---|---:|
| 11.5 | 12,507 | XDD:358 | 3 |
| 12 | 10,548 | XDD:351 | 3 |
| 14 | 11,021 | XDD:385 | 3 |
| Total | 34,076 | 3 donors | 9 libraries |

It contains 128 original numeric clusters. Libraries are not independent donors.
Source Region/Subregion/Subdivision/Tissue all equal Cerebellum for this cohort;
mixed Pons/Cereb dissections labeled Region=Pons are excluded by the exact region
criterion. No additional QC or minimum-cell threshold is applied. Clusters with
few selected cells should be interpreted cautiously; the full selection table
also records original clusters absent after selection.

This replaces the older project's exact stage set {9,12,13,15}, which would retain
only 10,548 cerebellar cells from one donor at source Age=12. All explicitly
labeled cerebellum would contain 142,362 cells from 12 donors, ages 6.6–14; that
broader cohort was not selected. No gestational-age offset or age rounding is used.
The [authors' repository](https://github.com/linnarsson-lab/developing-human-brain/blob/main/README.md)
identifies the HDF5 as their published processed reference, describes Expression
as spliced plus unspliced UMIs, Gene as gene names, and Accession as Ensembl IDs.
The [paper](https://pubmed.ncbi.nlm.nih.gov/37824650/) describes ages 5–14 PCW.
Using the authors' reference preserves their processing/annotations; it avoids
the project's subsequent filtering, Harmony, reclustering, and broad relabeling.

## How Science counts are read

`R/science_original_reference.R` uses existing system `h5dump` and `jsonlite`.
It needs neither hdf5r nor rhdf5 and never installs packages. Metadata checks
validate unique CellIDs, vector alignment, cluster lookup, CellClass/Class
agreement, NCells counts, and original gene identifiers. Duplicate original
symbols are reported with their accessions and source rows during planning;
they do not invalidate the cell-selection plan. If available, the configured
small Xenium panel gene-name RDS is used to report overlapping duplicates,
without loading the query object. Actual query features are checked before any
reference expression pass. Without panel-ID or explicit-map alignment, shared ambiguous symbols cause an informative stop;
duplicates outside the comparison gene set do not affect one-to-one alignment.
No source rows are silently selected, renamed, or combined.

The HPC plan confirms 116 repeated Science symbols across 1,175 source rows;
six overlap the configured panel: RGS5, SHOX, CSF2RA, IL3RA, CD99, IL9R.
Five have a corresponding _PAR_Y accession; RGS5 has two different stable
Ensembl records. Supply `--panel-json` with an original Xenium panel JSON.
The script validates its unique Ensembl IDs and symbols, checks its symbol set
against the configured small panel vector when available, and uses those IDs
for one-to-one Sepp and Science alignment. Only terminal numeric accession
versions are removed; `_PAR_Y` records remain distinct. Thus all six genes can
be retained using the specific canonical IDs reported by the panel, without
combining separate reference rows. Their scores represent those matched records;
separate Y records are not included in these gene totals. Aldinger retains its
original symbol alignment. Explicit per-reference feature maps/columns take
precedence over panel alignment.

`science_accession_alignment_audit.csv` records all original accessions, symbols,
and inclusion status; `science_panel_id_coverage.csv` records panel coverage.
Sepp writes `sepp_ensembl_alignment_audit.csv`. Original panel IDs, path, MD5,
and JSON are saved with the result provenance. Source references are unchanged.
The optional duplicate-symbol omission mode remains available only without
panel-ID alignment and is not used in the commands below; no omission was approved.

If required, `--science-feature-map=CSV` supplies an explicit one-to-one map,
with `reference_gene` containing exact original Accession values and `query_gene`
containing Xenium feature names. It must list all desired comparison genes,
not only ambiguous ones. Unknown or duplicated mapped accessions and repeated
query names are errors. No first-match or count-based selection is inferred.

The reader plans consecutive runs of retained cells, splits them into blocks of
at most 128 cells, and reads /shoji/Expression as little-endian uint16 binary
hyperslabs. No excluded cells are read to bridge gaps. At 59,480 source genes,
each default block contains at most about 14.5 MiB of binary counts (R objects
and summaries add memory). The source matrix is never loaded in full. Block
orientation and count semantics are checked by comparing every extracted cell's
full-gene sum with original TotalUMIs. Each retained cell must be read exactly
once. Source size/time changes during reading are errors.

Counts are summed within original clusters. Only compact expression summaries
on exactly matched query/reference genes are retained for comparison. Whole-
reference MeanExpression is not reused: it can include cells outside the approved
cohort. Genes outside Xenium's measured features cannot contribute to similarity.

## Methods and limits

1. Sum Xenium counts within sample × label. Normalize shared-gene counts to one
   million, log1p, and average profiles equally across samples containing a label.
   Reference counts are pooled within original identity and normalized over the
   same pairwise shared genes. Large reference donors can dominate pooled counts.
2. Report all-pairs Spearman correlations on shared genes varying across labels
   in both datasets. Report a second, specificity correlation after subtracting
   each gene's mean across labels separately within each dataset.
3. Rank by specificity, recording the best and runner-up, score difference,
   ordinary-correlation best label, and individual Xenium sample best matches.
4. With at least three contributing Xenium samples, resample those samples
   100 times (seed 42 by default), conditional on fixed reference profiles,
   fixed gene selection, and fixed background centering. Report 95% percentile
   intervals for the selected top score and its margin over the original
   runner-up, plus the fraction choosing the same best label. Bootstrap ties
   share weight; individual-sample ties are ranked alphabetically.

These are conditional stability summaries, not calibrated probabilities that a
cell type exists. They omit reference donor uncertainty and are not adjusted for
selecting the best match. Multiple libraries from one Xenium donor would also
require donor-aware aggregation before biological interpretation of resampling.
No cell-level p-values, automatic found thresholds, or new labels are produced.
A nearest match always exists, even when the true population is absent. Inspect
competing matches, shared-gene coverage, markers, age mismatch, and tiny reference
clusters. Pairwise scores using different gene sets are not calibrated across
references. These references contributed to annotation, so agreement is not
independent validation. Default Science grouping is the published broad CellClass.
Explicit --science-label=Clusters retains finer numeric IDs, with class names
added alongside IDs in plots; the 51 AnnotationName marker programs are not substituted
as mutually exclusive cell types, and no posterior threshold is introduced.

## Running on the HPC

Copy the driver, `R/science_original_reference.R`, and
`R/reference_similarity_plots.R` into their corresponding locations under the
HPC project's scripts directory. The optional CSV-only redraw command also
requires `xenium_redraw_reference_similarity.R`. Use the existing R 4.4.1
module and project renv. Commands below assume the project root.

```bash
Rscript --max-ppsize=500000 scripts/xenium_compare_reference_similarity.R \
  --dry-run --panel-json=data/FCXeniumProject/FB124_X_G/gene_panel.json
Rscript --max-ppsize=500000 scripts/xenium_compare_reference_similarity.R \
  --science-plan --panel-json=data/FCXeniumProject/FB124_X_G/gene_panel.json
```

The path/dependency dry run takes **under one minute** and reads the small panel JSON and configured panel-name vector, but no query/reference RDS or HDF5 data.
`--science-plan` is a separate **metadata-only** preflight, estimated **1–3 minutes**
based on the observed 67-second HPC metadata audit; gene-vector checks add work.
It does not load the Xenium/RDS references, read Expression, or write outputs.
Review its cohort counts (expected 34,076 cells, 3 donors, 9 libraries, 128
clusters) and block count before the production expression pass. It prints the
uncompressed selected data volume and maximum block size. A highly fragmented
plan exceeding 10,000 blocks stops for I/O review.

After both preflights pass, the full comparison command in an HPC allocation is:

```bash
Rscript --max-ppsize=500000 scripts/xenium_compare_reference_similarity.R \
  --references=aldinger,sepp,science \
  --panel-json=data/FCXeniumProject/FB124_X_G/gene_panel.json
```

Expected runtime: **1–6 hours or longer**, low confidence from source inspection,
original RDS loads, about 4 GiB of selected uncompressed Science count values,
HDF5 chunk decompression/open overhead, one query count aggregation with exact fine-to-broad pooling, and 100 bootstrap
iterations per query population. Measure the first expression chunks and collect
Slurm JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode to refine this estimate.
This is one process, not an array; queue wait is excluded. The estimate is not a
requested Slurm limit or a verified memory request. No full real-data run has
been performed locally. Do not launch the full analysis as a local validation test.

For a batch run, copy `run_xenium_compare_reference_similarity.slurm` to the HPC
scripts directory and submit from the project root:

```bash
sbatch scripts/run_xenium_compare_reference_similarity.slurm
```

The launcher requests one CPU, 256 GiB and a 12-hour limit on `short` (current
12-hour partition limit verified against UAB documentation on 2026-10-21).
Expected runtime remains **1–6 hours or longer**, excluding queue wait, with low
confidence. The memory request is a conservative first-run envelope: job
41005351 on the merged expression workflow used 22.86 GiB, while related VZ
job 40966278 used 152.09 GiB. Neither measured this complete reference analysis.
Original reference loading and validation can add allocations. The launcher
runs the lightweight preflight, preserves existing renv/module setup and output
protection, and logs peak resident memory with `/usr/bin/time -v`.

A pasted HPC attempt passed the 5,001-gene panel preflight and Science cohort
planning, then printed `Killed` after announcing query loading. This suggests an
external resource termination; it does not establish OOM without allocation
accounting. No reference expression pass or scores were reported. The output
directory is created later in the driver. Keep any existing results protected;
if needed select a new directory via `REFERENCE_SIMILARITY_OUT` at submission.
After completion collect `sacct -j JOBID --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode`
to refine requests. The unrelated `Unknown guide: edge_colourbar` traceback in
the same transcript belongs to the VZ resolution-stability job, not this script.

`--references=science` runs only the approved Science comparison; full default
references are aldinger,sepp,science. `--science-block-rows=N` (1–1024) controls
I/O memory and call count, not the selected cells or scientific parameters.
`--bootstrap=0` disables resampling (rough full-run estimate 30 minutes–4 hours,
low confidence; source reading may dominate). Existing output directories are
always protected, including partial runs. Choose a new `--out=PATH` for each run.

Original RDS metadata inventories use `--inspect --references=aldinger,sepp`
(5–60 minutes, low confidence; loads each full reference sequentially). Sepp
Ensembl IDs are aligned through the original panel JSON in the commands above.
Alternative explicit alignment can use an original feature field
(`--sepp-feature-column=NAME`) or CSV (`--sepp-feature-map=PATH`) with columns
reference_gene,query_gene. No mapping database is guessed or installed.
Equivalent options apply to Aldinger. Science also supports a complete exact
accession-based `--science-feature-map` instead of panel alignment.

The metadata inspector reads saved Seurat metadata/assay slots without invoking
Seurat S4 inheritance dispatch or its analysis dependencies. This supports legacy
objects whose serialized class points to package Seurat. The --inspect branch
does not require SeuratObject/Matrix to be loaded for Seurat inspection. SCE
inspection still uses SummarizedExperiment. Age/stage/tissue/region/donor/sample
value counts are printed (first 100 values per field) to identify suitable
original cohorts. The inspector now also prints original orig_cluster,
orig_sub_cluster, precisest_label, dev_state, and subtype values where present,
annotation completeness by original stage, and selected-cohort missing-label
counts. The supplied HPC inspection confirms 1,494 fetal cells are missing both
cell_type and precisest_label; orig_cluster and orig_sub_cluster are complete
for all 100,375 selected fetal cells. Across the full source, orig_cluster has
68 values and orig_sub_cluster has 662. The user selected original biological
names rather than numeric IDs. Sepp therefore keeps original cell_type as the
main comparison identity (for example VZ_neuroblast and Purkinje). Original
precisest_label, original IDs and available cell-type descriptions are recorded
in sepp_original_cell_annotations.csv. For a finer named comparison, explicitly
select --sepp-label=precisest_label with a fresh output directory; that uses
names such as VZ_neuroblast_1 and Purkinje_defined_EB.

As approved 2026-10-22, cells lacking the selected original published name are
excluded from matching, rather than placed in Unannotated comparison groups.
The earlier audit found 1,494 such fetal cells; that would leave 98,881 named
fetal cells if the inputs are unchanged. Actual counts are recorded in the new
run. `sepp_excluded_unannotated_cells.csv` preserves excluded CellIDs, original
labels/IDs, Stage and Tissue when present. The age-cohort audit still records
age eligibility before annotation filtering. The original source is unchanged.
Scores, backgrounds, competing matches and bootstrap summaries are recalculated
using the retained named groups; old Sepp CSVs are stale for this choice.
Named finer matching still requires --sepp-label=precisest_label.

If the configured original Xenium panel JSON is available, inspection prints
the original type/data annotation fields for the six ambiguous symbols and
examples of other gene targets. This can establish a direct Ensembl-to-query
mapping for Sepp and inform Science ambiguity resolution without downloading a
new annotation database. Field availability and ID cardinality must be verified.
The configured example JSON was unavailable in the HPC inspection. The new
`--panel-inspect` mode checks that path plus gene_panel.json at manifest-defined
input directories and one child-directory level, without recursive raw-data
searches or loading reference RDS files. It prints available paths and inspects
up to three JSON files, comparing their gene-name sets to the configured small
panel RDS and showing original identifier fields for the six ambiguous genes.
An explicit `--panel-json=PATH` is also supported. This read-only mode writes no
outputs and infers no mapping. Estimated runtime is under one minute, low
confidence from bounded file checks and small JSON reads; filesystem latency may
increase it. Syntax, path-only preflight and a synthetic explicit-panel test passed.
Inspection still loads the full saved object sequentially;
its 5–60-minute estimate remains low confidence. No source object is changed.

An HPC inspection attempt failed while autoloading Seurat because igraph.so
could not find libglpk.so.40. The revised metadata path avoids that dependency
for the inspection only; it does not repair the analysis environment. Before a
full run, inspect available GLPK modules and active project libraries. Use the
established project-root working directory/module/renv strategy, and do not
reinstall or upgrade packages merely on the basis of the R patch-version warning.

## Results and provenance

Default directory: `outputs/xenium/reference_similarity/`.
Each reference × query annotation level gets all_correlations and best_matches
CSV tables, a two-page heatmap PDF, a shared-gene audit, and individual sample
best-match tables. Per-sample files queryN correspond to row N of best_matches
(1-based, query labels sorted alphabetically).

Science additionally records selected CellIDs with original one-based source
rows, ages, donors, libraries, and clusters; the full original cluster-selection
lookup; count-block plan; cohort summary; shared-gene source indices/accessions;
and selected cluster cell counts. The full original duplicate-symbol audit is
saved as science_original_duplicate_symbols.csv. Original IDs/labels are unchanged. All other
original clusters are retained as zero-selected rows in the lookup.
For Aldinger/Sepp, reference_cohort_cells.csv and reference_cohort_counts.csv
record original cell IDs, stage values and selection without altering the source.

Input paths/sizes/times, reference settings, script/helper/config/manifest
snapshots, and methods/session information are saved. The output's source_genes
field records the full Science feature count; reference_genes in score summaries
counts the features retained in that reference's aggregate. Only COMPLETED.txt
marks a successful output set; a directory without it is incomplete.

## Validation and checkpoint

Both changed R files parse. Synthetic tests passed for approved region/age
selection, original/nonsequential cluster IDs, binary matrix orientation,
shared-gene order, chunked sums and row-total checks, changed-source rejection,
all six query/reference comparisons, known label recovery, split Seurat v5 query
count layers, direct original SCE counts, sample bootstrap output, and overwrite
protection. Local tests use R 4.5.1; production remains R 4.4.1/renv. The real-data
local dry run reports missing HPC-only input files; no large source was loaded.
The actual Science --science-plan and expression pass still need HPC validation.

The older metadata inspectors remain read-only diagnostic tools:
inspect_science_original_h5.R, inspect_science_shoji.sh, and
audit_science_shoji_metadata.R. The completed audit supplied by the user is
unchanged. The comparison cohort change affects comparison results only; original
reference files and Xenium annotations do not need regenerating. For rollback,
retain prior result directories and choose a new output path for each cohort.
Review changes and validation before a Git checkpoint; this local folder is not
currently recognized as a Git repository and no commit was created. No R package
or renv.lock changes were made, so no renv snapshot is needed for these code edits.

Panel-ID integration validation: both changed R scripts parsed (<10 seconds each),
a small-data end-to-end test passed all six reference/query-level comparisons,
including canonical/alternate RGS5 and separate PAR_Y rows, expected matches,
bootstrap behavior and output protection. This took under one minute locally;
it does not establish biological correctness or real-data HPC memory/runtime.

## Original Aldinger compatibility fix (job 41010260)

This job completed query aggregation but failed at Aldinger metadata access:
`no slot of name "images" for this object of class "Seurat"`. This is an
original-object/current-method compatibility error, not an OOM in that job.
Elapsed time was 5:16; sacct reported 33695788K MaxRSS (32.14 GiB), while
GNU time reported 37045212 KiB (35.33 GiB). Neither is a complete-run peak.
The 256-GiB request is unchanged pending a completed run across all references.

The driver now reads original Seurat metadata, assay inventories and traditional
Assay count slots directly, avoiding modern container methods. Approved cohort
selection is applied to count columns during aggregation, without subsetting,
updating or resaving the source Seurat container. It validates selected cell IDs
and count-layer coverage. Original SCE access and v5 split-layer support remain.
No labels, age windows, normalization or scoring settings were changed.

Validation: syntax parse (<10 seconds), supported lightweight dry run (<1 minute),
and a small all-reference integration test passed. The test includes an original
Seurat container missing `images`, a traditional count Assay, split v5 query
layers, Aldinger age exclusions and Sepp postnatal exclusions, plus expected
matches, panel-ID alignment, bootstrap and overwrite protection. This does not
validate the complete HPC inputs. Full expected runtime remains 1–6 hours or
longer with low confidence; exclude queue wait.

Copy the updated driver and rerun into a new output directory, preserving the
partial result directory from job 41010260:

```bash
REFERENCE_SIMILARITY_OUT=outputs/xenium/reference_similarity_legacy_fix \
  sbatch scripts/run_xenium_compare_reference_similarity.slurm
```

## 2026-10-22 revision: named Science groups, readable scales, faster scoring

The user approved published broad Science class names. Original numeric-cluster
count aggregates are summed by the original Class/CellClass mapping, before
normalization. Count pooling (not averaging cluster correlations) preserves
within-class cell weighting. Fine source IDs remain in audit tables. Class
pooling changes Science profiles, backgrounds, matches and uncertainty; all old
Science score tables and PDFs must be regenerated. Select --science-label=Clusters
only to explicitly recover the original fine-cluster comparison.

Every heatmap page now has a numeric Spearman-rho color legend. Raw-expression
plots use the full observed score range for contrast: blue means lower, and can
still be positive. White marks the scale midpoint, not necessarily zero. Broadly
shared expression can give positive raw correlations across most cell types;
a red page alone is not evidence that all populations match. Specificity plots
remain zero-centered and use symmetric observed limits. Per-plot scale limits
are explicit in legends, so colors must not be compared across PDFs as absolute
scores. Correlation values are not stretched or rescaled in the CSVs.

Performance changes retain the seed and use 100 bootstrap iterations (approved reduction from 200): count
layers are aggregated once by fine population/sample, and broad counts are
exact sums of those groups; unchanged reference ranks are cached, and bootstrap
profiles are scored together. Tests verified exact broad-count equality and
agreement with the original Spearman/bootstrap calculations within 1e-12.
A local tied-rank benchmark (5,000 genes, 12 query profiles, 60 reference groups,
five scoring passes) took 1.781 seconds with repeated cor calls and 0.046 seconds
with cached ranks, about 39x faster in that component. This is not a whole-job
speedup; input loading and HDF5 reads remain. Keep the low-confidence full-run
estimate of 1–6 hours or longer until measured; queue wait is excluded and the
launcher still requests 12 hours/256 GiB.

Run from the HPC scripts directory after copying the updated files:

```bash
REFERENCE_SIMILARITY_OUT=outputs/xenium/reference_similarity_20261022 \
  sbatch run_xenium_compare_reference_similarity.slurm
```

For later styling changes, redraw from a completed revised result directory
without loading any large expression objects (expected <1–2 minutes, based on
small CSV inputs; syntax check <10 seconds/script and dry run <1 minute):

```bash
Rscript --max-ppsize=500000 xenium_redraw_reference_similarity.R \
  --input=outputs/xenium/reference_similarity_20261022 \
  --out=outputs/xenium/reference_similarity_20261022_redraw --dry-run
Rscript --max-ppsize=500000 xenium_redraw_reference_similarity.R \
  --input=outputs/xenium/reference_similarity_20261022 \
  --out=outputs/xenium/reference_similarity_20261022_redraw
```

The redraw refuses old Science numeric-cluster scores or Sepp Unannotated groups:
those scientific changes require recomputing counts/profiles and correlations.
No existing analysis or plot directory is overwritten. The supplied local folder
contains only PDFs, so the actual revised analysis cannot be regenerated locally
from those files. Original PDFs were preserved. Small-data integration covered
legacy Aldinger slots, v5 split query layers, annotation/age exclusions, pooling
multiple Science clusters per published class, known matches, bootstrap stability
and output protection. Revised synthetic pages were rendered and visually checked.
After the logical update is reviewed and validated, consider a Git checkpoint
on the HPC checkout; no commit or package/environment change was made here.

On 2026-10-22 the user approved reducing bootstraps from 200 to 100. The driver
default and batch launcher now use 100. Point correlations and best matches
are unaffected; bootstrap intervals and stability estimates have more Monte
Carlo variability. Bootstrap work is approximately halved; expression loading
and Science count reading are unchanged. Full-run runtime remains a low-confidence
1–6 hours or longer until measured, excluding queue wait; the requested time
limit is unchanged. Use a new output directory to preserve prior summaries.

## Shared heatmap order (2026-10-22)

Both the full comparison and CSV-only redraw use the current color_palette.R
master_subcluster_order for comb_subcluster rows and celltype_order for broad
rows. The first configured label appears at the top. Only observed labels are
shown. Numbered Unknown groups follow the configured named groups; unexpected
unlisted labels fail explicitly rather than disappearing. Orders are read from
the literal vectors in the palette file without sourcing unrelated plot helpers.

This change affects display only: scoring order, reference columns, score values,
bootstrap random draws, sample-match filenames and the 100-bootstrap default
are unchanged. Each PDF has an accompanying *_heatmap_row_order.csv giving its
top-to-bottom order; the palette source is snapshotted alongside the plots.
The full driver validates query order coverage before count aggregation.

After copying the update to HPC scripts/, redraw completed revised results:

```bash
Rscript --max-ppsize=500000 xenium_redraw_reference_similarity.R \
  --input=outputs/xenium/reference_similarity_100boot \
  --out=outputs/xenium/reference_similarity_100boot_ordered --dry-run
Rscript --max-ppsize=500000 xenium_redraw_reference_similarity.R \
  --input=outputs/xenium/reference_similarity_100boot \
  --out=outputs/xenium/reference_similarity_100boot_ordered
```

The redraw should take under 1–2 minutes (small CSV inputs, no expression-object
loading); its path-only dry run should take under one minute. Use the actual
completed result directory if it differs. Existing outputs are protected.
Syntax checks (<10 seconds/script), supported full-driver and redraw preflights,
and a six-PDF synthetic redraw passed. Shared ordering, numbered Unknown handling
and top-to-bottom layout were checked, including rendered PDF inspection. No
full expensive analysis was used for validation and no original data changed.
