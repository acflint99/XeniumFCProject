# MetaNeighborUS against the original fetal cerebellum references

`xenium_compare_metaneighbor.R` is a separate validation stage. It compares the
existing Xenium `corrected_broad_label` and `comb_subcluster` populations against
each original reference separately. It does not replace the Spearman comparison,
change annotations, recluster, integrate datasets or write source RDS/HDF5 files.
Coarse published reference labels are allowed for both Xenium resolutions.
Several fine Xenium populations can legitimately match the same reference class.

## Original inputs and labels

The merged query defaults to
`outputs/xenium/vz_rl/nonspatial/03_merged/rds/Xenium_VZRL_merged.rds` (Xenium raw
counts, `orig.ident`). Its complete sample set must match `config/samples.csv`
(currently 41 entries), and the existing fine-to-broad map must match exactly.
Original reference roots come from config$inputs$published_data_root; processed
reference RDS paths in config$inputs$references are deliberately not used.

| Array task | Original source | Original labels | Approved cohort |
|---|---|---|---|
| 1 | AldingerFC/Aldinger_seurat.rds | figure_clusters | age 9–21 PCW inclusive |
| 2 | SeppFC/Sepp_hum_sce_final.rds | cell_type | fetal 7/8/9/11/17/20 wpc, reviewed cerebellar tissues |
| 3 | ScienceBraunFC/Science_HumanFetalBrainPool.h5 | CellClass | exact Cerebellum, numeric Age 9–15 inclusive |

The Sepp unannotated-cell exclusion and Science published broad Class/CellClass
choice are carried over from the user's approval. Source numeric cluster IDs
and all scalar original metadata remain available in audits. Non-scalar metadata
columns are not flattened into CSV; any such fields are listed explicitly and
remain unchanged in the source. The age-cohort audit distinguishes age eligibility
from original annotation availability. Missing Aldinger published labels are errors.

By default, the driver verifies that saved Xenium fine labels reproduce the
current approved broad-label mapping. `--legacy-query-labels` is an explicit
exception for the older saved merged object containing standalone `Endothelial`
and `Pericytes`. It preserves both saved query columns as-is, places those
legacy labels immediately before `Vascular` in heatmaps, and writes
`comparison_mode.csv` plus a LEGACY title/provenance note. It does not edit the
object, harmonize either label to Vascular, or make the results equivalent to
the current convention.

`--panel-json` defaults to the discovered original
`data/FCXeniumProject/FB124_X_G/gene_panel.json`. Panel Ensembl IDs align Sepp rows
and Science accessions one to one, removing terminal numeric accession versions
only. Separate `_PAR_Y` records and the alternate RGS5 accession are not combined.
Aldinger uses exact shared symbols. All original reference cells satisfying the
approved cohort and named-label requirement remain in the scoring population.
No full reference is saved as a panel-subset or reprocessed Seurat object.

## Method and interpretation

The implementation calls the installed package's public functions directly:

- `MetaNeighbor::variableGenes`: recurrent HVGs in both datasets, selected from
  shared measured panel genes without using cell-type labels. A seeded random
  sample of up to 10,000 cells per dataset is used only for feature selection;
  this is not cell downsampling for AUROC scoring. Selected source CellIDs and
  every shared gene's selection status are written to audit tables.
- `MetaNeighbor::MetaNeighborUS`: official published low-memory approximation,
  `fast_version=TRUE`, `node_degree_normalization=TRUE`, `one_vs_best=FALSE`,
  `symmetric_output=FALSE`. Raw-count within-cell ranks avoid needing matched
  expression units. The same selected HVGs are used for both query label levels.
- Directed cross-dataset AUROCs are retained separately. The displayed
  bidirectional summary is their arithmetic mean, matching the package's usual
  symmetrization. It is not a new classification probability.
- Reciprocal directional best matches are recorded without forcing a one-to-one
  taxonomy. Tied best matches are counted and retained in the all-pairs table;
  the summary's representative label uses lexical tie-breaking. A missing
  reciprocal match is not, by itself, an invalid annotation.

The package uses rank-normalized cell relationships, rather than correlations
between count pseudobulks. No enormous cell-by-cell similarity matrix is built.
A guard estimates the size of one dense selected-gene-by-cell matrix (default
32 GiB maximum); several such allocations may coexist, so this is not a job
memory guarantee. No cells or genes are silently discarded to satisfy the guard.
There must be at least five recurrent matched variable genes. Cells whose
selected-gene expression is constant cannot be ranked. They are excluded only
from MetaNeighbor scoring, never from a source object or annotation. Every
eligible cell remains in source-cell audits with an `included_in_scoring` flag,
and `rank_constant_cell_audit.csv` records its ID, selected-HVG count and reason.
The driver stops if the exclusion would remove every cell from any Xenium or
reference label.

References:
[Crow et al., Nature Communications 2018](https://www.nature.com/articles/s41467-018-03282-0),
[Fischer et al., Nature Protocols 2021](https://www.nature.com/articles/s41596-021-00575-5),
and the authors' [MetaNeighborUS source](https://github.com/gillislab/MetaNeighbor/blob/master/R/MetaNeighborUS.R).

This stage reports one-versus-rest recovery, not a within-lineage subtype test
or a one-versus-best analysis. Broad reference names support broad correspondence;
they do not establish a fine subtype name. Scores near 0.5 indicate chance-level
discrimination, but no universal 'found' threshold or p-value is applied. There
are no bootstrap intervals in this stage; the existing 100-bootstrap Spearman
stage remains separate. Donors/samples are not held out here, and the same
references contributed to annotation, so these results support transcriptional
consistency rather than independent annotation accuracy or full pipeline validity.

## Dependencies and preflight

Use the established R/4.4.1-foss-2023b module and project renv. No installation,
restore, R upgrade, package upgrade or lockfile update occurs in the analysis or
launcher. Required packages are here, yaml, Matrix, SeuratObject,
SummarizedExperiment, jsonlite and MetaNeighbor. Science also uses system h5dump.
Legacy Aldinger objects are read through original slots without requiring modern
Seurat container slots such as images. The new code checks the installed
MetaNeighbor API for the low-memory/directional options.

Copy the update archive into the HPC project's scripts directory. From scripts/:

```bash
tar -xzf xenium_metaneighbor_setup_20261022.tar.gz
Rscript --max-ppsize=500000 xenium_compare_metaneighbor.R --list
Rscript --max-ppsize=500000 xenium_compare_metaneighbor.R --reference=aldinger --dry-run
Rscript --max-ppsize=500000 xenium_compare_metaneighbor.R --reference=sepp --dry-run
Rscript --max-ppsize=500000 xenium_compare_metaneighbor.R --reference=science --dry-run
```

Each path/package preflight is expected to take under one minute; all three
should take under three minutes. They load no query/reference RDS, HDF5 contents
or expression matrix. The tiny engine smoke test uses synthetic data only:

```bash
Rscript --max-ppsize=500000 test_xenium_metaneighbor_engine.R
```

Estimate: 10–60 seconds, low confidence from 240 genes and 72 synthetic cells plus
package startup. It tests both query-label resolutions, directional naming and
known reciprocal matches using the actual installed MetaNeighbor engine.

If preflight reports MetaNeighbor missing, first review the project environment.
An explicit package setup outside analysis jobs can use the existing project
Bioconductor version; do not upgrade R/Bioconductor to the newest release:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
Rscript --max-ppsize=500000 -e 'renv::status()'
Rscript --max-ppsize=500000 -e 'renv::install("bioc::MetaNeighbor")'
```

These are suggested environment actions, not actions already performed. Status
may take 1–5 minutes on a large project; installation may take 5–30 minutes or
longer depending on missing dependencies and compilation (low confidence).
Review the package plan and rerun the preflights/smoke test afterward. If the
package environment changes, review renv::status() and snapshot renv.lock when
ready, with explicit authorization. Do not run renv::restore() in production jobs.
[renv Bioconductor guidance](https://pkgs.rstudio.com/renv/articles/bioconductor.html).

## Batch run and runtime

After preflights and the engine smoke test pass, submit from scripts/:

```bash
sbatch run_xenium_compare_metaneighbor.slurm
```

This array has three tasks, one per reference, at concurrency two. Each task
loads the query once and runs both query label levels. The launcher prints the
reference/task map and runs a supported preflight before expression loading.
It also runs the tiny synthetic engine test before loading the query.
It preserves project working directory, module setup, renv startup and job/task
log identifiers. Optional METANEIGHBOR_OUT selects a fresh output root; existing
per-reference directories (including partial runs) are protected.
For the explicit legacy query comparison, set
`METANEIGHBOR_LEGACY_QUERY_LABELS=1` at submission; the launcher passes the
corresponding driver flag to every array task.

Science reads the same selected cells and genes in bounded HDF5 blocks of 512
cells by default (`METANEIGHBOR_SCIENCE_BLOCK_ROWS`, permitted range 1–1024).
This is an I/O-only setting: it reduces `h5dump` process launches compared with
the former 128-cell blocks and does not change counts, selected cells, HVGs or
AUROCs. The larger maximum raw block is about 58 MiB before R-object overhead.

Expected runtime is **4–12 hours per task**, and **8–24 hours total** at
concurrency two, excluding queue wait. This is a low-confidence static estimate
from ~2.36 million query cells, original reference loads, feature selection,
Science bounded HDF5 reads and two low-memory scoring passes per reference.
The requested limit is **12 hours per task**, with **one CPU and 256 GiB** on
short. The [current Cheaha table](https://docs.rc.uab.edu/cheaha/hardware/) confirms
short's 12-hour limit and 3,072-GiB per-person memory QoS; two 256-GiB tasks fit
within that published ceiling, although fair share and node availability govern
when they start. Resources are provisional: related VZ history reached
152.09 GiB, and the earlier partial reference-similarity run reached ~35.3 GiB
before any complete MetaNeighbor computation. No completed MetaNeighbor sacct
history exists. Rank matrices add allocations; measure the first complete run.

## Outputs, validation and reproducibility

Default output root: outputs/xenium/metaneighbor_original. Each reference has its
own directory and COMPLETED.txt marker. Each of the two query levels produces:

- *_all_AUROCs.csv: both directional AUROCs, their mean, best-hit and reciprocal flags.
- *_best_matches.csv: best/runner-up original names, margins, ties and reciprocity.
- *_AUROC_heatmaps.pdf: three pages (query recovery, reference recovery, mean),
  fixed 0–1 numeric legends, white at 0.5, shared Xenium top-to-bottom order,
  and an asterisk for the highest score in each Xenium row (all tied maxima marked).
- Query and reference cell-count tables.

Audits include source and comparison CellIDs, original metadata, cohort/name
selection, gene alignment, feature-selection cell IDs, rank-constant cell
exclusions, selected HVGs, panel IDs
and MD5, input paths/size/mtime, source/config/manifest snapshots, session/package
versions and an output manifest. A completion marker is written only after the
required output set exists and files are nonempty. Source RDS/HDF5 files and
cell labels are never overwritten. Generated large data/logs should not be added
to Git. After completion collect JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
with sacct and check all three per-reference completion markers/manifests.

Local validation: all new R scripts parsed (<10 seconds each), Slurm bash syntax
checked (<5 seconds), and reference mapping/list mode verified against the active
41-sample manifest. Small tests passed for legacy Seurat/v5/SCE/HDF5 original
readers, approved cohorts, Sepp missing-name exclusion, canonical gene alignment,
unchanged input checksums, invalid-profile guards, directed AUROC formatting,
ties and reciprocal maxima. A synthetic three-page heatmap was rendered and
inspected. The actual MetaNeighbor package is absent locally: supported preflight
correctly stops with that missing dependency, and no official-engine numerical
or full HPC analysis is claimed as tested. The supplied smoke test addresses
engine validation on HPC. No package installation or renv.lock change was made.
Review git status/diff and these checks before creating a logical Git checkpoint
on the HPC checkout; no commit or push has been made here.
