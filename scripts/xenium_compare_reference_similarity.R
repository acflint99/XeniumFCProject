#!/usr/bin/env Rscript
# Read-only comparison of corrected Xenium populations with published labels.
# Expected HPC runtime: 1–6 h, low confidence (original full-reference sizes unknown).
# --dry-run checks paths only (<1 min); --inspect loads one object at a time.
args <- commandArgs(trailingOnly = TRUE)
flags <- c("--dry-run", "--list", "--inspect", "--science-plan", "--panel-inspect", "--help")
keys <- c("query", "out", "query-assay", "sample-column", "bootstrap", "seed",
          "aldinger", "sepp", "science", "aldinger-label", "sepp-label", "science-label",
          "aldinger-assay", "sepp-assay", "science-assay", "references",
          "aldinger-feature-column", "sepp-feature-column", "science-feature-column",
          "aldinger-feature-map", "sepp-feature-map", "science-feature-map", "science-block-rows", "science-duplicate-policy", "panel-json")
opt <- list()
for (a in args) {
  if (a %in% flags) next
  if (!grepl("^--[^=]+=.+$", a)) stop("Use --option=value: ", a)
  k <- sub("=.*$", "", substring(a, 3L))
  if (!k %in% keys || !is.null(opt[[k]])) stop("Unknown or repeated option: ", a)
  opt[[k]] <- sub("^[^=]+=", "", a)
}
if ("--help" %in% args) {
  cat("Compare corrected_broad_label and comb_subcluster against each reference.\n",
      "Modes: --dry-run / --list (paths only), --inspect (metadata inventory), or full run.\n",
      "Options use --key=value. Keys: ", paste(keys, collapse = ", "), "\n",
      "Original sources only: Aldinger_seurat.rds and Sepp_hum_sce_final.rds.\n",
      "Science: original /shoji/Expression, Region=Cerebellum, ages 9–15 inclusive, original CellClass (use --science-label=Clusters for numeric clusters).\n",
      "--science-plan reads Science metadata only (~1–3 min); no query RDS/expression/output read/write.\n",
      "--panel-inspect finds original gene_panel.json files using config/manifest directories; no reference RDS loaded.\n",
      "Use --references=aldinger,sepp to run the two established original sources.\n",
      "Seurat assay or SCE assay: --REF-assay=NAME (Sepp default umi).\n",
      "Gene alignment: --panel-json=PATH supplies original panel IDs for Sepp/Science; alternatively exact IDs, --REF-feature-column=NAME, or --REF-feature-map=CSV.\n",
      "Map CSV columns: reference_gene,query_gene; one-to-one shared-gene alignment only.\n",
      "Science duplicate symbols: default error; explicit --science-duplicate-policy=omit audits and omits shared ambiguous genes.\n",
      "Sepp missing named annotations are excluded and audited; source files remain unchanged.\n",
      "Default bootstrap=100; uncertainty conditions on fixed reference profiles.\n",
      "Use a NEW --out directory for each run; existing directories are protected.\n", sep = "")
  quit(save = "no")
}
for (p in c("here", "yaml")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing dependency in project renv: ", p)
}
here::i_am("scripts/xenium_compare_reference_similarity.R")
source(here::here("scripts", "R", "config.R"))
config <- load_pipeline_config()
samples <- load_sample_manifest(config)$sample_id
valid_text <- function(x) !anyNA(x) && all(nzchar(trimws(as.character(x))))
if (!length(samples) || anyDuplicated(samples) || !valid_text(samples)) stop("Invalid sample manifest.")
value <- function(k, default) if (is.null(opt[[k]])) default else opt[[k]]
path <- function(x) if (grepl("^(/|[A-Za-z]:)", x)) x else here::here(x)
if ("--panel-inspect" %in% args) {
  if (any(c("--dry-run", "--list", "--inspect", "--science-plan") %in% args)) stop("Run --panel-inspect separately from other modes.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Panel inspection needs existing jsonlite; no packages are installed.")
  manifest <- load_sample_manifest(config)
  input_dirs <- unique(vapply(seq_len(nrow(manifest)), function(i) resolve_sample_paths(manifest[i, , drop = FALSE], config)$input_dir, ""))
  if (!is.null(opt[["panel-json"]])) {
    candidates <- path(opt[["panel-json"]])
  } else {
    # Bounded search: configured example, manifest input directories and one
    # child directory level. Never recurse through raw data or the HPC filesystem.
    child_dirs <- unlist(lapply(input_dirs[dir.exists(input_dirs)], list.dirs, recursive = FALSE, full.names = TRUE), use.names = FALSE)
    candidates <- unique(c(resolve_config_path(config$inputs$xenium_panel_json_example, config),
                           file.path(c(input_dirs, child_dirs), "gene_panel.json")))
  }
  existing <- candidates[file.exists(candidates)]
  cat("Original panel discovery: checked", length(candidates), "candidate paths; found", length(existing), "files.\n")
  if (!length(existing)) {
    cat("No original panel JSON found. Manifest input directories (first 15):\n", paste(head(input_dirs, 15L), collapse = "\n"), "\n")
    stop("Supply --panel-json=PATH if the original panel metadata is stored elsewhere. No analysis was run.")
  }
  cat("Available files (first 15):\n", paste(head(existing, 15L), collapse = "\n"), "\n")
  panel_names_path <- resolve_config_path(config$inputs$xenium_panel_genes, config)
  expected <- if (file.exists(panel_names_path)) readRDS(panel_names_path) else NULL
  if (!is.null(expected) && (!is.character(expected) || !is.null(dim(expected)) || !valid_text(expected))) stop("Invalid configured gene-name vector.")
  for (panel_file in head(existing, 3L)) {
    panel <- jsonlite::fromJSON(panel_file, simplifyVector = FALSE)
    targets <- Filter(function(x) identical(x$type$descriptor, "gene"), panel$payload$targets)
    names_found <- vapply(targets, function(x) if (is.null(x$type$data$name)) "" else as.character(x$type$data$name), "")
    cat("\nPANEL:", panel_file, "\nGene targets:", length(targets), "\n")
    if (!is.null(expected)) cat("Gene-name set matches configured panel:", setequal(names_found, expected), "\n")
    cat("Original gene-annotation fields and examples:\n")
    for (x in head(targets, 2L)) utils::str(x$type$data, max.level = 1L)
    for (symbol in c("RGS5", "SHOX", "CSF2RA", "IL3RA", "CD99", "IL9R")) {
      selected <- targets[names_found == symbol]
      cat("\nOriginal panel annotation for", symbol, ":", length(selected), "targets\n")
      for (x in selected) utils::str(x$type$data, max.level = 1L)
    }
  }
  cat("\nInspection only: no mapping inferred, no reference object loaded, no expression read, no outputs written.\n")
  quit(save = "no")
}
query <- path(value("query", file.path(config$project$outputs_dir, "xenium", "vz_rl",
                                      "nonspatial", "03_merged", "rds", "Xenium_VZRL_merged.rds")))
out <- path(value("out", file.path(config$project$outputs_dir, "xenium", "reference_similarity")))
refs <- strsplit(value("references", "aldinger,sepp,science"), ",", fixed = TRUE)[[1L]]
if (!length(refs) || anyDuplicated(refs) || any(!refs %in% c("aldinger", "sepp", "science"))) stop("Invalid --references.")
# Use only paths upstream of project filtering, gene-panel subsetting, and relabeling.
# Science_clusters.rds has unverified provenance; deliberately do not guess it.
original_paths <- c(aldinger = file.path(config$inputs$published_data_root, "AldingerFC", "Aldinger_seurat.rds"),
                    sepp = file.path(config$inputs$published_data_root, "SeppFC", "Sepp_hum_sce_final.rds"),
                    science = file.path(config$inputs$published_data_root, "ScienceBraunFC", "Science_HumanFetalBrainPool.h5"))
paths <- setNames(vapply(refs, function(r) {
  p <- value(r, original_paths[[r]])
  if (nzchar(p)) path(p) else ""
}, ""), refs)
science_h5 <- "science" %in% refs && grepl("\\.(h5|hdf5)$", paths[["science"]], ignore.case = TRUE)
labels <- setNames(vapply(refs, function(r) value(paste0(r, "-label"),
             c(aldinger = "figure_clusters", sepp = "cell_type", science = if (science_h5) "CellClass" else "")[[r]]), ""), refs)
assays <- setNames(vapply(refs, function(r) value(paste0(r, "-assay"),
             c(aldinger = "RNA", sepp = "umi", science = if (science_h5) "Expression" else "originalexp")[[r]]), ""), refs)
feature_columns <- setNames(vapply(refs, function(r) value(paste0(r, "-feature-column"), ""), ""), refs)
feature_maps <- setNames(vapply(refs, function(r) {
  p <- value(paste0(r, "-feature-map"), "")
  if (nzchar(p)) path(p) else ""
}, ""), refs)
if (any(nzchar(feature_columns) & nzchar(feature_maps))) stop("Select a feature column OR feature map per reference.")
if (science_h5 && (!labels[["science"]] %in% c("Clusters", "CellClass") || assays[["science"]] != "Expression" ||
                   nzchar(feature_columns[["science"]]))) {
  stop("Audited Science HDF5 uses original CellClass (or Clusters)/Expression; use --science-feature-map for explicit accession alignment, not --science-feature-column.")
}
science_block_rows <- suppressWarnings(as.numeric(value("science-block-rows", "128")))
science_duplicate_policy <- value("science-duplicate-policy", "error")
if (!science_duplicate_policy %in% c("error", "omit")) stop("Invalid --science-duplicate-policy (error or omit).")
if (science_h5 && nzchar(feature_maps[["science"]]) && science_duplicate_policy != "error") stop("Choose a Science feature map OR duplicate omission, not both.")
science_plan <- NULL
source(here::here("scripts", "R", "science_original_reference.R"))
source(here::here("scripts", "R", "reference_similarity_plots.R"))
if (science_h5) {
  if (!is.finite(science_block_rows) || science_block_rows < 1 || science_block_rows > 1024 || science_block_rows != floor(science_block_rows)) stop("Invalid --science-block-rows (1–1024).")
}
# The original panel supplies the measured target's identity, without a mapping database.
panel_file <- if (!is.null(opt[["panel-json"]])) path(opt[["panel-json"]]) else ""
panel_map <- NULL
if (nzchar(panel_file)) {
  if (!file.exists(panel_file)) stop("Missing original --panel-json: ", panel_file)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Panel alignment requires jsonlite in project renv.")
  doc <- jsonlite::fromJSON(panel_file, simplifyVector = FALSE)
  targets <- Filter(function(x) identical(x$type$descriptor, "gene"), doc$payload$targets)
  field <- function(k) vapply(targets, function(x) {
    a <- x$type$data[[k]]
    if (!is.character(a) || length(a) != 1L || is.na(a)) stop("Invalid panel gene annotation: ", k)
    a
  }, "")
  panel_map <- data.frame(gene = field("name"), ensembl_id = field("id"), stringsAsFactors = FALSE)
  if (!nrow(panel_map) || !valid_text(panel_map$gene) || anyDuplicated(panel_map$gene) ||
      anyDuplicated(panel_map$ensembl_id) || any(!grepl("^ENSG[0-9]{11}$", panel_map$ensembl_id))) stop("Original panel must have unique gene names and unique unversioned human Ensembl IDs.")
  expected_path <- resolve_config_path(config$inputs$xenium_panel_genes, config)
  if (file.exists(expected_path) && !setequal(panel_map$gene, readRDS(expected_path))) stop("Original panel gene names differ from configured panel.")
  cat("Validated original panel ID map:", nrow(panel_map), "genes\n")
}
qassay <- value("query-assay", "Xenium")
sample_col <- value("sample-column", "orig.ident")
B <- suppressWarnings(as.numeric(value("bootstrap", "100")))
seed <- suppressWarnings(as.numeric(value("seed", "42")))
if (length(B) != 1L || !is.finite(B) || B < 0 || B != floor(B) || B > 100000) stop("Invalid bootstrap count.")
if (!is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) stop("Invalid seed.")
map_inputs <- if (any(nzchar(feature_maps))) {
  setNames(feature_maps[nzchar(feature_maps)], paste0(names(feature_maps)[nzchar(feature_maps)], "_feature_map"))
} else character()
inputs <- c(query = query, paths, map_inputs, if (nzchar(panel_file)) c(original_panel_json = panel_file))
print(data.frame(dataset = names(inputs), path = unname(inputs), exists = file.exists(inputs)), row.names = FALSE)
print(data.frame(reference = refs, assay = assays, original_label_column = labels), row.names = FALSE)
cat("Manifest samples:", length(samples), "\nOutput:", out, "\nBootstrap:", B, "seed:", seed, "\n")
if ("--science-plan" %in% args) {
  if (!science_h5) stop("--science-plan requires the original Science /shoji HDF5 source.")
  if (any(c("--dry-run", "--list", "--inspect") %in% args)) stop("Run --science-plan separately from other modes.")
  science_plan <- plan_science_original(paths[["science"]], science_block_rows)
  panel_path <- resolve_config_path(config$inputs$xenium_panel_genes, config)
  if (file.exists(panel_path)) {
    panel <- readRDS(panel_path) # Small configured gene-name vector, not a Seurat object.
    if (!is.character(panel) || !is.null(dim(panel)) || !valid_text(panel)) stop("Configured Xenium panel must be a nonblank character gene-name vector.")
    overlap <- science_plan$gene_duplicates[science_plan$gene_duplicates$Gene %in% panel, , drop = FALSE]
    cat("\nDUPLICATE SYMBOLS OVERLAPPING CONFIGURED XENIUM PANEL (not a query-object coverage check):\n")
    cat("Panel source:", panel_path, "\nAmbiguous panel symbols:", length(unique(overlap$Gene)), "\n")
    print(head(overlap, 100L), row.names = FALSE)
    if (nrow(overlap) > 100L) cat("First 100 overlapping duplicate rows shown.\n")
  } else {
    cat("Configured panel gene-name file is unavailable; overlap must be checked against actual query features before expression reading.\n")
  }
  if (!is.null(panel_map)) {
    a <- resolve_science_genes(science_plan, panel_map$gene, feature_maps[["science"]], science_duplicate_policy, panel_map)
    cat("Panel-ID-aligned Science genes:", nrow(a), "\n")
    print(a[a$gene %in% c("RGS5", "SHOX", "CSF2RA", "IL3RA", "CD99", "IL9R"), ], row.names = FALSE)
  }
  quit(save = "no")
}
if (any(c("--dry-run", "--list") %in% args)) {
  ok <- compact_dry_run("Reference similarity (paths and optional small panel JSON; no reference/query RDS loaded)", inputs,
                        outputs = out, checks = c(new_output_directory = !file.exists(out),
                                                  original_reference_paths_supplied = valid_text(paths),
                                                  published_label_columns_supplied = valid_text(labels),
                                                  science_reader_ready = !science_h5 || (nzchar(Sys.which("h5dump")) && requireNamespace("jsonlite", quietly = TRUE))))
  quit(save = "no", status = if (ok) 0L else 1L)
}
if ("science" %in% refs && !nzchar(paths[["science"]])) stop("Supply --science=ORIGINAL_RDS. Original Science provenance has not been established; no processed fallback is used.")
needed <- if ("--inspect" %in% args) inputs[names(inputs) != "query"] else inputs
if (any(!file.exists(needed))) stop("Missing input(s): ", paste(names(needed)[!file.exists(needed)], collapse = ", "))
# Inspect saved Seurat slots before S4 inheritance dispatch. Original objects
# can carry class package='Seurat', which otherwise autoloads analysis DLLs.
is_sce <- function(obj) {
  if (is.data.frame(attr(obj, "meta.data", exact = TRUE))) return(FALSE)
  inherits(obj, "SummarizedExperiment")
}
object_metadata <- function(obj) {
  md <- attr(obj, "meta.data", exact = TRUE)
  if (is.data.frame(md)) return(md)
  if (is_sce(obj)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Missing SummarizedExperiment in project renv.")
    return(as.data.frame(SummarizedExperiment::colData(obj)))
  }
  stop("Expected an original Seurat or SingleCellExperiment/SummarizedExperiment object.")
}
# Avoid current Seurat container methods on original objects missing newer slots.
object_cell_ids <- function(obj) {
  if (is.data.frame(attr(obj, "meta.data", exact = TRUE))) return(rownames(attr(obj, "meta.data", exact = TRUE)))
  colnames(obj)
}
object_assays <- function(obj) {
  if (is.data.frame(attr(obj, "meta.data", exact = TRUE))) return(names(attr(obj, "assays", exact = TRUE)))
  if (is_sce(obj)) return(SummarizedExperiment::assayNames(obj))
  stop("Unsupported object for assay inventory.")
}
reference_cohort_selection <- function(md, ref) {
  if (ref == "aldinger") {
    if (!"age" %in% names(md)) stop("Original Aldinger metadata lacks age.")
    stage <- as.character(md$age)
    if (anyNA(stage) || any(!grepl("^[0-9]+(\\.[0-9]+)? PCW$", stage))) stop("Unrecognized original Aldinger age values; expected explicit PCW strings.")
    age <- as.numeric(sub(" PCW$", "", stage))
    keep <- age >= 9 & age <= 21
    description <- "Original Aldinger age 9–21 PCW inclusive; original figure_clusters retained."
  } else if (ref == "sepp") {
    required <- c("Stage", "Tissue")
    if (!all(required %in% names(md))) stop("Original Sepp metadata lacks Stage/Tissue.")
    stage <- as.character(md$Stage)
    fetal <- c("7 wpc", "8 wpc", "9 wpc", "11 wpc", "17 wpc", "20 wpc")
    known <- c(fetal, "newborn", "infant", "toddler", "adult")
    if (anyNA(stage) || any(!stage %in% known)) stop("Unrecognized original Sepp Stage values; review cohort before analysis.")
    tissue <- as.character(md$Tissue)
    cerebellar <- c("cerebellar cx", "cerebellar dn", "cerebellar flocculonodular lobe", "cerebellar vermis", "cerebellum")
    if (anyNA(tissue) || any(!tissue %in% cerebellar)) stop("Original Sepp Tissue values differ from reviewed cerebellar inventory.")
    keep <- stage %in% fetal
    description <- "Original Sepp fetal cerebellar stages 7/8/9/11/17/20 wpc; postnatal/adult cells excluded."
  } else stop("No approved RDS cohort rule for reference: ", ref)
  if (!any(keep)) stop("Approved reference cohort contains no cells: ", ref)
  list(keep = keep, stage = stage, description = description)
}
reference_identity_groups <- function(md, ref, label_column) {
  if (!label_column %in% names(md)) stop(ref, " lacks original identity column: ", label_column)
  published <- as.character(md[[label_column]])
  missing <- is.na(published) | !nzchar(trimws(published))
  groups <- published
  if (any(missing)) {
    if (ref != "sepp" || !label_column %in% c("cell_type", "precisest_label")) stop(ref, ": missing/blank original identities in ", label_column)
    groups[missing] <- NA_character_
    message("Sepp: excluding ", sum(missing), " cells without original ", label_column,
            " from matching; source annotations are unchanged.")
  }
  list(groups = groups, published = published, unannotated = missing)
}
is_h5 <- function(p) grepl("\\.(h5|hdf5|h5ad)$", p, ignore.case = TRUE)
inspect_h5 <- function(p) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("HDF5 inspection needs hdf5r in the project renv; no package is installed automatically.")
  f <- hdf5r::H5File$new(p, mode = "r")
  on.exit(f$close_all())
  cat("\nHDF5 schema only (no expression matrix loaded):", p, "\n")
  print(f$ls(recursive = FALSE))
  for (g in c("obs", "var", "matrix", "matrix/features", "raw", "layers")) {
    if (grepl("/", g, fixed = TRUE)) {
      parent <- sub("/[^/]+$", "", g)
      if (!f$exists(parent) || !inherits(f[[parent]], "H5Group")) next
    }
    if (f$exists(g) && inherits(f[[g]], "H5Group")) {
      cat("\nGroup:", g, "\n"); print(f[[g]]$ls(recursive = FALSE))
    }
  }
  cat("Inspect original annotation names and raw-count storage before choosing a reader.\n")
}
if ("--inspect" %in% args) {
  for (nm in refs) {
    if (nm == "science" && science_h5) {
      science_plan <- plan_science_original(paths[[nm]], science_block_rows)
      next
    }
    if (is_h5(inputs[[nm]])) { inspect_h5(inputs[[nm]]); next }
    obj <- readRDS(inputs[[nm]]); md <- object_metadata(obj)
    cat("\n", nm, ": ", nrow(md), " cells\nAssays: ",
        paste(object_assays(obj), collapse = ", "), "\nMetadata: ",
        paste(colnames(md), collapse = ", "), "\n", sep = "")
    if (is_sce(obj)) {
      cat("Feature metadata:", paste(colnames(SummarizedExperiment::rowData(obj)), collapse = ", "), "\n")
      cat("Example feature IDs:", paste(head(rownames(obj)), collapse = ", "), "\n")
    } else {
      for (a in object_assays(obj)) {
        aa <- attr(obj, "assays", exact = TRUE)[[a]]
        fm <- attr(aa, "meta.features", exact = TRUE)
        if (is.null(fm)) fm <- attr(aa, "meta.data", exact = TRUE)
        cat("Assay", a, "feature metadata:", paste(colnames(fm), collapse = ", "), "\n")
      }
    }
    inventory_columns <- unique(c(labels[[nm]], grep("(age|pcw|gest|week|stage|region|tissue|donor|sample|patient|individual|subject|orig.ident|cell_type|figure_clusters|orig_cluster|orig_sub_cluster|precisest_label|dev_state|subtype|^type$|^sex$)", names(md), value = TRUE, ignore.case = TRUE)))
    for (k in intersect(inventory_columns, names(md))) {
      if (!is.atomic(md[[k]])) next
      tab <- table(md[[k]], useNA = "ifany")
      cat("\nMetadata values:", k, "(", length(tab), "distinct values; first 100 shown)\n")
      print(head(tab, 100L))
    }
    if (nm %in% c("aldinger", "sepp")) {
      cohort <- reference_cohort_selection(md, nm)
      cat("\nAPPROVED ORIGINAL COHORT:", cohort$description, "\n")
      print(table(stage = cohort$stage, selected = cohort$keep), max = 1000L)
      for (k in intersect(c("figure_clusters", "cell_type", "orig_cluster", "orig_sub_cluster", "precisest_label"), names(md))) {
        bad <- is.na(md[[k]]) | !nzchar(trimws(as.character(md[[k]])))
        cat("Selected cells missing", k, ":", sum(bad & cohort$keep), "\n")
      }
    }
    stage_column <- if ("Stage" %in% names(md)) "Stage" else if ("age" %in% names(md)) "age" else NULL
    annotation_columns <- intersect(c("figure_clusters", "cell_type", "orig_cluster", "orig_sub_cluster", "precisest_label"), names(md))
    if (!is.null(stage_column)) {
      cat("\nANNOTATION COMPLETENESS BY ORIGINAL STAGE (counts only; no selection applied):\n")
      completeness <- lapply(sort(unique(as.character(md[[stage_column]]))), function(stage) {
        ii <- which(as.character(md[[stage_column]]) == stage)
        row <- data.frame(stage = stage, cells = length(ii))
        for (k in annotation_columns) row[[paste0("missing_", k)]] <- sum(is.na(md[[k]][ii]) | !nzchar(trimws(as.character(md[[k]][ii]))))
        row
      })
      if (length(completeness)) print(do.call(rbind, completeness), row.names = FALSE)
    }
    rm(obj, md); gc()
  }
  panel_json <- if (nzchar(panel_file)) panel_file else resolve_config_path(config$inputs$xenium_panel_json_example, config)
  if (file.exists(panel_json) && requireNamespace("jsonlite", quietly = TRUE)) {
    cat("\nORIGINAL XENIUM PANEL GENE-IDENTIFIER FIELD PREVIEW:\nSource:", panel_json, "\n")
    panel_document <- jsonlite::fromJSON(panel_json, simplifyVector = FALSE)
    gene_targets <- Filter(function(x) identical(x$type$descriptor, "gene"), panel_document$payload$targets)
    for (x in head(gene_targets, 3L)) utils::str(x$type$data, max.level = 1L)
    ambiguous_panel_genes <- c("RGS5", "SHOX", "CSF2RA", "IL3RA", "CD99", "IL9R")
    for (symbol in ambiguous_panel_genes) {
      matching_targets <- Filter(function(x) identical(x$type$data$name, symbol), gene_targets)
      cat("\nOriginal panel annotation for", symbol, ":", length(matching_targets), "targets\n")
      for (x in matching_targets) utils::str(x$type$data, max.level = 1L)
    }
  } else {
    cat("\nOriginal panel JSON unavailable for gene-ID preview. Sepp Ensembl-to-query mapping remains to be supplied.\n")
  }
  quit(save = "no")
}
for (p in c("SeuratObject", "Matrix")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing dependency in project renv: ", p)
}
if (!valid_text(labels)) stop("Supply --science-label=PUBLISHED_COLUMN (use --inspect first). No relabeled annotation is substituted.")
if (any(vapply(paths[setdiff(refs, if (science_h5) "science" else character())], is_h5, logical(1)))) stop("Only the audited original Science /shoji HDF5 format is supported for HDF5 comparisons.")
# Never replace previous results, including a partial run; choose a new --out.
if (file.exists(out)) stop("Output already exists. Review it and choose a new --out directory: ", out)
set.seed(seed)
started <- Sys.time()
if (science_h5) science_plan <- plan_science_original(paths[["science"]], science_block_rows)

metadata <- function(obj, columns, dataset) {
  ids <- object_cell_ids(obj); md <- object_metadata(obj)
  if (!length(ids) || anyDuplicated(ids) || !identical(rownames(md), ids)) stop(dataset, ": invalid cell/metadata alignment.")
  missing <- setdiff(columns, names(md))
  if (length(missing)) stop(dataset, ": missing metadata: ", paste(missing, collapse = ", "), ". Use --inspect.")
  for (k in columns) if (!valid_text(md[[k]])) stop(dataset, ": missing/blank ", k)
  md
}
# Aggregate each count layer separately, so v5 split layers need not be joined.
# Every cell must occur in exactly one count layer. Dense cells x genes never formed.
aggregate_counts <- function(obj, assay, grouping, selected_cells = object_cell_ids(obj)) {
  source_cells <- object_cell_ids(obj)
  if (!length(selected_cells) || anyDuplicated(selected_cells) || anyNA(match(selected_cells, source_cells)) ||
      length(grouping) != length(selected_cells) || anyNA(grouping) || any(grouping < 1L)) stop("Invalid selected-cell aggregation mapping.")
  if (is_sce(obj)) {
    if (!assay %in% object_assays(obj)) stop("Missing original SCE counts assay: ", assay)
    m <- SummarizedExperiment::assay(obj, assay)
    if (!identical(colnames(m), source_cells) || anyDuplicated(rownames(m)) || !valid_text(rownames(m))) stop("Invalid original SCE matrix IDs.")
    if (!identical(selected_cells, colnames(m))) m <- m[, match(selected_cells, colnames(m)), drop = FALSE]
    vals <- if (inherits(m, "sparseMatrix")) m@x else as.vector(m)
    if (any(!is.finite(vals)) || any(vals < 0) || any(abs(vals - round(vals)) > 1e-6)) stop("Expected original raw integer counts.")
    design <- Matrix::sparseMatrix(i = seq_len(ncol(m)), j = grouping, x = 1,
                                  dims = c(ncol(m), max(grouping)))
    return(as.matrix(m %*% design))
  }
  if (!assay %in% object_assays(obj)) stop("Missing assay: ", assay)
  a <- attr(obj, "assays", exact = TRUE)[[assay]]; ids <- selected_cells
  # Original Assay counts are read directly; v5 split layers use assay methods.
  legacy_counts <- attr(a, "counts", exact = TRUE)
  layers <- if (!is.null(legacy_counts)) "counts" else SeuratObject::Layers(a, search = "^counts($|\\.)")
  if (!length(layers)) stop("No raw count layer in assay ", assay)
  features <- if (!is.null(legacy_counts)) rownames(legacy_counts) else rownames(a)
  if (!valid_text(features) || anyDuplicated(features)) stop("Invalid/duplicate feature IDs in ", assay)
  ng <- max(grouping)
  sums <- matrix(0, length(features), ng, dimnames = list(features, as.character(seq_len(ng))))
  seen <- integer(length(ids))
  for (layer in layers) {
    m <- if (!is.null(legacy_counts)) legacy_counts else SeuratObject::LayerData(a, layer = layer)
    if (anyDuplicated(colnames(m)) || anyNA(match(colnames(m), source_cells))) stop("Original count-layer cell IDs are invalid.")
    keep <- colnames(m) %in% selected_cells
    if (!any(keep)) next
    if (!all(keep)) m <- m[, keep, drop = FALSE]
    ci <- match(colnames(m), ids); gi <- match(rownames(m), features)
    if (anyNA(ci) || anyDuplicated(ci) || anyNA(gi) || anyDuplicated(gi)) stop("Count layer IDs do not align.")
    if (!setequal(rownames(m), features)) stop("Count layers have different feature sets; reconcile explicitly.")
    vals <- if (inherits(m, "sparseMatrix")) m@x else as.vector(m)
    if (any(!is.finite(vals)) || any(vals < 0) || any(abs(vals - round(vals)) > 1e-6)) stop("Expected finite, nonnegative raw integer counts.")
    seen[ci] <- seen[ci] + 1L
    design <- Matrix::sparseMatrix(i = seq_along(ci), j = grouping[ci], x = 1,
                                    dims = c(length(ci), ng))
    sums[gi, ] <- sums[gi, , drop = FALSE] + as.matrix(m %*% design)
  }
  if (any(seen != 1L)) stop("Each object cell must occur in exactly one counts layer.")
  sums
}
# Pairwise common genes define the normalization denominator in BOTH platforms.
# Equal-weight sample profiles avoid letting the largest Xenium sample dominate.
log_cpm <- function(m) {
  totals <- colSums(m)
  if (any(!is.finite(totals)) || any(totals <= 0)) stop("Zero/nonfinite pseudobulk library on common genes.")
  log1p(sweep(m, 2L, totals / 1e6, "/"))
}
score <- function(q, r) {
  ans <- stats::cor(q, r, method = "spearman")
  if (any(!is.finite(ans))) stop("Undefined correlation: constant profile on comparison genes.")
  ans
}
write_table <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE)
align_features <- function(obj, assay, counts, ref) {
  ids <- rownames(counts); comparison <- ids
  if (nzchar(feature_columns[[ref]])) {
    fm <- if (is_sce(obj)) as.data.frame(SummarizedExperiment::rowData(obj)) else {
      a <- attr(obj, "assays", exact = TRUE)[[assay]]
      saved <- attr(a, "meta.features", exact = TRUE)
      if (is.data.frame(saved)) saved else a[[]]
    }
    k <- feature_columns[[ref]]
    if (!k %in% names(fm) || !identical(rownames(fm), ids)) stop(ref, ": feature metadata column/alignment missing: ", k)
    comparison <- as.character(fm[[k]])
  } else if (nzchar(feature_maps[[ref]])) {
    m <- read.csv(feature_maps[[ref]], stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("reference_gene", "query_gene") %in% names(m))) stop("Feature map needs reference_gene,query_gene columns.")
    if (!valid_text(m$reference_gene) || !valid_text(m$query_gene) || anyDuplicated(m$reference_gene) || anyDuplicated(m$query_gene)) stop("Feature map must be nonblank and one-to-one; resolve ambiguous mappings explicitly.")
    comparison <- m$query_gene[match(ids, m$reference_gene)]
  } else if (ref == "sepp" && !is.null(panel_map)) {
    a <- panel_reference_alignment(ids, panel_map)
    comparison <- a$comparison_feature
    write_table(a, paste0(ref, "_ensembl_alignment_audit.csv"))
  }
  usable <- !is.na(comparison) & nzchar(trimws(comparison))
  if (anyDuplicated(comparison[usable])) stop(ref, ": duplicate comparison gene IDs; provide an explicit one-to-one feature map.")
  write_table(data.frame(original_feature = ids, comparison_feature = comparison, usable = usable), paste0(ref, "_feature_alignment.csv"))
  counts <- counts[usable, , drop = FALSE]; rownames(counts) <- comparison[usable]
  counts
}

message("Loading query; only compact count aggregates are retained.")
obj <- readRDS(query)
levels_q <- c("corrected_broad_label", "comb_subcluster")
md <- metadata(obj, c(sample_col, levels_q), "Xenium")
if (!setequal(as.character(md[[sample_col]]), samples)) stop("Xenium sample set differs from active manifest.")
# Enforce the existing approved map without changing labels or object.
source(here::here("scripts", "R", "corrected_broad_labels.R"))
if (!identical(as.character(md$corrected_broad_label), map_corrected_broad_labels(md$comb_subcluster))) {
  stop("Query corrected broad labels do not match the approved fine-to-broad map.")
}
for (level in levels_q) reference_similarity_order(unique(as.character(md[[level]])), level)
query_data <- list()
for (level in rev(levels_q)) {
  pairs <- unique(data.frame(sample = as.character(md[[sample_col]]), label = as.character(md[[level]])))
  # Integer grouping avoids delimiter collisions in biological labels or sample IDs.
  grouping <- integer(nrow(md))
  for (s in unique(pairs$sample)) {
    ii <- which(md[[sample_col]] == s); jj <- which(pairs$sample == s)
    grouping[ii] <- jj[match(as.character(md[[level]][ii]), pairs$label[jj])]
  }
  if (anyNA(grouping) || any(grouping < 1L)) stop("Invalid pseudobulk mapping.")
  pairs$cells <- tabulate(grouping, nbins = nrow(pairs))
  if (level == "comb_subcluster") {
    counts <- aggregate_counts(obj, qassay, grouping)
  } else {
    fine <- query_data[["comb_subcluster"]]
    fine_broad <- map_corrected_broad_labels(fine$pairs$label)
    mapping <- integer(nrow(fine$pairs))
    for (sample_id in unique(pairs$sample)) {
      ii <- which(fine$pairs$sample == sample_id); jj <- which(pairs$sample == sample_id)
      mapping[ii] <- jj[match(fine_broad[ii], pairs$label[jj])]
    }
    if (anyNA(mapping) || any(mapping < 1L)) stop("Fine-to-broad pseudobulk mapping failed.")
    design <- Matrix::sparseMatrix(i = seq_along(mapping), j = mapping, x = 1,
                                   dims = c(length(mapping), nrow(pairs)))
    counts <- as.matrix(fine$counts %*% design)
    if (!identical(as.numeric(rowsum(fine$pairs$cells, mapping, reorder = TRUE)), as.numeric(pairs$cells))) stop("Fine-to-broad cell counts differ.")
  }
  query_data[[level]] <- list(counts = counts, pairs = pairs)
}
query_data <- query_data[levels_q]
rm(obj, md, counts); gc()
if (science_h5) {
  # Reject ambiguous feature matches before any reference expression pass.
  science_alignment <- resolve_science_genes(science_plan, rownames(query_data[[1L]]$counts), feature_maps[["science"]], science_duplicate_policy, panel_map)
  rm(science_alignment)
}
if (!dir.create(out, recursive = TRUE)) stop("Cannot create output directory.")
writeLines("INCOMPLETE: successful runs contain COMPLETED.txt", file.path(out, "RUN_STATUS.txt"))
for (level in levels_q) write_table(query_data[[level]]$pairs, paste0(level, "_sample_counts.csv"))

for (ref in refs) {
  message("Loading reference: ", ref)
  if (ref == "science" && science_h5) {
    result <- aggregate_science_original(science_plan, rownames(query_data[[1L]]$counts), out, feature_maps[["science"]], science_duplicate_policy, panel_map)
    rc <- result$counts; rl <- colnames(rc)
    if (labels[["science"]] == "CellClass") {
      class_names <- science_plan$lookup$Class[match(rl, as.character(science_plan$lookup$ClusterID))]
      if (anyNA(class_names)) stop("Science original class mapping failed.")
      write_table(data.frame(original_cluster_id = rl, original_cell_class = class_names), "science_comparison_label_mapping.csv")
      rc <- t(rowsum(t(rc), group = class_names, reorder = TRUE)); rl <- colnames(rc)
      class_cells <- rowsum(result$cell_counts$cells, group = class_names, reorder = TRUE)
      write_table(data.frame(reference_label = rownames(class_cells), cells = as.numeric(class_cells)), "science_cell_counts.csv")
      reference_status <- rep("published_cell_class", length(rl))
    } else {
      reference_status <- rep("original_cluster_id", length(rl))
      write_table(result$cell_counts, "science_cell_counts.csv")
    }
    rm(result)
  } else {
    obj <- readRDS(paths[[ref]])
    selected_cells <- object_cell_ids(obj)
    if (ref %in% c("aldinger", "sepp")) {
      source_md <- object_metadata(obj)
      cohort <- reference_cohort_selection(source_md, ref)
      if (!valid_text(rownames(source_md)) || anyDuplicated(rownames(source_md))) stop("Invalid original reference metadata cell IDs.")
      write_table(data.frame(CellID = rownames(source_md), original_stage = cohort$stage,
                             selected = cohort$keep), paste0(ref, "_reference_cohort_cells.csv"))
      write_table(as.data.frame(table(original_stage = cohort$stage, selected = cohort$keep)), paste0(ref, "_reference_cohort_counts.csv"))
      selected_cells <- rownames(source_md)[cohort$keep]
      rm(source_md, cohort)
    }
    md <- metadata(obj, character(), ref)[selected_cells, , drop = FALSE]
    identities <- reference_identity_groups(md, ref, labels[[ref]])
    if (ref == "sepp") {
      excluded_columns <- intersect(c(labels[[ref]], "orig_cluster", "orig_sub_cluster", "Stage", "Tissue"), names(md))
      write_table(data.frame(CellID = rownames(md)[identities$unannotated],
                             reason = rep("Missing original named annotation", sum(identities$unannotated)),
                             md[identities$unannotated, excluded_columns, drop = FALSE], check.names = FALSE), "sepp_excluded_unannotated_cells.csv")
      retained <- !identities$unannotated
      if (!any(retained)) stop("No originally named Sepp cells remain.")
      md <- md[retained, , drop = FALSE]; selected_cells <- selected_cells[retained]
      identities <- lapply(identities, `[`, retained)
    }
    rl <- sort(unique(identities$groups))
    grouping <- match(identities$groups, rl)
    unannotated_groups <- unique(identities$groups[identities$unannotated])
    reference_status <- ifelse(rl %in% unannotated_groups, "unannotated_original_cluster", "published_annotation")
    original_columns <- intersect(c(labels[[ref]], "cell_type", "precisest_label", "orig_cluster", "orig_sub_cluster", "figure_clusters"), names(md))
    write_table(data.frame(CellID = rownames(md), comparison_group = identities$groups,
                           published_identity_missing = identities$unannotated,
                           md[, original_columns, drop = FALSE], check.names = FALSE),
                paste0(ref, "_original_cell_annotations.csv"))
    rc <- aggregate_counts(obj, assays[[ref]], grouping, selected_cells); colnames(rc) <- rl
    rc <- align_features(obj, assays[[ref]], rc, ref)
    write_table(data.frame(reference_label = rl, annotation_status = reference_status,
                           cells = tabulate(grouping, length(rl))), paste0(ref, "_cell_counts.csv"))
    rm(obj, md, identities)
  }
  if (length(rl) < 2L) stop(ref, " needs at least two labels for specificity comparisons.")
  gc()
  for (level in levels_q) {
    qd <- query_data[[level]]; common <- intersect(rownames(qd$counts), rownames(rc))
    if (length(common) < 3L) stop("Fewer than 3 matching genes for ", ref, ". Inspect original gene IDs; provide --", ref, "-feature-column or --", ref, "-feature-map for explicit gene alignment. Original objects are never renamed.")
    q <- log_cpm(qd$counts[common, , drop = FALSE]); r <- log_cpm(rc[common, , drop = FALSE])
    ql <- sort(unique(qd$pairs$label))
    qi <- lapply(ql, function(x) which(qd$pairs$label == x)); names(qi) <- ql
    qm <- vapply(qi, function(i) rowMeans(q[, i, drop = FALSE]), numeric(length(common)))
    # Remove genes constant across either dataset's label profiles; record every gene.
    keep <- apply(qm, 1L, stats::sd) > 0 & apply(r, 1L, stats::sd) > 0
    if (sum(keep) < 3L) stop("Fewer than 3 variable shared genes for ", ref, " / ", level)
    stem <- paste(ref, level, sep = "_")
    write_table(data.frame(gene = common, variable_in_both = keep), paste0(stem, "_genes.csv"))
    q <- q[keep, , drop = FALSE]; qm <- qm[keep, , drop = FALSE]; r <- r[keep, , drop = FALSE]
    qb <- rowMeans(qm); rb <- rowMeans(r)
    raw <- score(qm, r)
    fast_specific <- make_reference_spearman(r - rb)
    specific <- fast_specific(qm - qb)
    full <- expand.grid(query_label = ql, reference_label = rl, stringsAsFactors = FALSE)
    full$spearman <- as.vector(raw)
    full$specificity_spearman <- as.vector(specific)
    full$reference_annotation_status <- reference_status[match(full$reference_label, rl)]
    write_table(full, paste0(stem, "_all_correlations.csv"))
    summaries <- list()
    for (k in seq_along(ql)) {
      order_k <- order(-specific[k, ], rl)
      top <- order_k[1L]; second <- order_k[2L]; ii <- qi[[k]]
      n <- length(ii)
      per_sample <- fast_specific(q[, ii, drop = FALSE] - qb)
      sample_best <- apply(per_sample, 1L, function(x) rl[order(-x, rl)[1L]])
      write_table(data.frame(sample = qd$pairs$sample[ii], best_reference_label = sample_best,
                             specificity_to_group_best = per_sample[, top]),
                  paste0(stem, "_query", k, "_sample_matches.csv"))
      ci <- c(NA_real_, NA_real_); margin_ci <- ci; win <- NA_real_
      # Resample biological samples within the target population. Reference and
      # background centering stay fixed: conditional stability, not full uncertainty.
      if (n >= 3L && B > 0L) {
        # Batch query bootstraps; fixed reference ranks are computed only once.
        # Same sample draws/order and seed as the scalar implementation.
        boot_profiles <- replicate(B, rowMeans(q[, sample(ii, n, replace = TRUE), drop = FALSE]))
        boot <- t(fast_specific(boot_profiles - qb))
        ci <- stats::quantile(boot[top, ], c(.025, .975), names = FALSE)
        margin_ci <- stats::quantile(boot[top, ] - boot[second, ], c(.025, .975), names = FALSE)
        # Ties split their weight equally; do not fabricate a certain winner.
        win <- mean(apply(boot, 2L, function(z) {
          winners <- which(abs(z - max(z)) < 1e-12)
          if (top %in% winners) 1 / length(winners) else 0
        }))
      }
      summaries[[k]] <- data.frame(query_label = ql[k], n_samples = n,
        cells = sum(qd$pairs$cells[ii]), reference_label = rl[top], runner_up = rl[second],
        reference_annotation_status = reference_status[top],
        specificity_spearman = specific[k, top], spearman_at_specificity_best = raw[k, top],
        raw_best_label = rl[order(-raw[k, ], rl)[1L]],
        top2_margin = specific[k, top] - specific[k, second],
        conditional_ci_low = ci[1L], conditional_ci_high = ci[2L],
        conditional_margin_ci_low = margin_ci[1L], conditional_margin_ci_high = margin_ci[2L],
        bootstrap_best_fraction = win, sample_best_fraction = mean(sample_best == rl[top]),
        query_genes = nrow(qd$counts), reference_genes = nrow(rc), shared_genes = length(common),
        variable_shared_genes = sum(keep), uncertainty = if (n < 3L) "insufficient_samples" else if (B == 0L) "disabled" else "query_samples_only_fixed_reference_and_background")
    }
    write_table(do.call(rbind, summaries), paste0(stem, "_best_matches.csv"))
    display <- rl
    if (ref == "science" && labels[["science"]] == "Clusters") {
      display <- paste0(science_plan$lookup$Class[match(rl, as.character(science_plan$lookup$ClusterID))], " [", rl, "]")
    }
    plot_order <- write_reference_heatmaps(raw, specific, file.path(out, paste0(stem, "_heatmaps.pdf")),
                             paste(ref, level), display, query_level = level)
    write_table(data.frame(top_to_bottom = seq_along(plot_order), query_label = plot_order), paste0(stem, "_heatmap_row_order.csv"))
  }
  rm(rc); gc()
}
provenance <- data.frame(dataset = names(inputs), path = unname(inputs),
                         size_bytes = file.info(inputs)$size, mtime = as.character(file.info(inputs)$mtime))
write_table(provenance, "input_provenance.csv")
if (!is.null(panel_map)) {
  write_table(panel_map, "original_panel_gene_ids.csv")
  write_table(data.frame(path = panel_file, md5 = unname(tools::md5sum(panel_file))), "original_panel_provenance.csv")
  file.copy(panel_file, file.path(out, "original_panel_snapshot.json"))
}
write_table(data.frame(reference = refs, assay = assays, label_column = labels,
                       feature_column = feature_columns, feature_map = feature_maps), "reference_settings.csv")
writeLines(c(paste("Started:", started), paste("Finished:", Sys.time()),
             paste("Seed:", seed, "Bootstrap:", B, "Manifest samples:", length(samples)),
             paste("Query assay:", qassay, "sample column:", sample_col),
             paste("Science duplicate-symbol policy:", science_duplicate_policy),
             "Panel ID mode: Sepp/Science align by original panel Ensembl IDs, removing terminal numeric accession versions only; separate _PAR_Y and other unmatched records are not combined. Explicit reference feature maps take precedence.",
             "Gene matching: exact IDs or explicit one-to-one comparison alignment; pairwise common-gene CPM denominator. Original objects are unchanged.",
             "Query profiles: equally weighted log1p CPM sample-label count pseudobulks.",
             "References: pooled label count pseudobulks; cell weighted within reference label.",
             "Specificity: Spearman after subtracting per-gene mean across labels within each dataset.",
             "Bootstrap: query samples within label; fixed reference and fixed background; 95% percentile intervals.",
             "Heatmaps: raw correlations use observed-range colors (blue may be positive); specificity uses symmetric zero-centered colors; legends show actual Spearman rho.",
             "No p-values, automatic found threshold, label changes, object writes, or inferred ontology mappings.",
             "Correlations use only variable shared genes. Missing sample-label combinations are absent, not zero.",
             "References used in annotation: agreement is not independent validation.",
             "Original source IDs/labels are retained within each approved cohort. Science HDF5 uses Region=Cerebellum and source Age 9–15 inclusive, original CellClass (or explicit Clusters)/Gene, sum of spliced+unspliced UMIs from Expression.",
             "Approved broader fetal coverage: Aldinger age 9–21 PCW inclusive; Sepp original fetal stages 7,8,9,11,17,20 wpc in reviewed cerebellar tissues. Cell/cluster labels are not changed; Sepp cells lacking original named annotations are excluded from comparison and audited, as approved 2026-10-22.",
             "Science summaries are recomputed for selected cells in bounded blocks, never taken from whole-reference MeanExpression. Source objects and published labels are unchanged.",
             capture.output(sessionInfo())), file.path(out, "methods_and_session.txt"))
file.copy(here::here("scripts", "xenium_compare_reference_similarity.R"), file.path(out, "script_snapshot.R"))
file.copy(attr(config, "config_path"), file.path(out, "config_snapshot.yml"))
file.copy(resolve_config_path(config$manifests$samples, config), file.path(out, "sample_manifest_snapshot.csv"))
file.copy(here::here("scripts", "R", "science_original_reference.R"), file.path(out, "science_reader_snapshot.R"))
file.copy(here::here("scripts", "R", "reference_similarity_plots.R"), file.path(out, "plot_and_score_helpers_snapshot.R"))
file.copy(here::here("scripts", "color_palette.R"), file.path(out, "color_palette_snapshot.R"))
writeLines("Complete; review best_matches tables alongside all_correlations and gene coverage.", file.path(out, "COMPLETED.txt"))
message("Completed in ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " minutes: ", out)
