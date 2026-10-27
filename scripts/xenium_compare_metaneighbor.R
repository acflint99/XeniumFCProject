#!/usr/bin/env Rscript
# Original-reference MetaNeighborUS, broad + subcluster comparisons per reference.
# Expected HPC runtime: 0.5-4 h/reference, low confidence; no local full analysis.
# --dry-run / --list: <1 minute, paths/dependencies only; no large RDS/HDF5 reads.
args <- commandArgs(trailingOnly = TRUE)
flags <- c("--dry-run", "--list", "--legacy-query-labels", "--help")
keys <- c("reference", "query", "aldinger", "sepp", "science", "panel-json", "out", "seed", "feature-cells", "max-dense-gib", "science-block-rows")
opt <- list()
for (a in args) {
  if (a %in% flags) next
  if (!grepl("^--[^=]+=.+$", a)) stop("Use --option=value: ", a)
  k <- sub("=.*$", "", substring(a,3))
  if (!k %in% keys || !is.null(opt[[k]])) stop("Unknown/repeated option: ", a)
  opt[[k]] <- sub("^[^=]+=", "", a)
}
if ("--help" %in% args) {
  cat("MetaNeighborUS against original references; no integration or reference relabeling.\n",
      "--reference=aldinger|sepp|science (one original reference; both Xenium label levels)\n",
      "--list gives the three reference-array tasks; --dry-run checks paths/packages only.\n",
      "--panel-json=PATH --out=NEW_DIRECTORY; other keys: ", paste(keys, collapse=", "), "\n",
      "Original labels: figure_clusters, cell_type, CellClass. Original approved age/tissue cohorts.\n",
      "All eligible named cells are considered; rank-constant cells are excluded only from MetaNeighbor scoring and audited. feature-cells=10000 caps feature-selection sampling only.\n",
      "--legacy-query-labels preserves saved legacy Endothelial/Pericytes query labels and marks every output as legacy; it does not modify the object.\n",
      "seed=42; max-dense-gib=32 bounds a single dense rank matrix, not total job memory.\n",
      "science-block-rows=512 controls Science HDF5 I/O block size only (1–1024); it does not change selected cells or scores.\n",
      "Uses official fast_version=TRUE, one_vs_best=FALSE, symmetric_output=FALSE.\n",
      "No bootstrapping, automatic found threshold, label changes or package installation.\n", sep="")
  quit(save="no")
}
for (p in c("here", "yaml")) if (!requireNamespace(p,quietly=TRUE)) stop("Missing project dependency: ", p)
here::i_am("scripts/xenium_compare_metaneighbor.R")
source(here::here("scripts", "R", "config.R"))
config <- load_pipeline_config(); manifest <- load_sample_manifest(config)
value <- function(k,d) if (is.null(opt[[k]])) d else opt[[k]]
path <- function(p) if (grepl("^/",p)) p else here::here(p)
refs <- c("aldinger", "sepp", "science")
paths <- c(aldinger=file.path(config$inputs$published_data_root,"AldingerFC","Aldinger_seurat.rds"),
           sepp=file.path(config$inputs$published_data_root,"SeppFC","Sepp_hum_sce_final.rds"),
           science=file.path(config$inputs$published_data_root,"ScienceBraunFC","Science_HumanFetalBrainPool.h5"))
paths <- setNames(vapply(refs,function(r) path(value(r,paths[[r]])), ""),refs)
query <- path(value("query",file.path(config$project$outputs_dir,"xenium","vz_rl","nonspatial","03_merged","rds","Xenium_VZRL_merged.rds")))
panel_file <- path(value("panel-json","data/FCXeniumProject/FB124_X_G/gene_panel.json"))
ref <- value("reference", "aldinger")
if (!ref %in% refs) stop("Invalid --reference.")
out_root <- path(value("out",file.path(config$project$outputs_dir,"xenium","metaneighbor_original")))
out <- file.path(out_root, ref)
required <- c("Matrix","SeuratObject","SummarizedExperiment","jsonlite","MetaNeighbor")
available <- setNames(vapply(required,requireNamespace,TRUE,quietly=TRUE),required)
api <- available[["MetaNeighbor"]] && all(c("fast_version","symmetric_output","one_vs_best") %in% names(formals(MetaNeighbor::MetaNeighborUS))) &&
  all(c("downsampling_size","min_recurrence") %in% names(formals(MetaNeighbor::variableGenes)))
if ("--list" %in% args) {
  print(data.frame(task=seq_along(refs),reference=refs,path=unname(paths),label=c("figure_clusters","cell_type","CellClass")),row.names=FALSE)
  cat("Manifest samples:",nrow(manifest),"; reference array 1-3; both query levels per task.\n")
  quit(save="no")
}
inputs <- c(query=query, reference=paths[[ref]], original_panel=panel_file,
            palette=here::here("scripts","color_palette.R"), manifest=resolve_config_path(config$manifests$samples,config))
print(data.frame(package=required,available=available),row.names=FALSE)
print(data.frame(input=names(inputs),path=unname(inputs),exists=file.exists(inputs)),row.names=FALSE)
if ("--dry-run" %in% args) {
  ok <- compact_dry_run("MetaNeighborUS (paths/packages only; no expression or panel contents read)", inputs, out,
        c(dependencies=all(available), low_memory_API=api, new_output=!file.exists(out),
          science_h5dump=ref!="science" || nzchar(Sys.which("h5dump"))))
  if (!all(available)) cat("Missing packages:",paste(names(available)[!available],collapse=", "),". No packages are installed by this script.\n")
  quit(save="no",status=if(ok) 0L else 1L)
}
if (!all(available) || !api) stop("Required project packages / low-memory MetaNeighbor API unavailable. Run --dry-run; no package is installed automatically.")
if (any(!file.exists(inputs)) || file.exists(out)) stop("Missing input or protected output exists. Run --dry-run and select a new --out root.")
integer_option <- function(k,d) {
  x <- suppressWarnings(as.numeric(value(k,d)))
  if (!is.finite(x) || x < 1 || x != floor(x) || x > .Machine$integer.max) stop("Invalid --",k)
  as.integer(x)
}
seed <- integer_option("seed","42"); feature_cells <- integer_option("feature-cells","10000")
science_block_rows <- integer_option("science-block-rows","512")
if (science_block_rows > 1024L) stop("Invalid --science-block-rows (1–1024).")
max_dense <- suppressWarnings(as.numeric(value("max-dense-gib","32")))
if (!is.finite(max_dense) || max_dense <= 0) stop("Invalid dense-matrix memory guard.")
source(here::here("scripts","R","science_original_reference.R"))
source(here::here("scripts","R","metaneighbor_original_io.R"))
source(here::here("scripts","R","metaneighbor_results.R"))
source(here::here("scripts","R","reference_similarity_plots.R"))
source(here::here("scripts","R","corrected_broad_labels.R"))
set.seed(seed); started <- Sys.time()
legacy_query_labels <- "--legacy-query-labels" %in% args
input_info <- file.info(inputs)
expected_path <- resolve_config_path(config$inputs$xenium_panel_genes,config)
panel <- mn_panel(panel_file,if(file.exists(expected_path)) readRDS(expected_path) else NULL)
message("Loading original merged Xenium object.")
obj <- readRDS(query); qmd <- mn_metadata(obj)
levels_q <- c("corrected_broad_label","comb_subcluster")
if (!all(c("orig.ident",levels_q) %in% names(qmd)) || any(!vapply(qmd[,c("orig.ident",levels_q)],mn_text,TRUE))) stop("Required Xenium metadata is missing/blank.")
if (anyDuplicated(manifest$sample_id) || !mn_text(manifest$sample_id) || !setequal(qmd$orig.ident,manifest$sample_id)) stop("Query samples differ from authoritative manifest.")
if (!legacy_query_labels && !identical(as.character(qmd$corrected_broad_label),map_corrected_broad_labels(qmd$comb_subcluster))) stop("Query labels differ from approved fine-to-broad map. Use --legacy-query-labels only for the saved legacy Endothelial object.")
if (legacy_query_labels && !"Endothelial" %in% as.character(qmd$comb_subcluster)) stop("--legacy-query-labels requires the saved legacy Endothelial fine label; refusing an unnecessary legacy mode.")
for(level in levels_q) reference_similarity_order(unique(qmd[[level]]),level,legacy_query_labels=legacy_query_labels)
qx <- mn_extract_counts(obj,"Xenium",rownames(qmd),panel); qc <- qx$counts
rm(obj); gc()
if (!dir.create(out,recursive=TRUE)) stop("Cannot create protected output.")
heatmap_dir <- file.path(out,"heatmaps")
if (!dir.create(heatmap_dir,recursive=TRUE)) stop("Cannot create heatmap output directory.")
write <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE)
writeLines("INCOMPLETE; successful results contain COMPLETED.txt",file.path(out,"RUN_STATUS.txt"))
write(data.frame(comparison_mode=if (legacy_query_labels) "legacy_saved_Endothelial_Pericytes_query_labels" else "current_approved_query_labels",
                 legacy_query_labels=legacy_query_labels,
                 note=if (legacy_query_labels) "Saved query labels are compared as-is; standalone Endothelial and Pericytes are legacy and not harmonized to Vascular." else "Saved query labels match the approved fine-to-broad mapping."),"comparison_mode.csv")
write(qx$alignment,"xenium_feature_alignment.csv"); rm(qx)
message("Loading original reference: ",ref)
if (ref == "science") {
  plan <- plan_science_original(paths[[ref]],science_block_rows)
  rmd <- plan$cells; rownames(rmd) <- rmd$CellID
  rc <- mn_science_counts(plan,panel,out); original_label <- "CellClass"
  write(plan$summary,"science_cohort_summary.csv")
} else {
  obj <- readRDS(paths[[ref]]); original <- mn_metadata(obj)
  keep <- mn_cohort(original,ref)
  original_label <- if(ref=="aldinger") "figure_clusters" else "cell_type"
  if (!original_label %in% names(original)) stop("Missing original published label: ", original_label)
  named <- !is.na(original[[original_label]]) & nzchar(trimws(as.character(original[[original_label]])))
  if(ref=="aldinger" && any(keep & !named)) stop("Original Aldinger labels are missing.")
  write(data.frame(CellID=rownames(original),age_eligible=keep,named_annotation=named,selected=keep & named),"reference_cell_selection.csv")
  rmd <- original[keep & named,,drop=FALSE]
  rx <- mn_extract_counts(obj,if(ref=="aldinger") "RNA" else "umi",rownames(rmd),panel,ensembl=ref=="sepp")
  rc <- rx$counts; write(rx$alignment,"reference_feature_alignment.csv")
  rm(obj,original,rx); gc()
}
if (!identical(colnames(qc),rownames(qmd)) || !identical(colnames(rc),rownames(rmd))) stop("Expression/metadata IDs differ.")
if (!mn_text(rmd[[original_label]]) || any(grepl("|",rmd[[original_label]],fixed=TRUE)) || any(vapply(qmd[,levels_q],function(x) any(grepl("|",x,fixed=TRUE)),TRUE))) stop("Invalid published labels or reserved '|' delimiter in labels.")
flat_columns <- vapply(rmd,function(x) is.atomic(x) && is.null(dim(x)),TRUE)
if(any(!flat_columns)) write(data.frame(field=names(rmd)[!flat_columns],reason="Non-scalar metadata; original object unchanged"),"reference_metadata_fields_not_exported.csv")
shared <- intersect(rownames(qc),rownames(rc))
if(length(shared)<5L) stop("Fewer than five matched panel genes.")
qsample <- sort(sample.int(ncol(qc),min(feature_cells,ncol(qc))))
rsample <- sort(sample.int(ncol(rc),min(feature_cells,ncol(rc))))
features <- cbind(qc[shared,qsample,drop=FALSE],rc[shared,rsample,drop=FALSE])
colnames(features) <- paste0("feature_cell_",seq_len(ncol(features)))
feature_study <- c(rep("xenium",length(qsample)),rep(ref,length(rsample)))
se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=features))
genes <- MetaNeighbor::variableGenes(dat=se,exp_labels=feature_study,min_recurrence=2,downsampling_size=0)
if(length(genes)<5L || anyDuplicated(genes) || any(!genes %in% shared)) stop("Fewer than five valid recurrent MetaNeighbor variable genes; review panel coverage. No gene-set fallback is applied.")
write(data.frame(gene=shared,selected_HVG=shared %in% genes),"comparison_genes.csv")
rm(se,features); gc()
dense_gib <- as.double(length(genes))*(ncol(qc)+ncol(rc))*8/1024^3
cat("Shared genes:",length(shared),"HVGs:",length(genes),"cells:",ncol(qc),"+",ncol(rc),"single dense matrix GiB:",round(dense_gib,2),"\n")
if(dense_gib > max_dense) stop("Estimated dense rank matrix exceeds --max-dense-gib. Review resources or a streaming design; no cells/genes were silently discarded.")
qc <- qc[genes,,drop=FALSE]; rc <- rc[genes,,drop=FALSE]
q_rank <- mn_rankable_cells(qc,"Xenium")
r_rank <- mn_rankable_cells(rc,ref)
rank_audit <- rbind(q_rank,r_rank)
rank_audit$selected_HVGs <- length(genes)
rank_audit$included_in_scoring <- rank_audit$rankable
write(rank_audit,"rank_constant_cell_audit.csv")
q_keep <- q_rank$rankable; r_keep <- r_rank$rankable
if (!any(q_keep) || !any(r_keep)) stop("No rankable cells remain in one study after selecting MetaNeighbor genes.")
if (any(!q_keep) || any(!r_keep)) message("Excluding ",sum(!q_keep)," Xenium and ",sum(!r_keep)," ",ref," rank-constant cells from MetaNeighbor scoring; see rank_constant_cell_audit.csv.")
qmd_full <- qmd; rmd_full <- rmd
qc <- qc[,q_keep,drop=FALSE]; qmd <- qmd[q_keep,,drop=FALSE]
rc <- rc[,r_keep,drop=FALSE]; rmd <- rmd[r_keep,,drop=FALSE]
if (!setequal(unique(qmd_full$corrected_broad_label),unique(qmd$corrected_broad_label)) ||
    !setequal(unique(qmd_full$comb_subcluster),unique(qmd$comb_subcluster)) ||
    !setequal(unique(rmd_full[[original_label]]),unique(rmd[[original_label]]))) stop("A MetaNeighbor label lost every cell after rank-constant exclusion; inspect rank_constant_cell_audit.csv.")
q_rank$comparison_cell_id <- NA_character_; q_rank$comparison_cell_id[q_keep] <- paste0("xenium_cell_",seq_len(sum(q_keep)))
r_rank$comparison_cell_id <- NA_character_; r_rank$comparison_cell_id[r_keep] <- paste0(ref,"_cell_",seq_len(sum(r_keep)))
write(data.frame(CellID=q_rank$source_cell_id,comparison_cell_id=q_rank$comparison_cell_id,
                 included_in_scoring=q_rank$rankable,qmd_full[q_rank$source_cell_id,c("orig.ident",levels_q),drop=FALSE],check.names=FALSE),"xenium_cell_labels.csv")
write(data.frame(source_cell_id=r_rank$source_cell_id,comparison_cell_id=r_rank$comparison_cell_id,
                 included_in_scoring=r_rank$rankable,rmd_full[r_rank$source_cell_id,flat_columns,drop=FALSE],check.names=FALSE),"reference_original_metadata.csv")
write(data.frame(study=c(rep("xenium",length(qsample)),rep(ref,length(rsample))),
                 CellID=c(q_rank$source_cell_id[qsample],r_rank$source_cell_id[rsample]),
                 included_in_scoring=c(q_rank$rankable[qsample],r_rank$rankable[rsample])),"feature_selection_cells.csv")
colnames(qc) <- paste0("xenium_cell_",seq_len(ncol(qc))); colnames(rc) <- paste0(ref,"_cell_",seq_len(ncol(rc)))
dat <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=cbind(qc,rc)))
study <- c(rep("xenium",ncol(qc)),rep(ref,ncol(rc)))
rl <- sort(unique(as.character(rmd[[original_label]])))
if(length(rl)<2L) stop("Reference has fewer than two original named populations.")
for(level in levels_q) {
  ql <- reference_similarity_order(unique(qmd[[level]]),level,legacy_query_labels=legacy_query_labels)
  if(length(ql)<2L) stop("Fewer than two query populations.")
  message("MetaNeighborUS fast cell-level comparison: ",ref," / ",level)
  result <- MetaNeighbor::MetaNeighborUS(var_genes=genes,dat=dat,study_id=study,
                cell_type=c(as.character(qmd[[level]]),as.character(rmd[[original_label]])),
                fast_version=TRUE,node_degree_normalization=TRUE,one_vs_best=FALSE,symmetric_output=FALSE)
  scores <- mn_pair_scores(result,ql,rl,ref)
  stem <- paste(ref,level,sep="_")
  write(scores$full,paste0(stem,"_all_AUROCs.csv")); write(scores$best,paste0(stem,"_best_matches.csv"))
  write(data.frame(query_label=ql,cells=as.integer(table(factor(qmd[[level]],levels=ql)))),paste0(stem,"_query_cell_counts.csv"))
  write(data.frame(reference_label=rl,cells=as.integer(table(factor(rmd[[original_label]],levels=rl)))),paste0(stem,"_reference_cell_counts.csv"))
  mn_heatmaps(scores,file.path(heatmap_dir,paste0(stem,"_AUROC_heatmaps.pdf")),paste("MetaNeighborUS",ref,level,if (legacy_query_labels) "[LEGACY QUERY LABELS]" else ""))
  rm(result,scores); gc()
}
current_info <- file.info(inputs)
if(any(current_info$size != input_info$size) || any(as.numeric(current_info$mtime) != as.numeric(input_info$mtime))) stop("An input changed during analysis; no completion marker will be written.")
write(data.frame(input=names(inputs),path=unname(inputs),size_bytes=input_info$size,mtime=as.character(input_info$mtime)),"input_provenance.csv")
write(panel,"original_panel_ids.csv")
write(data.frame(path=panel_file,md5=unname(tools::md5sum(panel_file))),"panel_provenance.csv")
writeLines(c(paste("Reference:",ref,"original label:",original_label),paste("Seed:",seed,"feature selection cap per study:",feature_cells),
  paste("Started:",started,"Finished:",Sys.time()),paste("Shared panel genes:",length(shared),"recurrent HVGs:",length(genes),"Science HDF5 block rows:",science_block_rows),
  "Official MetaNeighborUS fast_version=TRUE (published low-memory approximation), node_degree_normalization=TRUE, one_vs_best=FALSE, symmetric_output=FALSE.",
  "Raw-count within-cell ranks. Feature selection uses random cells only. Cells with constant values across selected HVGs cannot be ranked; they alone are excluded from AUROC scoring and recorded in rank_constant_cell_audit.csv.",
  if (legacy_query_labels) "LEGACY QUERY LABEL MODE: saved standalone Endothelial and Pericytes labels are preserved as-is for this comparison; results are not under the current Vascular = Endothelial + Pericytes convention." else "CURRENT QUERY LABEL MODE: the saved labels were validated against the approved fine-to-broad map.",
  "Cross-dataset directional AUROCs and their arithmetic mean; chance discrimination=0.5. Ties in best matches are retained and counted.",
  "No resampling confidence intervals, accuracy claims, found threshold, cell relabeling, integration or source writes.",
  "No donor-held-out validation: same references contributed to annotation, so this is transcriptional consistency evidence.",
  "Broad reference labels support broad correspondence even when query labels are finer. Within-lineage validation is not implemented in this stage.",
  capture.output(sessionInfo())),file.path(out,"methods_and_session.txt"))
snapshot_sources <- c("xenium_compare_metaneighbor.R","R/metaneighbor_original_io.R","R/metaneighbor_results.R","R/science_original_reference.R","R/reference_similarity_plots.R","R/corrected_broad_labels.R","color_palette.R")
for(f in snapshot_sources) if(!file.copy(here::here("scripts",f),file.path(out,paste0(gsub("/","_",f),".snapshot")))) stop("Could not snapshot analysis source: ", f)
file.copy(attr(config,"config_path"),file.path(out,"config_snapshot.yml"))
file.copy(inputs[["manifest"]],file.path(out,"sample_manifest_snapshot.csv"))
stems <- paste(ref,levels_q,sep="_")
expected_outputs <- c("xenium_cell_labels.csv","reference_original_metadata.csv","comparison_genes.csv",
                      "feature_selection_cells.csv","rank_constant_cell_audit.csv","comparison_mode.csv","input_provenance.csv","methods_and_session.txt",
                      "config_snapshot.yml","sample_manifest_snapshot.csv",
                      paste0(gsub("/","_",snapshot_sources),".snapshot"),
                      unlist(lapply(stems,function(stem) c(paste0(stem,c("_all_AUROCs.csv","_best_matches.csv","_query_cell_counts.csv","_reference_cell_counts.csv")), file.path("heatmaps",paste0(stem,"_AUROC_heatmaps.pdf"))))))
if(any(!file.exists(file.path(out,expected_outputs)))) stop("Required output set is incomplete.")
output_files <- list.files(out,full.names=TRUE,recursive=TRUE)
if(any(!is.finite(file.info(output_files)$size)) || any(file.info(output_files)$size <= 0)) stop("Missing or empty analysis outputs.")
write(data.frame(file=basename(output_files),size_bytes=file.info(output_files)$size),"output_manifest.csv")
writeLines("Complete: two query levels against one original reference; inspect both AUROC directions and metadata audits.",file.path(out,"COMPLETED.txt"))
message("Complete: ",out,"; elapsed minutes: ",round(as.numeric(difftime(Sys.time(),started,units="mins")),1))
