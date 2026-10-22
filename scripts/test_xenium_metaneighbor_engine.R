#!/usr/bin/env Rscript
# Tiny official-engine smoke test: ~10-60 seconds, low confidence; no real inputs.
args <- commandArgs(trailingOnly=TRUE)
if(any(!args %in% "--dry-run")) stop("Only --dry-run is supported.")
needed <- c("here","MetaNeighbor","Matrix","SummarizedExperiment")
ok <- setNames(vapply(needed,requireNamespace,TRUE,quietly=TRUE),needed)
print(data.frame(package=needed,available=ok),row.names=FALSE)
if(!all(ok)) stop("Missing smoke-test dependencies; no installation is performed.")
if("--dry-run" %in% args) quit(save="no")
here::i_am("scripts/test_xenium_metaneighbor_engine.R")
source(here::here("scripts","R","metaneighbor_results.R"))
set.seed(42)
group <- rep(1:3,each=12)
make_counts <- function() {
  m <- matrix(rpois(240*length(group),2),240)
  for(j in seq_along(group)) {
    ii <- seq.int((group[j]-1)*60+1,group[j]*60)
    m[ii,j] <- m[ii,j]+rpois(60,60)
  }
  rownames(m) <- paste0("G",seq_len(240));m
}
q <- make_counts();r <- make_counts();m <- cbind(q,r)
colnames(m) <- paste0("synthetic_cell_",seq_len(ncol(m)))
study <- rep(c("xenium","reference"),each=length(group))
dat <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=m))
genes <- MetaNeighbor::variableGenes(dat=dat,exp_labels=study,min_recurrence=2,downsampling_size=0)
if(length(genes)<5L) stop("Smoke test produced too few HVGs.")
for(query_names in list(c("VZ","Purkinje","Glia"),c("VZPs","Maturing PCs","BG"))) {
  reference_names <- c("original_VZ","original_PC","original_glia")
  a <- MetaNeighbor::MetaNeighborUS(var_genes=genes,dat=dat,study_id=study,
        cell_type=c(query_names[group],reference_names[group]),fast_version=TRUE,
        node_degree_normalization=TRUE,one_vs_best=FALSE,symmetric_output=FALSE)
  scores <- mn_pair_scores(a,query_names,reference_names,"reference")
  stopifnot(identical(scores$best$best_reference_label,reference_names),
            all(scores$best$reciprocal_best),all(scores$best$query_recovery_AUROC > .9))
}
cat("PASS: official MetaNeighborUS fast engine recovered known synthetic matches for both query label levels, with correct directional names and reciprocal hits.\n")
