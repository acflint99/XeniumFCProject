# Clear the environment
rm(list = ls())

# load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(jsonlite)

panel <- fromJSON("/data/user/acflint/FCXeniumProject/FB124_X_G/gene_panel.json", simplifyVector = FALSE)

#extract genes
targets <- panel$payload$targets
length(targets)

# extract only entries where descriptor is "gene"
xenium_genes <- sapply(targets, function(x) {
  if (!is.null(x$type$descriptor) && x$type$descriptor == "gene") {
    x$type$data$name
  } else {
    NA
  }
})

# remove NAs
xenium_genes <- unique(na.omit(xenium_genes))

length(xenium_genes)  # should now be ~5000
head(xenium_genes)

saveRDS(xenium_genes, "/data/user/acflint/FC_published/Xenium/xenium_5k_genes.rds")
writeLines(xenium_genes, "/data/user/acflint/FC_published/Xenium/xenium_5k_genes.txt")

