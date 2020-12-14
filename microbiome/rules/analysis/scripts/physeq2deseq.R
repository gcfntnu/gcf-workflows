#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(phyloseq)

physeq.fn <- args[1]
sample.info.fn <- args[2]
out.fn <- args[3]

physeq2deseq <- function(physeq){
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(diagdds), 1, gm_mean)
    diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
    diagdds = DESeq(diagdds, fitType="local")
}

physeq <- readRDS(physeq.fn)
sample.info <- read.table(sample.info.fn, sep="\t", stringsAsFactors=FALSE)
rownames(sample.info) <- sample.info[,"Sample_ID"]
physeq.names <- sample_names(physeq)
sample.info <- sample.info[physeq.names,]

saveRDS(physeq, out.fn)
    
