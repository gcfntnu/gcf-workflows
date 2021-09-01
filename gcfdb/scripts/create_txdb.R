#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(AnnotationDbi)
library(GenomicFeatures)
library(stringr)


GTF = args[[1]]
## is symlink, get symlinked name
if (nzchar(Sys.readlink(GTF))){
    gtf.abspath <- normalizePath(GTF)
    GTF <- file.path(dirname(gtf.abspath), Sys.readlink(GTF))
}
ORG <- str_replace(str_to_title(args[[2]]), "_", " ")
DB <- paste("GCF", args[[3]], sep=".")
output.file <- args[[4]]

txdb <- GenomicFeatures::makeTxDbFromGFF(GTF, format="gtf", dataSource=DB, organism=ORG)
AnnotationDbi::saveDb(txdb, output.file)


