#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(Seurat))
library(sctransform)

if (!require(Seurat)) {
    stop("This script needs to have Seurat installed.")
}
if (!require(optparse)) {
    stop("This script needs to have optparse installed.")
}


##parser <- ArgumentParser()
##parser$add_argument("-i", "input", default="data/tmp/test_sample/outs/filtered_gene_bc_matrices", nargs="+",
##                    help="Directory containing the matrix.mtx.gz, genes.tsv.gz, and barcodes.tsv.gz files provided by 10X. Multiple dirs allowed.")
##parser$add_argument("-s", "--sample-info", dest="samples", default=NULL,
##                    help="Sample info. Tab delimited file, needs column with `sample_id`")
##parser$add_argument("-f", "--format", choices=['seurat', 'loom'], default='loom', help='output file format')
##parser$add_argument("-o", "--output", default="data.loom",
##                    help="Output filename (loom/rds file)")
##parser$add_argument("-v", "--verbose", action="store_true",
##                    default=FALSE, help="Print extra output")
##
##args <- parser$parse_args(args=commandArgs(TRUE))

option_list <- list(
    make_option(c("-i", "--input"),  action="append", type="character", help="input file(s)"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="verbose output [default]"),
    make_option(c("-s", "--sample-info"), default=NULL,
        help="sample info. Tab delimited file, needs column with `sample_id`[default %default]",
        metavar="samples"),
    make_option(c("-f", "--format"), default="loom", 
        help = "output format [default \"%default\"]"),
    make_option(c("-o", "--output"), default="data.loom", 
        help = "output filename [default \"%default\"]")
    )

args <- parse_args(OptionParser(option_list=option_list))

if (args$verbose == TRUE) options(echo=TRUE)

if (length(args$input)> 1){
    input.names <- NULL
    for (dir in args$input){
        parent.dir <- dirname(dir)
        if (basename(parent.dir) == "outs"){
            parent.dir <- dirname(parent.dir)
        }
        input.names <- c(input.names, basename(parent.dir))
    }
    names(args$input) <- input.names
}

data <- Read10X(args$input)
if (!is.null(args$samples)){
    meta.info <- read.delim(args$samples, sep="\t", row.names=1)
} else{
    meta.info <- NULL
}

seurat.obj  <- CreateSeuratObject(data, meta.data = meta.info, assay="RNA", names.delim="_")
data <- SubsetData(data, subset = nFeature_RNA > 300 & nFeature_RNA < 30000)
seurat.obj <- SCTransform(seurat.obj, verbose=args$verbose)

if (args$format == "seurat"){
    saveRDS(seurat.obj, args$output)
}


if (args$format == "loom"){
    library(loomR)
    seurat.obj.loom <- as.loom(seurat.obj, filename = args$output, verbose = FALSE)
    seurat.obj.loom$close_all()
}


