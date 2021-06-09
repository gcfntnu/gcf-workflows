#/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(PCAtools))
library(cowplot)
library(ggplotify)




parser <- ArgumentParser(description="tximport obj csv export")

parser$add_argument("input", help="tximport obj RDS file")

parser$add_argument("--txinfo", required=FALSE,
                    help="Transcript info (optional). Tab delimited file, needs columns with `gene_id`, `gene_name`")
parser$add_argument("--sample-info", required=FALSE,
                    help="Sample info (optional). Tab delimited file, needs columns with `sample_id`, `Sample_Group`")
parser$add_argument("-o", "--output", required=TRUE, help="Output tsv file")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose==TRUE){
    print(args)
    options(echo=TRUE)
}

exprs <- read.table(args$input, sep="\t", header=TRUE)
n.lim <- max(1, round(ncol(exprs) * 0.1))
keep.genes <- rowSums(exprs > 1) > n.lim
exprs <- exprs[keep.genes,]
                                     
if (!is.null(args$txinfo)){
    gene.info <- read.table(args$txinfo, sep="\t", header=TRUE)
    gene.info <- gene.info[rownames(exprs),]
}

if (!is.null(args[,"sample_info"])){
    sample.info <- read.table(args$txinfo, sep="\t", header=TRUE)
    sample.info <- samle.info[colnames(exprs),]
}

p <- PCAtools::pca(t(exprs), metadata=sample.info)
elbow <- findElbowPoint(p$variance)
g <- pairsplot(p,
               components = getComponents(p, c(1,2,3,4)),
               triangle = FALSE,
               hline = 0, vline = 0,
               pointSize = 0.8,
               gridlines.major = FALSE, gridlines.minor = FALSE,
               colby = 'Sample_Group',
               title = 'Pairs plot', plotaxes = TRUE,
               margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))
