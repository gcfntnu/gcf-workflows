#/usr/bin/env Rscript
##Filters for RNA-Seq gene count tables


if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}



parser <- ArgumentParser()

parser$add_argument("-i", "--input", help="Gene quantification file (.quant)", default="a")

parser$add_argument("-o", "--output", default="./featurecounts", dest="prefix",
                    help="Output prefix")

parser$add_argument("--lef", type="character", default="k_over_limit",
                    help="Low Expression Filter [default %(default)s]")

parser$add_argument("--lef-limit", default=0.5,
                    help="detection limit (in cpm) [default %(default)s]")

parser$add_argument("--lvf", type="character", default="none",
                    help="Low Variance Filter [default %(default)s]")

parser$add_argument("--gene-info", default="gene_info.txt",
                    help="Gene info. Tab delimited file, needs column with `gene_id`")


parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


if (!require(edgeR)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    library("edgeR")
}

if (!require(genefilter)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("genefilter")
    library("genefilter")
}


gene.info <- fread(gene.info)
gene.info <- subset(gene.info, gene_id != "")
if (sum(duplicated(gene_id)) > 0){
    stop('duplicated gene ids found!')
}



gene.info <- data.frame(gene.info)
