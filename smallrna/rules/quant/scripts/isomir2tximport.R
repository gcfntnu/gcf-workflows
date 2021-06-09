
#!/usr/bin/env Rscript

library(argparse)
suppressPackageStartupMessages(library(edgeR))

parser <- ArgumentParser()

parser$add_argument("-i", "--input", required=TRUE,
                    help="tabulated isomirs (.tsv). Rownames in first column, colnames in header")
parser$add_argument("-f", "--feature-info", required=TRUE, dest="features",
                    help="isomir to mirna mapping(.tsv). Rownames in first column, colnames: isomir_id, mirna_id")
parser$add_argument("-o", "--output", default="tximport.rds",
                    help="Output filename")
parser$add_argument("-t", "--type", type="character", default="unitas",
                    help="Data origins (quickmirseq, unitas, mirge, none)")

parser$add_argument("--min-isomirs", default=1, dest="minc",
                    help="minimum number of isomiR sequences to be included in output counts")
parser$add_argument("--min-samples", default=1, dest="mins", 
                    help="INT: minimum number of samples with number of sequences bigger than ‘minc’ counts. FLOAT: minimum fraction ...")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")
args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose==TRUE){
    options(echo=TRUE)
    print(args)
}

count2TPM <- function(counts, prior.count=0){
    tpm <- function(counts, lengths=1) {
        rate <- (counts + prior.count) / lengths
        rate / sum(rate) * 1e6
    }

    TPM <- apply(counts, 2, function(x) tpm(x, 1))
    colnames(TPM) <- colnames(counts)
    rownames(TPM) <- rownames(counts)
    return(TPM)
}

if (args$type == "unitas"){
    counts <- read.delim(args$input, sep="\t", row.names=1, check.names=FALSE)
    anno <- read.delim(args$features, sep="\t", row.names=1, check.names=FALSE)
    anno <- anno[colnames(counts),]
}

countsMat <- t(as.matrix(round(counts)))
keep.ids <- rowSums(countsMat >= args$minc) >= args$mins
countsMat <- countsMat[keep.ids,]
lengthMat <- matrix(1, nrow=nrow(countsMat), ncol=ncol(countsMat))
abundanceMat <- edgeR::cpm(countsMat)

txi.tx <- list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
               countsFromAbundance="no")

txi.tx$type <- args$type
print(head(anno))
txi.tx$tx2gene <- anno[,c("isomir_id", "mirna_id")]

saveRDS(txi.tx, args$output)
