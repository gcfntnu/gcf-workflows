#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DESeq2))

parser <- ArgumentParser(description="Transform miRNA counts with DESeq2 (features by samples)")
parser$add_argument("--input", required=TRUE, help="Count TSV with feature identifiers in the first column")
parser$add_argument("--output", required=TRUE, help="Transformed TSV retaining feature and sample order")
args <- parser$parse_args()

counts <- as.matrix(read.delim(args$input, row.names=1, check.names=FALSE))
if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("Counts must be finite non-negative numbers")
}
if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Counts must contain features and samples")
}

# Preserve the historical pseudocount and fit; DESeq2 expects features by samples.
counts <- round(counts + 1)
if (ncol(counts) > 50) {
    transformed <- varianceStabilizingTransformation(counts, fitType="mean")
} else {
    transformed <- rlogTransformation(counts, fitType="mean")
}
write.table(transformed, file=args$output, sep="\t", quote=FALSE, col.names=NA)
