#!/usr/bin/env Rscript
suppressPackageStartupMessages({
	library(argparse)
	library(TENxIO)
	library(scDblFinder)}
	)

parser <- ArgumentParser() 
parser$add_argument("-o", "--output", required = TRUE, help="Output file")
parser$add_argument("-i", "--input", required = TRUE, type = "character", help = "Path to the 10x filtered matrix.mtx")
args <- parser$parse_args()

con <- TENxMTX(args$input)
sce  <- import(con)
barcodes.fn <- file.path(dirname(args$input), "barcodes.tsv.gz")
stopifnot(file.exists(barcodes.fn))
obs.names <- import(TENxTSV(barcodes.fn))
colnames(sce) <- obs.names$barcode

sce <- as(sce, "SingleCellExperiment")
doublet_ratio <- ncol(sce)/1000*0.008
sce <- scDblFinder(sce, dbr=doublet_ratio)

barcodes <- rownames(colData(sce))
barcodes.numerical <- paste0(sapply(strsplit(barcodes, split="-"), function(c) c[[1]]), "-1")
results <- data.frame("Barcode"=barcodes.numerical,
                      "doublet"=sce$scDblFinder.class,
		      "doublet_score"=sce$scDblFinder.score)

write.table(results, file=args$output, sep="\t", row.names=FALSE, quote=FALSE)
