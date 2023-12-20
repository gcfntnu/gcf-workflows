#!/usr/bin/env Rscript
suppressPackageStartupMessages({
	library(argparse)
	library(TENxIO)
	library(scds)
	}
	)

parser <- ArgumentParser() 
parser$add_argument("-o", "--output", required = TRUE, help="Output file")
parser$add_argument("-i", "--input", required = TRUE, type = "character", help = "Path to the 10x filtered matrix.mtx")
args <- parser$parse_args()
barcodes.fn <- file.path(dirname(args$input), "barcodes.tsv.gz")
stopifnot(file.exists(barcodes.fn))
sce <- import(TENxMTX(args$input))
obs.names <- import(TENxTSV(barcodes.fn))
colnames(sce) <- obs.names$barcode
sce <- as(sce, "SingleCellExperiment")

sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)

results <- data.frame("Barcode"=rownames(colData(sce)),
		      "doublet" = ifelse(colData(sce)$hybrid_call, "doublet", "singlet"),
		      "doublet_score"=colData(sce)$hybrid_score
		      )
write.table(results, file=args$output, sep="\t", row.names=FALSE, quote=FALSE)


