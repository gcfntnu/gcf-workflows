#!/usr/bin/env Rscript

if (!require(optparse)) {
    stop("This script needs to have optparse installed.")
}


parser <- OptionParser()
option_list <- list(
    make_option(c("-i", "--input"), type="character", help="input file"),
    make_option(c("-f", "--input-format"), type="character", default="anndata", help="input file format [%default]"),
    make_option(c("-o", "--output"), type="character", help="output directory"),
    make_option(c("-F", "--output-format"), type="character", default="mtx", help="output file format [%default]"),
    make_option(c("-m", "--model"), type="character", default="binomial", help="model (binomial, multinomal, poisson, geometric)[%default]"),
    make_option(c("-t", "--type"), type="character", default="deviance", help="null residual (deviance, pearson), [%default]"),
    make_option(c("-v", "--verbose"), action="store_true", help="print verbose output")
    )

parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list,
                       description="Variance transform of single cell RNA data using null residual approximation to glm-pca.")

args <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0, convert_hyphens_to_underscores=TRUE)$options


suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(Seurat)
suppressPackageStartupMessages(require(Matrix))
source("scrna2019/util/functions.R")


if (args$input_format == "anndata"){
    obj <- Seurat::ReadH5AD(args$input, overwrite=FALSE)
}


x <- GetAssayData(object = obj, slot = "counts")
U <- null_residuals(x, mod=args$model, type=args$type)
SetAssayData(object=obj, slot="counts", new.data=U)

## this is from aron lun's dropletutils
write_mtx <- function(path, x, barcodes, gene.id, gene.symbol, gene.type="Gene Expression", version="3") {
    require(Matrix)
    require(R.utils)
    dir.create(path, showWarnings=FALSE)
    gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

    if (version=="3") {
        gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
        fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
        on.exit({
            close(bhandle)
            close(fhandle)
        })
    } else {
        mhandle <- file.path(path, "matrix.mtx")
        bhandle <- file.path(path, "barcodes.tsv")
        fhandle <- file.path(path, "genes.tsv")
    }

    writeMM(x, file=mhandle)
    write(barcodes, file=bhandle)
    write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    if (version=="3") {
        # Annoyingly, writeMM doesn't take connection objects.
        R.utils::gzip(mhandle)
    }

    return(NULL)
}



if (args$output_format == "mtx"){
    x <- GetAssayData(object = obj, slot = "counts")
    x <- as(object = x, Class = 'dgCMatrix')
    barcodes = colnames(x)
    gene.ids <- rownames(x)
    write_mtx(args$output, x, barcodes, gene.ids, gene.symbol=gene.ids)
}

##if (args$output_format == "loom"){
    require(loomR)
    lfile <- as.loom(x = pbmc_small, filename=file.path(args$output, "adata.loom"))
    lfile$close_all()
##}
