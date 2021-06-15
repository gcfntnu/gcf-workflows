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
    make_option(c("-d", "--denoise"), action="store_true", help="denoise by PCA smoothing"),
    make_option(c("-t", "--threads"), type="integer", default="1", help="Number of threads [%default]"),
    make_option(c("-v", "--verbose"), action="store_true", help="print verbose output")
    )

parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list,
                       description="Variance transform of single cell rna data using sctransform (Seurat).")

args <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0, convert_hyphens_to_underscores=TRUE)$options


suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(sctransform))
suppressPackageStartupMessages(require(Matrix))

## 0 read data

if (args$input_format == "anndata"){
    obj <- Seurat::ReadH5AD(args$input, overwrite=FALSE)
}

##FIXME: somehow this does not work, SCTransform is using all threads
future::plan(strategy = 'multicore', workers = args$threads)

obj <- SCTransform(object=obj, verbose=FALSE)

##if (args$denoise == TRUE){
    ## FIXME this needs to be implemented in SCTransform
 ##   stop("NotImplementedError")
  ##  umi <- GetAssayData(object = obj, slot = "counts")
   ## umi_smooth <- smooth_via_pca(umi)
   ## umi_corrected <- correct(vst_out, data = y_smooth)
   ## }

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
