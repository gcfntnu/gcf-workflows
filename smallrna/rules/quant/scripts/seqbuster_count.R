#!/usr/bin/env Rscript

library(argparse)
suppressPackageStartupMessages(library(isomiRs))

parser <- ArgumentParser()

parser$add_argument("input", nargs='+', help="Gene alignment seqbuster files (.mirna)")

parser$add_argument("-o", "--output-dir", default="data/test", dest="outdir",
                    help="Output prefix")

parser$add_argument("--design", nargs="?", dest="design", default="~1",
                    help="Design formula. Used if variance transformed data is asked for.")

parser$add_argument("--gene-info", nargs="?", dest="features", 
                    help="Gene info (optional). Tab delimited file, needs column with `gene_id`")

parser$add_argument("--sample-info", required=TRUE, dest="samples",
                    help="Sample info (required). Tab delimited file, needs column with `sample_id`.")

parser$add_argument("--min-isomirs", default=1, dest="minc",
                    help="minimum number of isomiR sequences to be included in output counts")

parser$add_argument("--min-samples", default=1, dest="mins", 
                    help="INT: minimum number of samples with number of sequences bigger than ‘minc’ counts. FLOAT: minimum fraction ...")

parser$add_argument("--output-expression", action="store_true",
                    default=FALSE,
                    help="Additional output of variance transformed data, (DESeq2: vst)")

parser$add_argument("--output-psi", action="store_true",
                    default=FALSE,
                    help="Additional output of percent spliced in fractions.")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
args$design <- as.formula(args$design)

if (args$verbose==TRUE){
    options(echo=TRUE)
    print(args)
}

##fn_list = list.files("data/tmp/smallrna/quant/seqbuster/", pattern="*.mirna$", full.names=TRUE)

df <- read.delim(args$samples, sep="\t", row.names=1)
sample.ids <- tools::file_path_sans_ext(basename(args$input))
df <- df[sample.ids,]

data <- IsomirDataSeqFromFiles(args$input, coldata=df)

mir.data <- isoCounts(data, all=TRUE, minc=args$minc, mins=args$minc)
mir.counts <- counts(mir.data)
mir.anno <- cbind(mirna_id=rownames(mir.counts), row.names=rownames(mir.counts))
mir.counts.fn <- file.path(args$outdir, "mir_counts.txt")
write.table(data.frame(mirna_id=rownames(mir.counts), mir.counts), file=mir.counts.fn, sep="\t", row.names=FALSE, quote=FALSE)

isomir.data <- isoCounts(data, minc=args$minc, mins=args$minc, ref=TRUE, iso5=TRUE, iso3=TRUE, add=TRUE, snv=TRUE, seed=FALSE, all=TRUE)
isomir.counts <- counts(isomir.data)
mirna_id <- sapply(strsplit(rownames(isomir.counts), "\\;"), function(x) x[[1]])
isomir.anno <- data.frame(isomirna_id=rownames(isomir.counts), mirna_id=mirna_id, row.names=rownames(isomir.counts))
isomir.counts.fn <- file.path(args$outdir, "isomir_counts.txt")
write.table(data.frame(mirna_id=rownames(isomir.counts), isomir.counts), file=isomir.counts.fn, sep="\t", row.names=FALSE, quote=FALSE)


isomir.a <- isoAnnotate(mir.data)
rownames(isomir.a) <- isomir.a[,"uid"]
anno.index <- c(1, 2, ncol(isomir.a))
isomir.psi <- isomir.a[,-anno.index]
isomir.psi.anno <- isomir.a[,anno.index]
isomir.psi <- isomir.psi[rownames(isomir.counts),]
isomir.psi[is.na(isomir.psi)] <- 0
isomir.psi.anno <- isomir.psi.anno[rownames(isomir.psi),]
isomir.anno <- cbind(isomir.anno, isomir.psi.anno)
isomir.anno <- isomir.anno[,-which(colnames(isomir.anno)=="uid")]

if (args$output_psi == TRUE){
    isomir.psi.fn <- file.path(args$outdir, "isomir_psi.txt")
    write.table(data.frame(isomirna_id=rownames(isomir.psi), isomir.psi), file=isomir.psi.fn, sep="\t", row.names=FALSE, quote=FALSE)
}

if (args$output_expression == TRUE){
    require(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData=mir.counts, colData = colData(mir.data), design=args$design)
    vst <- varianceStabilizingTransformation(dds, blind=FALSE)
    E.mir <- assay(vst)
    mir.exprs.fn <- file.path(args$outdir, "mir_vst.txt")
    E.tab <- cbind(mirna_id=rownames(E.mir), E.mir)
    write.table(E.tab, file=mir.exprs.fn, sep="\t", row.names=FALSE, quote=FALSE)

    dds.iso <- DESeqDataSetFromMatrix(countData=isomir.counts, colData = colData(isomir.data), design=args$design)
    vst.iso <- varianceStabilizingTransformation(dds.iso, blind=FALSE)
    E.isomir <- assay(vst.iso)
    isomir.exprs.fn <- file.path(args$outdir, "isomir_vst.txt")
    E.tab.iso <- cbind(isomirna_id=rownames(E.isomir), E.isomir)
    write.table(E.tab.iso, file=isomir.exprs.fn, sep="\t", row.names=FALSE, quote=FALSE)
}


mir.anno.fn <- file.path(args$outdir, "mir_anno.txt")
isomir.anno.fn <- file.path(args$outdir, "isomir_anno.txt")

if (!is.null(args$features)){
    df <- read.delim(args$features, sep="\t", row.names=1)
    mir.anno <- cbind(mir.anno, df[rownames(mir.anno),])
    isomir.anno <- cbind(isomir.anno, df[as.character(isomir.anno$mirna_id),])
    }
write.table(mir.anno, file=mir.anno.fn, sep="\t", row.names=FALSE, quote=FALSE)
write.table(isomir.anno, file=isomir.anno.fn, sep="\t", row.names=FALSE, quote=FALSE)
