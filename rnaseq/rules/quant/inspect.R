#/usr/bin/env Rscript

library(argparse)
library(INSPEcT)
library(BiocParallel)
library(AnnotationDbi)

parser <- ArgumentParser(description="RNA kinetics quantification by INSPEcT-")

parser$add_argument("input", nargs="+", help="bam files")
parser$add_argument("--txdb", required=TRUE,
                    help="Gene info. TxDB sqllite file")
parser$add_argument("--sample_info", required=TRUE,
                    help="Sample info. Tab delimited file, needs column named `Sample_ID`")
parser$add_argument("--sample_group", type="character", default="Sample_Group",
                    help="Column name in sample_info defining subgroups [optional]")

parser$add_argument("--timepoint", type="character", default="Timepoint", required=TRUE,
                    help="Column name in sample_info defining numerical sample timepoints")
parser$add_argument("--timepoint-mutiplier", default=1,
                    help="timepoint mutiplier. ex: hrs-> min , mutiplier=60")
parser$add_argument("--output", default="inspect.",
                    help="Output rds file")
parser$add_argument("--threads", default=1, type="int",
                    help="Number of threads")
parser$add_argument("--verbose", action="store_true", default=FALSE,
                    help="Print extra output")
parser$add_argument("--tzero", action="store_true", default=FALSE,
                    help="data contains common T0 reference points")
parser$add_argument("--paired_end", action="store_true", default=FALSE,
                    help="bam files from paired end fastq")
parser$add_argument("--read_orientation", type="character", default="unstranded",
                    help="read strandedness type wrt R1 [forward/reverse/unstranded]")
parser$add_argument("--organism", type="character", default="homo_sapiens",
                    help="organism")

args <- parser$parse_args(args=commandArgs(TRUE))

args$input <- normalizePath(args$input)

if (args$threads > 1){
    BPPARAM <- BiocParallel::MulticoreParam(workers=args$threads)
} else{
    BPPARAM <- BiocParallel::SerialParam()
}
args$threads <- BPPARAM

if (args$read_orientation == "reverse"){
    s = 2
} else if(args$read_orientation == "forward"){
    s = 1
} else{
    s = 0
}
args$read_orientation <- s

#if (args$organism == "homo_sapiens"){
#    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#    txdb <-  TxDb.Hsapiens.UCSC.hg38.knownGene
#} else if (args$organism == "mus_musculus"){
#    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#} else{
#    stop("only human/mouse is supported")
#}
#args$organism <- txdb

txdb <- AnnotationDbi::loadDb(args$txdb)


df <- read.delim(args$sample_info, sep="\t", row.names="Sample_ID")
sample.ids <- sapply(str_split(basename(args$input), "\\.sorted\\.bam"), function(x) x[[1]])
df <- df[sample.ids,]
args$samples <- sample.ids

if (!args$timepoint %in% colnames(df)){
    stop("arg `timpoints` not a column name of sample_info")
} else{
    args$timepoint <- as.numeric(df[,args$timepoint]) * args$timepoint_multiplier 
}

if (!is.null(args$sample_group)){
    if (!args$sample_group %in% colnames(df)){
        stop("arg `sample_group` not a column name of sample_info")
    } else{
        args$sample_group <- as.factor(df[,args$sample_group])
    }
}


inspectFromBAM2 <- function (args){
    res <- list()
    if (is.null(args$sample_group)){
        E <- quantifyExpressionsFromBAMs(txdb=txdb,
                                         BAMfiles=args$input,
                                         experimentalDesign=args$timepoint,
                                         strandSpecific=args$read_orientation,
                                         isPairedEnd=args$paired_end,
                                         BPPARAM=args$thread)
        tpts <- sort(unique(args$timepoint))
        obj <- newINSPEcT(tpts, matureExpressions=E, BPPARAM=args$thread)
        obj <- modelRates(obj, estimateRatesWith="der", useSigmoidFun=TRUE, BPPARAM=args$thread)
        res[["all_samples"]] <- obj
        saveRDS(obj, paste0(args$output, "all_samples_obj.rds"))
        saveRDS(E, paste0(args$output, "all_samples_E.rds"))
    } else{
        t0.index <- NULL
        if (args$tzero == TRUE){
            t0.index <- which(args$timepoint == 0)
            print(t0.index)
        }
        for (name in levels(args$sample_group)){
            print(name)
            keep.index <- which(args$sample_group == name)
            keep.index <- c(t0.index, keep.index)
            E <- quantifyExpressionsFromBAMs(txdb=txdb,
                                         BAMfiles=args$input[keep.index],
                                         experimentalDesign=args$timepoint[keep.index],
                                         strandSpecific=args$read_orientation,
                                         isPairedEnd=args$paired_end,
                                         BPPARAM=args$thread)
            tpts <- sort(unique(args$timepoint[keep.index]))
            obj <- newINSPEcT(tpts, matureExpressions=E)
            obj <- modelRates(obj, estimateRatesWith="der", useSigmoidFun=TRUE, BPPARAM=BPPARAM)
            saveRDS(obj, paste0(args$output, name, "_obj.rds"))
            saveRDS(E, paste0(args$output, name, "_E.rds"))
            res[[name]] <- obj
            print(paste(name, "done ..."))
        }
    }
    return(res)
}
