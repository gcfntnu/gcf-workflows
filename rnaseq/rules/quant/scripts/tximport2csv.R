#/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(DESeq2))


parser <- ArgumentParser(description="tximport obj csv export")

parser$add_argument("input", help="tximport obj RDS file")

parser$add_argument("--txinfo", required=FALSE,
                    help="Transcript info. Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("--geneinfo", required=FALSE,
                    help="Gene info (optional). Tab delimited file, needs column with `gene_id`")
parser$add_argument("-t", "--type", type="character", default="gene",
                    help="csv count output option (gene, gene_tpm, gene_tpm_scaled, gene_tpm_length_scaled, gene_rlog, gene_vst, tx, tx_tpm, tx_tpm_scaled, tx_rlog, tx_vst, gene_length, variances, gene_info, tx_info)")

parser$add_argument("-o", "--output", required=TRUE, help="Output tsv file")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose==TRUE){
    print(args)
    options(echo=TRUE)
}
txi.tx <- readRDS(args$input)

if (is.null(args$txinfo)){
    tx.info <- txi.tx$tx2gene
} else{
    tx.info <- read.table(args$txinfo, sep="\t", header=TRUE)
}


#small testdata may contain a lot of zero counts 
safeEstimateSizeFactors <- function(dds){
    cts <- counts(dds)
    geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
    out <- estimateSizeFactors(dds, geoMeans=geoMeans)
}


if (args$type == "gene"){
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene)$counts
} else if (args$type == "gene_tpm"){
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene)$abundance
} else if (args$type == "gene_tpm_scaled"){
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene, countsFromAbundance="scaledTPM")$counts
} else if (args$type == "gene_tpm_length_scaled"){
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene, countsFromAbundance="lengthScaledTPM")$counts
} else if (args$type == "gene_length"){
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene)$length
} else if (args$type == "tx"){
    out <- txi.tx$counts
} else if (args$type == "tx_tpm"){
    out <- txi.tx$abundance
} else if (args$type == "tx_tpm_scaled"){
    length4CFA <- tximport:::medianLengthOverIsoform(txi.tx$length, txi.tx$tx2gene)
    out <- tximport:::makeCountsFromAbundance(countsMat = txi.tx$counts, 
                                              abundanceMat = txi.tx$abundance,
                                              lengthMat=length4CFA,
                                              countsFromAbundance = "lengthScaledTPM") 
} else if (args$type == "tx_length"){
    out <- txi.tx$length
} else if (args$type == "variances"){
    if (! args$type %in% c("salmon", "kallisto")){
        stop("only salmon, kallisto estimate variances!")
        }
    out <- summarizeToGene(txi.tx, txi.tx$tx2gene, varReduce=TRUE)$variance
}  else if (args$type == "gene_vst"){
    counts <- summarizeToGene(txi.tx, txi.tx$tx2gene)$counts
    samples <- data.frame(row.names=colnames(counts))
    dds <-DESeqDataSetFromMatrix(round(counts), colData=samples, design = ~1)
    dds <- safeEstimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds)
    out <- assay(vsd)
} else if (args$type == "gene_rlog"){
    counts <- summarizeToGene(txi.tx, txi.tx$tx2gene)$counts
    samples <- data.frame(row.names=colnames(counts))
    dds <-DESeqDataSetFromMatrix(round(counts), colData=samples, design = ~1)
    dds <- safeEstimateSizeFactors(dds)
    R <- rlog(dds)
    out <- assay(R)
} else if (args$type == "tx_vst"){
    samples <- data.frame(row.names=colnames(txi.tx$counts))
    dds <-DESeqDataSetFromMatrix(round(txi.tx$counts),colData=samples, design = ~1)
    dds <- safeEstimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds)
    out <- assay(vsd)
} else if (args$type == "tx_rlog"){
    samples <- data.frame(row.names=colnames(txi.tx$counts))
    dds <-DESeqDataSetFromMatrix(round(txi.tx$counts), colData=samples, design = ~1)
    dds <- safeEstimateSizeFactors(dds)
    R <- rlog(dds)
    out <- assay(R)
} else if (args$type == "gene_info"){
    gene.info <- read.table(args$geneinfo, sep="\t", header=TRUE)
    rownames(gene.info) <- gene.info[,"gene_id"]
    counts <- summarizeToGene(txi.tx, txi.tx$tx2gene)$counts
    out <- gene.info[rownames(counts),]
    keep.cols <- colSums(!is.na(out)) > 0
    out <- out[,keep.cols]
} else if (args$type == "tx_info"){
    if ('terminus_id' %in% colnames(tx.info)){
        tx.info <- tx.info[!duplicated(tx.info[,"terminus_id"]),]
        rownames(tx.info) <- tx.info[,"terminus_id"]
    } else{
        rownames(tx.info) <- tx.info[,"transcript_id"]
    }
    out <- tx.info[rownames(txi.tx$counts),]
    keep.cols <- colSums(!is.na(out)) > 0
    out <- out[,keep.cols]
}



write.table(as.data.frame(out), file=args$output, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)






