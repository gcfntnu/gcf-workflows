#/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(DESeq2))


parser <- ArgumentParser(description="tximport obj csv export")

parser$add_argument("input", help="tximport obj RDS file")

parser$add_argument("--txinfo", required=TRUE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

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

tx.info <- read.table(args$txinfo, sep="\t", header=TRUE)
tx2gene <- tx.info[,c("transcript_id", "gene_id")]

txi.tx <- readRDS(args$input)

if (args$type == "gene"){
    out <- summarizeToGene(txi.tx, tx2gene)$counts
} else if (args$type == "gene_tpm"){
    out <- summarizeToGene(txi.tx, tx2gene)$abundance
} else if (args$type == "gene_tpm_scaled"){
    out <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance="scaledTPM")$counts
} else if (args$type == "gene_tpm_length_scaled"){
    out <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance="lengthScaledTPM")$counts
} else if (args$type == "gene_length"){
    out <- summarizeToGene(txi.tx, tx2gene)$length
} else if (args$type == "tx"){
    out <- txi.tx$counts
} else if (args$type == "tx_tpm"){
    out <- txi.tx$abundance
} else if (args$type == "tx_tpm_scaled"){
    length4CFA <- tximport:::medianLengthOverIsoform(txi.tx$length, tx2gene)
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
    out <- summarizeToGene(txi.tx, tx2gene, varReduce=TRUE)$variance
}  else if (args$type == "gene_vst"){
    counts <- summarizeToGene(txi.tx, tx2gene)$counts
    samples <- data.frame(row.names=colnames(counts))
    dds <-DESeqDataSetFromMatrix(round(counts), colData=samples, design = ~1)
    dds <- estimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds)
    out <- assay(vsd)
} else if (args$type == "gene_rlog"){
    counts <- summarizeToGene(txi.tx, tx2gene)$counts
    samples <- data.frame(row.names=colnames(counts))
    dds <-DESeqDataSetFromMatrix(round(counts), colData=samples, design = ~1)
    dds <- estimateSizeFactors(dds)
    R <- rlog(dds)
    out <- assay(R)
} else if (args$type == "tx_vst"){
    samples <- data.frame(row.names=colnames(txi.tx$counts))
    dds <-DESeqDataSetFromMatrix(round(txi.tx$counts),colData=samples, design = ~1)
    dds <- estimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds)
    out <- assay(vsd)
} else if (args$type == "tx_rlog"){
    samples <- data.frame(row.names=colnames(txi.tx$counts))
    dds <-DESeqDataSetFromMatrix(round(txi.tx$counts), colData=samples, design = ~1)
    dds <- estimateSizeFactors(dds)
    R <- rlog(dds)
    out <- assay(R)
} else if (args$type == "gene_info"){
    gene.info <- read.table(args$geneinfo, sep="\t", header=TRUE, row.names=1)
    counts <- summarizeToGene(txi.tx, tx2gene)$counts
    out <- gene.info[rownames(counts),]
    keep.cols <- colSums(!is.na(out)) > 0
    out <- out[,keep.cols]
} else if (args$type == "tx_info"){
    print(head(tx.info))
    rownames(tx.info) <- tx.info[,"transcript_id"]
    print(head(tx.info))
    print(head(txi.tx$counts))
    out <- tx.info[rownames(txi.tx$counts),]
    keep.cols <- colSums(!is.na(out)) > 0
    out <- out[,keep.cols]
}



write.table(as.data.frame(out), file=args$output, sep="\t", quote=FALSE)






