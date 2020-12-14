#/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(DESeq2))


parser <- ArgumentParser(description="tabulated counts csv export")

parser$add_argument("input", help="tab separated counts file with colnames/rownames")

parser$add_argument("--txinfo", required=FALSE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("--geneinfo", required=FALSE,
                    help="Gene info (optional). Tab delimited file, needs column with `gene_id`")
parser$add_argument("-t", "--type", type="character", default="gene",
                    help="csv count output option (gene, gene_tpm, gene_tpm, gene_rlog, gene_vst, tx, tx_tpm, tx_tpm_scaled, tx_rlog, tx_vst, gene_length, gene_info, tx_info)")

parser$add_argument("-o", "--output", required=TRUE, help="Output tsv file")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose==TRUE){
    print(args)
    options(echo=TRUE)
}


