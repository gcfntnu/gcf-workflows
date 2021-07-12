
suppressPackageStartupMessages({
    library(argparse)
    library(tximeta)
})


parser <- ArgumentParser(description="tximeta db for Salmon index")

parser$add_argument("-i", "--index", required=TRUE, help="Salmon index")

parser$add_argument("-t", "--transcriptome", required=TRUE, help="Transcriptome fasta file")

parser$add_argument("-m", "--organism", type="character", default="homo sapiens", help="Organsim")

parser$add_argument("-g", "--gtf", type="character", required=TRUE, help="gtf filename")

parser$add_argument("-s", "--source", type="character", help="reference database", default="ensembl")

parser$add_argument("-r", "--release", type="character", help="reference database release", default="100")

parser$add_argument("-a", "--assembly", type="character", help="reference database assembly", default="GRCh38")

parser$add_argument("-o", "--output", default="tximeta.json", help="Output json file")

parser$add_argument("-c", "--cachedir", default="/tmp/.cache", help="Cachedir location")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


if (args$verbose == TRUE){
    options(echo=TRUE)
}

if (!dir.exists(args$cachedir)){
    dir.create(args$cachedir, showWarnings=FALSE, recursive=TRUE)
    }

setTximetaBFC(args$cachedir, quiet=TRUE)

makeLinkedTxome(indexDir = dirname(args$index),
                source = args$source,
                organism = args$organism,
                release = args$release,
                genome = args$assembly,
                fasta = args$transcriptome,
                gtf = args$gtf,
                write = TRUE,
                jsonFile = args$output)
