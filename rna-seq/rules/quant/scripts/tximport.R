#/usr/bin/env Rscript

library(argparse)
library(tximport)
library(readr)


parser <- ArgumentParser(description="Transcript quantification from Kallisto/Salmon/Sailfish/RSEM")

parser$add_argument("input", nargs="+", help="Kallisto/Salmon/Sailfish/RSEM/Alevin/Stringtie/none files")

parser$add_argument("--txinfo", required=TRUE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("-t", "--type", type="character", default="salmon",
                    help="Data origins (kallisto, salmon, sailfish, rsem, stringtie, alevin, none)")

parser$add_argument("-o", "--output", default="tximport.rds",
                    help="Output rds file")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


## add sample names to input files
nn <- as.character(sapply(args$input, function(x) {a <- strsplit(x, "/")[[1]]; a[[length(a)-1]]}))
names(args$input) <- nn

if (args$verbose == TRUE){
    print(args)
    options(echo=TRUE)
}


tx.info <- read.delim(args$txinfo, sep="\t")
tx2gene <- tx.info[,c("transcript_id", "gene_id")]


if (args$type %in% c("kallisto", "salmon", "stringtie", "sailfish", "rsem")){
    txi.tx <- tximport(args$input, type=args$type, tx2gene=tx2gene, txOut=TRUE)
} else if (args$type == "none"){
    txi.tx <- tximport(args$input, type=args$type, tx2gene=tx2gene,
                       importer=importer, txIn=TRUE, txOut=TRUE)
} else{
    stop(cat("method: ", args$type, " not valid \n"))
}

txi.tx$type <- args$type
txi.tx$tx2gene <- tx2gene

saveRDS(txi.tx, args$output, compress=FALSE)
