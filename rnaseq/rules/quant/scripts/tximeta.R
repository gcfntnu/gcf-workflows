#/usr/bin/env Rscript

library(argparse)
library(tximport)
library(tximeta)
library(readr)


parser <- ArgumentParser(description="Transcript quantification from Kallisto/Salmon/Sailfish/RSEM")

parser$add_argument("input", nargs="+", help="Kallisto/Salmon/Sailfish/RSEM/Alevin/Stringtie/none files")

parser$add_argument("--txome", required=TRUE, help="linked txome file (.json)")

parser$add_argument("--sample-info", required=TRUE, default="sample_info.tsv", help="sample info file (.tsv)")

parser$add_argument("-t", "--method", type="character", default="salmon",
                    help="Data origins (kallisto, salmon, sailfish, rsem, stringtie, alevin, none)")

parser$add_argument("-o", "--output", default="tximeta.rds",
                    help="Output rds file")

parser$add_argument("-c", "--cachedir", default="/tmp/.cache",
                    help="Cachedir location")

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

setTximetaBFC(args$cachedir, quiet=TRUE)

tximeta::loadLinkedTxome(args$txome)

desc <- read.delim(args$sample_info, header=TRUE, as.is=TRUE, sep="\t")
rownames(desc) <- desc[,"Sample_ID"]
desc <- desc[nn,]
remove.cols <- intersect(c("R1", "R2", "Organism", "Project_ID", "paired_end", "Fragment_Length"), colnames(desc))
coldata = desc[,!colnames(desc) %in% remove.cols]
coldata$names <- coldata[,"Sample_ID"]
coldata <- cbind(coldata, files = args$input, stringsAsFactors = FALSE)

if (args$method %in% c("salmon", "sailfish", "alevin")){
    txi.tx <- tximeta::tximeta(coldata, type=args$method, skipMeta=TRUE)
} else{
    stop(cat("method: ", args$method, " not valid \n"))
}

txi.tx$type <- args$method
txi.tx$method <- args$method
#txi.tx$tx2gene <- summarizeToGene(txi.tx)

saveRDS(txi.tx, args$output, compress=FALSE)
