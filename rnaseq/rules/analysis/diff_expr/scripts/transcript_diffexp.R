#/usr/bin/env Rscript

if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}



parser <- ArgumentParser()

parser$add_argument("input", nargs="+", help="Kallisto/Salmon files")

parser$add_argument("--transcript-info", required=TRUE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("--sample-info",
                    help="Sample info. Tab delimited file, needs column with `sample_id`")

parser$add_argument("-c", "--condition", type="character", default="sample_group",
                    help="A column name in sample-info used to define condition (sample_group)")

parser$add_argument("--reference-level", type="character", default="ctrl",
                    help="Name of reference level in `condition`")

parser$add_argument("--organism", type="character", default="homo_sapiens",
                    help="Organism. Used to fetch annotations from biomart if no --transcript-info is present.")

parser$add_argument("-o", "--output", default="auto",
                    help="Output prefix [optional]")

parser$add_argument("--output-tpm", action="store_true",
                    default=FALSE,
                    help="Output count (.tpm) TPM")

parser$add_argument("--output-rds", action="store_true",
                    default=FALSE,
                    help="Additional output of the sleuth result obj. (RDS)")

parser$add_argument("--output-transcripts", action="store_true",
                    default=FALSE,
                    help="Additional output Counts/TPM per transcript.")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose == TRUE) options(echo=TRUE)

if (!require(data.table)) {
    install.packages("data.table", repos="http://cran.rstudio.com")
    library("data.table")
}

if (!require(sleuth)) {
    source("http://bioconductor.org/biocLite.R")
    if (!require(rhdf5)){
        biocLite("rhdf5")
    }
    install.packages("devtools", repos="http://cran.rstudio.com")
    library(devtools)
    devtools::install_github("pachterlab/sleuth")
    library("sleuth")
}

create.s2c <- function(sample.info, condition.name, ref.level, paths){
    ## sample-to-covariates dataframe
    sample.info <- read.delim(sample.info, sep="\t", check.names=FALSE)
    sample.id.col <- grep("sample.?id", colnames(sample.info), ignore.case=TRUE)[1]
    rownames(sample.info) <- as.character(sample.info[,sample.id.col])
    stopifnot(condition.name %in% colnames(sample.info))
    condition <- factor(sample.info[,condition.name])
    stopifnot(ref.level %in% levels(condition))
    condition <- relevel(condition, ref.level)
    s2c <- data.frame(sample=rownames(sample.info), condition=condition)
    rownames(s2c) <- rownames(sample.info)
    if (!all(names(paths) %in% s2c$sample)){
        print(paths)
        print(s2c[names(paths),])
        stop("sample-info is missing some samples used in input.")
    }
    if (!all(s2c$sample %in% names(paths))){
        message("sample-info contains more samples compared to input. Subsetting sample info.")
        s2c <- s2c[names(paths),]
    }
    s2c <- s2c[names(paths),]
    s2c[,"path"] <- paths
    s2c
}

create.t2g <- function(transcript.info){
    ## transcript-to-gene dataframe
    tx.info <- fread(transcript.info, select=c("transcript_id",
                                               "gene_id",
                                               "transcript_version",
                                               "transcript_biotype",
                                               "gene_name"))
    tx.info <- subset(tx.info, transcript_id != "")
    tx.info[,transcript_id:=paste(transcript_id, transcript_version, sep=".")]
    tx.info <- subset(tx.info, !duplicated(transcript_id))
    tx2gene <- tx.info[,c("transcript_id", "gene_id", "gene_name", "transcript_biotype")]
    colnames(tx2gene) <- c("target_id", "gene_id", "gene_name", "transcript_biotype")
    tx2gene
}

input2paths <- function(files){
    stopifnot(all(file.exists(files)))
    files <- normalizePath(files)
    kallisto.paths <- dirname(files)
    sample.ids <- sapply(strsplit(kallisto.paths, .Platform$file.sep), function(x) x[[length(x)]])
    names(kallisto.paths) <- as.character(sample.ids)
    kallisto.paths
}

main <- function(args){
    if (args$verbose == TRUE) print(args)
    
    paths <- input2paths(args$input)

    ## parse sample info
    stopifnot(file.exists(args$sample_info))
    s2c <- create.s2c(args$sample_info, args$condition, args$reference_level, paths)
    
    ## parse transcript info if available
    if (!is.null(args$transcript_info)){
        stopifnot(file.exists(args$transcript_info))
        t2g <- create.t2g(args$transcript_info)
    }else{
        mart <- biomaRt::useMart(biomart = "ensembl", dataset=args$organism)
        t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=mart)
        colnames(t2g) <- c("target_id", "ens_gene", "ext_name")
    }
    
    ## create sleuth object
    so <- sleuth_prep(s2c, ~condition, target_mapping=t2g)
    ## fit model, all against ref
    design <- model.matrix(~condition, s2c)
    null <-  model.matrix(~1, s2c)
    so <- sleuth_fit(so, design, "full")
    so <- sleuth_fit(so, null, "reduced") 
    ##so <- sleuth_lrt(so, "reduced", "full") 

    levs <- levels(s2c$condition)
    for (cond in tail(levs, length(levs)-1)){
       so <- sleuth_wt(so, which_beta=paste0("condition", cond))  
    }
    
    ## write output
    if (args$output == "auto"){
        prefix <- "kallisto.transcript"
    } else{
        prefix <- args$output
    }
    
    for (cond in tail(levs, length(levs)-1)){
        test.df <- sleuth_results(so, paste0("condition",cond), test_type="wald")
        test.df.fn <- paste0(prefix, ".", cond, ".testtable")
        rownames(test.df) <- test.df[,"target_id"]
        colnames(test.df)[1] <- "transcript_id"
        fwrite(test.df, test.df.fn, row.names=TRUE, sep="\t")
        if (args$verbose == TRUE){
            message(cat("wrote", test.df.fn))
        }
    }
    
    
    ## write RDS object that may be shared for live session
    if (args$output_rds == TRUE){
        rds.fn <- "/tmp/so.RDS"
        saveRDS(so, file=rds.fn)
    }
}


main(args)
