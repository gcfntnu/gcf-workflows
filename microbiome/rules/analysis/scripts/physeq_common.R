#/usr/bin/env Rscript
### library imports ###     
if (!require(phyloseq)) {
    stop('This script needs the phyloseq package to be available')
}
suppressWarnings(library(phyloseq))
library(argparse)



parser <- ArgumentParser()

parser$add_argument("input", nargs="+", help="Input files")

parser$add_argument("--sample-info",
                    help="Sample info. Tab delimited file, needs column with `sample_id`")

parser$add_argument("--feature-info",
                    help="Feature info. Tab delimited file, needs column with `feature_id`")

parser$add_argument("--condition", type="character",
                    help="A column name in sample-info used to define condition (sample_group)")

parser$add_argument("--reference-level", type="character", default="ctrl",
                    help="Name of reference level in CONDITION")

parser$add_argument("--subset", help="Column name in samplesheet for subsetting samples/ remove outliers")

parser$add_argument("--subset-keep", help="Keep level in SUBSET. If set to `ALL` then a separate analysis will be done on each subset")

parser$add_argument("--batch", help="Column name in samplesheet for batch effect removal")

parser$add_argument("--block", type="character",
                    help="A column name in sample-info used to define blocks in a blocked design")

parser$add_argument("--test", type="character",
                    help="Test method")

parser$add_argument("--output-dir", default="models", help="Output prefix [optional]")

parser$add_argument("--logfile", help="Log filename")

parser$add_argument("--threads", help="number of threads", default=1)

parser$add_argument("--save-workspace", action="store_true",
                    default=FALSE, help="Save workspace to [output_dir]/workspace.Rdata")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose == TRUE) {
    print(args)
    options(echo=TRUE)
}

if (!is.null(args$output_dir)){
    if (!dir.exists(args$output_dir)){
        dir.create(args$output, showWarnings=FALSE, recursive=TRUE)
    }
}

read.sample.info <- function(sample.info.fn, condition.name, ref.level=NULL, blocked.name=NULL, batch.name=NULL,
                             subset.name=NULL, subset.keep=NULL, verbose=TRUE){
    ## sample-info dataframe
    stopifnot(file.exists(args$sample_info))
    sample.info <- read.delim(sample.info.fn, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
    sample.id.col <- grep("sample.?id", colnames(sample.info), ignore.case=TRUE)[1]
    if (is.na(sample.id.col)){
        warning("Failed to find `sample_id` column. Using first column as rownames")
        sample.id.col <- 1
    }
    sample.id.name <- colnames(sample.info)[sample.id.col]
    sample.ids <- as.character(sample.info[,sample.id.name])
    if (any(duplicated(sample.ids))){
        stop(sprintf("Duplicated sample-ids found in col (%s)", sample.id.name))
    }
    rownames(sample.info) <- sample.ids
    if (verbose == TRUE){
        print(sample.ids)
    }
    
    if (!is.null(condition.name)){
        if (!condition.name %in% colnames(sample.info)){
            err.msg <- sprintf("%s does not match %s, ", condition.name, colnames(sample.info))
            stop(err.msg)
        }
        condition <- factor(sample.info[,condition.name])
        if (!is.null(ref.level)){
           stopifnot(ref.level %in% levels(condition))
           condition <- relevel(condition, ref.level)
        }
        if (verbose == TRUE){
            print(condition)
            }
        sample.info[,condition.name] <- condition 
    }

    if (!is.null(blocked.name)){
        stopifnot(blocked.name %in% colnames(sample.info))
        blocks <- factor(sample.info[,blocked.name])
        block.sizes <- unique(table(blocks))
        if (!length(block.sizes) == 1){
            warning(sprintf("BLOCKS (%s) column has unbalanced experiment levels, %", blocked.name))
        }
        if (exists(condition) & (!all(table(paste(blocks, condition, sep=".")) >= 1))){
            warning(sprintf("BLOCKS (%s) column has levels do not match up with CONDITION (%s).", blocked.name, condition.name))
        }
        sample.info[,blocked.name] <- blocks
    }
    
    if (!is.null(batch.name)){
        stopifnot(batch.name %in% colnames(sample.info))
        batch <- factor(sample.info[,batch.name])
        sample.info[,batch.name] <- batch 
    }
    if (verbose == TRUE){
        print(head(sample.info))
    }
    
    sample.info
}



physeq <- readRDS(args$input)


if (!is.null(args$sample_info)){
    ## the general idea is that if sample_info is a subset of raw data, then subset raw data.
    ## If raw data (samples) is a subset of sample info, then fail
    sdata <- read.sample.info(args$sample_info, args$condition, ref.level=args$ref_level, blocked.name=args$block,
                              batch.name=args$batch, subset.name=args$subset, subset.keep=args$subset_keep,
                              verbose=args$verbose)
    
    ## update physeq obj
    sample.names.physeq <- phyloseq::sample_names(physeq)
    sample.names.sdata <- rownames(sdata)
    if (!all(sample.names.sdata %in% sample.names.physeq)){
        missing.samples <- sample.names.sdata[!sample.names.sdata %in% sample.names.physeq]
        stop(cat("Missing samples in raw data: ", missing.samples))
    }
    if (!all(sample.names.physeq %in% sample.names.sdata)){
        missing.samples <- sample.names.physeq[!sample.names.physeq %in% sample.names.sdata]
        physeq <- subset_samples(!missing.samples, physeq)
        sample.names.physeq <- phyloseq::sample_names(physeq)
        n.missing <- sum(missing.samples)
        cat(sprintf("%d samples removed from raw data to match sample info", n.missing), "\n")
        ## FIME: log sample names
        
    }
    sdata <- sdata[sample.names.physeq,]
    
    if (is.null(sample_data(physeq))) {
        sample_data(physeq) <- sample_data(sdata)
    }
    else {
        sdata <- phyloseq::sample_data(sdata)
        sdata.physeq <- as(sample_data(physeq), "data.frame")
        uniq.cols <- ! colnames(sdata) %in% colnames(sdata.physeq)
        sdata <- cbind(sdata.physeq, sdata[,uniq.cols])
        sample_data(physeq) <- sample_data(sdata)
    }
}

if (!is.null(args$subset) && !is.null(args$subset_keep)){
    ## use subset from sample info
    sdata <- sample_data(physeq)
    stopifnot(args$subset%in% colnames(sdata))
    subset.col <- sdata[,args$subset]
    keep.samples <- subset.col == args$subset_keep
    physeq <- phyloseq::subset_samples(physeq, keep.samples)
    sdata <- sdata[keep.samples,]
}

