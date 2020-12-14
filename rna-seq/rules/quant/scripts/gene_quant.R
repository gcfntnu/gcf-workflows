#/usr/bin/env Rscript

if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}



parser <- ArgumentParser()

parser$add_argument("input", nargs='+', help="Gene quantification (featurecounts, star, htseq) files.")

parser$add_argument("-t", "--type", type="character", default="featurecounts",
                    help="Data origins (featurecounts, htseq, star)[default %(default)s]")

parser$add_argument("-o", "--output", default="./featurecounts", dest="prefix",
                    help="Output prefix")

parser$add_argument("--gene-info", nargs="?",
                    help="Gene info (optional). Tab delimited file, needs column with `gene_id`")

parser$add_argument("--sample-info", nargs="?",
                    help="Sample info (optional). Tab delimited file, needs column with `sample_id`. If `fragment-length` is a column then this will override --fragment-length in the TPM estimation.")

parser$add_argument("--gene-length", nargs="?",
                    help="Gene length (optional). Tab delimited file of (effective) gene lengths per sample. This will override any gene length derived from input files.")

parser$add_argument("--fragment-length", default=200,
                    help="Mean fragment length. Used for TPM estimation.")

parser$add_argument("--output-tpm", action="store_true",
                    default=FALSE,
                    help="Additional output of TPM, Transcript Per Million (.quant) file")

parser$add_argument("--output-genelength", action="store_true",
                    default=FALSE,
                    help="Additional output of genelength")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


if (!require(data.table)) {
    install.packages("data.table", repos="http://cran.rstudio.com")
    library("data.table")
}


read.featurecounts.multisample <- function(filename){
    tab <- fread(filename, skip=1)
    fdata <- as.data.frame(tab[,1:6])
    counts <- as.data.frame(tab[,-seq(6)])
    rownames(fdata) <- fdata[,"Geneid"]
    rownames(counts) <- rownames(fdata)
    bams <- sapply(colnames(counts), basename)
    sample.ids <- sapply(strsplit(bams, "\\."), function(x) x[[1]])
    colnames(counts) <- sample.ids
    length <- counts
    for (i in seq(ncol(counts))){
        length[,i] <- fdata[,"Length"]
    }
    return(list(counts=counts, length=length))
}

read.gene.count <- function(filename, indices=c(1, 2, NULL)){
    gene.id.idx <- indices[1]
    counts.idx <- indices[2]
    length.idx <- indices[3]
    tab <- as.data.frame(fread(filename, skip=1))
    ids <- as.character(tab[ ,gene.id.idx])
    counts <- tab[ ,counts.idx]
    names(counts) <- ids
    if (!(is.null(length.idx) | is.na(length.idx))){
        length <- tab[, length.idx]
        names(length) <- ids
    } else{
        length <- NULL
    }
    return(list(counts=counts, "gene.ids"=ids, length=length))
}

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
    ## https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44

    ## note: mean fragment length per library
    ## af: added support for one mean fragment length across all samples
    if (length(meanFragmentLength) == 1){
        meanFragmentLength <- rep(meanFragmentLength, ncol(counts))
    }

    stopifnot(length(featureLength) == nrow(counts))
    stopifnot(length(meanFragmentLength) == ncol(counts))
    
  
    effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        featureLength - meanFragmentLength[i] + 1
    }))
    
    ## Exclude genes with length less than the mean fragment length.
    idx <- apply(effLen, 1, function(x) min(x) > 1)
    counts <- counts[idx,]
    effLen <- effLen[idx,]
    featureLength <- featureLength[idx]
    
    tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        rate = log(counts[,i]) - log(effLen[,i])
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }))
    
    colnames(tpm) <- colnames(counts)
    rownames(tpm) <- rownames(counts)
    return(tpm)
}


main <- function(args){
    stopifnot(all(file.exists(args$input)))
    type <- match.arg(args$type, c("featurecounts", "star"))
    length <- NULL
    if (length(args$input) == 1){
        if (args$type == "featurecounts"){
            out <- read.featurecounts.multisample(args$input)
            counts <- out$counts
            length <- out$length
        }
    } else{
        indices <- c(1, 2, NULL)
        if (type == "featurecounts") indices <- c(1, 7, 6)
        if (type == "star") indices <- c(1, 3, NULL)
        
        for (i in seq_along(args$input)){
                out <- read.gene.count(args$input[i], indices=indices)

                if (i == 1){
                    ids.1 <- out$gene.ids
                    n.genes <- length(ids.1)
                    n.samples <- length(args$input)
                    counts <- matrix(nrow=n.genes, ncol=n.samples)
                    rownames(counts) <- ids.1
                    if (!is.null(out$length)){
                        length <- matrix(nrow=n.genes, ncol=n.samples)
                    } else{
                        length <- NULL
                    }
                }
                
                ids <- out$gene.ids
                stopifnot(all(ids == ids.1)) #fixme: replace with outer join?
                counts[,i] <- out$counts
                if (!is.null(length)) length[,i] <- out$length
        }
        ## assume filenames are on the form '{sample-id}.*'
        ##basenames <- sapply(args$input, basename)
        ##sample.ids <- sapply(strsplit(basenames, "\\."), function(x) x[[1]])

        ## sample-id is  parent dir
        sample.ids <- sapply(sapply(args$input, dirname), basename)

        colnames(counts) <- sample.ids
        if (!is.null(length)){
            colnames(length) <- sample.ids
        }
    }

    
    gene.quant.fn <- paste0(args$prefix, ".gene.quant")
    write.table(counts, file=gene.quant.fn, sep="\t", row.names=TRUE)
    if (args$verbose == TRUE){
        message(cat("wrote", gene.quant.fn))
    }
    
    if (args$output_tpm == TRUE){
        if (is.null(length)){
            stop("no length calculations available from input. Use --gene-info with a column named `gene_length`")
        }
        length <- rowMeans(length) #fixme
        tpm <- counts_to_tpm(counts, length, args$fragment_length)
        if (args$verbose == TRUE){
            filtered.genes <- setdiff(rownames(counts), rownames(tpm))
            n.filtered <- length(filtered.genes)
            ##fixme: filter count table / gene_info also? 
            message(cat(n.filtered, "genes were filtered out in tpm table due to gene_length < fragment_length"))
        }
        tpm.fn <- paste0(args$prefix, ".gene.tpm")
        write.table(tpm, file=tpm.fn, sep="\t", quote=FALSE)
        if (args$verbose == TRUE){
            message(cat("wrote", tpm.fn))
        }
    }
    if (!is.null(args$gene_info)){
        gene.info <- as.data.frame(fread(args$gene_info))
        rownames(gene.info) <- gene.info[,"gene_id"]
        gene.info <- gene.info[rownames(counts),]
        gene.info.fn <- paste0(args$prefix, ".gene_info.tsv")
        write.table(gene.info, file=gene.info.fn, sep="\t")
        if (args$verbose == TRUE){
            message(cat("wrote", gene.info.fn))
        }
    }
    
    if (args$output_genelength == TRUE){
        if (is.null(length)){
            stop("no length calculations available from input. Use --gene-info with a column named `gene_length`")
        }
        length.fn <- paste0(args$prefix, ".gene.length")
        write.table(as.data.frame(length), length.fn, sep="\t")
        if (args$verbose == TRUE){
            message(cat("wrote", length.fn))
        }
    }
}


if (args$verbose == TRUE){
    print(args)
    options(echo=TRUE)
    #debug(main)
}


main(args)
