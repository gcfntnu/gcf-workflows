#/usr/bin/env Rscript

if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}


parser <- ArgumentParser()

parser$add_argument("input", nargs="+", help="Tab delimited file with gene counts.")

parser$add_argument("-m", "--method", default="voom",
                    help="Differential expression method. (default: voom)")

parser$add_argument("-g", "--gene-info",
                    help="Gene annotations. Tab delimited file, needs columns with `gene_id`")

parser$add_argument("-s", "--sample-info",
                    help="Sample info. Tab delimited file, needs column with `sample_id`")

parser$add_argument("-c", "--condition", type="character", default="condition",
                    help="A column name in sample-info used to define condition (sample_group)")

parser$add_argument("-p", "--pair", type="character", default="pair",
                    help="A column name in sample-info used to define pairs in a paired t-test. Only used if --test is set to `paired`")

parser$add_argument("-r", "--reference-level", type="character", default="ctrl",
                    help="Name of reference level in `condition`")

parser$add_argument("--normalization", type="character", default="TMM",
                    help="Between sample normalization")

parser$add_argument("--organism", type="character", default="homo_sapiens",
                    help="Organism. Used to fetch annotations from biomart if no --gene-info is present.")

parser$add_argument("--test", type="character", default="ALLvsREF",
                    help="Statistical test. (ALLvsREF, ALLvsALL, paired)")

parser$add_argument("-o", "--output", default="auto",
                    help="Output prefix [optional]")

parser$add_argument("--save-workspace", action="store_true",
                    default=FALSE, help="Save workspace (saveRDS) ")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

if (args$verbose == TRUE) options(echo=TRUE)

if (!require(data.table)) {
    install.packages("data.table", repos="http://cran.rstudio.com")
    library("data.table")
}

if (!require(readr)) {
    install.packages("readr", repos="http://cran.rstudio.com")
    library("readr")
}

if (!require(limma)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("limma")
    library("limma")
}

if (!require(edgeR)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("edgeR")
    library("edgeR")
}


create.s2c <- function(sample.info, condition.name, ref.level, paired.name=NULL){
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
    if (!is.null(paired.name)){
        stopifnot(paired.name %in% colnames(sample.info))
        pairs <- factor(sample.info[,paired.name])
        stopifnot(all(table(pairs) == 2))
        stopifnot(all(table(paste(pairs, condition, sep=".")) == 1))
        s2c[,"pairs"] <- pairs
    }
    s2c
}

##FIXME: this is specific for gtf dervied feature info, need to generalize so to can work for other scenarios
create.g2a <- function(gene.info){
    ## gene-to-annotation dataframe
    gene.info <- fread(gene.info, select=c("gene_id",
                                           "gene_name",
                                           "gene_biotype"))
    gene.info <- subset(gene.info, gene_id != "")
    gene.info <- subset(gene.info, !duplicated(gene_id))
    gene.info <- data.frame(gene.info)
    rownames(gene.info) <- gene.info[,"gene_id"]
    gene.info
}

create.g2a <- function(gene.info){
    ## gene-to-annotation dataframe
    gene.info <- read.delim(gene.info, sep="\t")
    gene.info <- gene.info[,-1]
    rownames(gene.info) <- gene.info[,1]
    gene.info
}

input.check <- function(filename){
    stopifnot(file.exists(filename))
    counts <- read.delim(filename, sep="\t", check.names=FALSE, row.names=1)
    sample.ids <- colnames(counts)
    if (length(sample.ids) != length(unique(sample.ids))){
        stop("Sample ids not unique!")
    }
    counts <- round(as.matrix(counts))
    counts
}


fun.voom <- function(counts, s2c, g2a, args){
    require(edgeR)
    require(limma)
    ## create voom object
    dge <- edgeR::DGEList(counts, group=s2c[,"condition"])
    dge <- edgeR::calcNormFactors(dge, args$normalization)
    if (args$test == "paired"){
        if (is.null(args$pair)) stop("missing pair info for a paired test.")
        design <- model.matrix(~pairs + condition, data=s2c)
    } else{
        design <- model.matrix(~condition, data=s2c)
    }

    ## filter non-expressed genes
    keep <- rowSums(cpm(dge)>1) >= 3
    n.filtered <- sum(keep)
    print(cat("Filtered ", n.filtered, "out of", length(keep),"\n"))
    dge <- edgeR::DGEList(counts[keep,], group=s2c[,"condition"])
    dge <- edgeR::calcNormFactors(dge, args$normalization)
    g2a <- g2a[rownames(dge),]
    
    v <- voom(dge, design)
    
    fit <- eBayes(lmFit(v, design))

    if (args$verbose == TRUE){
        print(summary(decideTests(fit)))
    }

    ## write output
    if (args$output == "auto"){
        prefix <- "voom"
    } else{
        prefix <- args$output
    }
    
    levs <- levels(s2c$condition)
    for (cond in tail(levs, length(levs)-1)){
        coef <- paste0("condition", cond)
        tab <- topTable(fit, coef=coef, adjust="BH", n=Inf, sort.by="P", genelist=g2a)
        test.df.fn <- paste0(prefix, ".", cond, ".testtable")
        fwrite(tab, test.df.fn, row.names=FALSE, sep="\t")
        if (args$verbose == TRUE){
            message(cat("wrote", test.df.fn, "\n"))
        }
    }
    
}


fun.deseq2 <- function(counts, s2c, g2a, args){
    if (!require(DESeq2)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("DESeq2")
        library("DESeq2")
    }

    ## create deseq object
    dge <- edgeR::DGEList(counts, group=s2c[,"condition"])
    dge <- edgeR::calcNormFactors(dge, args$normalization)
    if (args$test == "paired"){
        if (is.null(args$pair)) stop("missing pair info for a paired test.")
        design <- model.matrix(~pairs + condition, data=s2c)
    } else{
        design <- model.matrix(~condition, data=s2c)
    }

    ## filter non-expressed genes
    keep <- rowSums(cpm(dge)>1) >= 3
    n.filtered <- sum(keep)
    print(cat("Filtered ", n.filtered, "out of", length(keep),"\n"))
    dge <- edgeR::DGEList(counts[keep,], group=s2c[,"condition"])
    dge <- edgeR::calcNormFactors(dge, args$normalization)
    g2a <- g2a[rownames(dge),]
    
    v <- voom(dge, design)
    
    fit <- eBayes(lmFit(v, design))

    if (args$verbose == TRUE){
        print(summary(decideTests(fit)))
    }

    ## write output
    if (args$output == "auto"){
        prefix <- "voom"
    } else{
        prefix <- args$output
    }
    
    levs <- levels(s2c$condition)
    for (cond in tail(levs, length(levs)-1)){
        coef <- paste0("condition", cond)
        tab <- topTable(fit, coef=coef, adjust="BH", n=Inf, sort.by="P", genelist=g2a)
        test.df.fn <- paste0(prefix, ".", cond, ".testtable")
        fwrite(tab, test.df.fn, row.names=FALSE, sep="\t")
        if (args$verbose == TRUE){
            message(cat("wrote", test.df.fn, "\n"))
        }
    }
    
}


main <- function(args){
    if (args$verbose == TRUE) print(args)
    
    counts <- input.check(args$input)

    ## parse sample info
    stopifnot(file.exists(args$sample_info))
    s2c <- create.s2c(args$sample_info, args$condition, args$reference_level)
    stopifnot(all(rownames(s2c) %in% colnames(counts)))
    counts <- counts[,rownames(s2c)]
    
    ## parse gene info if available
    if (!is.null(args$gene_info)){
        stopifnot(file.exists(args$gene_info))
        g2a <- create.g2a(args$gene_info)
    }else{
        mart <- biomaRt::useMart(biomart = "ensembl", dataset=args$organism)
        g2a <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=mart)
    }
    g2a <- g2a[rownames(counts),]

    ## run stat method
    if (args$method == "voom"){
        fun.voom(counts, s2c, g2a, args)
    }
    
    ## write RDS object 
    if (args$save_workspace== TRUE){
        rds.fn <- "workspace.Rdata"
        save.image(file=rds.fn)
    }
    
}



main(args)
