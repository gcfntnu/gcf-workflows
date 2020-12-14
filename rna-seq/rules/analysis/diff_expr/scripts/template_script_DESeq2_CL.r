################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### March 20th, 2018
### designed to be executed with SARTools 1.6.6
### run "Rscript template_script_DESeq2_CL.r --help" to get some help
################################################################################
rm(list=ls())                                        # remove all the objects from the R session
library(optparse)                                    # to run the script in command lines

# options list with associated default value.
option_list <- list( 
make_option(c("-P", "--projectName"),
            default=basename(getwd()),
            dest="projectName",
            help="name of the project used for the report [default: name of the current directory]."),

make_option(c("-A", "--author"),
            default=Sys.info()[7],
            dest="author",
            help="name of the report author [default: %default]."),

make_option(c("-t", "--targetFile"),
            default="target.txt",
            dest="targetFile",
            help="path to the design/target file [default: %default]."),

make_option(c("-s", "--subset"),
            default=NULL,
            dest="subset",
            help="subset samples, use comma separated list of samples with no spaces  [default: %default]."),

make_option(c("-r", "--countsFile"),
            default="data/tmp/rnaseq/quant/salmon/tximport/tx_salmon.rds",
            dest="countsFile",
            help="path to the count matrix file [default: %default]."),

make_option(c("-m", "--metaFile"),
            default="data/tmp/rnaseq/quant/salmon/tximport/gene_info.tsv",
            dest="metaFile",
            help="path to the features info file [default: %default]."),

make_option(c("-R", "--templateFile"),
            default="src/rna-seq/rules/analysis/diff_expr/scripts/GCF_DESeq2.rmd",
            dest="templateFile",
            help="path to the directory R markdown template [default: %default]."),

make_option(c("-F", "--featuresToRemove"),
            default="alignment_not_unique,ambiguous,no_feature,not_aligned,too_low_aQual",
            dest="FTR",
            help="names of the features to be removed, more than once can be specified [default: %default]"),
			
make_option(c("-v", "--varInt"),
            default="group",
            dest="varInt", 
            help="factor of interest [default: %default]"),

make_option(c("-c", "--condRef"),
            default="WT",
            dest="condRef",
            help="reference biological condition [default: %default]"),

make_option(c("-b", "--batch"),
            default=NULL,
            dest="batch",
            help="blocking factor [default: %default] or \"batch\" for example"),

make_option(c("-f", "--fitType"),
            default="parametric",
            dest="fitType", 
            help="mean-variance relationship: [default: %default],local or mean"),

make_option(c("-o", "--cooksCutoff"),
            default=TRUE,
            dest="cooksCutoff", 
            help="perform the outliers detection (default is TRUE)"),

make_option(c("-i", "--independentFiltering"),
            default=TRUE,
            dest="independentFiltering",
            help="perform independent filtering (default is TRUE)"),

make_option(c("-a", "--alpha"),
            default=0.05,
            dest="alpha", 
            help="threshold of statistical significance [default: %default]"),

make_option(c("-p", "--pAdjustMethod"),
            default="BH",
            dest="pAdjustMethod", 
            help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),

make_option(c("-T", "--typeTrans"),
            default="VST",
            dest="typeTrans", 
            help="transformation for PCA/clustering: \"VST\" ou \"rlog\" [default: %default]"),

make_option(c("-l", "--locfunc"),
            default="median",
            dest="locfunc", 
            help="median or shorth to estimate the size factors [default: %default]"),

make_option(c("-C", "--colors"),
            default="dodgerblue,firebrick1,MediumVioletRed,SpringGreen,chartreuse,cyan,darkorchid,darkorange",
            dest="cols",
            help="colors of each biological condition on the plots\n\t\t\"col1,col2,col3,col4\"\n\t\t[default: %default]"),

make_option(c("-O", "--output"),
            default="data/tmp/rnaseq/quant/salmon/sartools",
            dest="output",
            help="output directory [default: %default]"),

make_option(c("--forceCairoGraph"),
            action="store_true",
            default=FALSE,
            dest="forceCairoGraph",
            help="activate cairo type")

)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with DESeq2.",
					   epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
projectName <- opt$projectName                       # name of the project
author <- opt$author	                             # author of the statistical analysis/report
targetFile <- opt$targetFile                         # path to the design/target file
rawDir <- "."       						 # path to the directory containing raw counts files
countsFile <- opt$countsFile
output <- opt$output
featuresToRemove <- unlist(strsplit(opt$FTR, ","))   # names of the features to be removed (specific HTSeq-count information and rRNA for example)
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
fitType <- opt$fitType                               # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- opt$cooksCutoff                       # outliers detection threshold (NULL to let DESeq2 choosing it)
independentFiltering <- opt$independentFiltering     # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
typeTrans <- opt$typeTrans                           # transformation for PCA/clustering: "VST" ou "rlog"
locfunc <- opt$locfunc                               # "median" (default) or "shorth" to estimate the size factors
colors <- unlist(strsplit(opt$cols, ","))            # vector of colors of each biologicial condition on the plots
forceCairoGraph <- opt$forceCairoGraph				 # force cairo as plotting device if enabled

if (!is.null(opt$subset)){
    subset <- unlist(strsplit(opt$subset, ","))
} else{
    subset <- opt$subset
}

print(paste("workDir", workDir))
print(paste("subset", subset))
options(echo=TRUE)

# print(paste("projectName", projectName))
# print(paste("author", author))
# print(paste("targetFile", targetFile))
# print(paste("rawDir", rawDir))
# print(paste("varInt", varInt))
# print(paste("condRef", condRef))
# print(paste("batch", batch))
# print(paste("fitType", fitType))
# print(paste("cooksCutoff", cooksCutoff))
# print(paste("independentFiltering", independentFiltering))
# print(paste("alpha", alpha))
# print(paste("pAdjustMethod", pAdjustMethod))
# print(paste("typeTrans", typeTrans))
# print(paste("locfunc", locfunc))
# print(paste("featuresToRemove", featuresToRemove))
# print(paste("colors", colors))

################################################################################
###                             running script                               ###
################################################################################

library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
problem <- checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)
if (problem) quit(save="yes")
					   
## loading target file
## update original func to accept subset as argument
loadTargetFile <- function(targetFile, varInt, condRef, batch, subset=NULL){
    target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="", check.names=FALSE)
    rownames(target) <- as.character(target[,1])
    if (!is.null(subset)){
        target <- target[subset,]
        factor_cols <- vapply(target, is.factor, logical(1))
        target[factor_cols] <- lapply(target[factor_cols], factor)
    } 
    if (!I(varInt %in% names(target))) stop(paste("The factor of interest", varInt, "is not in the target file"))
    if (!is.null(batch)){
        if (!I(batch %in% names(target))){
            stop(paste("The batch effect", batch, "is not in the target file"))
        } else{
            target[,batch] <- as.factor(target[,batch])
        }
    }
    target[,varInt] <- as.factor(target[,varInt])
    if (!I(condRef %in% as.character(target[,varInt]))) stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
    target[,varInt] <- relevel(target[,varInt],ref=condRef)
    ##target <- target[order(target[,varInt]),]
    rownames(target) <- as.character(target[,1])
    ## check if varInt contains replicates
    if (min(table(target[,varInt]))<2){
        stop(paste("The factor of interest", varInt, "has a level without replicates"))
    }
    ## check if NA in the target
    if (any(is.na(cbind(target[,c(varInt, batch)], target[,1:2])))) stop("NA are present in the target file")
    ## warning message if batch is numeric
    ##if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
    ##if (any(grepl("[[:punct:]]", as.character(target[,varInt])))) stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
    cat("Target file:\n")
    print(target)
    return(target)
}

target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch, subset=subset)


# loading counts
##counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
counts <- read.delim(countsFile, sep="\t", check.names=FALSE, as.is=TRUE, row.names=1)
counts <- as.matrix(round(counts))
counts <- counts[,rownames(target)]

## load features meta info
info <- read.delim(opt$metaFile, sep="\t", check.names=FALSE, as.is=TRUE)
i <- grep("gene_id.?", colnames(info))[1]
rownames(info) <- info[,i]
keep.cols <- intersect(c("gene_id", "gene_name", "gene_biotype", "seqname"), colnames(info))
info <- info[,keep.cols]
info <- info[rownames(counts),]

#setwd(output)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

## add feature data if available
if (!is.null(info)){
    ids <- rownames(out.DESeq2$dds)
    info <- info[ids,]
    info[,"Id"] = ids
}



# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering, 
                                          cooksCutoff=cooksCutoff, alpha=alpha)



## lets patch up the export
exportResults.DESeq2 <- function(out.DESeq2, group, alpha=0.05, info=NULL, export=TRUE){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(counts(dds, normalized=TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(Id=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="Id", all=TRUE)
  tmp <- base[,paste("norm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }
  
  complete <- list()
  for (name in names(results)){
    complete.name <- base

    # ajout d'elements depuis results
    res.name <- data.frame(Id=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(info, complete.name, by="Id", all.y=TRUE)
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(Id=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name
    
    if (export){
        ## s?lection des up et down
        up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]
        up.name <- up.name[order(up.name$padj),]
        
        down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]
        down.name <- down.name[order(down.name$padj),]
      
        ## exports
        name <- gsub("_","",name)
        write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
        write.table(up.name, file=paste0("tables/", name,".up.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
        write.table(down.name, file=paste0("tables/", name,".down.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
    }
  }

  return(complete)
}

summaryResults$complete <- exportResults.DESeq2(out.DESeq2, group=target[,varInt], alpha=0.05, info, export=TRUE)

# save image of the R session
#save.image(file=paste0(projectName, ".RData"))


writeReport.DESeq2 <- function(target, counts, out.DESeq2, summaryResults, majSequences,
                               workDir, projectName, author, targetFile, rawDir,
                               featuresToRemove, varInt, condRef, batch, fitType,
                               cooksCutoff, independentFiltering, alpha, pAdjustMethod,
                               typeTrans, locfunc, colors, info=NULL){
  rmarkdown::render(input=opt$templateFile,
                    output_file=paste0(projectName, "_report.html"),
                    output_dir=workDir,
                    intermediates_dir=workDir,
                    knit_root_dir=workDir,
                    run_pandoc=TRUE,
                    quiet=TRUE,
                    clean=TRUE)
  cat("HTML report created\n")
}

## generating HTML report
viz.counts <- data.frame(Id=rownames(counts), counts)
viz.counts <- merge(info[, c("Id", "gene_name", "gene_biotype")], viz.counts, by='Id')
writeReport.DESeq2(target=target, counts=viz.counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)


## mv to output folder
if (!dir.exists(output)){
    dir.create(output, showWarnings=TRUE, recursive=TRUE)
}

file.rename("tables", file.path(output, "tables"))
file.rename("figures", file.path(output, "figures"))
report.name <- paste(projectName, "report.html", sep="_")
file.rename(report.name, file.path(output, report.name))

#unlink("tables", recursive = TRUE)
#unlink("figures", recursive = TRUE)
#unlink(report.name)
