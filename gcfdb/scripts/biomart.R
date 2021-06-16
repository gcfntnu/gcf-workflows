suppressPackageStartupMessages({
    library(argparse)
    library(biomaRt)
    library(GRanges)
    library(Biostrings)
    library(stringr)
    library(data.table)
})


parser <- ArgumentParser(description="biomart queries")

parser$add_argument("-O", "--organism", type="character", default="homo sapiens", help="Organsim")

parser$add_argument("-g", "--gene-info", type="character", required=TRUE, help="gene info tsv filename")

parser$add_argument("-s", "--source", type="character", help="reference database", default="ensembl")

parser$add_argument("-r", "--release", type="character", help="reference database release", default="103")

parser$add_argument("-a", "--assembly", type="character", help="reference database assembly", default="GRCh38")

parser$add_argument("-o", "--output", default="biomart.txt", help="Output file")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


if (args$verbose == TRUE){
    options(echo=TRUE)
}

## format organism name
args$organism <- str_replace(args$organism, " ", "_")
dummy <- unlist(str_split(args$organsim, "_"))
org <- paste0(str_sub(dummy[1], 1, 1), dummy[2])

## read gtf (tsv type)
genes <- fread(args.gene_info)
gene <- genes[genes$type=="gene",]

if (args$assembly=="GRCh37"){
    args$release == "GRCh37"
}
ens.archives <- listEnsemblArchives()
if (!args$release %in% ens.archives$version){
    stop(paste("failed to find ensembl release:". args$release))
}
if (args$release == "current"){
    args$release <- ens.archives[ens.archives$current_release == "*", "version"]
}
ens.info <- ens.archives[ens.archives==args$release,]


ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host=ens.info$url)
ds <- listDatasets(ensembl)[,"dataset"]
ds <- grep(paste0("^", org), ds, value=TRUE)
if(length(ds) == 0)
    stop(paste("Mart not found for:", org))
else if(length(ds) > 1)
{
    message("Found several marts")
    sapply(ds, function(d)
        message(paste(which(ds==d), d, sep=": ")))
    n <- readline(paste0("Choose mart (1-", length(ds),") : "))
    ds <- ds[as.integer(n)]
}
ensembl <- useDataset(ds, mart=ensembl)

message( paste0( "Downloading sequence",
                ifelse(length(id) > 1, "s", ""), " ..."))
if(length(id) > 100) message("This may take a few minutes ...")

## download sequence
## (1) get exon coordinates
attrs <- c(id.type, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
coords <- getBM(filters=id.type, attributes=attrs, values=id, mart=ensembl)
id <- unique(coords[,id.type])
coords <- GRangesList(sapply(id,function(i) {i.coords <- coords[coords[,1]== i, 3:5] g <- GRanges(i.coords[,1], IRanges(i.coords[,2],i.coords[,3])) return(g)}), compress=FALSE)
coords <- reduce(coords)
len <- sum(width(coords))

## (2) get genes and sequences
sel <- c(id.type, "start_position", "end_position")
gene.pos <- getBM(attributes = sel, filters=id.type, values=id, mart=ensembl)
gene.seqs <- getSequence(id=id, type=id.type, seqType="gene_exon_intron", mart=ensembl)

## (3) get exonic sequences and correspondig GC content
gc.cont <- sapply(id,
                  function(i)
                  {
                      ## exon coordinates, gene position & sequence for current id i
                      ecoords <- coords[[i]]
                      gpos <- gene.pos[gene.pos[,id.type] == i,
                                       c("start_position", "end_position")]
                      gseq <- DNAString(
                          gene.seqs[gene.seqs[,id.type] == i, "gene_exon_intron"])
                      
                      ## exon coordinates relative to gene position
                      start <- start(ranges(ecoords)) - gpos[1,1] + 1
                      end <- end(ranges(ecoords)) - gpos[1,1] + 1
                      eseq <- gseq[IRanges(start, end)]
                      gc.cont <- sum(alphabetFrequency(eseq, as.prob=TRUE)[c("C","G")])
                      return(gc.cont)
                  }
                  )
}

res <- cbind(len, gc.cont)
colnames(res) <- c("length", "gc")
rownames(res) <- id
