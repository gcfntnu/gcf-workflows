library(Seurat)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", default="data/tmp/test_sample/outs/filtered_gene_bc_matrices", nargs="+",
                    help="Alevin output dir. Assumes files: alevin/quants_mat_rows.txt, alevin/quants_mat_cols.txt, alevin/quants_mat.csv")
parser$add_argument("-s", "--sample-info", dest="samples", default=NULL,
                    help="Sample info. Tab delimited file, needs column with `sample_id`")
parser$add_argument("-o", "--output", default="data/processed/test_sample/seurat_obj.rds",
                    help="Output filename (rds file)")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)


# from: https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/
# Parts of the function is taken from Seurat's Read10x parsing function
ReadAlevin <- function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }

    barcode.loc <- file.path( base.path, "alevin", "quants_mat_rows.txt" )
    gene.loc <- file.path( base.path, "alevin", "quants_mat_cols.txt" )
    matrix.loc <- file.path( base.path, "alevin", "quants_mat.csv" )
    if (!file.exists( barcode.loc )){
      print(barcode.loc)
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc) ){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc )){
      stop("Expression matrix file missing")
    }
    matrix <- as.matrix(read.csv( matrix.loc, header=FALSE))
    matrix <- t(matrix[,1:ncol(matrix)-1])

    cell.names <- readLines( barcode.loc )
    gene.names <- readLines( gene.loc )

    colnames(matrix) <- cell.names
    rownames(matrix) <- gene.names
    matrix[is.na(matrix)] <- 0
    return(matrix)
}

alv.data <- ReadAlevin(args$input)
sdata <- CreateSeuratObject(alv.data, min.cells = 3, project="10X_GCF")

saveRDS(sdata, args$output)
