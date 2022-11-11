library(argparse)
library(fishpond)
library(EmptyNN)

parser <- ArgumentParser()
parser$add_argument("-i", "--input", default="data/tmp/test_sample/outs/filtered_gene_bc_matrices",
                    help="Alevin output dir. Assumes files: alevin/quants_mat_rows.txt, alevin/quants_mat_cols.txt, alevin/quants_mat.csv")
parser$add_argument("-o", "--output", default="data/processed/test_sample/",
                    help="Output filename (rds file)")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)


sce <- fishpond::loadFry(args$input, outputFormat="velocity")

empty <- EmptyKNN::emptynn(assay(sce))

#A <- as.data.frame(empty$knn)
#write.table(A, file=args$output, sep="\t", row.names=FALSE)

saveRDS(empty, file=args$output)
