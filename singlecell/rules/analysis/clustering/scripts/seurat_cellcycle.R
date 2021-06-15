
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(argparse))


parser <- ArgumentParser()
parser$add_argument("-i", "--input", default="data/processed/singlecell/analysis/preprocess/test_sample/test_sample.rds",
                    help="Preprocessed single cell data (Serialized Seurat obj)", required=TRUE)
parser$add_argument("-o", "--output", default="data/processed/singlecell/analysis/clustering/test_sample/test_sample.rds",
                    help="Output filename (rds file)")
parser$add_argument("-m", "--markers", default="data/ext/regev_lab_cell_cycle_genes.txt",
                    help="Output filename (rds file)")
parser$add_argument"--write-scores", action="store_true",
                    default=FALSE, help="Write the cellcycle scores to a csv file")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)

cc.genes <- readLines(con=args$markers)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

res <- CellCycleScoring(object=data, s.genes=s.genes, g2m.genes=g2m.genes, 
                        set.ident = TRUE)
saveRDS(res, args$output)

if (args$write_scores == TRUE){
    scores <- res@meta.data[,c("S.score", "G2m.Score", "Phase")]
    write.table(scores, file=paste0(args$output, ".csv"), sep="\t", quote=FALSE)
}
