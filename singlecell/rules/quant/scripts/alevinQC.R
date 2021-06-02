
library(alevinQC)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--input",
                    help="Alevin output dir. Assumes files: alevin/quants_mat_rows.txt, alevin/quants_mat_cols.txt, alevin/quants_mat.csv")
parser$add_argument("-o", "--output", default="data/processed/alevin/test_sample/alevinqc/report.html",
                    help="Output filename (html file)")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)

if (! dir.exists(args$input)){
      stop("Input directory provided does not exist")
}



report <- alevinQCReport(
    baseDir=args$input,
    sampleId=dirname(args$input),
    outputFile=basename(args$output),
    outputDir=dirname(args$output),
    outputFormat='html_document',
    forceOverwrite=TRUE,
    knitrProgress=FALSE,
    quiet=TRUE
)

