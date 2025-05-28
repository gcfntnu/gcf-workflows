#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scDblFinder)
  library(BiocParallel)
})


normalize_barcodes <- function(bc) sub("-\\d+$", "-1", bc)

setup_threads <- function(n) {
  if (.Platform$OS.type == "windows") {
    register(SnowParam(n, progressbar = TRUE), default = TRUE)
  } else {
    register(MulticoreParam(n, progressbar = TRUE), default = TRUE)
  }
  cat(sprintf("[INFO] Using %d thread(s)\n", n), file = stderr())
}

main <- function() {
  # Argument parsing
  parser <- ArgumentParser()
  parser$add_argument("-o", "--output", required = TRUE, help = "Output file")
  parser$add_argument("-i", "--input", required = TRUE, type = "character", help = "Path to the 10x filtered matrix.mtx")
  parser$add_argument("-t", "--threads", type = "integer", default = NULL, help = "Number of threads to use")
  args <- parser$parse_args()

  # Threads
  n_threads <- if (!is.null(args$threads)) args$threads else as.integer(Sys.getenv("THREADS", unset = "1"))
  if (is.na(n_threads) || n_threads < 1) stop("[ERROR] Invalid thread count")
  setup_threads(n_threads)

  # I/O
  if (!file.exists(args$input)) stop("[ERROR] Input file does not exist")
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

  barcodes.fn <- file.path(dirname(args$input), "barcodes.tsv.gz")
  if (!file.exists(barcodes.fn)) {
    barcodes.fn <- file.path(dirname(args$input), "barcodes.tsv")
  }
  stopifnot(file.exists(barcodes.fn))

  cat("[INFO] Reading input matrix...\n", file = stderr())
  con <- TENxMTX(args$input)
  sce <- import(con)
  obs.names <- import(TENxTSV(barcodes.fn))
  colnames(sce) <- obs.names$barcode
  sce <- as(sce, "SingleCellExperiment")

  # Run
  doublet_ratio <- ncol(sce) / 1000 * 0.008
  cat("[INFO] Running scDblFinder...\n", file = stderr())
  sce <- scDblFinder(sce, dbr = doublet_ratio)

  # Output
  barcodes <- normalize_barcodes(rownames(colData(sce)))
  results <- data.frame(
    Barcode = barcodes,
    doublet = sce$scDblFinder.class,
    doublet_score = sce$scDblFinder.score
  )
  write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[INFO] Done.\n", file = stderr())
}

main()
