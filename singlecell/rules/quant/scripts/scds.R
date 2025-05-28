#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scds)
  library(SingleCellExperiment)
  library(BiocParallel)
})

# Utils
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
  parser$add_argument("--threshold", type = "double", default = 0.5, help = "Score threshold to classify doublets")
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

  # Run scds
  cat("[INFO] Running scds (cxds_bcds_hybrid)...\n", file = stderr())
  sce <- cxds_bcds_hybrid(sce, estNdbl = TRUE)

  barcodes <- normalize_barcodes(rownames(colData(sce)))
  results <- data.frame(
    Barcode = barcodes,
    doublet = ifelse(colData(sce)$hybrid_score > args$threshold, "doublet", "singlet"),
    doublet_score = colData(sce)$hybrid_score
  )

  write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[INFO] Done.\n", file = stderr())
}

main()
