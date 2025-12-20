#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scds)
  library(SingleCellExperiment)
  library(BiocParallel)
})

setup_threads <- function(n) {
  if (.Platform$OS.type == "windows") {
    BiocParallel::register(BiocParallel::SnowParam(n, progressbar = TRUE), default = TRUE)
  } else {
    BiocParallel::register(BiocParallel::MulticoreParam(n, progressbar = TRUE), default = TRUE)
  }
  cat(sprintf("[INFO] Using %d thread(s)\n", n), file = stderr())
}

read_barcodes <- function(input_mtx) {
  # Expect barcodes.tsv(.gz) next to matrix.mtx(.gz)
  d <- dirname(input_mtx)
  candidates <- c(file.path(d, "barcodes.tsv.gz"), file.path(d, "barcodes.tsv"))
  barcodes_fn <- candidates[file.exists(candidates)][1]
  if (is.na(barcodes_fn)) stop("[ERROR] barcodes.tsv(.gz) not found next to input matrix")

  df <- import(TENxTSV(barcodes_fn))
  if (!is.data.frame(df) || ncol(df) < 1) stop("[ERROR] Failed to read barcodes.tsv(.gz)")
  bc <- df[[1]]
  if (!is.character(bc)) bc <- as.character(bc)

  # Trim whitespace and drop empty
  bc <- trimws(bc)
  bc <- bc[nzchar(bc)]
  if (length(bc) == 0) stop("[ERROR] No barcodes read from barcodes.tsv(.gz)")

  bc
}

main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-o", "--output", required = TRUE, help = "Output TSV")
  parser$add_argument("-i", "--input", required = TRUE, type = "character",
                      help = "Path to the 10x filtered matrix.mtx(.gz)")
  parser$add_argument("-t", "--threads", type = "integer", default = NULL,
                      help = "Number of threads to use")
  parser$add_argument("--threshold", type = "double", default = 0.5,
                      help = "hybrid_score threshold to classify doublets")
  args <- parser$parse_args()

  n_threads <- if (!is.null(args$threads)) args$threads else as.integer(Sys.getenv("THREADS", unset = "1"))
  if (!is.finite(n_threads) || is.na(n_threads) || n_threads < 1) stop("[ERROR] Invalid thread count")
  setup_threads(as.integer(n_threads))

  if (!file.exists(args$input)) stop("[ERROR] Input file does not exist")
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

  cat("[INFO] Reading input matrix...\n", file = stderr())
  con <- TENxMTX(args$input)
  sce <- import(con)
  sce <- as(sce, "SingleCellExperiment")

  bc <- read_barcodes(args$input)
  if (ncol(sce) != length(bc)) {
    stop(sprintf("[ERROR] Matrix has %d cells but barcodes.tsv has %d rows",
                 ncol(sce), length(bc)))
  }

  # IMPORTANT: do NOT normalize "-2" -> "-1" etc here. That can create duplicates.
  if (anyDuplicated(bc)) {
    stop("[ERROR] Duplicate barcodes found in barcodes.tsv(.gz). Refuse to merge lanes silently.")
  }
  colnames(sce) <- bc
  # Keep full barcode list/order
  all_bc <- colnames(sce)

  # Cheap QC stats
  umis  <- Matrix::colSums(counts(sce) > 0)  # genes detected (nonzero)
  cnts  <- Matrix::colSums(counts(sce))      # total UMIs/counts

  # Define "valid" cells for scoring (tune thresholds as needed)
  valid <- (cnts >= 10) & (umis >= 5)
  if (sum(valid) < 50) {
    # still allow running, but expect instability; you can raise/lower this gate
    cat(sprintf("[WARN] Only %d/%d cells pass minimal gate; scores may be unstable.\n",
                sum(valid), length(valid)), file = stderr())
  }

  # Subset for scds computation
  sce_sub <- sce[, valid, drop = FALSE]

  cat("[INFO] Running scds (cxds_bcds_hybrid)...\n", file = stderr())
  sce_sub <- cxds_bcds_hybrid(sce_sub, estNdbl = TRUE)
  
  if (!"hybrid_score" %in% colnames(colData(sce_sub))) {
   stop("[ERROR] scds did not produce colData(sce_sub)$hybrid_score")
  }

  score_sub <- as.numeric(colData(sce_sub)$hybrid_score)
  if (all(is.na(score_sub))) {
    stop("[ERROR] hybrid_score is all NA. Input is degenerate (too few cells / zero counts / bad matrix).")
  }
  score_full <- rep(NA_real_, length(all_bc))
  score_full[valid] <- score_sub

  # Deterministic handling of NA scores: classify as singlet, keep NA score in output.
  call_full <- ifelse(
    is.na(score_full), "singlet",
    ifelse(score_full > args$threshold, "doublet", "singlet")
  )
  results <- data.frame(
    Barcode = all_bc,
    doublet = call_full,
    doublet_score = score_full,
    stringsAsFactors = FALSE
  )

  write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[INFO] Done.\n", file = stderr())
}

main()
