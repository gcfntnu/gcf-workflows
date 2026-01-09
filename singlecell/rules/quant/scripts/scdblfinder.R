#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(Matrix)
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
  d <- dirname(input_mtx)
  candidates <- c(file.path(d, "barcodes.tsv.gz"), file.path(d, "barcodes.tsv"))
  fn <- candidates[file.exists(candidates)][1]
  if (is.na(fn)) stop("[ERROR] barcodes.tsv(.gz) not found next to input matrix")

  df <- import(TENxTSV(fn))
  if (!is.data.frame(df) || ncol(df) < 1) stop("[ERROR] Failed to read barcodes.tsv(.gz)")
  bc <- as.character(df[[1]])
  bc <- trimws(bc)
  bc <- bc[nzchar(bc)]
  if (length(bc) == 0) stop("[ERROR] No barcodes read from barcodes.tsv(.gz)")
  if (anyDuplicated(bc)) stop("[ERROR] Duplicate barcodes in barcodes.tsv(.gz)")
  bc
}

clamp <- function(x, lo, hi) max(lo, min(hi, x))

main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-o", "--output", required = TRUE, help = "Output file")
  parser$add_argument("-i", "--input", required = TRUE, type = "character",
                      help = "Path to the 10x matrix.mtx(.gz)")
  parser$add_argument("-t", "--threads", type = "integer", default = NULL,
                      help = "Number of threads to use")
  parser$add_argument("--min-umis", type = "integer", default = 250,
                      help = "Minimum total counts to include a cell in scDblFinder fit")
  parser$add_argument("--min-genes", type = "integer", default = 0,
                      help = "Minimum detected genes to include a cell in scDblFinder fit")
  parser$add_argument("--dbr", type = "double", default = NULL,
                      help = "Doublet rate")			
  parser$add_argument("--dbr-sd", type = "double", default = NULL,
                      help = "Doublet rate standard deviation")		      
  parser$add_argument("--gene-min-cells", type="integer", default=0,
                      help="Keep genes detected in at least this many cells (after cell gating)")
  parser$add_argument("--gene-max-frac", type="double", default=1,
                      help="Drop genes detected in more than this fraction of cells (after cell gating)")
  parser$add_argument("--gene-top-k", type="integer", default=10000,
                      help="Cap genes to top-K by total counts (after prevalence filtering)")

  args <- parser$parse_args()

  n_threads <- if (!is.null(args$threads)) args$threads else as.integer(Sys.getenv("THREADS", unset = "1"))
  if (is.na(n_threads) || !is.finite(n_threads) || n_threads < 1) stop("[ERROR] Invalid thread count")
  setup_threads(as.integer(n_threads))

  if (!file.exists(args$input)) stop("[ERROR] Input file does not exist")
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

  cat("[INFO] Reading input matrix...\n", file = stderr())
  con <- TENxMTX(args$input)
  sce <- import(con)
  sce <- as(sce, "SingleCellExperiment")

  bc <- read_barcodes(args$input)
  if (ncol(sce) != length(bc)) {
    stop(sprintf("[ERROR] Matrix has %d cells but barcodes.tsv has %d rows", ncol(sce), length(bc)))
  }
  colnames(sce) <- bc
  all_bc <- colnames(sce)
  if (is.null(rownames(sce)) || length(rownames(sce)) != nrow(sce)) {
     rownames(sce) <- sprintf("g%05d", seq_len(nrow(sce)))
     }

  # Identify valid cells for fitting
  cnts  <- Matrix::colSums(counts(sce))
  genes <- Matrix::colSums(counts(sce) > 0)
  valid <- (cnts >= args$min_umis) & (genes >= args$min_genes)

  cat(sprintf("[INFO] Valid cells for fitting: %d/%d (min_umis=%d, min_genes=%d)\n",
              sum(valid), length(valid), args$min_umis, args$min_genes), file = stderr())

  # Prepare full outputs (same length/order as input)
  out_class <- rep("singlet", length(all_bc))
  out_score <- rep(NA_real_, length(all_bc))

  # If nothing valid, still write output (pipeline invariant)
  if (sum(valid) < 2) {
    cat("[WARN] Too few valid cells for scDblFinder fit; writing all singlets with NA scores.\n", file = stderr())
    results <- data.frame(Barcode = all_bc, doublet = out_class, doublet_score = out_score,
                          stringsAsFactors = FALSE)
    write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("[INFO] Done.\n", file = stderr())
    return(invisible(NULL))
  }

  sce_sub <- sce[, valid, drop = FALSE]

  cts_sub <- counts(sce_sub)
  n_cells_sub <- ncol(sce_sub)

  det <- Matrix::rowSums(cts_sub > 0)
  min_cells <- as.integer(args$gene_min_cells)
  max_cells <- as.integer(floor(args$gene_max_frac * n_cells_sub))
  
  keep <- det >= min_cells & det <= max_cells
  idx <- which(keep)

  if (length(idx) >= 50) {
    gene_sum <- Matrix::rowSums(cts_sub)
    idx <- idx[order(gene_sum[idx], decreasing=TRUE)]
    k <- min(as.integer(args$gene_top_k), length(idx))
    idx <- idx[seq_len(k)]
    sce_sub <- sce_sub[idx, , drop=FALSE]
    cat(sprintf("[INFO] scDblFinder gene filter kept %d genes\n", nrow(sce_sub)), file=stderr())
    } else {
    cat(sprintf("[WARN] Only %d genes pass prevalence filter; proceeding without gene cap.\n", length(idx)),
    			  file=stderr())
    }


  sce_sub <- scDblFinder(sce_sub, dbr = args$dbr, dbr.sd = args$dbr_sd)

  # Re-expand to full barcode set
  out_class[valid] <- as.character(sce_sub$scDblFinder.class)
  out_score[valid] <- as.numeric(sce_sub$scDblFinder.score)

  results <- data.frame(
    Barcode = all_bc,
    doublet = out_class,
    doublet_score = out_score,
    stringsAsFactors = FALSE
  )
  write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[INFO] Done.\n", file = stderr())
}

main()
