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
  bc <- trimws(as.character(df[[1]]))
  bc <- bc[nzchar(bc)]
  if (length(bc) == 0) stop("[ERROR] No barcodes read from barcodes.tsv(.gz)")
  if (anyDuplicated(bc)) stop("[ERROR] Duplicate barcodes in barcodes.tsv(.gz)")
  bc
}

resolve_expected_doublet_rate <- function(explicit_rate, derive_10x, n_cells) {
  if (!is.null(explicit_rate)) {
    if (!is.finite(explicit_rate) || explicit_rate <= 0 || explicit_rate >= 1) {
      stop("[ERROR] --expected-doublet-rate must be between 0 and 1")
    }
    return(as.numeric(explicit_rate))
  }

  if (derive_10x) {
    if (n_cells < 1) stop("[ERROR] Cannot derive 10x doublet rate with zero cells")
    return(min((n_cells / 1000.0) * 0.008, 0.40))
  }

  NULL
}

main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-o", "--output", required = TRUE, help = "Output TSV")
  parser$add_argument("-i", "--input", required = TRUE, type = "character",
                      help = "Path to matrix.mtx(.gz)")
  parser$add_argument("-t", "--threads", type = "integer", default = NULL,
                      help = "Number of threads")
  parser$add_argument("--expected-doublet-rate", type = "double", default = NULL,
                      help = "Explicit biological expected doublet rate")
  parser$add_argument("--derive-10x-doublet-rate", action = "store_true",
                      help = "Derive expected rate from input cell count using 0.8%% per 1000 cells")
  parser$add_argument("--min-umis", type = "integer", default = 250,
                      help = "Minimum total counts to include a cell in scDblFinder fit")
  parser$add_argument("--min-genes", type = "integer", default = 0,
                      help = "Minimum detected genes to include a cell in scDblFinder fit")
  parser$add_argument("--gene-min-cells", type = "integer", default = 0,
                      help = "Keep genes detected in at least this many cells")
  parser$add_argument("--gene-max-frac", type = "double", default = 1,
                      help = "Drop genes detected in more than this fraction of cells")
  parser$add_argument("--gene-top-k", type = "integer", default = 10000,
                      help = "Cap genes to top-K by total counts")
  args <- parser$parse_args()

  if (!is.null(args$expected_doublet_rate) && args$derive_10x_doublet_rate) {
    stop("[ERROR] --expected-doublet-rate and --derive-10x-doublet-rate are mutually exclusive")
  }

  n_threads <- if (!is.null(args$threads)) args$threads else as.integer(Sys.getenv("THREADS", unset = "1"))
  if (is.na(n_threads) || !is.finite(n_threads) || n_threads < 1) stop("[ERROR] Invalid thread count")
  if (!is.finite(args$gene_max_frac) || args$gene_max_frac <= 0 || args$gene_max_frac > 1) {
    stop("[ERROR] --gene-max-frac must be in (0, 1]")
  }
  setup_threads(as.integer(n_threads))

  if (!file.exists(args$input)) stop("[ERROR] Input file does not exist")

  cat("[INFO] Reading input matrix...\n", file = stderr())
  sce <- as(import(TENxMTX(args$input)), "SingleCellExperiment")

  bc <- read_barcodes(args$input)
  if (ncol(sce) != length(bc)) {
    stop(sprintf("[ERROR] Matrix has %d cells but barcodes.tsv has %d rows", ncol(sce), length(bc)))
  }
  colnames(sce) <- bc
  all_bc <- colnames(sce)
  if (is.null(rownames(sce)) || length(rownames(sce)) != nrow(sce)) {
    rownames(sce) <- sprintf("g%05d", seq_len(nrow(sce)))
  }

  expected_rate <- resolve_expected_doublet_rate(
    args$expected_doublet_rate,
    args$derive_10x_doublet_rate,
    ncol(sce)
  )

  cnts <- Matrix::colSums(counts(sce))
  genes <- Matrix::colSums(counts(sce) > 0)
  valid <- (cnts >= args$min_umis) & (genes >= args$min_genes)

  out_class <- rep("singlet", length(all_bc))
  out_score <- rep(NA_real_, length(all_bc))

  if (sum(valid) < 2) {
    results <- data.frame(Barcode = all_bc, doublet = out_class, doublet_score = out_score,
                          stringsAsFactors = FALSE)
    write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
    return(invisible(NULL))
  }

  sce_sub <- sce[, valid, drop = FALSE]
  cts_sub <- counts(sce_sub)
  det <- Matrix::rowSums(cts_sub > 0)
  max_cells <- as.integer(floor(args$gene_max_frac * ncol(sce_sub)))
  idx <- which(det >= args$gene_min_cells & det <= max_cells)

  if (length(idx) >= 50) {
    gene_sum <- Matrix::rowSums(cts_sub)
    idx <- idx[order(gene_sum[idx], decreasing = TRUE)]
    idx <- idx[seq_len(min(as.integer(args$gene_top_k), length(idx)))]
    sce_sub <- sce_sub[idx, , drop = FALSE]
  }

  if (is.null(expected_rate)) {
    # Rate-free mode: do not constrain threshold by an expected proportion.
    sce_sub <- scDblFinder(sce_sub, dbr = NULL, dbr.sd = 1)
  } else {
    sce_sub <- scDblFinder(sce_sub, dbr = expected_rate)
  }

  out_class[valid] <- as.character(sce_sub$scDblFinder.class)
  out_score[valid] <- as.numeric(sce_sub$scDblFinder.score)

  results <- data.frame(
    Barcode = all_bc,
    doublet = out_class,
    doublet_score = out_score,
    stringsAsFactors = FALSE
  )
  write.table(results, file = args$output, sep = "\t", row.names = FALSE, quote = FALSE)
}

main()
