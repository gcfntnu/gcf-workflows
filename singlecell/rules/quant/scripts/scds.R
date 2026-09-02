#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scds)
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

  if (is.na(fn)) {
    stop("[ERROR] barcodes.tsv(.gz) not found next to input matrix")
  }

  df <- import(TENxTSV(fn))

  if (!is.data.frame(df) || ncol(df) < 1) {
    stop("[ERROR] Failed to read barcodes.tsv(.gz)")
  }

  bc <- trimws(as.character(df[[1]]))
  bc <- bc[nzchar(bc)]

  if (length(bc) == 0) {
    stop("[ERROR] No barcodes read from barcodes.tsv(.gz)")
  }

  if (anyDuplicated(bc)) {
    stop("[ERROR] Duplicate barcodes in barcodes.tsv(.gz)")
  }

  bc
}


read_feature_ids <- function(input_mtx) {
  d <- dirname(input_mtx)

  candidates <- c(
    file.path(d, "features.tsv.gz"),
    file.path(d, "features.tsv"),
    file.path(d, "genes.tsv.gz"),
    file.path(d, "genes.tsv")
  )

  fn <- candidates[file.exists(candidates)][1]

  if (is.na(fn)) {
    return(NULL)
  }

  df <- import(TENxTSV(fn))

  if (!is.data.frame(df) || ncol(df) < 1) {
    return(NULL)
  }

  fid <- trimws(as.character(df[[1]]))
  fid[!nzchar(fid) | is.na(fid)] <- ""

  fid
}


validate_rate <- function(rate, name) {
  if (!is.finite(rate) || rate <= 0 || rate >= 1) {
    stop(sprintf("[ERROR] %s must be between 0 and 1", name))
  }

  as.numeric(rate)
}


derive_10x_doublet_rate <- function(n_cells) {
  if (n_cells < 1) {
    stop("[ERROR] Cannot derive 10x doublet rate with zero cells")
  }

  min((n_cells / 1000.0) * 0.008, 0.40)
}


write_all_singlets <- function(output, barcodes) {
  results <- data.frame(
    Barcode = barcodes,
    doublet = rep("singlet", length(barcodes)),
    doublet_score = rep(NA_real_, length(barcodes)),
    stringsAsFactors = FALSE
  )

  write.table(results, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
}


main <- function() {
  parser <- ArgumentParser()

  parser$add_argument(
    "-o", "--output",
    required = TRUE,
    help = "Output TSV"
  )

  parser$add_argument(
    "-i", "--input",
    required = TRUE,
    type = "character",
    help = "Path to matrix.mtx(.gz)"
  )

  parser$add_argument(
    "-t", "--threads",
    type = "integer",
    default = NULL,
    help = "Number of threads"
  )

  parser$add_argument(
    "--expected-doublet-rate",
    type = "double",
    default = NULL,
    help = paste(
      "Explicit expected doublet rate for QC reporting.",
      "Does not override native scds classification."
    )
  )

  parser$add_argument(
    "--derive-10x-doublet-rate",
    action = "store_true",
    help = paste(
      "Derive the empirical 10x expected doublet rate from the input cell count",
      "using 0.8%% per 1000 cells, capped at 40%%.",
      "Used for QC reporting only and does not override native scds classification."
    )
  )

  parser$add_argument(
    "--threshold",
    type = "double",
    default = NULL,
    help = "Manual hybrid-score threshold; overrides native scds classification"
  )

  parser$add_argument(
    "--min-umis",
    type = "integer",
    default = 0,
    help = "Minimum total counts to include a cell in scds scoring"
  )

  parser$add_argument(
    "--min-genes",
    type = "integer",
    default = 0,
    help = "Minimum detected genes to include a cell in scds scoring"
  )

  parser$add_argument(
    "--min-valid-cells",
    type = "integer",
    default = 0,
    help = "Skip scds if fewer cells pass gating"
  )

  parser$add_argument(
    "--gene-min-cells",
    type = "integer",
    default = 0,
    help = "Keep genes detected in at least this many gated cells"
  )

  parser$add_argument(
    "--gene-max-frac",
    type = "double",
    default = 1,
    help = "Drop genes detected in more than this fraction of gated cells"
  )

  parser$add_argument(
    "--gene-top-k",
    type = "integer",
    default = 25000,
    help = "Cap genes to top-K by total counts"
  )

  args <- parser$parse_args()


  if (!is.null(args$expected_doublet_rate) && args$derive_10x_doublet_rate) {
    stop("[ERROR] --expected-doublet-rate and --derive-10x-doublet-rate are mutually exclusive")
  }

  if (!is.null(args$expected_doublet_rate)) {
    args$expected_doublet_rate <- validate_rate(
      args$expected_doublet_rate,
      "--expected-doublet-rate"
    )
  }

  if (!is.null(args$threshold) && !is.finite(args$threshold)) {
    stop("[ERROR] --threshold must be finite")
  }

  if (!is.finite(args$gene_max_frac) || args$gene_max_frac <= 0 || args$gene_max_frac > 1) {
    stop("[ERROR] --gene-max-frac must be in (0, 1]")
  }

  if (args$gene_min_cells < 0) {
    stop("[ERROR] --gene-min-cells must be >= 0")
  }

  if (args$gene_top_k < 1) {
    stop("[ERROR] --gene-top-k must be >= 1")
  }

  if (args$min_umis < 0) {
    stop("[ERROR] --min-umis must be >= 0")
  }

  if (args$min_genes < 0) {
    stop("[ERROR] --min-genes must be >= 0")
  }

  if (args$min_valid_cells < 0) {
    stop("[ERROR] --min-valid-cells must be >= 0")
  }


  n_threads <- if (!is.null(args$threads)) {
    args$threads
  } else {
    as.integer(Sys.getenv("THREADS", unset = "1"))
  }

  if (is.na(n_threads) || !is.finite(n_threads) || n_threads < 1) {
    stop("[ERROR] Invalid thread count")
  }

  setup_threads(as.integer(n_threads))


  if (!file.exists(args$input)) {
    stop("[ERROR] Input file does not exist")
  }

  sce <- as(import(TENxMTX(args$input)), "SingleCellExperiment")


  bc <- read_barcodes(args$input)

  if (ncol(sce) != length(bc)) {
    stop(sprintf(
      "[ERROR] Matrix has %d cells but barcodes.tsv has %d rows",
      ncol(sce),
      length(bc)
    ))
  }

  colnames(sce) <- bc
  all_bc <- colnames(sce)


  fid <- read_feature_ids(args$input)

  if (!is.null(fid) && length(fid) == nrow(sce)) {
    fid[fid == ""] <- sprintf("g%05d", which(fid == ""))
    rownames(sce) <- make.unique(fid, sep = "__dup")
  } else {
    rownames(sce) <- sprintf("g%05d", seq_len(nrow(sce)))
  }


  n_input_cells <- ncol(sce)
  expected_rate <- NULL
  expected_rate_source <- NULL

  if (!is.null(args$expected_doublet_rate)) {
    expected_rate <- args$expected_doublet_rate
    expected_rate_source <- "explicit"

    cat(sprintf(
      "[INFO] Explicit expected doublet rate for QC: %.4f\n",
      expected_rate
    ), file = stderr())
  }

  if (args$derive_10x_doublet_rate) {
    expected_rate <- derive_10x_doublet_rate(n_input_cells)
    expected_rate_source <- "derived_10x"

    cat(sprintf(
      "[INFO] Derived 10x expected doublet rate from %d input cells: %.4f\n",
      n_input_cells,
      expected_rate
    ), file = stderr())
  }


  cts <- counts(sce)
  n_genes <- Matrix::colSums(cts > 0)
  n_umis <- Matrix::colSums(cts)

  valid <- (
    n_umis >= args$min_umis
    & n_genes >= args$min_genes
    & is.finite(n_umis)
    & is.finite(n_genes)
  )

  n_valid <- sum(valid)

  cat(sprintf(
    "[INFO] %d/%d cells pass scds input gating\n",
    n_valid,
    n_input_cells
  ), file = stderr())


  if (n_valid < 2 || n_valid < args$min_valid_cells) {
    cat(
      "[WARN] Too few valid cells for scds; writing all cells as singlets\n",
      file = stderr()
    )

    write_all_singlets(args$output, all_bc)
    return(invisible(NULL))
  }


  sce_sub <- sce[, valid, drop = FALSE]

  cts_sub <- counts(sce_sub)
  det <- Matrix::rowSums(cts_sub > 0)
  max_cells <- as.integer(floor(args$gene_max_frac * ncol(sce_sub)))

  idx <- which(
    det >= args$gene_min_cells
    & det <= max_cells
  )


  if (length(idx) >= 50) {
    gene_sum <- Matrix::rowSums(cts_sub)
    idx <- idx[order(gene_sum[idx], decreasing = TRUE)]
    idx <- idx[seq_len(min(as.integer(args$gene_top_k), length(idx)))]

    sce_sub <- sce_sub[idx, , drop = FALSE]
  }


  cat(sprintf(
    "[INFO] Running scds on %d cells and %d genes\n",
    ncol(sce_sub),
    nrow(sce_sub)
  ), file = stderr())

  cat(
    "[INFO] Running scds cxds_bcds_hybrid(estNdbl=TRUE)...\n",
    file = stderr()
  )


  sce_sub <- cxds_bcds_hybrid(
    sce_sub,
    estNdbl = TRUE,
    verb = TRUE
  )


  required_columns <- c("hybrid_score", "hybrid_call")
  missing_columns <- setdiff(required_columns, colnames(colData(sce_sub)))

  if (length(missing_columns) > 0) {
    stop(sprintf(
      "[ERROR] scds did not return required column(s): %s",
      paste(missing_columns, collapse = ", ")
    ))
  }


  score <- as.numeric(colData(sce_sub)$hybrid_score)
  native_calls <- as.logical(colData(sce_sub)$hybrid_call)
  finite <- is.finite(score)


  if (!is.null(args$threshold)) {
    calls <- score >= args$threshold
    call_source <- "manual_threshold"

    cat(sprintf(
      "[INFO] Using manual hybrid-score threshold: %.6g\n",
      args$threshold
    ), file = stderr())
  } else {
    calls <- native_calls
    call_source <- "scds_native"
  }


  out_class <- rep("singlet", length(all_bc))
  out_score <- rep(NA_real_, length(all_bc))

  out_score[valid] <- score

  sub_class <- ifelse(calls, "doublet", "singlet")
  sub_class[!finite | is.na(calls)] <- "unassigned"

  out_class[valid] <- sub_class


  native_rate <- mean(native_calls[finite], na.rm = TRUE)
  output_rate <- mean(sub_class == "doublet")

  cat(sprintf(
    "[INFO] Native scds hybrid_call rate among finite-score cells: %.4f\n",
    native_rate
  ), file = stderr())

  cat(sprintf(
    "[INFO] scds call source: %s; output doublet rate among valid cells: %.4f\n",
    call_source,
    output_rate
  ), file = stderr())


  if (!is.null(expected_rate)) {
    cat(sprintf(
      "[INFO] QC comparison: scds native rate %.4f; %s expected rate %.4f\n",
      native_rate,
      expected_rate_source,
      expected_rate
    ), file = stderr())
  }


  results <- data.frame(
    Barcode = all_bc,
    doublet = out_class,
    doublet_score = out_score,
    stringsAsFactors = FALSE
  )

  write.table(
    results,
    file = args$output,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


main()