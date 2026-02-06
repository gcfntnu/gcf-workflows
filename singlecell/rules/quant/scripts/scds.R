#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TENxIO)
  library(scds)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(Matrix)
})

clamp <- function(x, lo, hi) {
  max(lo, min(hi, x))
}

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

read_features_ids <- function(input_mtx) {
  d <- dirname(input_mtx)
  candidates <- c(
    file.path(d, "features.tsv.gz"), file.path(d, "features.tsv"),
    file.path(d, "genes.tsv.gz"),    file.path(d, "genes.tsv")
  )
  fn <- candidates[file.exists(candidates)][1]
  if (is.na(fn)) return(NULL)

  df <- import(TENxTSV(fn))
  if (!is.data.frame(df) || ncol(df) < 1) return(NULL)

  # re-read without dropping to preserve row count
  df2 <- import(TENxTSV(fn))
  fid2 <- as.character(df2[[1]])
  fid2 <- trimws(fid2)
  fid2[!nzchar(fid2) | is.na(fid2)] <- ""
  fid2
}

write_all_singlets <- function(output, barcodes) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  results <- data.frame(
    Barcode = barcodes,
    doublet = rep("singlet", length(barcodes)),
    doublet_score = rep(NA_real_, length(barcodes)),
    stringsAsFactors = FALSE
  )
  write.table(results, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
}

derive_threshold_from_rate <- function(scores, rate) {
  x <- scores[is.finite(scores)]
  if (length(x) < 2) return(Inf)
  rate <- clamp(rate, 1e-6, 1 - 1e-6)
  q <- 1 - rate
  as.numeric(stats::quantile(x, probs = q, names = FALSE, type = 7))
}

main <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-o", "--output", required = TRUE, help = "Output TSV")
  parser$add_argument("-i", "--input", required = TRUE, type = "character",
                      help = "Path to the 10x matrix.mtx(.gz)")
  parser$add_argument("-t", "--threads", type = "integer", default = NULL,
                      help = "Number of threads to use")

  # Prior for expected doublet rate (guides ROC threshold via rel_loss; fallback uses quantile)
  parser$add_argument("--expected-doublet-rate", type = "double", default = NULL,
                      help = "Expected doublet rate prior. If NULL, use 10x-style fallback: (n_cells/1000)*0.008")
  parser$add_argument("--dbr-min", type = "double", default = 0.005,
                      help = "Lower clamp for expected doublet rate")
  parser$add_argument("--dbr-max", type = "double", default = 0.40,
                      help = "Upper clamp for expected doublet rate")

  # Manual override threshold
  parser$add_argument("--threshold", type = "double", default = NULL,
                      help = "Manual score threshold override. If set, used directly; otherwise derived from expected doublet rate.")

  # Cell gating for scds computation
  parser$add_argument("--min-umis", type = "integer", default = 0,
                      help = "Minimum total counts to include a cell in scds scoring")
  parser$add_argument("--min-genes", type = "integer", default = 0,
                      help = "Minimum detected genes to include a cell in scds scoring")
  parser$add_argument("--min-valid-cells", type = "integer", default = 0,
                      help = "If fewer than this many cells pass gating, skip scds and write all singlets (NA scores)")

  # Gene filtering to prevent bcds dense coercion
  parser$add_argument("--gene-min-cells", type = "integer", default = 0,
                      help = "Keep genes detected in at least this many gated cells")
  parser$add_argument("--gene-max-frac", type = "double", default = 1,
                      help = "Drop genes detected in more than this fraction of gated cells")
  parser$add_argument("--gene-top-k", type = "integer", default = 25000,
                      help = "Cap genes to top-K by total counts after prevalence filtering")

  # Numerical sanity for variance checks
  parser$add_argument("--sd-eps", type = "double", default = 1e-8,
                      help = "Minimum SD required to treat a score vector as non-degenerate")

  args <- parser$parse_args()

  n_threads <- if (!is.null(args$threads)) args$threads else as.integer(Sys.getenv("THREADS", unset = "1"))
  if (is.na(n_threads) || !is.finite(n_threads) || n_threads < 1) stop("[ERROR] Invalid thread count")
  setup_threads(as.integer(n_threads))

  if (!file.exists(args$input)) stop("[ERROR] Input file does not exist")
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

  if (!is.finite(args$gene_max_frac) || is.na(args$gene_max_frac) || args$gene_max_frac <= 0 || args$gene_max_frac > 1) {
    stop("[ERROR] --gene-max-frac must be in (0, 1].")
  }
  if (!is.finite(args$sd_eps) || is.na(args$sd_eps) || args$sd_eps <= 0) stop("[ERROR] --sd-eps must be > 0")
  if (!is.finite(args$dbr_min) || !is.finite(args$dbr_max) ||
      args$dbr_min <= 0 || args$dbr_min >= 1 || args$dbr_max <= 0 || args$dbr_max >= 1 ||
      args$dbr_min > args$dbr_max) {
    stop("[ERROR] --dbr-min/--dbr-max must be in (0,1) and min<=max")
  }

  cat("[INFO] Reading input matrix...\n", file = stderr())
  con <- TENxMTX(args$input)
  sce <- import(con)
  sce <- as(sce, "SingleCellExperiment")

  # Barcodes
  bc <- read_barcodes(args$input)
  if (ncol(sce) != length(bc)) {
    stop(sprintf("[ERROR] Matrix has %d cells but barcodes.tsv has %d rows", ncol(sce), length(bc)))
  }
  colnames(sce) <- bc
  all_bc <- colnames(sce)

  # Features: prefer features.tsv/genes.tsv if present; else dummy IDs
  fid <- read_features_ids(args$input)
  if (!is.null(fid)) {
    if (length(fid) != nrow(sce)) {
      cat(sprintf("[WARN] features.tsv/genes.tsv has %d rows but matrix has %d genes; using dummy IDs.\n",
                  length(fid), nrow(sce)), file = stderr())
      rownames(sce) <- sprintf("g%05d", seq_len(nrow(sce)))
    } else {
      fid[fid == ""] <- sprintf("g%05d", which(fid == ""))
      rownames(sce) <- make.unique(fid, sep = "__dup")
    }
  } else {
    rownames(sce) <- sprintf("g%05d", seq_len(nrow(sce)))
  }

  # QC stats
  cts <- counts(sce)
  n_genes <- Matrix::colSums(cts > 0)
  n_umis  <- Matrix::colSums(cts)

  valid <- (n_umis >= args$min_umis) & (n_genes >= args$min_genes) & is.finite(n_umis) & is.finite(n_genes)

  cat(sprintf("[INFO] Valid cells for scds scoring: %d/%d (min_umis=%d, min_genes=%d)\n",
              sum(valid), length(valid), args$min_umis, args$min_genes), file = stderr())

  # Always maintain pipeline invariant output
  if (sum(valid) < 2) {
    cat("[WARN] Too few valid cells; writing all singlets with NA scores.\n", file = stderr())
    write_all_singlets(args$output, all_bc)
    cat("[INFO] Done.\n", file = stderr())
    return(invisible(NULL))
  }
  if (sum(valid) < args$min_valid_cells) {
    cat(sprintf("[WARN] Only %d valid cells (<%d); skipping scds for stability. Writing all singlets with NA scores.\n",
                sum(valid), args$min_valid_cells), file = stderr())
    write_all_singlets(args$output, all_bc)
    cat("[INFO] Done.\n", file = stderr())
    return(invisible(NULL))
  }

  # Subset cells for scds computation
  sce_sub <- sce[, valid, drop = FALSE]

  # Gene filtering for bcds memory safety
  cts_sub <- counts(sce_sub)
  n_cells_sub <- ncol(sce_sub)

  det <- Matrix::rowSums(cts_sub > 0)
  min_cells <- as.integer(args$gene_min_cells)
  max_cells <- as.integer(floor(args$gene_max_frac * n_cells_sub))

  keep <- det >= min_cells & det <= max_cells
  idx <- which(keep)

  if (length(idx) >= 50) {
    gene_sum <- Matrix::rowSums(cts_sub)
    idx <- idx[order(gene_sum[idx], decreasing = TRUE)]
    k <- min(as.integer(args$gene_top_k), length(idx))
    idx <- idx[seq_len(k)]
    sce_sub <- sce_sub[idx, , drop = FALSE]
  } else {
    cat(sprintf("[WARN] Only %d genes pass prevalence filter; proceeding without gene cap.\n", length(idx)),
        file = stderr())
  }

  cat(sprintf("[INFO] Running scds on %d cells and %d genes (gene_min_cells=%d, gene_max_frac=%.2f, gene_top_k=%d)\n",
              ncol(sce_sub), nrow(sce_sub), min_cells, args$gene_max_frac, as.integer(args$gene_top_k)),
      file = stderr())

  # Scores + simulated scores (needed for ROC option)
  cat("[INFO] Running scds (cxds_bcds_hybrid, estNdbl=TRUE)...\n", file = stderr())
  sce_sub <- cxds_bcds_hybrid(sce_sub, estNdbl = TRUE, verb = TRUE)

  # Pull scores
  cx <- as.numeric(colData(sce_sub)$cxds_score)
  bc_score <- as.numeric(colData(sce_sub)$bcds_score)
  hy <- as.numeric(colData(sce_sub)$hybrid_score)

  # Degeneracy handling: keep score definition consistent with calling
  cx_sd <- sd(cx, na.rm = TRUE)
  if (!is.finite(cx_sd) || is.na(cx_sd) || cx_sd <= args$sd_eps) {
    cat("[WARN] cxds_score is degenerate; using bcds_score instead of hybrid_score.\n", file = stderr())
    hy <- bc_score
  }
  hy_sd <- sd(hy, na.rm = TRUE)
  if (!is.finite(hy_sd) || is.na(hy_sd) || hy_sd <= args$sd_eps) {
    cat("[WARN] hybrid_score is degenerate; using bcds_score.\n", file = stderr())
    hy <- bc_score
  }

  if (all(is.na(hy))) {
    cat("[WARN] score vector all NA after fallbacks; writing all singlets with NA scores.\n", file = stderr())
    write_all_singlets(args$output, all_bc)
    cat("[INFO] Done.\n", file = stderr())
    return(invisible(NULL))
  }

  # Prior rate (user or 10x fallback), then clamp
  n_valid <- length(hy)
  if (is.null(args$`expected_doublet_rate`) || !is.finite(args$`expected_doublet_rate`)) {
    prior_rate <- (n_valid / 1000.0) * 0.008
    prior_src <- "fallback_10x"
  } else {
    prior_rate <- as.numeric(args$`expected_doublet_rate`)
    prior_src <- "user"
  }
  prior_rate <- clamp(prior_rate, as.numeric(args$dbr_min), as.numeric(args$dbr_max))

  # Threshold selection per your aim:
  # 1) user threshold override
  # 2) else ROC guided by prior (if sim_scores available)
  # 3) else quantile threshold from prior_rate
  thresh <- NA_real_
  thresh_src <- NA_character_
  rel_loss <- NA_real_

  if (!is.null(args$threshold) && is.finite(args$threshold)) {
    thresh <- as.numeric(args$threshold)
    thresh_src <- "user"
  } else {
    # Try ROC if we have sim_scores
    scrs_sim <- NULL
    if (!is.null(metadata(sce_sub)$hybrid) && !is.null(metadata(sce_sub)$hybrid$sim_scores)) {
      scrs_sim <- metadata(sce_sub)$hybrid$sim_scores
    }

    if (!is.null(scrs_sim) && length(scrs_sim) >= 10 &&
        sum(is.finite(scrs_sim)) >= 10 && sum(is.finite(hy)) >= 10) {

      p_sim <- length(scrs_sim) / (length(scrs_sim) + length(hy))
      p_sim <- clamp(p_sim, 1e-6, 1 - 1e-6)

      sim_odds <- p_sim / (1 - p_sim)
      prior_odds <- prior_rate / (1 - prior_rate)
      rel_loss <- sim_odds / prior_odds

      est <- tryCatch(
        scds:::get_dblCalls_ROC(hy, scrs_sim, rel_loss = rel_loss),
        error = function(e) {
          cat(sprintf("[WARN] get_dblCalls_ROC failed: %s\n", conditionMessage(e)), file = stderr())
          NULL
        }
      )

      if (!is.null(est)) {
        t <- as.numeric(est["threshold"])
        if (is.finite(t) && !is.na(t)) {
          thresh <- t
          thresh_src <- "roc"
        }
      }
    } else {
      cat("[WARN] sim_scores missing/too small; cannot compute ROC threshold.\n", file = stderr())
    }

    # ROC failed → quantile from prior
    if (!is.finite(thresh) || is.na(thresh)) {
      thresh <- derive_threshold_from_rate(hy, prior_rate)
      thresh_src <- "prior_quantile"
    }
  }

  call_rate <- mean(hy >= thresh, na.rm = TRUE)

  if (thresh_src == "roc") {
    cat(sprintf("[INFO] Using ROC threshold=%.6f (prior_rate=%.4f, prior_src=%s, rel_loss=%.4f, call_rate=%.4f)\n",
                thresh, prior_rate, prior_src, rel_loss, call_rate), file = stderr())
  } else if (thresh_src == "user") {
    cat(sprintf("[INFO] Using user threshold=%.6f (prior_rate=%.4f, prior_src=%s, call_rate=%.4f)\n",
                thresh, prior_rate, prior_src, call_rate), file = stderr())
  } else {
    cat(sprintf("[WARN] Using prior-quantile threshold=%.6f (prior_rate=%.4f, prior_src=%s, call_rate=%.4f)\n",
                thresh, prior_rate, prior_src, call_rate), file = stderr())
  }

  # Re-expand scores to full barcode set (same order as input)
  score_full <- rep(NA_real_, length(all_bc))
  score_full[valid] <- hy

  # Calls using chosen threshold (NA => unassigned)
  call_full <- ifelse(
    is.na(score_full), "unassigned",
    ifelse(score_full >= thresh, "doublet", "singlet")
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
