#!/usr/bin/env Rscript
# Title: Convert gene identifiers between species using ortholog mapping
# Description:
#   Reads gene metadata from an AnnData (.h5ad) file or a
#   Matrix Market (.mtx) file (with accompanying features.tsv or features.tsv.gz),
#   converts identifiers from a source species to orthologs in
#   a destination species using orthogene, and writes results
#   to a TSV file. Optionally caches mappings for reuse.
#
# Usage:
#   run_orthogene.R --input INPUT --output OUTPUT --src SRC --dst DST [--no-cache]
#
# Dependencies:
#   orthogene, schard, argparse, Matrix

suppressPackageStartupMessages({
    if (!requireNamespace("orthogene", quietly = TRUE)) stop("Missing package: orthogene")
    if (!requireNamespace("schard", quietly = TRUE)) stop("Missing package: schard")
    if (!requireNamespace("argparse", quietly = TRUE)) stop("Missing package: argparse")
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix")
})

library(orthogene)
library(schard)
library(argparse)
library(Matrix)

# Argument parser
tparser <- ArgumentParser(description = "Convert gene identifiers using ortholog mapping")
tparser$add_argument("--input",   dest = "input",   required = TRUE,
                     help = "Path to .h5ad or .mtx input file")
tparser$add_argument("--output",  dest = "output",  required = TRUE,
                     help = "Path to output TSV file")
tparser$add_argument("--src",     dest = "src",     required = TRUE,
                     help = "Source species (e.g. 'human')")
tparser$add_argument("--dst",     dest = "dst",     required = TRUE,
                     help = "Destination species (e.g. 'mouse')")
tparser$add_argument("--no-cache", dest = "no_cache", action = "store_true", default = FALSE,
                     help = "Disable caching of ortholog mapping (default: FALSE)")
args <- tparser$parse_args()

# Determine default cache path
use_cache <- !args$no_cache
if (use_cache) {
    tmpdir <- Sys.getenv("TMPDIR", unset = tempdir())
    cache_dir <- file.path(tmpdir, "orthogene")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(cache_dir, paste0(args$src, "_", args$dst, "_map.rds"))
}

# Load gene metadata
message("[1/6] Reading input file: ", args$input)
if (grepl("\\.mtx$", args$input, ignore.case = TRUE)) {
    dir <- dirname(args$input)
    legacy_tsv_path =  file.path(dir, "genes.tsv")
    tsv_path <- file.path(dir, "features.tsv")
    gz_path  <- file.path(dir, "features.tsv.gz")
    if (file.exists(gz_path)) {
        feat_file <- gz_path
    } else if (file.exists(tsv_path)) {
        feat_file <- tsv_path
    } else if (file.exists(legacy_tsv_path)) {
        feat_file <- legacy_tsv_path
    } else {
        stop("No features.tsv or features.tsv.gz (or genes.tsv) found in ", dir)
    }
    message("[2/6] Loading gene metadata from: ", feat_file)
    features <- read.delim(feat_file, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(features) >= 2) {
        var <- data.frame(
            gene_id     = features[[1]],
            gene_symbol = features[[2]],
            stringsAsFactors = FALSE
        )
    } else {
        var <- data.frame(
            gene_id     = features[[1]],
            gene_symbol = features[[1]],
            stringsAsFactors = FALSE
        )
    }
    rownames(var) <- var$gene_id
} else {
    message("[2/6] Loading gene metadata from AnnData (.h5ad)")
    var <- schard::h5ad2data.frame(args$input, 'var')
    # Fallbacks for gene_id and gene_symbol
    if (!"gene_id" %in% colnames(var)) {
        var$gene_id <- rownames(var)
    }
    if (!"gene_symbol" %in% colnames(var)) {
        var$gene_symbol <- rownames(var)
    }
}
src_ids <- var$gene_id

# Handle mapping
if (args$src == args$dst) {
    message("[3/6] No mapping needed. Copying gene_id and gene_symbol.")
    dst_df <- var[, c("gene_id", "gene_symbol")]
} else {
    message("[3/6] Mapping genes from ", args$src, " to ", args$dst)

    # Load from cache or run conversion
    if (use_cache && file.exists(cache_file)) {
        message("[4/6] Using cached mapping: ", cache_file)
        dst_ids <- readRDS(cache_file)
    } else {
        message("[4/6] Running ortholog conversion...")
        dst_ids <- convert_orthologs(
            gene_df        = src_ids,
            input_species  = args$src,
            output_species = args$dst,
            method         = "gprofiler",
            mthreshold     = 1
        )
        if (use_cache) {
            message("[5/6] Saving cache to: ", cache_file)
            saveRDS(dst_ids, file = cache_file)
        }
    }
    
  message("[4.5/6] Map to dst organism...")

  # Map to obtain destination symbols (KEEP NAs, then align back to dst_ids)
  # dst_map: KEEP NAs
dst_map <- map_genes(
  genes      = rownames(dst_ids),
  species    = args$dst,
  mthreshold = 1,
  drop_na    = FALSE
)

sym <- dst_map$target
names(sym) <- dst_map$input
aligned_sym <- unname(sym[rownames(dst_ids)])

dst_df <- data.frame(
  gene_id      = dst_ids[, "input_gene", drop = TRUE],
  dst_gene_id  = rownames(dst_ids),
  dst_symbol   = aligned_sym,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(dst_df)[names(dst_df) == "dst_gene_id"] <- paste0(args$dst, "_gene_id")
names(dst_df)[names(dst_df) == "dst_symbol"]  <- paste0(args$dst, "_gene_symbol")

}

# Write output
message("[6/6] Writing result to: ", args$output)
write.table(
    dst_df,
    file      = args$output,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)
