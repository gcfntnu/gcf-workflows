#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    if (!requireNamespace("orthogene", quietly = TRUE)) stop("Missing package: orthogene")
    if (!requireNamespace("schard", quietly = TRUE)) stop("Missing package: schard")
    if (!requireNamespace("argparse", quietly = TRUE)) stop("Missing package: argparse")

    library(orthogene)
    library(schard)
    library(argparse)
})

# Argument parser
parser <- ArgumentParser(description = "Convert gene identifiers using ortholog mapping")
parser$add_argument("input", help = "Path to .h5ad input file")
parser$add_argument("output", help = "Path to output TSV file")
parser$add_argument("src", help = "Source species (e.g. 'human')")
parser$add_argument("dst", help = "Destination species (e.g. 'mouse')")
parser$add_argument("--no-cache", action="store_true", default=FALSE,
                    help = "Disable caching of ortholog mapping (default: FALSE)")

args <- parser$parse_args()

# Determine default cache path
use_cache <- !args$no_cache
if (use_cache) {
    tmpdir <- Sys.getenv("TMPDIR", unset = tempdir())
    cache_dir <- file.path(tmpdir, "orthogene")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(cache_dir, paste0(args$src, "_", args$dst, "_map.rds"))
}

# Load var from h5ad
message("[1/5] Reading h5ad file: ", args$input)
var <- schard::h5ad2data.frame(args$input, 'var')

# Fallbacks for gene_id and gene_symbol
if (!"gene_id" %in% colnames(var)) {
    var$gene_id <- rownames(var)
}
if (!"gene_symbol" %in% colnames(var)) {
    var$gene_symbol <- rownames(var)
}
src_ids <- var$gene_id

# Handle mapping
if (args$src == args$dst) {
    message("[2/5] No mapping needed. Copying gene_id and gene_symbol.")
    dst_df <- var[, c("gene_id", "gene_symbol")]
} else {
    message("[2/5] Mapping genes from ", args$src, " to ", args$dst)

    # Load from cache if available
    if (use_cache && file.exists(cache_file)) {
        message("[3/5] Using cached mapping: ", cache_file)
        dst_ids <- readRDS(cache_file)
    } else {
        message("[3/5] Running ortholog conversion...")
        dst_ids <- convert_orthologs(
            gene_df = src_ids,
            input_species = args$src,
            output_species = args$dst,
            method = "gprofiler",
            mthreshold = 1
        )

        if (use_cache) {
            message("[4/5] Saving cache to: ", cache_file)
            saveRDS(dst_ids, file=cache_file)
        }
    }

    dst_df <- map_genes(rownames(dst_ids), species=args$dst, mthreshold=1, drop_na=TRUE)
    dst_df <- dst_df[, c("target", "input")]
    colnames(dst_df) <- paste(args$dst, c("gene_id", "gene_symbol"), sep="_")
    dst_df <- cbind(gene_id=dst_ids[,"input_gene"], dst_df)
}

# Write output
message("[5/5] Writing result to: ", args$output)
write.table(dst_df, file=args$output, sep="\t", quote=FALSE, row.names=FALSE)
