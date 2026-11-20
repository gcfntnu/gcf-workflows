#!/usr/bin/env Rscript
# Differential expression on the small-RNA benchmark's ratiometric pools (#414): for each row (the oracle and each
# method), pool A against pool B with edgeR's exact test and limma-voom, on the counts, library sizes and design that
# `benchmark.py evaluate` writes (benchmark/de/input/). Library sizes are the workflow's, without TMM. edgeR runs as in
# the McCall-lab benchmark's best configuration (common dispersion, exact test; two groups, no covariates); voom takes
# the design's covariates. A sequence without reads in a row is written as NA for both tests.
#
# Usage: de.R --counts TSV --library-sizes TSV --design TSV --results TSV --record TSV
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
})

options <- list(make_option('--counts'), make_option('--library-sizes', dest = 'library_sizes'), make_option('--design'),
                make_option('--results'), make_option('--record'))
opt <- parse_args(OptionParser(option_list = options))
for (name in c('counts', 'library_sizes', 'design', 'results', 'record')) {
  if (is.null(opt[[name]])) stop('--', gsub('_', '-', name), ' is required')
}
read <- function(path) read.delim(path, check.names = FALSE, quote = '', colClasses = 'character')

design <- read(opt$design)
counts <- read(opt$counts)
sizes <- read(opt$library_sizes)
samples <- design$sample
if (!identical(colnames(design)[1:2], c('sample', 'group'))) stop(opt$design, ': the first columns must be sample and group')
if (!identical(sort(unique(design$group)), c('A', 'B'))) stop(opt$design, ': groups must be A and B')
if (!identical(colnames(counts), c('row', 'pool_id', samples))) stop(opt$counts, ': columns must be row, pool_id and the design samples in order')
if (!identical(colnames(sizes), c('row', samples))) stop(opt$library_sizes, ': columns must be row and the design samples in order')

covariates <- setdiff(colnames(design), c('sample', 'group'))
frame <- data.frame(group = factor(design$group, levels = c('B', 'A')), row.names = samples)
for (covariate in covariates) frame[[covariate]] <- factor(design[[covariate]])
formula <- as.formula(paste('~', paste(c(covariates, 'group'), collapse = ' + ')))
model <- model.matrix(formula, data = frame)
if (qr(model)$rank < ncol(model)) stop('the design ', deparse(formula), ' is not full rank: a covariate is confounded with the pools')

results <- list()
for (row in unique(counts$row)) {
  part <- counts[counts$row == row, , drop = FALSE]
  if (anyDuplicated(part$pool_id)) stop(opt$counts, ': ', row, ' lists a sequence twice')
  values <- suppressWarnings(matrix(as.numeric(as.matrix(part[samples])), ncol = length(samples), dimnames = list(part$pool_id, samples)))
  if (anyNA(values) || any(values < 0) || any(values != round(values))) stop(opt$counts, ': ', row, ' has counts that are not whole numbers of 0 or more')
  size <- suppressWarnings(as.numeric(unlist(sizes[sizes$row == row, samples])))
  if (length(size) != length(samples) || anyNA(size) || any(size <= 0)) stop(opt$library_sizes, ': ', row, ' needs one positive library size per sample')
  keep <- rowSums(values) > 0
  tested <- values[keep, , drop = FALSE]

  y <- DGEList(counts = tested, lib.size = size, group = frame$group)
  y <- estimateCommonDisp(y)
  exact <- exactTest(y, pair = c('B', 'A'), dispersion = 'common')$table
  edger <- data.frame(pool_id = rownames(exact), log2fc = exact$logFC, pvalue = exact$PValue, padj = p.adjust(exact$PValue, 'BH'))

  fit <- eBayes(lmFit(voom(tested, model, lib.size = size), model))
  top <- topTable(fit, coef = 'groupA', number = Inf, sort.by = 'none', adjust.method = 'BH')
  voomed <- data.frame(pool_id = rownames(top), log2fc = top$logFC, pvalue = top$P.Value, padj = top$adj.P.Val)

  none <- rep(NA_real_, sum(!keep))
  untested <- data.frame(pool_id = part$pool_id[!keep], log2fc = none, pvalue = none, padj = none)
  for (test in c('edger', 'voom')) {
    table <- rbind(if (test == 'edger') edger else voomed, untested)
    table <- table[match(part$pool_id, table$pool_id), ]
    results[[length(results) + 1]] <- data.frame(row = row, test = test, table, check.names = FALSE)
  }
}
write.table(do.call(rbind, results), opt$results, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')

script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)))
md5 <- function(path) unname(tools::md5sum(path))
record <- data.frame(
  field = c('md5_counts', 'md5_library_sizes', 'md5_design', 'md5_script', 'R', 'edgeR', 'limma', 'formula', 'edger', 'voom'),
  value = c(md5(opt$counts), md5(opt$library_sizes), md5(opt$design), md5(script), R.version.string,
            as.character(packageVersion('edgeR')), as.character(packageVersion('limma')), deparse(formula),
            'DGEList with the library sizes and no normalisation factors, estimateCommonDisp, exactTest of A against B',
            'voom with the library sizes, lmFit, eBayes, coefficient groupA'))
write.table(record, opt$record, sep = '\t', quote = FALSE, row.names = FALSE)
