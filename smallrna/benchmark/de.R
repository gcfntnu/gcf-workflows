#!/usr/bin/env Rscript
# Differential expression in the small-RNA benchmark (#414, #447): for each row (a reference row such as the oracle or
# the truth, and each method), one group against the other with edgeR's exact test, edgeR's GLM likelihood-ratio test
# and limma-voom, on the counts, library sizes, design and model that `benchmark.py evaluate` writes
# (benchmark/de/input/). The model names the contrast (`numerator` against `denominator`, so a positive log2 fold
# change is more in the numerator), the covariates that are numbers (`numeric`; the others are factors) and the least
# mean count of a feature a row tests (`min_mean`), and optionally `normalisation: tmm`, `rle` or `tmmwsp` (#453, #527):
# edgeR's TMM, RLE or TMMwsp factors (calcNormFactors) on the given library sizes, used by all three tests (voom on the
# normalised library sizes; without it, the library sizes as given, no normalisation factors); and optionally
# `filter: filterByExpr` with `min_count` (#527): a row tests only the features edgeR's filterByExpr keeps on its counts,
# the given library sizes and the set's design (besides `min_mean`). The exact test runs
# as in the McCall-lab benchmark's best configuration (common dispersion, exact test; two groups, no covariates). The
# GLM test is its counterpart on the design with the covariates: a common dispersion estimated under the design and a
# likelihood-ratio test of the group. Without covariates it calls nearly what the exact test calls, so the two differ
# mainly by the covariates. voom takes the same design. A feature without reads in a row, or with a mean count under
# `min_mean` there, is written as NA for every test, and so is every feature of a row with fewer than two features to
# test: one feature gives no common dispersion and no mean-variance trend; so is every feature of a row whose factors
# are not finite (RLE without a tested feature that has reads in every sample). Each row's two common dispersions go to
# the log. Group labels, covariate values and names are read as text (a group may be called NA), and covariate names
# are quoted in the formula (one may be an R keyword such as `repeat`).
#
# The inputs and this script are copied when it starts; the fits read the copies and the record names their MD5s, so
# it describes the bytes the fits used even if `evaluate` rewrites an input meanwhile. The record also names the MD5
# of the results table written.
#
# Usage: de.R --counts TSV --library-sizes TSV --design TSV --model TSV --results TSV --record TSV
suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
})
tests <- c('edger', 'edger_glm', 'voom')  # as de.py's TESTS

options <- list(make_option('--counts'), make_option('--library-sizes', dest = 'library_sizes'), make_option('--design'),
                make_option('--model'), make_option('--results'), make_option('--record'),
                make_option('--filter-library-sizes', dest = 'filter_library_sizes'))  # optional: the sizes filterByExpr uses
opt <- parse_args(OptionParser(option_list = options))
for (name in c('counts', 'library_sizes', 'design', 'model', 'results', 'record')) {
  if (is.null(opt[[name]])) stop('--', gsub('_', '-', name), ' is required')
}
snapshot <- function(path) {
  copy <- tempfile(fileext = '.tsv')
  if (!file.copy(path, copy)) stop('cannot read ', path)
  copy
}
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)))
copies <- list(counts = snapshot(opt$counts), library_sizes = snapshot(opt$library_sizes), design = snapshot(opt$design),
               model = snapshot(opt$model), script = snapshot(script),
               filter_library_sizes = if (is.null(opt$filter_library_sizes)) NULL else snapshot(opt$filter_library_sizes))
md5 <- function(path) unname(tools::md5sum(path))
read <- function(path) read.delim(path, check.names = FALSE, quote = '', colClasses = 'character', na.strings = character())  # NA is text

design <- read(copies$design)
counts <- read(copies$counts)
sizes <- read(copies$library_sizes)
filter_sizes <- if (is.null(copies$filter_library_sizes)) sizes else read(copies$filter_library_sizes)
settings <- read(copies$model)
samples <- design$sample
if (!identical(colnames(settings), c('field', 'value'))) stop(opt$model, ': columns must be field and value')
setting <- function(field) settings$value[settings$field == field]
numerator <- setting('numerator')
denominator <- setting('denominator')
if (length(numerator) != 1 || length(denominator) != 1 || numerator == denominator) stop(opt$model, ': needs one numerator and one other denominator')
min_mean <- if (length(setting('min_mean'))) suppressWarnings(as.numeric(setting('min_mean'))) else 0
if (length(min_mean) != 1 || is.na(min_mean) || min_mean < 0) stop(opt$model, ': min_mean must be one number of 0 or more')
numeric <- setting('numeric')
normalisation <- setting('normalisation')
methods <- c(tmm = 'TMM', rle = 'RLE', tmmwsp = 'TMMwsp')  # de.py's NORMALISATIONS: edgeR's calcNormFactors methods
if (length(normalisation) > 1 || (length(normalisation) == 1 && !normalisation %in% names(methods)))
  stop(opt$model, ': normalisation must be tmm, rle or tmmwsp, or absent')
method <- if (length(normalisation) == 1) methods[[normalisation]] else NA
filter <- setting('filter')  # #527: filterByExpr on the row's counts, library sizes and the set's design, with min_count
if (length(filter) > 1 || (length(filter) == 1 && filter != 'filterByExpr')) stop(opt$model, ': filter must be filterByExpr, or absent')
min_count <- if (length(setting('min_count'))) suppressWarnings(as.numeric(setting('min_count'))) else NA
if (length(filter) == 1 && (length(min_count) != 1 || !is.finite(min_count) || min_count <= 0)) stop(opt$model, ': filterByExpr needs one finite min_count above 0')
if (length(filter) == 0 && length(setting('min_count'))) stop(opt$model, ': min_count goes with filter: filterByExpr')
unknown <- setdiff(settings$field, c('numerator', 'denominator', 'min_mean', 'numeric', 'normalisation', 'filter', 'min_count'))
if (length(unknown)) stop(opt$model, ': unknown field ', paste(unknown, collapse = ', '))
if (!identical(colnames(design)[1:2], c('sample', 'group'))) stop(opt$design, ': the first columns must be sample and group')
if (any(!nzchar(design$sample)) || any(!nzchar(design$group))) stop(opt$design, ': every sample needs a name and a group')
if (!identical(sort(unique(design$group)), sort(c(numerator, denominator)))) stop(opt$design, ': groups must be ', numerator, ' and ', denominator)
if (!identical(colnames(counts), c('row', 'feature', samples))) stop(opt$counts, ': columns must be row, feature and the design samples in order')
if (!identical(colnames(sizes), c('row', samples))) stop(opt$library_sizes, ': columns must be row and the design samples in order')
if (!identical(colnames(filter_sizes), c('row', samples))) stop(opt$filter_library_sizes, ': columns must be row and the design samples in order')
if (!is.null(copies$filter_library_sizes) && length(filter) == 0) stop('--filter-library-sizes goes with filter: filterByExpr')

covariates <- setdiff(colnames(design), c('sample', 'group'))
if (length(setdiff(numeric, covariates))) stop(opt$model, ': numeric names ', paste(setdiff(numeric, covariates), collapse = ', '), ', not a covariate of the design')
frame <- data.frame(group = factor(design$group, levels = c(denominator, numerator)), row.names = samples)
for (covariate in covariates) {
  if (covariate %in% numeric) {
    frame[[covariate]] <- suppressWarnings(as.numeric(design[[covariate]]))
    if (anyNA(frame[[covariate]])) stop(opt$design, ': ', covariate, ' is numeric in the model but has values that are not numbers')
  } else {
    frame[[covariate]] <- factor(design[[covariate]])
  }
}
formula <- as.formula(paste('~', paste(sprintf('`%s`', c(covariates, 'group')), collapse = ' + ')))  # quoted: a name may be a keyword
formula_text <- paste(deparse(formula, width.cutoff = 500L), collapse = ' ')  # one string, however long the formula
model <- model.matrix(formula, data = frame)
coefficient <- ncol(model)  # the group is the formula's last term, with two levels: its one column, group<numerator>
if (qr(model)$rank < ncol(model)) stop('the design ', formula_text, ' is not full rank: a covariate is confounded with the groups')
if (nrow(model) - ncol(model) < 1) stop('the design ', formula_text, ' leaves no residual degrees of freedom for ', nrow(model), ' samples')

results <- list()
add <- function(row, test, table) results[[length(results) + 1]] <<- data.frame(row = row, test = test, table, check.names = FALSE)
for (row in unique(counts$row)) {
  part <- counts[counts$row == row, , drop = FALSE]
  if (anyDuplicated(part$feature)) stop(opt$counts, ': ', row, ' lists a feature twice')
  values <- suppressWarnings(matrix(as.numeric(as.matrix(part[samples])), ncol = length(samples), dimnames = list(part$feature, samples)))
  if (anyNA(values) || any(values < 0) || any(values != round(values))) stop(opt$counts, ': ', row, ' has counts that are not whole numbers of 0 or more')
  size <- suppressWarnings(as.numeric(unlist(sizes[sizes$row == row, samples])))
  if (length(size) != length(samples) || anyNA(size) || any(size <= 0)) stop(opt$library_sizes, ': ', row, ' needs one positive library size per sample')
  keep <- rowSums(values) > 0 & rowMeans(values) >= min_mean
  if (length(filter) == 1) {  # on the counts as given, before any normalisation factors (#527)
    # on the filter's library sizes: the fit's, or the workflow's for an arm with other sizes, so every arm tests the same features
    filter_size <- suppressWarnings(as.numeric(unlist(filter_sizes[filter_sizes$row == row, samples])))
    if (length(filter_size) != length(samples) || anyNA(filter_size) || any(filter_size <= 0)) stop('filter library sizes: ', row, ' needs one positive size per sample')
    keep <- keep & filterByExpr(DGEList(counts = values, lib.size = filter_size), design = model, min.count = min_count)
  }
  if (sum(keep) == 1) {  # one feature: no common dispersion and no mean-variance trend (voom needs two), so no test
    message(sprintf('%s: 1 feature to test (%s), at least 2 are needed; none tested', row, part$feature[keep]))
    keep[] <- FALSE
  }
  none <- rep(NA_real_, sum(!keep))
  untested <- data.frame(feature = part$feature[!keep], log2fc = none, pvalue = none, padj = none)
  if (!any(keep)) {  # no feature to test in this row: nothing to fit
    for (test in tests) add(row, test, untested)
    next
  }
  tested <- values[keep, , drop = FALSE]

  y <- DGEList(counts = tested, lib.size = size, group = frame$group)
  if (!is.na(method)) {
    y <- calcNormFactors(y, method = method)
    if (!all(is.finite(y$samples$norm.factors) & y$samples$norm.factors > 0)) {
      # RLE needs a feature with reads in every sample (#528 review): without one the row has no factors, so no test
      message(sprintf('%s: no finite %s factors (no tested feature with reads in every sample); none tested', row, method))
      all_na <- rep(NA_real_, nrow(part))
      for (test in tests) add(row, test, data.frame(feature = part$feature, log2fc = all_na, pvalue = all_na, padj = all_na))
      next
    }
  }
  y <- estimateCommonDisp(y)
  exact <- exactTest(y, pair = c(denominator, numerator), dispersion = 'common')$table
  edger <- data.frame(feature = rownames(exact), log2fc = exact$logFC, pvalue = exact$PValue, padj = p.adjust(exact$PValue, 'BH'))

  common <- estimateGLMCommonDisp(y, model)$common.dispersion  # under the design, so without the covariates' variance
  lrt <- glmLRT(glmFit(y, model, dispersion = common), coef = coefficient)$table
  edger_glm <- data.frame(feature = rownames(lrt), log2fc = lrt$logFC, pvalue = lrt$PValue, padj = p.adjust(lrt$PValue, 'BH'))
  message(sprintf('%s: %d features tested; common dispersion %.4g for the exact test, %.4g for the GLM test', row, nrow(tested),
                  y$common.dispersion, common))

  fit <- eBayes(lmFit(voom(tested, model, lib.size = size * y$samples$norm.factors), model))  # factors are 1 without a method
  top <- topTable(fit, coef = coefficient, number = Inf, sort.by = 'none', adjust.method = 'BH')
  voomed <- data.frame(feature = rownames(top), log2fc = top$logFC, pvalue = top$P.Value, padj = top$adj.P.Val)

  fits <- list(edger = edger, edger_glm = edger_glm, voom = voomed)
  for (test in tests) {
    table <- rbind(fits[[test]], untested)
    add(row, test, table[match(part$feature, table$feature), ])
  }
}
table <- if (length(results)) do.call(rbind, results) else  # no feature in any row: the header alone
  data.frame(row = character(), test = character(), feature = character(), log2fc = numeric(), pvalue = numeric(), padj = numeric())
write.table(table, opt$results, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')

record <- data.frame(
  field = c('md5_counts', 'md5_library_sizes', 'md5_design', 'md5_model', 'md5_script', 'md5_results', 'R', 'edgeR', 'limma', 'formula',
            'contrast', 'min_mean', 'filter', 'md5_filter_library_sizes', 'normalisation', 'edger', 'edger_glm', 'voom'),
  value = c(md5(copies$counts), md5(copies$library_sizes), md5(copies$design), md5(copies$model), md5(copies$script), md5(opt$results),
            R.version.string, as.character(packageVersion('edgeR')), as.character(packageVersion('limma')), formula_text,
            paste(numerator, 'against', denominator), as.character(min_mean),
            if (length(filter) == 1) sprintf('filterByExpr(design = model, min.count = %g) on the %s', min_count,
                                             if (is.null(copies$filter_library_sizes)) 'library sizes' else 'filter library sizes') else 'none',
            if (is.null(copies$filter_library_sizes)) 'none' else md5(copies$filter_library_sizes),
            if (!is.na(method)) sprintf('%s factors (calcNormFactors) on the library sizes', method) else 'none: the library sizes as given',
            sprintf('DGEList with the library sizes and %s, estimateCommonDisp, exactTest of %s against %s',
                    if (!is.na(method)) sprintf('%s normalisation factors', method) else 'no normalisation factors', numerator, denominator),
            sprintf('the same DGEList, estimateGLMCommonDisp under the design, glmFit with that common dispersion, glmLRT of coefficient %s',
                    colnames(model)[coefficient]),
            sprintf('voom with the library sizes%s, lmFit, eBayes, coefficient %s', if (!is.na(method)) sprintf(' times the %s factors', method) else '',
                    colnames(model)[coefficient])))
write.table(record, opt$record, sep = '\t', quote = FALSE, row.names = FALSE)
