"""Differential expression on the ratiometric pools (#414).

Each row (the oracle and each method) is tested pool A against pool B over the scored exact sequences with edgeR's exact
test and limma-voom (de.R, which de.smk runs in the `sartools` image), and the calls are scored against the design.

Normalisation is the workflow's: a sample's library size is every miRNA read the method reported in it (all its rows of
common.tsv), as smallrna's report figures scale counts per million. There is no TMM and no centring, so an offset
between the pools stays in the fold changes, and `offset` reports it. The oracle's library size is its reads of the exact
and shifted pool sequences: what a method counting every pool miRNA read exactly would report. edgeR runs as in the
McCall-lab benchmark's best configuration (library sizes as given, a common dispersion, the exact test; two groups, so
no covariates); voom takes the covariates the dataset declares (`de: {covariates: [...]}`).

`benchmark.py evaluate` writes the inputs (input_tables) to benchmark/de/input/; de.smk writes benchmark/de/results.tsv
and record.tsv; `evaluate --de` checks the record against the current inputs and de.R, and scores the results.
"""
import hashlib
import math
import statistics
from pathlib import Path

import pool

TESTS = ('edger', 'voom')
TEST_NAMES = {'edger': "edgeR's exact test", 'voom': 'limma-voom'}
FDR = 0.05
KEYS = ('called', 'false_calls', 'fdr', 'fc_bias', 'offset')
COVARIATES = ('lab',)  # design columns of samples.tsv a declaration may name
INPUTS = ('counts.tsv', 'library_sizes.tsv', 'design.tsv')
RESULT_HEADER = ['row', 'test', 'pool_id', 'log2fc', 'pvalue', 'padj']
SCRIPT = Path(__file__).resolve().with_name('de.R')


def metric_keys():
    """The DE metrics' names in metrics.tsv: `<test>_<metric>`."""
    return ['{}_{}'.format(test, key) for test in TESTS for key in KEYS]


def designed(entry):
    """The designed log2(A / B)."""
    return math.log2(entry['amount_a'] / entry['amount_b'])


def fold(entry):
    """The designed fold change, at least 1."""
    return max(entry['amount_a'], entry['amount_b']) / min(entry['amount_a'], entry['amount_b'])


def library_sizes(reads, oracle, entries, rows, samples):
    """{(row, sample): library size}: for a method its every miRNA read in the sample (method_counts' 'all'), for the
    oracle its reads of the exact and shifted pool sequences."""
    sequences = [e['sequence'] for e in entries if e['group'] in ('exact', 'shifted')]
    return {(row, s): sum(oracle[s].get(q, 0) for q in sequences) if row == 'oracle' else reads.get((row, s), {}).get('all', 0)
            for row in rows for s in samples}


def input_tables(counts, oracle, entries, rows, design, sizes):
    """{file name: text} of the DE inputs; `design`: [(sample, group, {covariate: value})] in the order of samples.tsv.

    counts.tsv: row, pool_id, then one column per sample, over the scored exact sequences. Counts are rounded to whole
    numbers (half to even): a read a method splits between names gives fractional counts, which edgeR and voom do not
    model. library_sizes.tsv: row, then one column per sample, at full precision. design.tsv: sample, group, then the
    covariates."""
    samples = [s for s, _, _ in design]
    covariates = list(design[0][2]) if design else []
    count_lines = ['\t'.join(['row', 'pool_id'] + samples)]
    for row in rows:
        for e in pool.scored(entries):
            values = [oracle[s].get(e['sequence'], 0) if row == 'oracle' else counts.get((row, s, e['pool_id']), 0) for s in samples]
            count_lines.append('\t'.join([row, e['pool_id']] + [str(round(v)) for v in values]))
    size_lines = ['\t'.join(['row'] + samples)] + ['\t'.join([row] + [repr(sizes[(row, s)]) for s in samples]) for row in rows]
    design_lines = ['\t'.join(['sample', 'group'] + covariates)] + ['\t'.join([s, g] + [str(c[k]) for k in covariates])
                                                                   for s, g, c in design]
    return {name: '\n'.join(lines) + '\n' for name, lines in zip(INPUTS, (count_lines, size_lines, design_lines))}


def write_inputs(folder, texts):
    """Writes each input whose text changed: an unchanged file keeps its time, so de.smk does not run again."""
    folder.mkdir(parents=True, exist_ok=True)
    for name, text in texts.items():
        path = folder / name
        if not path.is_file() or path.read_text() != text:
            path.write_text(text)


def md5(data):
    return hashlib.md5(data if isinstance(data, bytes) else data.encode()).hexdigest()


def read_record(path):
    """de.R's record.tsv: {field: value}."""
    with open(path) as fh:
        lines = [line.rstrip('\n').split('\t', 1) for line in fh if line.strip()]
    return {fields[0]: fields[1] if len(fields) > 1 else '' for fields in lines[1:]}


def stale(record, texts):
    """The inputs, and `script` for de.R, whose MD5 in de.R's record differs from the current one ([] when all agree)."""
    current = {Path(name).stem: md5(text) for name, text in texts.items()}
    current['script'] = md5(SCRIPT.read_bytes())
    return [name for name, value in current.items() if record.get('md5_' + name) != value]


def read_results(path, rows, pool_ids):
    """{(row, test, pool_id): (log2fc, pvalue, padj)} from de.R's results.tsv, None for NA; ValueError unless there is
    exactly one line for every row, test and scored exact sequence."""
    table = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header != RESULT_HEADER:
            raise ValueError('{}: columns {}, expected {}'.format(path, ', '.join(header), ', '.join(RESULT_HEADER)))
        for number, line in enumerate(fh, 2):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != len(RESULT_HEADER):
                raise ValueError('{}: line {} has {} fields, expected {}'.format(path, number, len(fields), len(RESULT_HEADER)))
            key = tuple(fields[:3])
            if key in table:
                raise ValueError('{}: line {}: {} repeated'.format(path, number, ' '.join(key)))
            values = [None if v == 'NA' else float(v) for v in fields[3:]]
            table[key] = tuple(None if v is None or math.isnan(v) else v for v in values)
    expected = {(row, test, pool_id) for row in rows for test in TESTS for pool_id in pool_ids}
    missing, extra = sorted(expected - set(table)), sorted(set(table) - expected)
    if missing or extra:
        raise ValueError('{}: {} expected lines missing and {} unexpected ones (e.g. {}): run de.smk again on this '
                         'evaluation\'s inputs'.format(path, len(missing), len(extra), ' '.join((missing or extra)[0])))
    return table


def significant(fit):
    return fit[2] is not None and fit[2] <= FDR


def share(flags):
    return round(sum(flags) / len(flags), 4) if flags else float('nan')


def scores(results, exact, rows):
    """{row: {<test>_<metric>: value}} over the scored exact sequences given (a bootstrap resample repeats some), with the
    fits fixed; NaN where a metric has no sequences, never ValueError.

    `called`: share of sequences designed at least 2-fold apart that are significant (BH-adjusted p <= FDR) in the
    designed direction; an untested sequence (NA) is not called. `false_calls`: share of the 1:1 sequences called; a call
    on a sequence designed 1.5-fold apart in its direction is a true one. `fdr`: share of all calls that are on 1:1
    sequences or in the wrong direction. `fc_bias`: median of sign(d) * (log2fc - d) over the tested sequences with a
    designed change d, negative when fold changes are compressed. `offset`: median log2fc of the 1:1 sequences."""
    out = {}
    for row in rows:
        values = {}
        for test in TESTS:
            fits = [(results[(row, test, e['pool_id'])], designed(e)) for e in exact]
            right = [significant(f) and f[0] * d > 0 for f, d in fits if abs(d) >= 1]
            null = [f for f, d in fits if d == 0]
            calls = [(f, d) for f, d in fits if significant(f)]
            bias = [math.copysign(1, d) * (f[0] - d) for f, d in fits if d != 0 and f[0] is not None]
            offsets = [f[0] for f in null if f[0] is not None]
            values.update({test + '_called': share(right),
                           test + '_false_calls': share([significant(f) for f in null]),
                           test + '_fdr': share([d == 0 or not f[0] * d > 0 for f, d in calls]),
                           test + '_fc_bias': round(statistics.median(bias), 4) if bias else float('nan'),
                           test + '_offset': round(statistics.median(offsets), 4) if offsets else float('nan')})
        out[row] = values
    return out


def report_lines(results, exact, rows, record, results_sha256):
    """report.md's DE section: for each test, the share of sequences called in the designed direction per designed fold
    change (at 1:1 the share called at all) and the share tested, then the tools' versions and settings."""
    folds = sorted({round(fold(e), 2) for e in exact})
    lines = ['', '## Differential expression', '',
             'Pool A against pool B over the {} scored `exact` sequences (#414). Library sizes are all the miRNA reads a method '
             'reported in a sample, as the workflow scales counts per million (the oracle: its pool reads); no TMM, no '
             'centring. A call is a BH-adjusted p of at most {}. Per designed fold change: the share called in the designed '
             'direction, at 1:1 the share called at all; `tested`: the share of sequences with reads, so with a test.'.format(
                 len(exact), FDR)]
    for test in TESTS:
        lines += ['', '**{}**'.format(TEST_NAMES[test]), '', '| fold | n | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 2) + '|']
        for value in folds:
            group = [e for e in exact if round(fold(e), 2) == value]
            cells = []
            for row in rows:
                fits = [(results[(row, test, e['pool_id'])], designed(e)) for e in group]
                called = [significant(f) for f, _ in fits] if value == 1 else [significant(f) and f[0] * d > 0 for f, d in fits]
                cells.append('{:.2f}'.format(share(called)))
            lines.append('| {:g}× | {} | '.format(value, len(group)) + ' | '.join(cells) + ' |')
        tested = [share([results[(row, test, e['pool_id'])][1] is not None for e in exact]) for row in rows]
        lines.append('| tested | {} | '.format(len(exact)) + ' | '.join('{:.2f}'.format(t) for t in tested) + ' |')
    lines += ['', 'edgeR {}: {}. limma {}: {}, design `{}`. Both: BH within each row and test. {}. `de/results.tsv` SHA-256 '
              '`{}`.'.format(record.get('edgeR', '?'), record.get('edger', '?'), record.get('limma', '?'), record.get('voom', '?'),
                             record.get('formula', '?'), record.get('R', '?'), results_sha256)]
    return lines
