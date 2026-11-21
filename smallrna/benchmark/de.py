"""Differential expression in the benchmark (#414, #444, #447): one step for the ratiometric pools, the genome-wide
synthetic dataset and biological datasets.

Each row (a reference row, the pools' oracle or the synthetic truth, and each method) is tested one group against the
other with edgeR's exact test, edgeR's GLM likelihood-ratio test and limma-voom (de.R, which de.smk runs in the
`sartools` image). Where the changes are known (pools: the designed amounts; synthetic: the designed fold changes), the
calls are scored against them, and each method's against the reference row's fits on the same samples
(reference_keys). A DE set is one design: the whole dataset's, and others in folders of their own whose metrics carry a
suffix (the pools' `de: {within: lab}` tests each lab's runs on their own, suffix `.<level>`, #444; the synthetic
dataset adds the true nuisance scores as a covariate, suffix `.adjusted`).

Normalisation is the workflow's: a sample's library size is every miRNA read the method reported in it (all its rows of
common.tsv), as smallrna's report figures scale counts per million. There is no TMM and no centring, so an offset
between the groups stays in the fold changes, and `offset` reports it. Two normalisation arms (#453) test the pools'
and the synthetic dataset's whole design again, side by side: `fixed` (suffix `.fixed`, folder `arm/fixed`), each
sample's filtered reads as its library size for every row, so a route's miRNA calls no longer change it; and `tmm`
(suffix `.tmm`, folder `arm/tmm`), edgeR's TMM factors on the workflow's library sizes. The main metrics stay the
workflow's. The oracle's library size is its reads of the
exact and shifted pool sequences, the truth's its miRNA reads: what a method counting every miRNA read exactly would
report. The exact test runs as in the McCall-lab benchmark's best configuration (library sizes as given, a common
dispersion; two groups, so no covariates). The GLM test and voom take the set's covariates. The GLM test keeps the exact
test's model otherwise (a common dispersion, estimated under that design). Without covariates it calls nearly what the
exact test calls, so the two differ mainly by the covariates.

`benchmark.py evaluate` writes a set's inputs (input_tables) to benchmark/de/<folder>/input/ and the sets besides the
whole dataset's to benchmark/de/levels.tsv; de.smk writes results.tsv and record.tsv beside each input folder;
`evaluate --de` checks each record against the current inputs and de.R, and scores the results.
"""
import collections
import hashlib
import math
import random
import re
import statistics
from pathlib import Path

TESTS = ('edger', 'edger_glm', 'voom')  # as de.R's tests
TEST_NAMES = {'edger': "edgeR's exact test", 'edger_glm': "edgeR's GLM likelihood-ratio test", 'voom': 'limma-voom'}
FDR = 0.05
DESIGN_KEYS = ('called', 'false_calls', 'fdr', 'fc_bias', 'offset', 'false_calls_against')  # against the known changes
ONE_WAY = 2 / 3  # false_calls_against needs at least this share of the known changes in one direction (#453)
NAME = re.compile(r'[A-Za-z][A-Za-z0-9_]*')  # a covariate: a column of the design and a term of de.R's formula
LEVEL_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9._-]*')  # a group, or a level of `within`: it names a folder and a metric suffix
LEVELS = 'levels.tsv'  # benchmark/de/levels.tsv: the sets de.smk tests besides the whole dataset's, in <column>/<level>
INPUTS = ('counts.tsv', 'library_sizes.tsv', 'design.tsv', 'model.tsv')
ARMS = ('fixed', 'tmm')  # normalisation arms besides the workflow's (#453): the DE sets' `arm`, folder arm/<arm>, suffix .<arm>
RESULT_HEADER = ['row', 'test', 'feature', 'log2fc', 'pvalue', 'padj']
SCRIPT = Path(__file__).resolve().with_name('de.R')


def reference_keys(reference='oracle'):
    """The metrics of a method's fits against the reference row's (the pools' oracle, the synthetic truth)."""
    return ('own_false_calls', 'discordant_vs_' + reference, 'fc_shift_vs_' + reference, 'fc_error_vs_' + reference)


ORACLE_KEYS = reference_keys()
KEYS = DESIGN_KEYS + ORACLE_KEYS
LEVEL_KEYS = ('called', 'false_calls')  # with the reference keys, the report's table of each set besides the whole dataset's


def metric_keys(suffixes=('',), reference='oracle'):
    """The DE metrics' names in metrics.tsv: `<test>_<metric>` for the whole dataset's set (suffix ''),
    `<test>_<metric><suffix>` for another set (e.g. `.<level>`)."""
    return ['{}_{}{}'.format(test, key, suffix) for suffix in suffixes for test in TESTS for key in DESIGN_KEYS + reference_keys(reference)]


def pool_features(entries):
    """[(pool_id, designed log2(A / B))] of pool entries: the features a pool's DE set scores."""
    return [(e['pool_id'], math.log2(e['amount_a'] / e['amount_b'])) for e in entries]


def pool_count(counts, oracle, entries):
    """count(row, sample, pool_id) for input_tables: the oracle's reads of the entry's sequence, a method's count of it."""
    sequence = {e['pool_id']: e['sequence'] for e in entries}
    return lambda row, sample, feature: (oracle[sample].get(sequence[feature], 0) if row == 'oracle'
                                         else counts.get((row, sample, feature), 0))


def library_sizes(reads, oracle, entries, rows, samples):
    """{(row, sample): library size} of a pool: for a method its every miRNA read in the sample (method_counts' 'all'),
    for the oracle its reads of the exact and shifted pool sequences."""
    sequences = [e['sequence'] for e in entries if e['group'] in ('exact', 'shifted')]
    return {(row, s): sum(oracle[s].get(q, 0) for q in sequences) if row == 'oracle' else reads.get((row, s), {}).get('all', 0)
            for row in rows for s in samples}


def model_text(numerator, denominator, numeric=(), min_mean=0, normalisation=None):
    """model.tsv: the contrast (a positive log2 fold change is more in `numerator`), the covariates that are numbers
    (the others are factors), the least mean count of a feature a row tests and, for the `tmm` arm (#453),
    `normalisation: tmm` (edgeR's TMM factors on the library sizes; without the field, none)."""
    fields = [('numerator', numerator), ('denominator', denominator), ('min_mean', '{:g}'.format(min_mean))] + [('numeric', c) for c in numeric]
    if normalisation is not None:
        if normalisation != 'tmm':
            raise ValueError('normalisation must be tmm or None, got {!r}'.format(normalisation))
        fields.append(('normalisation', normalisation))
    return 'field\tvalue\n' + ''.join('{}\t{}\n'.format(*field) for field in fields)


def input_tables(count, rows, features, design, sizes, model):
    """{file name: text} of a DE set's inputs. `count(row, sample, feature)`: a row's count; `design`: [(sample, group,
    {covariate: value})] in the order of samples.tsv; `sizes`: {(row, sample): library size}; `model`: model_text.

    counts.tsv: row, feature, then one column per sample. Counts are rounded to whole numbers (half to even): a read a
    method splits between names gives fractional counts, which edgeR and voom do not model. library_sizes.tsv: row, then
    one column per sample, at full precision. design.tsv: sample, group, then the covariates."""
    samples = [s for s, _, _ in design]
    covariates = list(design[0][2]) if design else []
    count_lines = ['\t'.join(['row', 'feature'] + samples)]
    for row in rows:
        for feature in features:
            count_lines.append('\t'.join([row, feature] + [str(round(count(row, s, feature))) for s in samples]))
    size_lines = ['\t'.join(['row'] + samples)] + ['\t'.join([row] + [repr(sizes[(row, s)]) for s in samples]) for row in rows]
    design_lines = ['\t'.join(['sample', 'group'] + covariates)] + ['\t'.join([s, g] + [str(c[k]) for k in covariates])
                                                                   for s, g, c in design]
    texts = ['\n'.join(lines) + '\n' for lines in (count_lines, size_lines, design_lines)] + [model]
    return dict(zip(INPUTS, texts))


def fractional_share(count, row, features, samples):
    """The share of a row's features with reads whose count in some sample is not a whole number, so that input_tables
    rounds it (a method that divides a read between names)."""
    values = [[count(row, s, feature) for s in samples] for feature in features]
    values = [v for v in values if any(v)]
    return round(sum(any(x != round(x) for x in v) for v in values) / len(values), 4) if values else float('nan')


def write_changed(path, text):
    """Writes text unless the file holds it already: an unchanged file keeps its time, so de.smk does not run again."""
    if not path.is_file() or path.read_text() != text:
        path.write_text(text)


def write_inputs(folder, texts):
    """Writes each input whose text changed."""
    folder.mkdir(parents=True, exist_ok=True)
    for name, text in texts.items():
        write_changed(folder / name, text)


def levels_text(levels):
    """benchmark/de/levels.tsv for [(column, level)]: the sets de.smk tests besides the whole dataset's, each in
    benchmark/de/<column>/<level>/; only the header when there is none."""
    return 'column\tlevel\n' + ''.join('{}\t{}\n'.format(column, level) for column, level in levels)


def md5(data):
    return hashlib.md5(data if isinstance(data, bytes) else data.encode()).hexdigest()


def read_record(path):
    """de.R's record.tsv: {field: value}."""
    with open(path) as fh:
        lines = [line.rstrip('\n').split('\t', 1) for line in fh if line.strip()]
    return {fields[0]: fields[1] if len(fields) > 1 else '' for fields in lines[1:]}


def stale(record, texts):
    """The inputs, and `script` for de.R, whose MD5 in de.R's record differs from the current one ([] when all agree);
    the record also names the results table de.R wrote (`md5_results`), checked by the caller on the bytes it parses."""
    current = {Path(name).stem: md5(text) for name, text in texts.items()}
    current['script'] = md5(SCRIPT.read_bytes())
    return [name for name, value in current.items() if record.get('md5_' + name) != value]


def read_results(text, rows, features, path='results.tsv'):
    """{(row, test, feature): (log2fc, pvalue, padj)} from the text of de.R's results.tsv (`path` names it in messages),
    None for NA; ValueError unless there is exactly one line for every row, test and feature."""
    lines = text.splitlines()
    header = lines[0].split('\t') if lines else []
    if header != RESULT_HEADER:
        raise ValueError('{}: columns {}, expected {}'.format(path, ', '.join(header), ', '.join(RESULT_HEADER)))
    table = {}
    for number, line in enumerate(lines[1:], 2):
        fields = line.split('\t')
        if len(fields) != len(RESULT_HEADER):
            raise ValueError('{}: line {} has {} fields, expected {}'.format(path, number, len(fields), len(RESULT_HEADER)))
        key = tuple(fields[:3])
        if key in table:
            raise ValueError('{}: line {}: {} repeated'.format(path, number, ' '.join(key)))
        values = [None if v == 'NA' else float(v) for v in fields[3:]]
        table[key] = tuple(None if v is None or math.isnan(v) else v for v in values)
    expected = {(row, test, feature) for row in rows for test in TESTS for feature in features}
    missing, extra = sorted(expected - set(table)), sorted(set(table) - expected)
    if missing or extra:
        raise ValueError('{}: {} expected lines missing and {} unexpected ones (e.g. {}): run de.smk again on this '
                         'evaluation\'s inputs'.format(path, len(missing), len(extra), ' '.join((missing or extra)[0])))
    return table


def significant(fit):
    return fit[2] is not None and fit[2] <= FDR


def share(flags):
    return round(sum(flags) / len(flags), 4) if flags else float('nan')


def call(fit):
    """A fit's call: 1 up, -1 down, 0 none (not significant, or untested)."""
    return 0 if not significant(fit) else 1 if fit[0] > 0 else -1


def median(values):
    return round(statistics.median(values), 4) if values else float('nan')


def direction(features):
    """+1 or -1 when at least ONE_WAY of the known changes in `features` ([(feature, log2 fold change)]) go that way,
    else 0: the direction `false_calls_against` measures against (#453)."""
    up = [d > 0 for _, d in features if d != 0]
    if not up:
        return 0
    return 1 if sum(up) >= ONE_WAY * len(up) else -1 if len(up) - sum(up) >= ONE_WAY * len(up) else 0


def scores(results, features, rows, suffix='', reference='oracle', way=None):
    """{row: {<test>_<metric><suffix>: value}} of one set's fits over `features` ([(feature, known log2 fold change)],
    0 for a feature known not to change; a bootstrap resample repeats some), with the fits fixed; NaN where a metric has
    no features, never ValueError.

    Against the known changes: `called`, the share of features at least 2-fold apart that are significant (BH-adjusted
    p <= FDR) in the known direction; an untested feature (NA) is not called. `false_calls`: share of the unchanged
    features called; a call on a feature 1.5-fold apart in its direction is a true one. `fdr`: share of all calls that
    are on unchanged features or in the wrong direction. `fc_bias`: median of sign(d) * (log2fc - d) over the tested
    features with a known change d, negative when fold changes are compressed. `offset`: median log2fc of the unchanged
    features. `false_calls_against` (#453): share of the unchanged features called in the direction opposite to the
    known changes, when at least ONE_WAY of them go one way (NaN otherwise): the shift that library sizes made of all
    miRNA reads give the unchanged features when most changes go one way.

    Against the reference row's fits on the same samples, for each method when `rows` holds `reference`
    (reference_keys, #444): `own_false_calls`, share of the unchanged features the method calls and the reference does
    not; `discordant_vs_<reference>`, share of the features whose call (up, down or none) differs from the reference's.
    Over the features both test, the differences of the log2 fold changes (method - reference):
    `fc_shift_vs_<reference>`, their median, a shift common to all features that comes from the library sizes (a
    method's every miRNA read, the reference's own reads, in other ratios between the groups);
    `fc_error_vs_<reference>`, their median absolute deviation from that shift, the disagreement of the counts
    themselves. The shift also moves calls near the cutoff, so the two call metrics include it."""
    own, discordant, fc_shift, fc_error = reference_keys(reference)
    way = direction(features) if way is None else way  # a bootstrap passes the full design's, so every draw counts one event
    out = {}
    for row in rows:
        values = {}
        for test in TESTS:
            fits = [(results[(row, test, feature)], d) for feature, d in features]
            right = [significant(f) and f[0] * d > 0 for f, d in fits if abs(d) >= 1]
            null = [f for f, d in fits if d == 0]
            calls = [(f, d) for f, d in fits if significant(f)]
            bias = [math.copysign(1, d) * (f[0] - d) for f, d in fits if d != 0 and f[0] is not None]
            offsets = [f[0] for f in null if f[0] is not None]
            values.update({test + '_called': share(right),
                           test + '_false_calls': share([significant(f) for f in null]),
                           test + '_fdr': share([d == 0 or not f[0] * d > 0 for f, d in calls]),
                           test + '_fc_bias': median(bias),
                           test + '_offset': median(offsets),
                           # unchanged features called against the changes' direction, when most go one way (#453):
                           # the composition shift the workflow's library sizes turn into false calls
                           test + '_false_calls_against': share([significant(f) and f[0] * way < 0 for f in null]) if way else float('nan')})
            if row != reference and reference in rows:
                pairs = [(f, results[(reference, test, feature)], d) for (f, d), (feature, _) in zip(fits, features)]
                differences = [f[0] - o[0] for f, o, _ in pairs if f[0] is not None and o[0] is not None]
                shift = statistics.median(differences) if differences else 0.0
                values.update({'{}_{}'.format(test, own): share([significant(f) and not significant(o) for f, o, d in pairs if d == 0]),
                               '{}_{}'.format(test, discordant): share([call(f) != call(o) for f, o, _ in pairs]),
                               '{}_{}'.format(test, fc_shift): median(differences),
                               '{}_{}'.format(test, fc_error): median([abs(x - shift) for x in differences])})
        out[row] = {key + suffix: value for key, value in values.items()}
    return out


POOL_INTRO = (
    'Pool A against pool B over the {} scored `exact` sequences (#414). Library sizes are all the miRNA reads a method '
    'reported in a sample, as the workflow scales counts per million (the oracle: its pool reads); no TMM, no '
    'centring. A call is a BH-adjusted p of at most {}. Per designed fold change: the share called in the designed '
    'direction, at 1:1 the share called at all; `tested`: the share of sequences with reads, so with a test. Against '
    'the oracle\'s fits on the same runs (the metrics table, #444): `own_false_calls`, the 1:1 sequences a method calls '
    'and the oracle does not; `discordant_vs_oracle`, the sequences whose call (up, down or none) differs from the '
    'oracle\'s; over the sequences both test, `fc_shift_vs_oracle`, the median difference of the log2 fold changes '
    '(method - oracle), a shift that comes from the library sizes (all of a method\'s miRNA reads, the oracle\'s pool '
    'reads) and also moves calls near the cutoff; `fc_error_vs_oracle`, the median absolute deviation of those '
    'differences from the shift, the disagreement of the counts themselves.')
SYNTHETIC_INTRO = (
    'The {numerator} samples against the {denominator} samples over {features} scored miRNAs (#447): the {changed} designed '
    'to change with at least {min_mean:g} truth reads on average, and every other miRNA of the truth or reported by a '
    'method that is designed not to change. The tests are those of the ratiometric pools (de.R). Library sizes are all '
    'the miRNA reads a method reported in a sample, as the workflow scales counts per million (the truth: its miRNA '
    'reads); no TMM, no centring. Each row tests the miRNAs with at least {min_mean:g} reads on average in its own '
    'counts, so reads spilled onto a relative of a changing miRNA can make a false call. A call is a BH-adjusted p of at '
    'most {fdr}. Per designed fold change: the share called in the designed direction, at 1× the share called at all; '
    '`tested`: the share of the scored miRNAs a row tests. Against the truth\'s fits on the same samples (the metrics '
    'table): `own_false_calls`, the unchanged miRNAs a method calls and the truth does not; `discordant_vs_truth`, the '
    'miRNAs whose call (up, down or none) differs from the truth\'s; over the miRNAs both test, `fc_shift_vs_truth`, '
    'the median difference of the log2 fold changes (method - truth), a shift that comes from the library sizes; '
    '`fc_error_vs_truth`, the median absolute deviation of those differences from the shift.')


def report_lines(results, features, rows, record, results_sha256, intro):
    """report.md's DE section: `intro`, then for each test the share of features called in the known direction per known
    fold change (at 1× the share called at all) and the share tested, then the tools' versions and settings."""
    fold = lambda d: round(2 ** abs(d), 2)
    folds = sorted({fold(d) for _, d in features})
    lines = ['', '## Differential expression', '', intro]
    for test in TESTS:
        lines += ['', '**{}**'.format(TEST_NAMES[test]), '', '| fold | n | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 2) + '|']
        for value in folds:
            group = [(feature, d) for feature, d in features if fold(d) == value]
            cells = []
            for row in rows:
                fits = [(results[(row, test, feature)], d) for feature, d in group]
                called = [significant(f) for f, _ in fits] if value == 1 else [significant(f) and f[0] * d > 0 for f, d in fits]
                cells.append('{:.2f}'.format(share(called)))
            lines.append('| {:g}× | {} | '.format(value, len(group)) + ' | '.join(cells) + ' |')
        tested = [share([results[(row, test, feature)][1] is not None for feature, _ in features]) for row in rows]
        lines.append('| tested | {} | '.format(len(features)) + ' | '.join('{:.2f}'.format(t) for t in tested) + ' |')
    lines += ['', 'edgeR {}, exact test: {}; GLM test: {}, design `{}`. limma {}: {}, design `{}`. BH within each row and test. '
              '{}. `de/results.tsv` SHA-256 `{}`.'.format(record.get('edgeR', '?'), record.get('edger', '?'), record.get('edger_glm', '?'),
                                                          record.get('formula', '?'), record.get('limma', '?'), record.get('voom', '?'),
                                                          record.get('formula', '?'), record.get('R', '?'), results_sha256)]
    return lines


def text(value):
    """A metric value for a report cell: empty when missing or NaN."""
    return '' if value is None or (isinstance(value, float) and math.isnan(value)) else str(value)


def set_lines(title, intro, label, sets, values, rows, reference='oracle', keys=LEVEL_KEYS):
    """report.md's section on DE sets: for each test, `keys` (default LEVEL_KEYS) and the reference keys of every set
    and row from `values` ({row: {metric: value}}, the sets' metrics with their suffix), then each set's samples, the
    design and tool versions de.R recorded, and the SHA-256 of its results. `sets`: [(name, folder, samples,
    results_sha256, record)], `name` the metric suffix without its dot ('' for the whole dataset's set, shown as
    `workflow`), `samples` describing the set's samples; `label` heads the table's first column."""
    if not sets:
        return []
    keys = tuple(keys) + reference_keys(reference)
    lines = ['', '### ' + title, '', intro]
    for test in TESTS:
        lines += ['', '**{}**'.format(TEST_NAMES[test]), '', '| {} | metric | '.format(label) + ' | '.join(rows) + ' |',
                  '|---' * (len(rows) + 2) + '|']
        for name, _, _, _, _ in sets:
            lines += ['| {} | {} | '.format((name or 'workflow') if i == 0 else '', key) +
                      ' | '.join(text(values.get(row, {}).get('{}_{}{}'.format(test, key, '.' + name if name else ''))) for row in rows) + ' |'
                      for i, key in enumerate(keys)]
    lines += [''] + ['{}: {}, design `{}`{} (edgeR {}, limma {}); `{}` SHA-256 `{}`.'.format(
        name or 'workflow', samples, record.get('formula', '?'),
        ', TMM factors' if record.get('normalisation', '').startswith('TMM') else '',
        record.get('edgeR', '?'), record.get('limma', '?'), '/'.join(p for p in ('de', folder, 'results.tsv') if p), sha)
        for name, folder, samples, sha, record in sets]
    return lines


def arm_lines(sets, values, rows, reference='oracle'):
    """report.md's section on the normalisation arms (#453): the whole dataset's set (`workflow`) and each arm, with
    DESIGN_KEYS and the reference keys. `sets` as for set_lines, the whole dataset's first."""
    if len(sets) < 2:
        return []
    intro = ('The whole design tested again with other library sizes (#453): `workflow` is the main result above (a sample\'s '
             'library size is every miRNA read the route reported); `fixed` gives every row the same library size, the sample\'s '
             'filtered reads (the cutadapt route\'s kept reads), so what a route counts as miRNA no longer changes it; `tmm` '
             'applies edgeR\'s TMM factors to the workflow\'s library sizes (all three tests, voom on the normalised sizes). '
             'A shift in the unchanged features (`offset`) that differs between arms is the normalisation\'s, not the route\'s. '
             'In metrics.tsv the metrics carry `.fixed` and `.tmm`; their intervals and the differences between methods are in '
             'metrics_ci.tsv and differences.tsv. The metrics against the {} are empty for the {}.'.format(reference, reference))
    return set_lines('Normalisation arms', intro, 'arm', sets, values, rows, reference, DESIGN_KEYS)


def level_lines(column, levels, values, rows):
    """set_lines for the levels of a pool's `within` (#444); `levels`: [(level, runs, results_sha256, record)]."""
    if not levels:
        return []
    intro = ("The runs of each {0} tested on their own, pool A against pool B with the design below, and scored as above. In "
             "metrics.tsv the metrics carry the {0}, e.g. `edger_called.{1}`; their intervals and the differences between "
             "methods are in metrics_ci.tsv and differences.tsv, not in the uncertainty section below. The metrics against "
             "the oracle are empty for the oracle.".format(column, levels[0][0]))
    return set_lines('Within each {}'.format(column), intro, column,
                     [(level, '{}/{}'.format(column, level), runs, sha, record) for level, runs, sha, record in levels], values, rows)


def call_lines(sets, features, rows, numerator, denominator, library='all the miRNA reads a method reported in a sample'):
    """report.md's DE section of a dataset without known changes (biological samples): for each DE set and test, the
    features each row tests and those it calls up and down, then each set's design and results. Descriptive: no call
    is scored. `sets`: [(name, folder, samples, results, record, results_sha256)], the whole contrast's first (name
    ''), then each level of `within`."""
    lines = ['', '## Differential expression', '',
             '{} against {} over the {} miRNAs any method reports (identical sequences as one), with the tests of the '
             'ratiometric pools (de.R, #447). Library sizes are {}; no TMM. A '
             'row tests the miRNAs with at least the mean count below in its own counts; a call is a BH-adjusted p of at most '
             '{}. Nothing is known about the true changes, so the calls are not scored: `up` is more in {}.'.format(
                 numerator, denominator, len(features), library, FDR, numerator)]
    for name, _, samples, results, _, _ in sets:
        lines += ['', '**{}** ({})'.format(name or 'All samples of the contrast', samples), '',
                  '| test | count | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 2) + '|']
        for test in TESTS:
            fits = {row: [results[(row, test, feature)] for feature in features] for row in rows}
            for i, (label, counted) in enumerate((('tested', lambda f: f[1] is not None), ('up', lambda f: call(f) > 0),
                                                  ('down', lambda f: call(f) < 0))):
                lines.append('| {} | {} | '.format(TEST_NAMES[test] if i == 0 else '', label) +
                             ' | '.join(str(sum(1 for f in fits[row] if counted(f))) for row in rows) + ' |')
    lines += [''] + ['{}: least mean count {}; design `{}`. edgeR {}, limma {}; {}. `de/{}results.tsv` SHA-256 `{}`.'.format(
        name or 'All samples', record.get('min_mean', '?'), record.get('formula', '?'), record.get('edgeR', '?'), record.get('limma', '?'),
        record.get('R', '?'), folder + '/' if folder else '', sha) for name, folder, _, _, record, sha in sets]
    return lines


# Null comparisons (#449): the samples of one group split in two, so every call is false.
NULL = 'null'  # their sets' folders: benchmark/de/null/<group>-<k>
NULL_KEYS = ('null_called', 'null_any', 'null_p05')
NULL_MIN = 4  # samples a group needs: two in each half
NULL_TOP = 10  # the miRNAs called most often, in the report


def null_keys():
    """The null comparisons' metrics in metrics.tsv: `<test>_null_called`, `_null_any` and `_null_p05`."""
    return ['{}_{}'.format(test, key) for test in TESTS for key in NULL_KEYS]


def folder_name(group):
    """A group's name as a folder name: characters other than letters, digits, `.`, `_` and `-` become `_`."""
    name = re.sub(r'[^A-Za-z0-9._-]', '_', group)
    return name if LEVEL_NAME.fullmatch(name) else 'g' + name


def null_splits(units, splits, seed, columns=()):
    """{group: [(x samples, y samples)]} of the null comparisons (#449): within each group of `units` ([(sample, group,
    {column: value})]) that has at least NULL_MIN samples, up to `splits` distinct random halvings drawn from
    random.Random(seed). Each stratum (the samples with the same values of `columns`, the design's columns) is halved, a
    stratum's odd sample going to the halves in turn, so the halves are balanced on the design. A split and its mirror
    are one split; when fewer distinct splits exist than asked for, the drawing stops after 100 draws per split asked
    for, with all or nearly all of them. Groups and samples are taken in sorted order, so the splits depend on the seed
    and the units only."""
    rng = random.Random(seed)
    groups = collections.defaultdict(list)
    for sample, group, fields in units:
        groups[group].append((sample, fields))
    out = {}
    for group in sorted(groups):
        members = sorted(groups[group], key=lambda member: member[0])
        if len(members) < NULL_MIN:
            continue
        strata = collections.defaultdict(list)
        for sample, fields in members:
            strata[tuple(str(fields.get(c)) for c in columns)].append(sample)
        everyone, seen, found = [sample for sample, _ in members], set(), []
        for _ in range(100 * splits):
            x, odd = [], 0
            for key in sorted(strata):
                stratum = list(strata[key])
                rng.shuffle(stratum)
                x += stratum[:len(stratum) // 2 + (odd % 2 if len(stratum) % 2 else 0)]
                odd += len(stratum) % 2
            y = [sample for sample in everyone if sample not in x]
            split = tuple(sorted([tuple(sorted(x)), tuple(y)]))
            if split not in seen:
                seen.add(split)
                found.append((sorted(x), y))
                if len(found) == splits:
                    break
        out[group] = found
    return out


def null_scores(splits, rows):
    """({row: {<test>_<metric>: value}}, {row: {test: (lowest, highest) share called}}, {test: {(row, feature): splits
    called}}) of the null comparisons, `splits` [(results, features)] of every split of every group. Per row and test:
    `null_called`, the median over the splits of the share of tested features called (BH-adjusted p <= FDR), and its
    range; `null_any`, the share of splits with a call; `null_p05`, the share of tested p-values of at most 0.05, pooled
    over the splits (0.05 if they are uniform)."""
    values, spread, called = {}, {}, {test: collections.Counter() for test in TESTS}
    for row in rows:
        values[row], spread[row] = {}, {}
        for test in TESTS:
            shares, any_call, low_p, tested_total = [], [], 0, 0
            for results, features in splits:
                fits = [(feature, results[(row, test, feature)]) for feature in features]
                tested = [fit for _, fit in fits if fit[1] is not None]
                calls = [feature for feature, fit in fits if significant(fit)]
                called[test].update((row, feature) for feature in calls)
                shares.append(len(calls) / len(tested) if tested else float('nan'))
                any_call.append(bool(calls))
                low_p += sum(1 for fit in tested if fit[1] <= 0.05)
                tested_total += len(tested)
            finite = [v for v in shares if not math.isnan(v)]
            values[row].update({test + '_null_called': median(finite), test + '_null_any': share(any_call),
                                test + '_null_p05': round(low_p / tested_total, 4) if tested_total else float('nan')})
            spread[row][test] = (round(min(finite), 4), round(max(finite), 4)) if finite else None
    return values, spread, called


def null_lines(values, spread, called, rows, groups, label, mean):
    """report.md's section on the null comparisons (#449): per test, each row's metrics with the range of `null_called`
    over the splits; the splits per group; the features called in the most splits. `groups`: [(group, samples, splits)];
    `label(feature)`: its name; `mean(row, feature)`: the row's mean count of it."""
    lines = ['', '## Null comparisons', '',
             'The samples of one group split in two at random, balanced on the design, and tested against each other with the tests '
             'of the ratiometric pools (#449): every call is false. Per route and test, over the splits: `null_called`, the median '
             'share of the tested miRNAs called (BH-adjusted p of at most {}), with its range; `null_any`, the share of splits '
             'with a call; `null_p05`, the share of tested p-values of at most 0.05 (5% when they are uniform). Splits: {}.'.format(
                 FDR, '; '.join('{} ({} samples) {}'.format(group, samples, n) for group, samples, n in groups))]
    for test in TESTS:
        lines += ['', '**{}**'.format(TEST_NAMES[test]), '', '| metric | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 1) + '|',
                  '| null_called (range) | ' + ' | '.join('{} ({}–{})'.format(text(values[r][test + '_null_called']), *spread[r][test])
                                                           if spread[r][test] else '' for r in rows) + ' |']
        lines += ['| {} | '.format(key) + ' | '.join(text(values[r]['{}_{}'.format(test, key)]) for r in rows) + ' |'
                  for key in NULL_KEYS[1:]]
    totals = collections.Counter()
    for test in TESTS:
        for (_, feature), n in called[test].items():
            totals[feature] += n
    if not totals:
        return lines + ['', 'No route called any miRNA in any split.']
    lines += ['', 'The miRNAs called in the most splits: per route, the splits in which edgeR\'s exact test, its GLM test and voom '
              'call it, and the route\'s mean count.', '', '| miRNA | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 1) + '|']
    for feature, _ in sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))[:NULL_TOP]:
        lines.append('| {} | '.format(label(feature)) + ' | '.join('{} ({:.0f})'.format(
            '/'.join(str(called[test][(r, feature)]) for test in TESTS), mean(r, feature)) for r in rows) + ' |')
    return lines
