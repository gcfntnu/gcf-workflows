"""Sequencing-depth (saturation) curves of an evaluated biological project (#409).

Each method's count table is thinned binomially to a series of depths: every read is kept with probability
depth / the sample's sequenced reads, the same probability for every method. Depths are thinned from the deepest down,
so each is a subset of the one above. For one method's detection this is equivalent to subsampling the reads when a
read counts for one miRNA; for methods that split a multimapping read into fractions (tiny-count, mirflowz) counts are
rounded first, so it is approximate there. Only depths every sample reaches are used, and the full depth (each
sample's own) ends the series.

Each table is thinned on its own (each method's, and the isomiR output apart from the common table), while real
subsampling removes the same reads for every method: the agreement between methods is measured under independent
thinning, and typically falls faster with depth than on real shallow data.

Per method and depth: median miRNA reads per sample, median miRNAs with at least 1, 5 and 10 reads per sample,
miRNAs with at least 10 reads in every sample, the median replicate CV within the groups, and, from the isomiR
output, the median distinct read sequences per sample; per pair of methods the biological report's agreement
(Spearman, median |log2 ratio|, top-100 overlap). The depth at which detection reaches 90% and 95% of its full-depth
value is interpolated on log depth between depths every sample was measured at (the full depth only when all samples
have the same); it is relative to the dataset's own depth, so the rise per doubling of depth over the last step tells
whether the curve is still rising there.

Differential expression by depth (#414): a dataset that declares `de:` has the evaluation's DE sets (the contrast, each
level of `within`, the null comparisons of #449) tested on each depth's thinned tables, the same draw as the detection
metrics, and on the full depth, with de.R (de.smk). Each route's calls at a depth are scored against its own calls at
the full depth (recall and precision), and the null comparisons' false calls are reported per depth. A mixture's
contrasts between fractions (#451) are scored against each route's expected changes instead (#488).
"""
import collections
import math
import statistics

import numpy as np

import biological
import de

DEPTHS = '500k,1M,2M,4M,8M,16M,32M'
LEVELS = (1, 5, 10)            # miRNAs with at least this many reads in a sample
SATURATION = (0.9, 0.95)       # shares of full-depth detection whose depth is reported
SATURATION_KEY = 'detected_10'
AGREEMENT_KEYS = ('spearman_either', 'median_abs_log2_ratio', 'top{}_overlap'.format(biological.TOP))
FULL = 'full'
BEYOND = 'beyond'  # a saturation share not reached by the last depth every sample was measured at


def parse_depths(text):
    """Read depths from '500k,1M,2M' (k and M suffixes; plain numbers are reads), ascending, each above 0."""
    depths = []
    for part in str(text).split(','):
        part = part.strip()
        scale = {'k': 1e3, 'K': 1e3, 'M': 1e6}.get(part[-1:], 1)
        try:
            value = float(part[:-1] if scale != 1 else part) * scale
        except ValueError:
            raise ValueError('depths: {!r} is not a number of reads (e.g. 500k, 2M)'.format(part))
        if value <= 0 or value != int(value):
            raise ValueError('depths: {!r} must be a whole number of reads above 0'.format(part))
        depths.append(int(value))
    if len(set(depths)) != len(depths):
        raise ValueError('depths: {} lists a depth twice'.format(text))
    return sorted(depths)


def thin(rows, fraction, rng):
    """Rows (..., sample at index 1, ..., count last) with each count thinned binomially, every read kept with
    probability fraction[sample]; counts are rounded first; rows left without reads are dropped."""
    if not rows:
        return []
    counts = np.array([int(round(float(row[-1]))) for row in rows], dtype=np.int64)
    kept = rng.binomial(counts, np.array([fraction[row[1]] for row in rows]))
    return [tuple(row[:-1]) + (int(k),) for row, k in zip(rows, kept) if k > 0]


def replicate_groups(groups):
    """The groups of two samples or more with a name: those the replicate CV uses."""
    return {g: members for g, members in groups.items() if str(g).strip() and len(members) > 1}


def cv_sets(common, methods, groups):
    """{(method, group): the miRNAs with at least MIN_READS reads in every sample of the group}, for replicate_groups():
    the replicate CV's miRNAs."""
    by = biological.counts_by(common)
    return {(m, g): set.intersection(*({i for i, v in by.get((m, s), {}).items() if v >= biological.MIN_READS} for s in members))
            for m in methods for g, members in replicate_groups(groups).items()}


def metrics(common, methods, samples, groups, isomir_rows=None, fixed=None, with_isomirs=None):
    """Per method: median miRNA reads, median miRNAs with at least 1, 5 and 10 reads per sample, miRNAs with at least
    10 in every sample, the median over groups (two samples or more) of the median CPM coefficient of variation of
    the miRNAs in `fixed` (cv_sets(); by default those of these tables), and with `isomir_rows` the median distinct
    nonempty read sequences per sample ('' for a method not in `with_isomirs`, by default those with rows here)."""
    fixed = cv_sets(common, methods, groups) if fixed is None else fixed
    by = biological.counts_by(common)
    sequences = collections.defaultdict(set)
    for row in isomir_rows or []:
        if float(row[-1]) > 0 and row[4]:
            sequences[(row[0], row[1])].add(row[4])
    with_isomirs = {method for method, _ in sequences} if with_isomirs is None else with_isomirs
    out = {}
    for m in methods:
        per = {s: by.get((m, s), {}) for s in samples}
        values = {'mirna_reads': round(statistics.median(sum(c.values()) for c in per.values()))}
        for level in LEVELS:
            values['detected_{}'.format(level)] = statistics.median(sum(1 for v in c.values() if v >= level) for c in per.values())
        detected = [{i for i, v in c.items() if v >= biological.MIN_READS} for c in per.values()]
        values['detected_all_samples'] = len(set.intersection(*detected)) if detected else 0
        cvs = []
        for g, members in replicate_groups(groups).items():
            levels = [biological.cpm(per[s]) for s in members]
            values_of = lambda i: [c.get(i, 0.0) for c in levels]
            cv = [statistics.stdev(values_of(i)) / statistics.mean(values_of(i)) for i in fixed[(m, g)] if any(values_of(i))]
            if cv:
                cvs.append(statistics.median(cv))
        values['replicate_cv'] = round(statistics.median(cvs), 4) if cvs else ''
        if isomir_rows is not None:
            values['sequences'] = statistics.median(len(sequences[(m, s)]) for s in samples) if m in with_isomirs else ''
        out[m] = values
    return out


def usable_depths(depths, samples, reads_in):
    """The depths below every sample's own depth (reads_in), ascending as given."""
    return [d for d in depths if all(reads_in[s] > d for s in samples)]


def thinned(common, samples, reads_in, depths, seed, isomir_rows=None):
    """({depth or FULL: (common rows, isomir rows or None)}, [the depths used, ascending]): the tables thinned from the
    full depth down, each depth from the one above (so a subset of it), with one random generator seeded with `seed`.
    Only depths below every sample's own are used: a depth reaching it is the full depth itself. The detection metrics
    (series_points) and the differential expression (#414) of a depth read the same thinned tables."""
    usable = usable_depths(depths, samples, reads_in)
    rng = np.random.default_rng(seed)
    tables = {FULL: (common, isomir_rows)}
    above = None
    for depth in sorted(usable, reverse=True):
        fraction = {s: depth / (reads_in[s] if above is None else above) for s in samples}
        rows, iso = tables[FULL if above is None else above]
        tables[depth] = (thin(rows, fraction, rng), thin(iso, fraction, rng) if iso is not None else None)
        above = depth
    return tables, usable


def series(common, methods, samples, reads_in, groups, depths, seed, isomir_rows=None):
    """The depth series: (series_points() of the thinned() tables, the depths left out because a sample has fewer
    reads)."""
    tables, usable = thinned(common, samples, reads_in, depths, seed, isomir_rows)
    return series_points(tables, usable, methods, samples, groups, isomir_rows), [d for d in depths if d not in usable]


def series_points(tables, usable, methods, samples, groups, isomir_rows=None):
    """[(depth or FULL, metrics(), biological.agreement())] of thinned() tables, from the lowest depth up; the replicate
    CV of a method and group uses the miRNAs of the lowest depth where it has any (cv_sets()) all along."""
    order = usable + [FULL]
    fixed = {}
    for depth in order:  # from the lowest depth up, the first where a method and group have any
        for key, kept in cv_sets(tables[depth][0], methods, groups).items():
            if kept and not fixed.get(key):
                fixed[key] = kept
    fixed = {key: fixed.get(key, set()) for key in cv_sets([], methods, groups)}
    with_isomirs = {row[0] for row in isomir_rows or [] if float(row[-1]) > 0 and row[4]}
    return [(depth, metrics(tables[depth][0], methods, samples, groups, tables[depth][1], fixed, with_isomirs),
             biological.agreement(tables[depth][0], methods, samples)) for depth in order]


def rise(points, methods, full_depth, key=SATURATION_KEY):
    """{method: the rise of `key` per doubling of depth over the series' last step (the last thinned depth to the full
    depth), as a share of its full-depth value}: near 0 when the curve has levelled off; '' without a thinned depth,
    a full depth no deeper than it, or a full-depth value of 0."""
    if len(points) < 2 or full_depth <= points[-2][0]:
        return dict.fromkeys(methods, '')
    doublings = math.log2(full_depth / points[-2][0])
    out = {}
    for m in methods:
        last, full = points[-2][1][m][key], points[-1][1][m][key]
        out[m] = round((full - last) / doublings / full, 4) if full else ''
    return out


def saturation(points, methods, full_depth, key=SATURATION_KEY, shares=SATURATION, equal=True):
    """{method: {share: reads}}: the depth at which `key` first reaches `share` of its full-depth value, interpolated
    on log depth between the depths around it. The full depth (`full_depth`) is a point of the curve only when every
    sample has it (`equal`); otherwise a share not reached by the last thinned depth is BEYOND. '' when the full-depth
    value is 0, and the smallest depth when the first point already reaches it."""
    out = {}
    for m in methods:
        top = points[-1][1][m][key]
        curve = [(full_depth if depth == FULL else depth, values[m][key]) for depth, values, _ in points if depth != FULL or equal]
        out[m] = {}
        for share in shares:
            target = share * top
            reached = ''
            for i, (x, y) in enumerate(curve):
                if top and y >= target:
                    if i == 0 or y == curve[i - 1][1]:
                        reached = x
                    else:
                        (x0, y0) = curve[i - 1]
                        t = (target - y0) / (y - y0)
                        reached = 2 ** (math.log2(x0) + t * (math.log2(x) - math.log2(x0)))
                    break
            if reached == '' and top and not equal:
                reached = BEYOND
            out[m][share] = round(reached) if reached not in ('', BEYOND) else reached
    return out


DE_KEYS = ('calls', 'recall', 'precision')
WHOLE = 'all'  # the folder of a depth's DE set of the whole contrast


def de_set_name(folder):
    """The name of a DE set (#414) by depth, from its folder in benchmark/de: `all` for the whole contrast's (folder
    ''), the folder with `.` for `/` for the others (`lab.L1`, `null.case-01`)."""
    return folder.replace('/', '.') if folder else WHOLE


def de_folder(depth, folder):
    """The folder of a DE set at a depth under the depth DE root, `<depth>/<de_set_name()>`: one column and level of
    de.smk's levels.tsv."""
    return '{}/{}'.format(label(depth), de_set_name(folder))


def de_scores(results, full, features, rows):
    """{row: {<test>_<key>: value}} of a DE set's fits at a depth against the same row's fits at the full depth, over
    the set's `features` (the same at every depth); `results` and `full` as de.read_results, `full` None for the full
    depth itself (`calls` only). `calls`: the features called (BH-adjusted p of at most de.FDR); `recall`: the share of
    the full-depth calls called in the same direction; `precision`: the share of the calls that the full depth makes
    in the same direction; NaN without full-depth calls, or without calls."""
    out = {}
    for row in rows:
        out[row] = {}
        for test in de.TESTS:
            at = [de.call(results[(row, test, f)]) for f in features]
            out[row][test + '_calls'] = sum(1 for c in at if c)
            if full is not None:
                pairs = list(zip(at, (de.call(full[(row, test, f)]) for f in features)))
                out[row].update({test + '_recall': de.share([c == top for c, top in pairs if top]),
                                 test + '_precision': de.share([c == top for c, top in pairs if c])})
    return out


BETWEEN = 'between'  # the set name of a mixture's contrasts between fractions in the tables (#488)


def de_tables(contrasts, nulls, methods, between=()):
    """depth_de.tsv's rows (method, depth, set, metric, value): `contrasts` [(set, title, [(depth, de_scores())])], the
    set its folder name (WHOLE, `lab.L1`); `nulls` [(depth, de.null_scores' values)], the set `null`; `between` [(depth,
    {method: values})] of a mixture's contrasts between fractions, the set `between`."""
    rows = [[m, label(depth), name, key, value] for name, _, points in contrasts for depth, values in points
            for m in methods for key, value in values[m].items()]
    rows += [[m, label(depth), de.NULL, key, value] for depth, values in nulls for m in methods for key, value in values[m].items()]
    return rows + [[m, label(depth), BETWEEN, key, value] for depth, values in between for m in methods for key, value in values[m].items()]


def signed(value):
    return '-' if value is None or (isinstance(value, float) and math.isnan(value)) else '{:+.2f}'.format(value)


def between_lines(labelled, header, methods, md_table, where):
    """Lines on a mixture's contrasts between fractions (#488): `labelled` [(label cells, {method: values})] with
    benchmark.mixture_between_scores' values, `header` the label columns, `where` what varies (each depth, ...)."""
    lines = ['', '## Between the mixtures, against the expected changes', '',
             "A mixture's contrasts between two fractions (#451) at {}, scored as the evaluation scores them against each route's "
             'expected changes, which come from its unthinned pure samples and content ratios and are the same throughout: '
             '`called` (in brackets `called_shared`, over the miRNAs every route expects at least 2-fold apart, an untested one not '
             'called) / `wrong_direction` / `fc_bias`. `called` is a share of the miRNAs a route still tests; `scored` counts those '
             'it expects at least 2-fold apart, over the contrasts.'.format(where)]
    for test in de.TESTS:
        lines += ['', '**{}**'.format(de.TEST_NAMES[test]), '']
        lines += md_table(header + methods, [cells + ['{} ({}) / {} / {}'.format(
            share_cell(values[m][test + '_called']), share_cell(values[m][test + '_called_shared']),
            share_cell(values[m][test + '_wrong_direction']), signed(values[m][test + '_fc_bias'])) for m in methods] for cells, values in labelled])
    return lines + ['', '**scored**', ''] + md_table(header + methods, [cells + [str(values[m]['scored']) for m in methods] for cells, values in labelled])


def share_cell(value):
    return '-' if value is None or (isinstance(value, float) and math.isnan(value)) else '{:.2f}'.format(value)


def de_lines(contrasts, nulls, methods, versions, md_table, between=()):
    """depth.md's section on differential expression by depth (#414): per DE set of the contrast and test, each
    route's recall / precision (calls) per depth; per test, the null comparisons' median share called (share of splits
    with a call) per depth; a mixture's contrasts between fractions (`between`, as for de_tables) against the expected
    changes. `versions`: a line on the tools de.R recorded."""
    lines = ['', '## Differential expression by depth', '',
             "The evaluation's DE sets tested on each depth's thinned tables, the draw of the detection metrics above, and on the "
             "full depth (de.R: edgeR's exact test, its GLM test and limma-voom, #447), each over the same miRNAs at every depth: "
             "those any route reports at the full depth in the set's samples. A row tests the miRNAs with at least {} reads on "
             'average; library sizes are every miRNA read the route reported at that depth, as the workflow scales counts per '
             'million. Each route against its own calls at the full depth (BH-adjusted p of at most {}, in the same direction): '
             '`recall / precision (calls)`; `full` shows the calls. These measure how depth changes a route\'s calls, not '
             'whether the calls are right: the full-depth calls are no truth, and a route that calls more there has more to '
             'recall. They are optimistic twice over: each depth is a subset of the full-depth reads, so its calls share their '
             'noise, and the full depth is what the project holds (a read cap where one was set), not what was sequenced.'.format(
                 biological.MIN_READS, de.FDR)]
    for _, title, points in contrasts:
        for test in de.TESTS:
            lines += ['', '**{}**, {}'.format(title, de.TEST_NAMES[test]), '']
            lines += md_table(['depth'] + methods, [[label(depth)] + [
                '{} / {} ({})'.format(share_cell(values[m][test + '_recall']), share_cell(values[m][test + '_precision']), values[m][test + '_calls'])
                if test + '_recall' in values[m] else str(values[m][test + '_calls']) for m in methods] for depth, values in points])
    if nulls:
        lines += ['', '**Null comparisons** (#449): the same splits of each group at every depth, so every call is false: the median '
                  'over the splits of the share of tested miRNAs called (`null_called`), and the share of splits with a call (`null_any`) '
                  'in brackets.']
        for test in de.TESTS:
            lines += ['', '**Null comparisons**, {}'.format(de.TEST_NAMES[test]), '']
            lines += md_table(['depth'] + methods, [[label(depth)] + ['{} ({})'.format(
                share_cell(values[m][test + '_null_called']), share_cell(values[m][test + '_null_any'])) for m in methods]
                for depth, values in nulls])
    if between:
        lines += between_lines([([label(depth)], values) for depth, values in between], ['depth'], methods, md_table, 'each depth')
    return lines + ['', versions]


def label(depth):
    if depth == FULL:
        return FULL
    return '{:g}M'.format(depth / 1e6) if depth >= 1e6 else '{:g}k'.format(depth / 1e3) if depth >= 1e3 else str(depth)


def tables(points, methods):
    """(depth.tsv rows: method, depth, metric, value; depth_agreement.tsv rows: method_a, method_b, depth, metric, value)."""
    rows = [[m, label(depth), key, value] for depth, values, _ in points for m in methods for key, value in values[m].items()]
    pairs = [[a, b, label(depth), key, value[key]] for depth, _, agreement in points for (a, b), value in agreement.items()
             for key in AGREEMENT_KEYS]
    return rows, pairs


def saturation_cell(value, last):
    if value == BEYOND:
        return 'beyond {}'.format(label(last))
    return '{:,}'.format(value) if value != '' else ''


def report(name, points, methods, skipped, reached, rising, reads_in, seed, md_table):
    """depth.md's lines."""
    keys = list(next(iter(points[0][1].values())))
    lines = ['# Sequencing depth: {}'.format(name), '',
             'Each method\'s count table thinned binomially to each depth (#409): every read kept with probability depth / the '
             'sample\'s sequenced reads (seed {}), from the full depth down, each depth a subset of the one above; fractional '
             'counts (tiny-count, mirflowz) rounded first. Each table is thinned on its own (each method\'s; the isomiR output '
             'apart from the common table), whereas real subsampling removes the same reads for every method, so the agreement '
             'between methods is measured under independent thinning and typically falls faster than on real shallow data. '
             '`full`: each sample\'s own depth, {:,} to '
             '{:,} reads sequenced (median {:,}). The replicate CV uses the miRNAs with at least {} reads in every sample of a '
             'group at the lowest depth where it has any, the same set at every depth.'.format(
                 seed, min(reads_in.values()), max(reads_in.values()), round(statistics.median(reads_in.values())), biological.MIN_READS)]
    if skipped:
        lines += ['', 'Left out, above the shallowest sample: {}.'.format(', '.join(label(d) for d in skipped))]
    lines += ['', '**Depth reaching {} of full-depth detection** (miRNAs with at least {} reads, median per sample; reads '
              'sequenced per sample). It is relative to this dataset\'s depth: `rise per doubling` is the share of full-depth '
              'detection gained per doubling of depth over the last step, near 0 when the curve has levelled off. With samples of '
              'different depths the full depth is no common measurement: a share the last common depth misses is `beyond` it.'.format(
                  ' and '.join('{:.0%}'.format(s) for s in SATURATION), biological.MIN_READS), '']
    lines += md_table(['method'] + ['{:.0%}'.format(s) for s in SATURATION] + ['rise per doubling'],
                      [[m] + [saturation_cell(reached[m][s], points[-2][0] if len(points) > 1 else FULL) for s in SATURATION] +
                       ['{:.1%}'.format(rising[m]) if rising[m] != '' else ''] for m in methods])
    for key in keys:
        lines += ['', '**{}**'.format(key), ''] + md_table(['depth'] + methods, [[label(depth)] + [values[m][key] for m in methods]
                                                                                for depth, values, _ in points])
    pairs = list(points[0][2])
    if pairs:
        lines += ['', '**Agreement between methods** (Spearman over miRNAs either detects, median per sample), under independent '
                  'thinning below full depth:', '']
        lines += md_table(['pair'] + [label(depth) for depth, _, _ in points],
                          [['{} / {}'.format(a, b)] + [agreement[(a, b)]['spearman_either'] for _, _, agreement in points] for a, b in pairs])
        lines += ['', 'Median |log2 ratio| and top-{} overlap per depth are in `depth_agreement.tsv`.'.format(biological.TOP)]
    return lines


def plot(points, methods, full_depth, path):
    """depth.png: miRNAs detected (at least 10 reads, median per sample) and replicate CV against depth, log scale."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x = [full_depth if depth == FULL else depth for depth, _, _ in points]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key, title in ((axes[0], SATURATION_KEY, 'miRNAs with at least {} reads'.format(biological.MIN_READS)),
                           (axes[1], 'replicate_cv', 'replicate CV (median)')):
        for m in methods:
            ys = [values[m][key] for _, values, _ in points]
            if all(y != '' for y in ys):
                ax.plot(x, ys, marker='o', label=m)
        ax.set_xscale('log')
        ax.set_xlabel('reads sequenced per sample')
        ax.set_title(title)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
