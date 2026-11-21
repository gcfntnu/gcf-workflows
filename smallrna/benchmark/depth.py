"""Sequencing-depth (saturation) curves of an evaluated biological project (#409).

Each method's count table is thinned binomially to a series of depths: every read is kept with probability
depth / the sample's sequenced reads, the same probability for every method. Depths are thinned from the deepest down,
so each is a subset of the one above. For one method's detection this is equivalent to subsampling the reads when a
read counts for one miRNA; for methods that split a multimapping read into fractions (tiny-count, mirflowz) counts are
rounded first, so it is approximate there. Only depths every sample reaches are used, and the full depth (each
sample's own) ends the series.

Each table is thinned on its own (each method's, and the isomiR output apart from the common table), while real
subsampling removes the same reads for every method: the agreement between methods therefore falls faster with depth
than on real shallow data, and is a lower bound.

Per method and depth: median miRNA reads per sample, median miRNAs with at least 1, 5 and 10 reads per sample,
miRNAs with at least 10 reads in every sample, the median replicate CV within the groups, and, from the isomiR
output, the median distinct read sequences per sample; per pair of methods the biological report's agreement
(Spearman, median |log2 ratio|, top-100 overlap). The depth at which detection reaches 90% and 95% of its full-depth
value is interpolated on log depth; it is relative to the dataset's own depth, so the rise per doubling of depth over
the last step tells whether the curve is still rising there.
"""
import collections
import math
import statistics

import numpy as np

import biological

DEPTHS = '500k,1M,2M,4M,8M,16M,32M'
LEVELS = (1, 5, 10)            # miRNAs with at least this many reads in a sample
SATURATION = (0.9, 0.95)       # shares of full-depth detection whose depth is reported
SATURATION_KEY = 'detected_10'
AGREEMENT_KEYS = ('spearman_either', 'median_abs_log2_ratio', 'top{}_overlap'.format(biological.TOP))
FULL = 'full'


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


def cv_sets(common, methods, groups):
    """{(method, group): the miRNAs with at least MIN_READS reads in every sample of the group}, for groups of two
    samples or more: the replicate CV's miRNAs, fixed over a series at its lowest depth."""
    by = biological.counts_by(common)
    return {(m, g): set.intersection(*({i for i, v in by.get((m, s), {}).items() if v >= biological.MIN_READS} for s in members))
            for m in methods for g, members in groups.items() if len(members) > 1}


def metrics(common, methods, samples, groups, isomir_rows=None, fixed=None):
    """Per method: median miRNA reads, median miRNAs with at least 1, 5 and 10 reads per sample, miRNAs with at least
    10 in every sample, the median over groups (two samples or more) of the median CPM coefficient of variation of
    the miRNAs in `fixed` (cv_sets(); by default those of these tables), and with `isomir_rows` the median distinct read
    sequences per sample ('' for a method without isomiR output)."""
    fixed = cv_sets(common, methods, groups) if fixed is None else fixed
    by = biological.counts_by(common)
    sequences = collections.defaultdict(set)
    for row in isomir_rows or []:
        if float(row[-1]) > 0:
            sequences[(row[0], row[1])].add(row[4])
    with_isomirs = {method for method, _ in sequences}
    out = {}
    for m in methods:
        per = {s: by.get((m, s), {}) for s in samples}
        values = {'mirna_reads': round(statistics.median(sum(c.values()) for c in per.values()))}
        for level in LEVELS:
            values['detected_{}'.format(level)] = statistics.median(sum(1 for v in c.values() if v >= level) for c in per.values())
        detected = [{i for i, v in c.items() if v >= biological.MIN_READS} for c in per.values()]
        values['detected_all_samples'] = len(set.intersection(*detected)) if detected else 0
        cvs = []
        for g, members in groups.items():
            if len(members) < 2:
                continue
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


def series(common, methods, samples, reads_in, groups, depths, seed, isomir_rows=None):
    """The depth series: ([(depth or FULL, metrics(), biological.agreement())] from the lowest depth up, the depths
    left out because a sample has fewer reads). The tables are thinned from the full depth down, each depth from the
    one above (so a subset of it), with one random generator seeded with `seed`; the replicate CV uses the miRNAs of
    the lowest depth (cv_sets()) all along."""
    usable = [d for d in depths if all(reads_in[s] >= d for s in samples)]
    rng = np.random.default_rng(seed)
    tables = {FULL: (common, isomir_rows)}
    above = None
    for depth in sorted(usable, reverse=True):
        fraction = {s: depth / (reads_in[s] if above is None else above) for s in samples}
        rows, iso = tables[FULL if above is None else above]
        tables[depth] = (thin(rows, fraction, rng), thin(iso, fraction, rng) if iso is not None else None)
        above = depth
    order = usable + [FULL]
    fixed = cv_sets(tables[order[0]][0], methods, groups)
    out = [(depth, metrics(tables[depth][0], methods, samples, groups, tables[depth][1], fixed),
            biological.agreement(tables[depth][0], methods, samples)) for depth in order]
    return out, [d for d in depths if d not in usable]


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


def saturation(points, methods, full_depth, key=SATURATION_KEY, shares=SATURATION):
    """{method: {share: reads}}: the depth at which `key` first reaches `share` of its full-depth value, interpolated
    on log depth between the depths around it; the full depth (the median of the samples' reads) is the last point.
    '' when the full-depth value is 0, and the smallest depth when the first point already reaches it."""
    out = {}
    for m in methods:
        curve = [(full_depth if depth == FULL else depth, values[m][key]) for depth, values, _ in points]
        top = curve[-1][1]
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
            out[m][share] = round(reached) if reached != '' else ''
    return out


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


def report(name, points, methods, skipped, reached, rising, reads_in, seed, md_table):
    """depth.md's lines."""
    keys = list(next(iter(points[0][1].values())))
    lines = ['# Sequencing depth: {}'.format(name), '',
             'Each method\'s count table thinned binomially to each depth (#409): every read kept with probability depth / the '
             'sample\'s sequenced reads (seed {}), from the full depth down, each depth a subset of the one above; fractional '
             'counts (tiny-count, mirflowz) rounded first. Each table is thinned on its own (each method\'s; the isomiR output '
             'apart from the common table), whereas real subsampling removes the same reads for every method, so the agreement '
             'between methods falls faster than on real shallow data: a lower bound. `full`: each sample\'s own depth, {:,} to '
             '{:,} reads sequenced (median {:,}). The replicate CV uses the miRNAs with at least {} reads in every sample of a '
             'group at the lowest depth, the same set at every depth.'.format(
                 seed, min(reads_in.values()), max(reads_in.values()), round(statistics.median(reads_in.values())), biological.MIN_READS)]
    if skipped:
        lines += ['', 'Left out, above the shallowest sample: {}.'.format(', '.join(label(d) for d in skipped))]
    lines += ['', '**Depth reaching {} of full-depth detection** (miRNAs with at least {} reads, median per sample; reads '
              'sequenced per sample). It is relative to this dataset\'s depth: `rise per doubling` is the share of full-depth '
              'detection gained per doubling of depth over the last step, near 0 when the curve has levelled off.'.format(
                  ' and '.join('{:.0%}'.format(s) for s in SATURATION), biological.MIN_READS), '']
    lines += md_table(['method'] + ['{:.0%}'.format(s) for s in SATURATION] + ['rise per doubling'],
                      [[m] + ['{:,}'.format(reached[m][s]) if reached[m][s] != '' else '' for s in SATURATION] +
                       ['{:.1%}'.format(rising[m]) if rising[m] != '' else ''] for m in methods])
    for key in keys:
        lines += ['', '**{}**'.format(key), ''] + md_table(['depth'] + methods, [[label(depth)] + [values[m][key] for m in methods]
                                                                                for depth, values, _ in points])
    pairs = list(points[0][2])
    if pairs:
        lines += ['', '**Agreement between methods** (Spearman over miRNAs either detects, median per sample; a lower bound below '
                  'full depth, as each method is thinned on its own):', '']
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
