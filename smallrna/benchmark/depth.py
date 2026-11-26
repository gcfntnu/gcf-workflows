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
import pool

DEPTHS = '500k,1M,2M,4M,8M,16M,32M'
LEVELS = (1, 5, 10)            # miRNAs with at least this many reads in a sample
SATURATION = (0.9, 0.95)       # shares of full-depth detection whose depth is reported
SATURATION_KEY = 'detected_10'
CPM = 0.1                      # detection that scales with depth: reads per million, about 3-5 reads at 30-50M (#409)
CPM_KEY = 'detected_cpm'
SEEN = 3                       # reads that count a full-depth miRNA (feature) as still seen at a depth
SEEN_KEY = 'full_cpm_seen'     # share of the full depth's miRNAs at CPM or more seen with SEEN reads at a depth
SATURATION_KEYS = (SATURATION_KEY, SEEN_KEY)  # the saturation table's two definitions of detected
AGREEMENT_KEYS = ('spearman_either', 'median_abs_log2_ratio', 'top{}_overlap'.format(biological.TOP))
STRATA = (1, 10, 100)          # abundance strata of the full depth's miRNAs: per million from `cpm` up to 1, 10, 100, and above
STRATUM_PREFIX = 'seen_cpm_'   # their metrics: the share of a stratum still seen with SEEN reads
GAIN_KEYS = ('detected_1', 'detected_5', 'detected_10')   # new miRNAs per step of depth, at a ladder of read counts
STRICT_FDR = 0.01              # depth DE: recall also of the full depth's calls at this BH-adjusted p (#414)
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


def cpm_sets(by, keys, cpm=CPM, totals=None):
    """{key: the features with at least `cpm` reads per million of the key's reads (and at least 1 read)}, for {key:
    {feature: reads}} `by`: a full depth's reference set (SEEN_KEY). `totals` {key: reads}: the per-million denominator
    instead of the key's own reads (e.g. the sample's filtered reads, the same for every route)."""
    out = {}
    for key in keys:
        counts = by.get(key, {})
        total = totals[key] if totals else sum(counts.values())
        out[key] = {f for f, v in counts.items() if v >= 1 and v * 1e6 >= cpm * total}
    return out


def stratum_names(cpm=CPM):
    """The abundance strata's metric names, from `cpm` per million up: seen_cpm_0.1-1, seen_cpm_1-10, ..., seen_cpm_100+."""
    edges = [cpm] + [e for e in STRATA if e > cpm]
    return ['{}{:g}-{:g}'.format(STRATUM_PREFIX, a, b) for a, b in zip(edges, edges[1:])] + ['{}{:g}+'.format(STRATUM_PREFIX, edges[-1])]


def strata_sets(by, keys, cpm=CPM, totals=None):
    """{key: {stratum name: features}}: the cpm_sets() of `by` (a full depth's counts) split by reads per million at
    the STRATA edges, each feature in the stratum its full-depth abundance falls in."""
    edges, names = [e for e in STRATA if e > cpm], stratum_names(cpm)
    out = {}
    for key, kept in cpm_sets(by, keys, cpm, totals).items():
        counts = by.get(key, {})
        total = totals[key] if totals else sum(counts.values())
        out[key] = {n: set() for n in names}
        for f in kept:
            out[key][names[sum(1 for e in edges if counts[f] * 1e6 >= e * total)]].add(f)
    return out


def seen_share(counts, reference):
    """The share of `reference` with at least SEEN reads in {feature: reads}; '' for an empty reference."""
    return sum(1 for f in reference if counts.get(f, 0) >= SEEN) / len(reference) if reference else ''


def median_of(values):
    values = [v for v in values if v != '']
    return round(statistics.median(values), 4) if values else ''


def metrics(common, methods, samples, groups, isomir_rows=None, fixed=None, with_isomirs=None, cpm=CPM, reference=None,
            strata=None, library=None):
    """Per method: median miRNA reads, median miRNAs with at least 1, 5 and 10 reads per sample, median miRNAs with
    at least `cpm` reads per million of the sample's miRNA reads, or of `library` {sample: reads} when given
    (detected_cpm; at least 1 read, so below the depth
    where `cpm` is a read it counts every miRNA with a read), with `reference` ({(method, sample): the full depth's
    cpm_sets()}) the median share of a sample's reference miRNAs seen with at least SEEN reads (full_cpm_seen), with
    `strata` ({(method, sample): strata_sets()}) the same per abundance stratum (stratum_names()),
    miRNAs with at least 10 in every sample, the median over groups (two samples or more) of the median CPM coefficient of variation of
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
        total = (lambda s: library[s]) if library else (lambda s: sum(per[s].values()))
        values[CPM_KEY] = statistics.median(sum(1 for v in per[s].values() if v >= 1 and v * 1e6 >= cpm * total(s)) for s in samples)
        if reference is not None:
            values[SEEN_KEY] = median_of(seen_share(per[s], reference.get((m, s), set())) for s in samples)
        if strata is not None:
            for name in stratum_names(cpm):
                values[name] = median_of(seen_share(per[s], strata.get((m, s), {}).get(name, set())) for s in samples)
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


def series(common, methods, samples, reads_in, groups, depths, seed, isomir_rows=None, cpm=CPM):
    """The depth series: (series_points() of the thinned() tables, the depths left out because a sample has fewer
    reads)."""
    tables, usable = thinned(common, samples, reads_in, depths, seed, isomir_rows)
    return series_points(tables, usable, methods, samples, groups, isomir_rows, cpm), [d for d in depths if d not in usable]


def full_strata(tables, methods, samples, cpm=CPM, library=None):
    """(the full depth's cpm_sets() and strata_sets() per (method, sample)), per million of each route's miRNA reads in
    the sample, or of `library` {FULL: {sample: reads}} when given."""
    keys = [(m, s) for m in methods for s in samples]
    totals = {(m, s): library[FULL][s] for m, s in keys} if library else None
    by = biological.counts_by(tables[FULL][0])
    return cpm_sets(by, keys, cpm, totals), strata_sets(by, keys, cpm, totals)


def strata_sizes(strata, methods, samples, cpm=CPM):
    """{method: {stratum name: median miRNAs per sample in it}} of full_strata()'s strata."""
    return {m: {n: whole(statistics.median(len(strata[(m, s)][n]) for s in samples)) for n in stratum_names(cpm)} for m in methods}


def library_at(library, reads_in, usable):
    """{depth or FULL: {sample: reads}}: a per-million denominator (`library` {sample: reads at full depth}, e.g. the
    filtered reads) at each thinned depth, scaled by depth / the sample's sequenced reads."""
    out = {FULL: dict(library)}
    for depth in usable:
        out[depth] = {s: library[s] * depth / reads_in[s] for s in library}
    return out


def series_points(tables, usable, methods, samples, groups, isomir_rows=None, cpm=CPM, library=None):
    """[(depth or FULL, metrics(), biological.agreement())] of thinned() tables, from the lowest depth up; the replicate
    CV of a method and group uses the miRNAs of the lowest depth where it has any (cv_sets()) all along. `library`
    (library_at()): the per-million denominator at each depth instead of each route's miRNA reads."""
    order = usable + [FULL]
    fixed = {}
    for depth in order:  # from the lowest depth up, the first where a method and group have any
        for key, kept in cv_sets(tables[depth][0], methods, groups).items():
            if kept and not fixed.get(key):
                fixed[key] = kept
    fixed = {key: fixed.get(key, set()) for key in cv_sets([], methods, groups)}
    with_isomirs = {row[0] for row in isomir_rows or [] if float(row[-1]) > 0 and row[4]}
    reference, strata = full_strata(tables, methods, samples, cpm, library)
    return [(depth, metrics(tables[depth][0], methods, samples, groups, tables[depth][1], fixed, with_isomirs, cpm, reference,
                            strata, library[depth] if library else None),
             biological.agreement(tables[depth][0], methods, samples)) for depth in order]


def gains(points, methods, full_depth, keys=GAIN_KEYS):
    """[(from depth, to depth, {method: {key: (new per doubling of depth, new per extra million reads)}})] over each
    step of the series, the full depth at `full_depth` (the median sample's; with unequal samples that last step is
    approximate): the marginal gain of detection, the published form of saturation."""
    xs = [full_depth if depth == FULL else depth for depth, _, _ in points]
    out = []
    for i in range(1, len(points)):
        if xs[i] <= xs[i - 1]:
            continue
        doublings, millions = math.log2(xs[i] / xs[i - 1]), (xs[i] - xs[i - 1]) / 1e6
        step = {m: {} for m in methods}
        for m in methods:
            for key in keys:
                new = points[i][1][m][key] - points[i - 1][1][m][key]
                step[m][key] = (round(new / doublings, 1), round(new / millions, 2))
        out.append((points[i - 1][0], points[i][0], step))
    return out


def gain_rows(steps, methods):
    """depth.tsv's rows of gains(): method, the step's upper depth, gain_doubling_<key> or gain_per_million_<key>, value."""
    return [[m, label(high), '{}_{}'.format(kind, key), values[m][key][i]] for _, high, values in steps for m in methods
            for key in values[m] for i, kind in enumerate(('gain_doubling', 'gain_per_million'))]


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
    in the same direction; `recall_strict`: `recall` of the full-depth calls at a BH-adjusted p of at most STRICT_FDR, a
    reference less made of borderline calls; `lfc_spearman`: Spearman of the log2 fold changes with the full depth's
    over the features both test; NaN without full-depth calls, or without calls (or fewer than 3 features both test)."""
    out = {}
    for row in rows:
        out[row] = {}
        for test in de.TESTS:
            at = [de.call(results[(row, test, f)]) for f in features]
            out[row][test + '_calls'] = sum(1 for c in at if c)
            if full is not None:
                fits = [full[(row, test, f)] for f in features]
                pairs = list(zip(at, (de.call(fit) for fit in fits)))
                strict = [fit[2] is not None and fit[2] <= STRICT_FDR for fit in fits]
                both = [(results[(row, test, f)][0], fit[0]) for f, fit in zip(features, fits)
                        if results[(row, test, f)][0] is not None and fit[0] is not None]
                out[row].update({test + '_recall': de.share([c == top for c, top in pairs if top]),
                                 test + '_precision': de.share([c == top for c, top in pairs if c]),
                                 test + '_recall_strict': de.share([c == top for (c, top), s in zip(pairs, strict) if s and top]),
                                 test + '_lfc_spearman': round(pool.spearman(*zip(*both)), 4) if len(both) >= 3 else float('nan')})
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


def de_lines(contrasts, nulls, methods, versions, md_table, between=(), library='every miRNA read the route reported'):
    """depth.md's section on differential expression by depth (#414): per DE set of the contrast and test, each
    route's recall / precision (calls) per depth; per test, the null comparisons' median share called (share of splits
    with a call) per depth; a mixture's contrasts between fractions (`between`, as for de_tables) against the expected
    changes. `versions`: a line on the tools de.R recorded; `library`: what a library size counts (a mixture's differ)."""
    lines = ['', '## Differential expression by depth', '',
             "The evaluation's DE sets tested on each depth's thinned tables, the draw of the detection metrics above, and on the "
             "full depth (de.R: edgeR's exact test, its GLM test and limma-voom, #447), each over the same miRNAs at every depth: "
             "those any route reports at the full depth in the set's samples. A row tests the miRNAs with at least {} reads on "
             'average; library sizes are {} at that depth. Each route against its own calls at the full depth (BH-adjusted p of at most {}, in the same direction): '
             '`recall / precision (calls)`; `full` shows the calls. These measure how depth changes a route\'s calls, not '
             'whether the calls are right: the full-depth calls are no truth, and a route that calls more there has more to '
             'recall. They are optimistic twice over: each depth is a subset of the full-depth reads, so its calls share their '
             'noise, and the full depth is what the project holds (a read cap where one was set), not what was sequenced. '
             'In brackets, the recall of the full-depth calls at a BH-adjusted p of at most {} (`recall_strict`): borderline '
             'full-depth calls are partly noise a shallower depth cannot recall, so this is the fairer recall; precision stays '
             'against the {} calls, and the null comparisons below measure false calls.'.format(
                 biological.MIN_READS, library, de.FDR, STRICT_FDR, de.FDR)]
    for _, title, points in contrasts:
        for test in de.TESTS:
            lines += ['', '**{}**, {}'.format(title, de.TEST_NAMES[test]), '']
            lines += md_table(['depth'] + methods, [[label(depth)] + [
                '{} [{}] / {} ({})'.format(share_cell(values[m][test + '_recall']), share_cell(values[m][test + '_recall_strict']),
                                           share_cell(values[m][test + '_precision']), values[m][test + '_calls'])
                if test + '_recall' in values[m] else str(values[m][test + '_calls']) for m in methods] for depth, values in points])
        lines += ['', '**{}**, log2 fold changes against the full depth\'s (Spearman over the miRNAs both test; {})'.format(
            title, ' / '.join(de.TEST_NAMES[t] for t in de.TESTS)), '']
        lines += md_table(['depth'] + methods, [[label(depth)] + [' / '.join(share_cell(values[m][t + '_lfc_spearman']) for t in de.TESTS)
                                                                  for m in methods] for depth, values in points if depth != FULL])
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


def detected_titles(cpm):
    """{saturation key: its column title}: the two definitions of detected."""
    return {SATURATION_KEY: 'at least {} reads'.format(biological.MIN_READS),
            SEEN_KEY: 'full depth\'s {:g} per million seen'.format(cpm)}


def saturation_table(names, points, reached, rising, titles, md_table, first='method', extra=None):
    """The saturation table's lines: per name (method or class), for each definition of detected in `titles`, the
    depths reaching SATURATION of full-depth detection and the rise per doubling; `extra` {name: [cells]} after the name."""
    last = points[-2][0] if len(points) > 1 else FULL
    header = [first] + (list(extra['header']) if extra else [])
    for key, title in titles.items():
        header += ['{}: {:.0%}'.format(title, s) for s in SATURATION] + ['{}: rise per doubling'.format(title)]
    rows = []
    for n in names:
        row = [n] + (list(extra[n]) if extra else [])
        for key in titles:
            cells = reached[key].get(n)
            row += [saturation_cell(cells[s], last) if cells else '' for s in SATURATION]
            rate = rising[key].get(n, '')
            row.append('{:.1%}'.format(rate) if rate != '' else '')
        rows.append(row)
    return md_table(header, rows)


def report(name, points, methods, skipped, reached, rising, reads_in, seed, md_table, cpm=CPM, steps=None, sizes=None, library=None):
    """depth.md's lines; `reached` and `rising` by saturation key (SATURATION_KEYS); `steps` gains(), `sizes`
    strata_sizes(); `library` {sample: reads}: the per-million denominator when not each route's miRNA reads."""
    of = 'of its filtered reads, the same for every route' if library else 'of its miRNA reads'
    typical = statistics.median(library.values()) if library else statistics.median(points[-1][1][m]['mirna_reads'] for m in methods)
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
    lines += ['', '**Depth reaching {} of full-depth detection** (median per sample; reads sequenced per sample), with two '
              'definitions of detected: at least {} reads, a fixed count that keeps rising with depth as rarer reads (permissive '
              'matches among them) pass it; and `full_cpm_seen`, the share of the miRNAs a sample has at {:g} reads per '
              'million or more at its full depth ({}; {:g} per million is {:.1f} reads at this dataset\'s median) that a '
              'depth still sees with at least {} reads, which rises to its full-depth value and so shows where the miRNAs above that abundance are seen. `detected_cpm` counts the miRNAs at {:g} per million or more at '
              'each depth itself; below the depth where that is a read it counts every miRNA with a read, so it overshoots '
              'there. It is relative to this dataset\'s depth: '
              '`rise per doubling` is the share of full-depth detection gained per doubling of depth over the last step, near 0 '
              'when the curve has levelled off. With samples of different depths the full depth is no common measurement: a '
              'share the last common depth misses is `beyond` it.'.format(
                  ' and '.join('{:.0%}'.format(s) for s in SATURATION), biological.MIN_READS, cpm, of, cpm,
                  cpm * typical / 1e6, SEEN, cpm), '']
    lines += saturation_table(methods, points, reached, rising, detected_titles(cpm), md_table)
    if steps:
        lines += ['', '**New miRNAs per step of depth** (median per sample): per doubling of depth, and in brackets per extra '
                  'million reads sequenced, at a ladder of read counts; the published measure of saturation, which does not '
                  'depend on how deep the dataset happens to be. The last step ends at the median sample\'s full depth{}.'.format(
                      '' if len(set(reads_in.values())) == 1 else ', approximate with samples of different depths')]
        for key in GAIN_KEYS:
            lines += ['', '`{}`'.format(key), ''] + md_table(['step'] + methods, [['{} to {}'.format(label(low), label(high))] + [
                '{:+g} ({:+g})'.format(*values[m][key]) for m in methods] for low, high, values in steps])
    if sizes:
        names = stratum_names(cpm)
        lines += ['', '**By abundance**: the full depth\'s miRNAs at {:g} per million or more ({}) split by their full-depth '
                  'abundance; `seen_cpm_<low>-<high>` below is the share of a stratum a depth still sees with at least {} reads '
                  '(median per sample). MiRNAs per stratum at full depth (median per sample):'.format(cpm, of, SEEN), '']
        lines += md_table(['method'] + [n[len(STRATUM_PREFIX):] for n in names], [[m] + [sizes[m][n] for n in names] for m in methods])
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


def plot(points, methods, full_depth, path, cpm=CPM, steps=None):
    """depth.png: miRNAs detected (the full depth's at `cpm` per million seen, and at least 10 reads; median per sample)
    and replicate CV against depth, log scale; with `steps` (gains()) also the new miRNAs with at least 10 reads per
    doubling of depth, at each step's upper depth."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x = [full_depth if depth == FULL else depth for depth, _, _ in points]
    panels = [(SEEN_KEY, 'full depth\'s miRNAs at {:g} per million, seen ({}+ reads)'.format(cpm, SEEN)),
              (SATURATION_KEY, 'miRNAs with at least {} reads'.format(biological.MIN_READS)),
              ('replicate_cv', 'replicate CV (median)')]
    fig, axes = plt.subplots(1, len(panels) + bool(steps), figsize=(5 * (len(panels) + bool(steps)), 4.2))
    for ax, (key, title) in zip(axes, panels):
        for m in methods:
            ys = [values[m][key] for _, values, _ in points]
            if all(y != '' for y in ys):
                ax.plot(x, ys, marker='o', label=m)
        ax.set_title(title)
    if steps:
        ax = axes[-1]
        xs = [full_depth if high == FULL else high for _, high, _ in steps]
        for m in methods:
            ax.plot(xs, [values[m][SATURATION_KEY][0] for _, _, values in steps], marker='o', label=m)
        ax.axhline(0, color='grey', linewidth=0.8)
        ax.set_title('new miRNAs ({}+ reads) per doubling of depth'.format(biological.MIN_READS))
    for i, ax in enumerate(axes):
        ax.set_xscale('log')
        ax.set_xlabel('reads sequenced per sample')
        if i < len(panels):
            ax.set_ylim(bottom=0)  # from 0: an autoscaled CV axis turns differences in the third decimal into a curve
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def strata_plot(points, methods, full_depth, path, cpm=CPM):
    """depth_strata.png: per abundance stratum of the full depth's miRNAs, the share still seen with SEEN reads against
    depth (median per sample), one line per method."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = stratum_names(cpm)
    x = [full_depth if depth == FULL else depth for depth, _, _ in points]
    fig, axes = plt.subplots(1, len(names), figsize=(4 * len(names), 3.8), sharey=True)
    for ax, name in zip(axes, names):
        for m in methods:
            ys = [values[m].get(name, '') for _, values, _ in points]
            if all(y != '' for y in ys):
                ax.plot(x, ys, marker='o', label=m)
        ax.set_xscale('log')
        ax.set_ylim(0, 1.02)
        ax.set_xlabel('reads sequenced per sample')
        ax.set_title('{} per million at full depth'.format(name[len(STRATUM_PREFIX):]))
    axes[0].set_ylabel('share seen ({}+ reads)'.format(SEEN))
    if any(ax.get_lines() for ax in axes):
        next(ax for ax in axes if ax.get_lines()).legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# --- Other RNA classes by depth (#409): UNITAS's read sequences thinned, each sequence's reads split over its features
# as UNITAS splits them (classes.unitas_sequences), so a depth is a real subsample of UNITAS's reads.

CLASS_KEYS = ('reads', 'share', 'features_1', 'features_10', 'features_cpm', 'full_cpm_seen', 'features_all_samples')
CLASS_SATURATION_KEY = 'features_10'
CLASS_CPM_KEY = 'features_cpm'
CLASS_SATURATION_KEYS = (CLASS_SATURATION_KEY, SEEN_KEY)
UNASSIGNED = 'unassigned'  # reads without annotation: a share, no features


def class_thinning(sequences, reads_in, depths, seed):
    """{depth or FULL: {(class, feature): reads}} of one sample: `sequences` (reads, {(class, feature): share}) thinned
    binomially from the full depth down, each depth from the one above (probability depth / the depth above, the full
    depth being the sample's sequenced reads `reads_in`), with a generator seeded with `seed`. Each sequence's whole
    reads are thinned, then divided by its shares."""
    counts, seq_of, feature_of, weight, keys = [], [], [], [], {}
    for i, (reads, share) in enumerate(sequences):
        counts.append(reads)
        for key, w in share.items():
            seq_of.append(i)
            feature_of.append(keys.setdefault(key, len(keys)))
            weight.append(w)
    counts = np.array(counts, dtype=np.int64)
    seq_of, feature_of, weight = np.array(seq_of, dtype=np.int64), np.array(feature_of, dtype=np.int64), np.array(weight)
    names = list(keys)

    def by_feature(kept):
        totals = np.bincount(feature_of, weights=kept[seq_of] * weight, minlength=len(names)) if len(names) else np.zeros(0)
        return {names[j]: float(v) for j in np.nonzero(totals > 1e-9)[0] for v in (totals[j],)}

    rng = np.random.default_rng(seed)
    out, kept, above = {FULL: by_feature(counts)}, counts, reads_in
    for depth in sorted(depths, reverse=True):
        kept = rng.binomial(kept, depth / above)
        above = depth
        out[depth] = by_feature(kept)
    return out


def class_points(per_sample, usable, samples, cpm=CPM):
    """[(depth or FULL, {class: {CLASS_KEYS: value}})] from the lowest depth up, of each sample's class_thinning():
    per class the median over samples of its reads, its share of the sample's reads (unassigned included), its features
    with at least 1 and 10 reads and with at least `cpm` reads per million of the sample's reads (at least 1 read), the
    share of the features the sample has at `cpm` per million or more at its full depth seen with at least SEEN reads
    (full_cpm_seen), and the features with at least 10 reads in every sample. `unassigned` has no features."""
    classes = sorted({cls for s in samples for table in per_sample[s].values() for cls, _ in table})
    reference = {}
    for s in samples:
        full_table = per_sample[s][FULL]
        cutoff = cpm * sum(full_table.values()) / 1e6
        for (c, feature), v in full_table.items():
            if v >= 1 and v >= cutoff:
                reference.setdefault((c, s), set()).add(feature)
    points = []
    for depth in usable + [FULL]:
        values = {}
        for cls in classes:
            reads, shares, f1, f10, fcpm, seen, everywhere = [], [], [], [], [], [], None
            for s in samples:
                table = per_sample[s][depth]
                total = sum(table.values())
                mine = {feature: v for (c, feature), v in table.items() if c == cls}
                reads.append(sum(mine.values()))
                shares.append(sum(mine.values()) / total if total else 0.0)
                f1.append(sum(1 for v in mine.values() if v >= 1))
                f10.append(sum(1 for v in mine.values() if v >= biological.MIN_READS))
                fcpm.append(sum(1 for v in mine.values() if v >= 1 and v * 1e6 >= cpm * total))
                seen.append(seen_share(mine, reference.get((cls, s), set())))
                level = {feature for feature, v in mine.items() if v >= biological.MIN_READS}
                everywhere = level if everywhere is None else everywhere & level
            none = cls == UNASSIGNED
            values[cls] = {'reads': round(statistics.median(reads)), 'share': round(statistics.median(shares), 5),
                           'features_1': '' if none else whole(statistics.median(f1)),
                           'features_10': '' if none else whole(statistics.median(f10)),
                           'features_cpm': '' if none else whole(statistics.median(fcpm)),
                           'full_cpm_seen': '' if none else median_of(seen),
                           'features_all_samples': '' if none else len(everywhere or ())}
        points.append((depth, values))
    return points


def whole(value):
    """A median as an int when it is whole (1, not 1.0), else as it is (1.5)."""
    return int(value) if float(value).is_integer() else value


def class_rows(points):
    """depth_classes.tsv rows: class, depth, metric, value."""
    return [[cls, label(depth), key, value] for depth, values in points for cls in values for key, value in values[cls].items()]


def class_lines(points, reached, rising, seed, md_table, cpm=CPM):
    """depth.md's section on the other RNA classes; `reached` and `rising` by CLASS_SATURATION_KEYS."""
    classes = sorted(points[-1][1], key=lambda c: -points[-1][1][c]['reads'])
    lines = ['', '## RNA classes by depth (UNITAS)', '',
             'UNITAS\'s read sequences (its full annotation matrix) thinned binomially per sample (seed {}, one generator per '
             'sample, apart from the miRNA tables\' draw above), from the full depth down, each sequence\'s reads then divided '
             'over its features as UNITAS divides them: a real subsample of UNITAS\'s reads. A feature is an Ensembl gene '
             '(records naming one), a tRNA gene (tRNA fragments), an rRNA subunit (SILVA copies) or a piRNA cluster; miRNAs by '
             'name. `share` is of the sample\'s reads UNITAS read, unassigned included. Medians per sample, except features in '
             'all samples.'.format(seed), '',
             '**Depth reaching {} of full-depth detection** (features, median per sample): at least {} reads, or the share '
             'of the features at {:g} per million of the sample\'s reads or more at full depth seen with at least {} reads '
             '(`full_cpm_seen`):'.format(' and '.join('{:.0%}'.format(s) for s in SATURATION), biological.MIN_READS, cpm, SEEN), '']
    titles = {CLASS_SATURATION_KEY: 'at least {} reads'.format(biological.MIN_READS), SEEN_KEY: 'full depth\'s {:g} per million seen'.format(cpm)}
    full = points[-1][1]
    extra = {'header': ['share at full depth', 'features at full depth ({:g} per million)'.format(cpm)]}
    extra.update({c: ['{:.2%}'.format(full[c]['share']), full[c][CLASS_CPM_KEY]] for c in classes})
    lines += saturation_table(classes, points, reached, rising, titles, md_table, first='class', extra=extra)
    for key, title in (('share', 'Share of the reads'),
                       (SEEN_KEY, 'Share of the full depth\'s features at {:g} per million seen with at least {} reads'.format(cpm, SEEN)),
                       (CLASS_CPM_KEY, 'Features at {:g} per million or more'.format(cpm)),
                       ('features_10', 'Features with at least {} reads'.format(biological.MIN_READS))):
        lines += ['', '**{}**'.format(title), '', '| depth | ' + ' | '.join(classes) + ' |', '|---|' + '---|' * len(classes)]
        for depth, values in points:
            cells = [values[c][key] for c in classes]
            lines.append('| {} | {} |'.format(label(depth), ' | '.join(
                '{:.2%}'.format(v) if key == 'share' and v != '' else '{:.1%}'.format(v) if key == SEEN_KEY and v != '' else str(v)
                for v in cells)))
    lines += ['', 'Reads, features with at least 1 read and features in all samples per depth are in `depth_classes.tsv`.']
    return lines


def class_plot(points, full_depth, path, top=7, cpm=CPM):
    """depth_classes.png: per class, the share of the full depth's features at `cpm` per million or more seen with SEEN
    reads, and the
    class's share of the reads, against depth (log scale); the `top` classes by reads at full depth (assigned ones)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import classes as class_names
    full = points[-1][1]
    shown = [c for c in sorted(full, key=lambda c: -full[c]['reads']) if c != UNASSIGNED and full[c][SEEN_KEY] != ''][:top]
    x = [full_depth if depth == FULL else depth for depth, _ in points]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, cls in enumerate(shown):
        colour = class_names.SERIES[i % len(class_names.SERIES)]
        axes[0].plot(x, [values[cls][SEEN_KEY] if values[cls][SEEN_KEY] != '' else float('nan') for _, values in points], marker='o',
                     color=colour, label=cls)
        axes[1].plot(x, [values[cls]['share'] for _, values in points], marker='o', color=colour, label=cls)
    axes[0].set_title('full depth\'s features at {:g} per million, seen ({}+ reads)'.format(cpm, SEEN))
    axes[1].set_title('share of the reads (median)')
    for ax in axes:
        ax.set_xscale('log')
        ax.set_xlabel('reads sequenced per sample')
        ax.set_ylim(bottom=0)
    axes[1].set_yscale('symlog', linthresh=0.01)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
