"""The summary's overview (#338): one row per route and one column per criterion across the evaluated projects, each
cell binned as good, intermediate or poor by fixed cut-offs (the form of Soneson & Robinson 2018, Fig. 5), with the
raw value printed. There is no overall score: the criteria rest on different kinds of truth.

Each column takes its projects by what they declare and hold, never by directory name: the genome-wide synthetic
sets with designed changes, the equimolar and ratiometric pools per library kit, the mixture, a biological project
with knockouts, one with null comparisons, the synthetic set with decoy probes, and the biological project with the
most samples for the cost. A column without a project is left out; a route the projects declare but did not run
shows its state (`not run`, `stopped`). Settings variants (#417) are left out.
"""
import csv
import math
import statistics
from pathlib import Path

import yaml

import classes
import de
import evaluation
import resources

GOOD, INTERMEDIATE, POOR = 'good', 'intermediate', 'poor'
UP, DOWN = 'up', 'down'
NULL_LEVEL = 0.05  # the nominal level of the null comparisons: only calls above it count against a route


def bin_of(value, direction, good, poor):
    """GOOD at or past `good` (inclusive), POOR strictly past `poor`, otherwise INTERMEDIATE; `direction` UP when
    higher is better, DOWN when lower is."""
    if direction == UP:
        return GOOD if value >= good else POOR if value < poor else INTERMEDIATE
    return GOOD if value <= good else POOR if value > poor else INTERMEDIATE


def kit_label(kit):
    """A short name of a library kit for a column header."""
    for key, label in (('TruSeq', 'TruSeq'), ('NEXTflex', 'NEXTflex'), ('NEBNext', 'NEBNext'), ('QIAseq', 'QIAseq')):
        if key in kit:
            version = kit.split()[-1] if kit.split()[-1].startswith('v') else ''
            return ' '.join(x for x in (label, version) if x)
    return kit


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def number(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


class Project:
    """What the overview reads of one evaluated project."""

    def __init__(self, path):
        self.path = Path(path)
        bench = self.path / 'benchmark'
        self.spec = yaml.safe_load((bench / 'dataset.yaml').read_text())
        table, self.record, intervals, routes = evaluation.read_evaluation(bench)
        self.columns = table[0][1:]
        self.metrics = {(row[0], c): number(v) for row in table[1:] for c, v in zip(self.columns, row[1:])}
        # bootstrap intervals (#413) as (value, low, high), those with both ends only
        self.intervals = {(ci['metric'], ci['row']): (number(ci['value']), number(ci['low']), number(ci['high'])) for ci in intervals}
        self.intervals = {k: v for k, v in self.intervals.items() if None not in v}
        self.status = {row['route']: row['state'] for row in routes if row['state'] != 'evaluated'}
        self.kind = self.spec.get('kind')
        self.pool = self.spec.get('pool')
        self.kit = self.spec.get('libprepkit', '')
        self.genome = (self.spec.get('synthetic') or {}).get('mode') == 'genome'
        # a synthetic set with designed changes: not a null scenario (#453), whose structure.group.fraction is 0
        self.designed = (((self.spec.get('synthetic') or {}).get('structure') or {}).get('group') or {}).get('fraction', 1) != 0
        self.preset = de.analysis(bench) == 'preset'  # DE metrics from #527's preset, not an earlier analysis (#531 review)
        self.bench = bench

        self.variant = bool(self.spec.get('variant'))
        self.seed = self.record.get('synthetic_seed', (self.spec.get('synthetic') or {}).get('seed'))
        samples = bench / 'samples.tsv'
        rows = read_tsv(samples) if samples.is_file() else []
        self.lanes = any(row.get('lane', 'merged') != 'merged' for row in rows)
        self.samples = len(rows)
        self.sample_names = [row['sample'] for row in rows]

    def imbalance_nulls(self):
        """{(arm, test): {route: null_p05}} of the null comparisons unthinned in each arm, from benchmark/imbalance's
        imbalance.tsv (`imbalance --ratios none`, #530); {} without it."""
        path = self.bench / 'imbalance' / 'imbalance.tsv'
        out = {}
        if path.is_file():
            for row in read_tsv(path):
                if row['design'] == 'none' and row['set'] == 'null' and row['metric'].endswith('_null_p05'):
                    out.setdefault((row['arm'], row['metric'][:-len('_null_p05')]), {})[row['method']] = number(row['value'])
        return out

    def value(self, metric, route):
        return self.metrics.get((metric, route))

    def has(self, metric):
        return any(v is not None for (m, _), v in self.metrics.items() if m == metric)


def mean(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None


def mean_metric(projects, metric, relative=False):
    """{route: mean over `projects` of `metric`}, each divided by the project's `truth` value when `relative`."""
    def one(p, r):
        v = p.value(metric, r)
        if v is None or not relative:
            return v
        t = p.value(metric, 'truth')
        return v / t if t else None
    return lambda routes: {r: mean(one(p, r) for p in projects) for r in routes}


def right_class(projects):
    """{route: mean over the decoy sources of the % of a source's reads in its right class}, from class_truth.tsv,
    over the sources whose right class the route reports (class_states.tsv `reported` or `totals only`); a source
    whose class the route does not report or lumps with others is left out, and a route without any is None."""
    def get(routes):
        out = {}
        for r in routes:
            values = []
            for p in projects:
                bench = p.path / 'benchmark'
                states = {(row['method'], row['class']): row['state'] for row in read_tsv(bench / 'class_states.tsv')}
                percent, scored = {}, set()
                for row in read_tsv(bench / 'class_truth.tsv'):
                    if row['method'] == r:
                        percent[(row['source'], row['class'])] = float(row['percent'])
                        scored.add(row['source'])
                for source, right in classes.DECOY_CLASSES.items():  # a source the route was scored on; its right class reported
                    if source in scored and any(states.get((r, c), '').startswith(('reported', 'totals only')) for c in right):
                        values.append(sum(percent.get((source, c), 0.0) for c in right))
            out[r] = mean(values)
        return out
    return get


def cost(project):
    """({route: wall minutes per sample}, {route: peak memory GB}) of a project's resources.tsv: every quantification
    step's wall time, per sample or per project, summed and divided by the project's samples; the peak resident memory
    of any step. A route is left out unless every quantification job resources.STEPS requires has a row, as its
    report's runtime section requires (a partial rerun otherwise shows a partial cost)."""
    rows = read_tsv(project.path / 'benchmark' / 'resources.tsv')
    samples = project.samples or len({row['sample'] for row in rows if row['sample']}) or 1
    wall, peak, seen = {}, {}, {}
    for row in rows:
        if row['kind'] != 'quant':
            continue
        seen.setdefault(row['method'], set()).add((row['step'], row['sample']))
        wall[row['method']] = wall.get(row['method'], 0.0) + (number(row['wall_s']) or 0.0)
        rss = number(row['max_rss_mb'])
        if rss is not None:
            peak[row['method']] = max(peak.get(row['method'], 0.0), rss)

    def complete(method):
        steps = resources.STEPS.get(method)
        return steps is None or all((step, sample) in seen[method] for step, per_sample, kind in steps if kind == 'quant'
                                    for sample in (project.sample_names if per_sample else ['']))
    measured = [m for m in wall if complete(m)]
    return ({m: wall[m] / samples / 60 for m in measured}, {m: peak[m] / 1024 for m in measured if m in peak})


def columns(projects):
    """[(group, label, direction, good, poor, format, {route: value}-getter, projects used)] for the projects at hand."""
    usable = [p for p in projects if not p.variant]
    out, preset = [], []  # preset: the DE columns of #527's preset, one group beside the knockouts
    synthetic = [p for p in usable if p.kind == 'synthetic' and p.genome and p.designed and p.has('recovered')]
    if synthetic:
        group = 'Synthetic ({} sets)'.format(len(synthetic))
        out += [(group, 'reads on the right miRNA', UP, 0.97, 0.90, '{:.3f}', mean_metric(synthetic, 'recovered'), synthetic),
                (group, 'counted reads that are right', UP, 0.99, 0.98, '{:.3f}', mean_metric(synthetic, 'precision'), synthetic),
                (group, 'counts exactly right', UP, 0.5, 0.2, '{:.2f}', mean_metric(synthetic, 'exact'), synthetic),
                ]
        called = [p for p in synthetic if p.preset and p.value('voom_called', 'truth')]   # the preset's test (#527), a usable truth
        if called:
            preset.append((PRESET_GROUP, 'synthetic: real changes found vs truth', UP, 0.95, 0.85, '{:.2f}',
                           mean_metric(called, 'voom_called', True), called))
    for pool, metric, label, direction, good, poor, fmt, relative in (
            ('equimolar', 'within_2fold', 'within 2-fold', UP, 0.95, 0.80, '{:.2f}', False),
            ('ratiometric', 'median_abs_error', 'ratio error vs exact match', DOWN, 1.05, 1.20, '{:.2f}', True)):
        for kit in sorted({p.kit for p in usable if p.kind == 'pool' and p.pool == pool}):
            chosen = [p for p in usable if p.kind == 'pool' and p.pool == pool and p.kit == kit and p.has(metric)]
            if not chosen:
                continue
            get = (lambda ps, m: lambda routes: {r: mean(p.value(m, r) / p.value(m, 'oracle') if p.value(m, r) is not None
                                                          and p.value(m, 'oracle') else None for p in ps) for r in routes})(chosen, metric) \
                if relative else mean_metric(chosen, metric)
            out.append(({'equimolar': 'Equal-amount pools', 'ratiometric': 'Set-ratio pools'}[pool], '{} {}'.format(label, kit_label(kit)), direction, good, poor, fmt, get, chosen))
    mixture = [p for p in usable if p.kind == 'mixture' and p.has('linearity_error')]
    if mixture:
        out += [('Mixture', 'off the mixing proportions', DOWN, 0.25, 0.30, '{:.3f}', mean_metric(mixture, 'linearity_error'), mixture),
                ('Mixture', 'within 1.5-fold of predicted', UP, 0.80, 0.70, '{:.2f}', mean_metric(mixture, 'within_1_5fold'), mixture)]
    knockouts = [p for p in usable if p.kind == 'biological' and p.has('knockout_surviving_unexpected')]
    if knockouts:
        out.append(('Knockout', 'knocked-out miRNAs counted', DOWN, 35, 80, '{:.0f}', mean_metric(knockouts, 'knockout_surviving_unexpected'), knockouts))
    nulls = [p for p in usable if p.kind == 'biological' and not p.lanes and p.preset and p.has('voom_null_p05')]
    if nulls:  # the preset's test (#527): limma-voom, TMMwsp factors, filterByExpr
        get = (lambda ps: lambda routes: {r: (lambda v: None if v is None else max(0.0, v - NULL_LEVEL))(mean(p.value('voom_null_p05', r) for p in ps))
                                          for r in routes})(nulls)
        preset.append((PRESET_GROUP, 'real data: false changes above 5%', DOWN, 0.01, 0.03, '{:.3f}', get, nulls))
    out += preset
    decoy = [p for p in usable if p.kind == 'synthetic' and p.has('decoy_as_mirna')]
    typed = [p for p in decoy if (p.path / 'benchmark' / 'class_truth.tsv').is_file() and (p.path / 'benchmark' / 'class_states.tsv').is_file()]
    if decoy:
        out.append(('Other RNA (decoys)', 'counted as miRNA', DOWN, 0.005, 0.02, '{:.3f}', mean_metric(decoy, 'decoy_as_mirna'), decoy))
    if typed:  # the right class needs the class tables, which an evaluation without class scoring lacks
        out.append(('Other RNA (decoys)', 'right RNA type (%)', UP, 90, 70, '{:.0f}', right_class(typed), typed))
    measured = [p for p in usable if p.kind == 'biological' and not p.lanes and (p.path / 'benchmark' / 'resources.tsv').is_file()]
    if measured:
        biggest = max(measured, key=lambda p: (p.samples, str(p.path)))
        minutes, gb = cost(biggest)
        group = 'Cost ({})'.format(biggest.spec['name'])
        out += [(group, 'min per sample', None, None, None, '{:.0f}', lambda routes, v=minutes: {r: v.get(r) for r in routes}, [biggest]),
                (group, 'peak GB', None, None, None, '{:.1f}', lambda routes, v=gb: {r: v.get(r) for r in routes}, [biggest])]
    return out


def build(paths):
    """(routes, [(group, label, direction, good, poor)], {(route, column index): (text, bin or state or None)}, the projects
    each column used), or None without any column."""
    projects = [Project(p) for p in paths]
    cols = columns(projects)
    if not cols:
        return None
    routes = list(dict.fromkeys(c for p in projects for c in p.columns + list(p.status) if c not in ('oracle', 'truth')))
    cells = {}
    for j, (group, label, direction, good, poor, fmt, get, used) in enumerate(cols):
        values = get(routes)
        for r in routes:
            state = next((p.status[r] for p in used if r in p.status), None)
            v = values.get(r)
            if v is None:
                cells[(r, j)] = (state or 'n/a', state)
            else:
                text = fmt.format(v)  # binned as shown, so a value printed at a cut-off is binned at it
                cells[(r, j)] = (text, bin_of(float(text), direction, good, poor) if direction else None)
    used = [sorted(p.path.name for p in c[7]) for c in cols]
    return routes, [c[:6] for c in cols], cells, used


ARROW = {UP: '↑', DOWN: '↓', None: ''}
TESTS = ('edger', 'edger_glm', 'voom')
TEST_LABELS = {'edger': 'edgeR exact', 'edger_glm': 'edgeR GLM', 'voom': 'limma-voom'}
PRESET_GROUP = 'DE under the preset'  # #527: limma-voom, TMMwsp factors, filterByExpr

NOTE = ('One row per route, one column per criterion; ↑ higher is better, ↓ lower is better. Each cell is the raw value, '
        'binned as good, intermediate or poor by the cut-offs below (good at or past the first, poor past the second). '
        'Synthetic: reads simulated from known miRNAs, with changes designed between two groups, averaged over the sets. Equal-amount and set-ratio pools: synthetic miRNAs '
        'mixed in equal amounts or at set ratios; the reads that exactly match a pool sequence are the reference a route should '
        'reproduce, and the ratio error is divided by that reference\'s own. Mixture: two tissues mixed at known proportions. '
        'Knockout: miRNAs counted in cells where they were knocked out. DE under the preset (#527: limma-voom on edgeR\'s TMMwsp '
        'factors, the miRNAs filterByExpr keeps), for projects evaluated with it: in the synthetic sets, the real changes found, '
        'divided by the same analysis of the true counts; in real data, samples of one group split in two halves, where only '
        'about 5% of miRNAs should pass p < 0.05 by chance, and only the share above 5% counts. '
        'Other RNA: reads from known non-miRNA RNAs added to the synthetic data; right RNA type over the types a route reports. '
        'Cost is not binned. There is no overall score.')


def cutoffs(cols):
    """The printed cut-offs: '<label> ≥good / <poor' per binned column."""
    def one(c):
        group, label, direction, good, poor, fmt = c
        return '{} {} {}{} / {}{}'.format(group, label, '≥' if direction == UP else '≤', fmt.format(good),
                                           '<' if direction == UP else '>', fmt.format(poor))
    return '; '.join(one(c) for c in cols if c[2])


NORMALISATIONS = ('workflow', 'fixed', 'tmm', 'rle', 'tmmwsp')  # the analysis overview's rows per test (#527)
ANALYSIS_NOTE = ('The downstream layer (#527): each test with each normalisation, over the projects that hold it; a cell is '
                 'the median over the routes, averaged over the projects (the projects in the last row). `workflow`: the library '
                 'sizes without factors; `fixed`: each sample\'s filtered reads; `tmm`, `rle`, `tmmwsp`: edgeR\'s factors on the '
                 'workflow\'s sizes. The preset is limma-voom with `tmmwsp` on the miRNAs filterByExpr keeps (the synthetic and '
                 'biological datasets\' main set, and their arms\' features); the set-ratio pools\' main set is `workflow`. Synthetic: '
                 'the designed changes called and the FDR, the false calls where nothing changes (the null scenario) and, where most '
                 'changes go one way, the unchanged miRNAs called against them. Set-ratio pools: the median log2 fold change of the 1:1 '
                 'sequences (`offset`, 0 when the normalisation is right). Real data: the share of p-values of at most 0.05 in the null '
                 'comparisons (0.05 when calibrated), from `imbalance --ratios none` where a project has it, else the preset alone.')


def analysis(paths):
    """(header, rows, sources) of the summary's analysis overview (#527): per test and normalisation, the medians over
    the routes of the downstream metrics, averaged over the projects at `paths`; None without any."""
    usable = [p for p in (Project(path) for path in paths) if not p.variant]
    columns = []   # (label, [(project, getter(project, test, normalisation) -> {route: value} or None)])

    def routes_of(p):
        return [c for c in p.columns if c not in ('truth', 'oracle')]

    def suffixed(main):
        return lambda norm: '' if norm == main else '.' + norm

    def from_metrics(metric, main):
        sfx = suffixed(main)
        return lambda p, test, norm: {r: p.value('{}_{}{}'.format(test, metric, sfx(norm)), r) for r in routes_of(p)}

    synthetic = [p for p in usable if p.kind == 'synthetic' and p.genome and p.preset]
    designed = [p for p in synthetic if p.designed]
    one_way = [p for p in designed if p.has('voom_false_calls_against')]
    null_scenarios = [p for p in synthetic if not p.designed]
    pools = [p for p in usable if p.kind == 'pool' and p.pool == 'ratiometric' and p.has('voom_offset')]
    real = [p for p in usable if p.kind == 'biological' and not p.lanes and p.preset and p.has('voom_null_p05')]
    columns += [('synthetic: changes called', designed, from_metrics('called', 'tmmwsp')),
                ('synthetic: FDR', designed, from_metrics('fdr', 'tmmwsp')),
                ('synthetic null: false calls', null_scenarios, from_metrics('false_calls', 'tmmwsp')),
                ('one-way changes: calls against them', one_way, from_metrics('false_calls_against', 'tmmwsp')),
                ('set-ratio pools: offset', pools, from_metrics('offset', 'workflow'))]

    def real_null(p, test, norm):
        found = p.imbalance_nulls().get((norm, test))
        if found is not None:
            return found
        return {r: p.value(test + '_null_p05', r) for r in routes_of(p)} if norm == 'tmmwsp' else None
    columns.append(('real data: null p ≤ 0.05', real, real_null))
    columns = [c for c in columns if c[1]]
    if not columns:
        return None
    rows = []
    for test in TESTS:
        for norm in NORMALISATIONS:
            cells = []
            for _, chosen, get in columns:
                per_project = []
                for p in chosen:
                    values = get(p, test, norm)
                    finite = [v for v in (values or {}).values() if v is not None]
                    if finite:
                        per_project.append(statistics.median(finite))
                cells.append('{:.3f}'.format(statistics.mean(per_project)) if per_project else '')
            if any(cells):
                rows.append([TEST_LABELS[test], norm] + cells)
    sources = [', '.join(p.path.name for p in chosen) for _, chosen, _ in columns]
    return ['test', 'normalisation'] + [c[0] for c in columns], rows, sources
