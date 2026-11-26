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
import evaluation

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
        self.variant = bool(self.spec.get('variant'))
        self.seed = self.record.get('synthetic_seed', (self.spec.get('synthetic') or {}).get('seed'))
        samples = bench / 'samples.tsv'
        rows = read_tsv(samples) if samples.is_file() else []
        self.lanes = any(row.get('lane', 'merged') != 'merged' for row in rows)
        self.samples = len(rows)

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
                percent = {}
                for row in read_tsv(bench / 'class_truth.tsv'):
                    if row['method'] == r:
                        percent[(row['source'], row['class'])] = float(row['percent'])
                for source, right in classes.DECOY_CLASSES.items():
                    if any(states.get((r, c), '').startswith(('reported', 'totals only')) for c in right):
                        values.append(sum(percent.get((source, c), 0.0) for c in right))
            out[r] = mean(values)
        return out
    return get


def cost(project):
    """({route: wall minutes per sample}, {route: peak memory GB}) of a project's resources.tsv: every quantification
    step's wall time, per sample or per project, summed and divided by the samples; the peak resident memory of any step."""
    rows = read_tsv(project.path / 'benchmark' / 'resources.tsv')
    samples = len({row['sample'] for row in rows if row['sample']}) or 1
    wall, peak = {}, {}
    for row in rows:
        if row['kind'] != 'quant':
            continue
        wall[row['method']] = wall.get(row['method'], 0.0) + (number(row['wall_s']) or 0.0)
        peak[row['method']] = max(peak.get(row['method'], 0.0), number(row['max_rss_mb']) or 0.0)
    return ({m: v / samples / 60 for m, v in wall.items()}, {m: v / 1024 for m, v in peak.items()})


def columns(projects):
    """[(group, label, direction, good, poor, format, {route: value}-getter, projects used)] for the projects at hand."""
    usable = [p for p in projects if not p.variant]
    out = []
    synthetic = [p for p in usable if p.kind == 'synthetic' and p.genome and p.value('de_sensitivity', 'truth')]
    if synthetic:
        group = 'Synthetic ({} sets)'.format(len(synthetic))
        out += [(group, 'reads on the right miRNA', UP, 0.97, 0.90, '{:.3f}', mean_metric(synthetic, 'recovered'), synthetic),
                (group, 'counted reads that are right', UP, 0.99, 0.98, '{:.3f}', mean_metric(synthetic, 'precision'), synthetic),
                (group, 'counts exactly right', UP, 0.5, 0.2, '{:.2f}', mean_metric(synthetic, 'exact'), synthetic),
                (group, 'real changes found vs truth', UP, 0.95, 0.85, '{:.2f}', mean_metric(synthetic, 'de_sensitivity', True), synthetic)]
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
    nulls = [p for p in usable if p.kind == 'biological' and not p.lanes and p.has('edger_null_p05')]
    if nulls:
        get = (lambda ps: lambda routes: {r: (lambda v: None if v is None else max(0.0, v - NULL_LEVEL))(mean(p.value('edger_null_p05', r) for p in ps))
                                          for r in routes})(nulls)
        out.append(('No change', 'false changes above 5%', DOWN, 0.01, 0.03, '{:.3f}', get, nulls))
    decoy = [p for p in usable if p.kind == 'synthetic' and p.has('decoy_as_mirna')
             and (p.path / 'benchmark' / 'class_truth.tsv').is_file() and (p.path / 'benchmark' / 'class_states.tsv').is_file()]
    if decoy:
        out += [('Other RNA (decoys)', 'counted as miRNA', DOWN, 0.005, 0.02, '{:.3f}', mean_metric(decoy, 'decoy_as_mirna'), decoy),
                ('Other RNA (decoys)', 'right RNA type (%)', UP, 90, 70, '{:.0f}', right_class(decoy), decoy)]
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
                cells[(r, j)] = (fmt.format(v), bin_of(v, direction, good, poor) if direction else None)
    used = [sorted(p.path.name for p in c[7]) for c in cols]
    return routes, [c[:6] for c in cols], cells, used


ARROW = {UP: '↑', DOWN: '↓', None: ''}

NOTE = ('One row per route, one column per criterion; ↑ higher is better, ↓ lower is better. Each cell is the raw value, '
        'binned as good, intermediate or poor by the cut-offs below (good at or past the first, poor past the second). '
        'Synthetic: reads simulated from known miRNAs, with changes designed between two groups, averaged over the sets; '
        'real changes found are divided by the same test on the true counts. Equal-amount and set-ratio pools: synthetic miRNAs '
        'mixed in equal amounts or at set ratios; the reads that exactly match a pool sequence are the reference a route should '
        'reproduce, and the ratio error is divided by that reference\'s own. Mixture: two tissues mixed at known proportions. '
        'Knockout: miRNAs counted in cells where they were knocked out. No change: samples of one group split in two halves, '
        'where only about 5% of miRNAs should pass p < 0.05 by chance; only the share above 5% counts (edgeR\'s exact test). '
        'Other RNA: reads from known non-miRNA RNAs added to the synthetic data; right RNA type over the types a route reports. '
        'Cost is not binned. There is no overall score.')


def cutoffs(cols):
    """The printed cut-offs: '<label> ≥good / <poor' per binned column."""
    def one(c):
        group, label, direction, good, poor, fmt = c
        return '{} {} {}{} / {}{}'.format(group, label, '≥' if direction == UP else '≤', fmt.format(good),
                                           '<' if direction == UP else '>', fmt.format(poor))
    return '; '.join(one(c) for c in cols if c[2])
