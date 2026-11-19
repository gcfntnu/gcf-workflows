"""Bootstrap confidence intervals for an evaluation's metrics and for differences between methods (#413).

A metric is computed over scoring units (truth miRNAs, pool sequences) and, where a dataset has several, samples. One
resample draws the units and the samples with replacement (samples within their group, e.g. the ratiometric pools A
and B) and recomputes every row's metrics on that same draw, so the difference between two methods is paired: its
interval reflects how the two move together, not each one's spread alone.

`intervals` returns, per metric and row, the point value and the percentile interval (2.5% and 97.5% by default) of
the resampled values, and per metric and pair of methods (the oracle left out) the point difference (first minus
second) and its interval. A resample that cannot be scored (the metric function raises ValueError, e.g. a ratiometric
draw without a sequence designed 1:1) is left out and counted. NaN values (e.g. no detected entry in a draw) are left
out of a metric's interval; an interval from fewer than half of the resamples is not given.
"""
import math
import random

RESAMPLES = 1000
SEED = 1
LEVEL = 0.95


def resample(rng, items):
    """A draw of len(items) items with replacement."""
    return [items[rng.randrange(len(items))] for _ in items]


def resample_groups(rng, pairs):
    """(item, group) pairs drawn with replacement within each group, each group keeping its size."""
    groups = {}
    for item, group in pairs:
        groups.setdefault(group, []).append(item)
    return [(item, group) for group, items in groups.items() for item in resample(rng, items)]


def percentile(values, q):
    """The q-quantile (0-1) of values, linear between order statistics."""
    values = sorted(values)
    position = q * (len(values) - 1)
    low = math.floor(position)
    high = min(low + 1, len(values) - 1)
    return values[low] + (values[high] - values[low]) * (position - low)


def tidy(value):
    """A value rounded to 4 decimals, an integer when it is whole (e.g. read counts)."""
    value = round(value, 4)
    return int(value) if value == int(value) else value


def number(value):
    """A metric value as a float, or None when it is empty, not a number or NaN."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def intervals(point, compute, draw, rows, keys, resamples=RESAMPLES, seed=SEED, level=LEVEL):
    """(ci_rows, difference_rows, skipped) for `point` ({row: {metric: value}}, the evaluation's values).

    compute(*draw(rng)) gives {row: {metric: value}} for one resample. ci_rows: [metric, row, value, low, high];
    difference_rows: [metric, row_a, row_b, difference, low, high], for pairs of rows other than 'oracle'."""
    rng = random.Random(seed)
    draws, skipped = [], 0
    for _ in range(resamples):
        try:
            draws.append(compute(*draw(rng)))
        except ValueError:
            skipped += 1
    lo, hi = (1 - level) / 2, 1 - (1 - level) / 2

    def bounds(values):
        values = [v for v in values if v is not None]
        if len(values) < max(1, resamples // 2):
            return '', ''
        return tidy(percentile(values, lo)), tidy(percentile(values, hi))

    ci_rows = []
    for key in keys:
        for row in rows:
            value = point.get(row, {}).get(key, '')
            if number(value) is None:
                continue
            ci_rows.append([key, row, value, *bounds([number(d.get(row, {}).get(key)) for d in draws])])
    methods = [row for row in rows if row != 'oracle']
    difference_rows = []
    for key in keys:
        for i, a in enumerate(methods):
            for b in methods[i + 1:]:
                pa, pb = number(point.get(a, {}).get(key)), number(point.get(b, {}).get(key))
                if pa is None or pb is None:
                    continue
                values = []
                for d in draws:
                    va, vb = number(d.get(a, {}).get(key)), number(d.get(b, {}).get(key))
                    values.append(None if va is None or vb is None else va - vb)
                difference_rows.append([key, a, b, tidy(pa - pb), *bounds(values)])
    return ci_rows, difference_rows, skipped


CI_HEADER = ['metric', 'row', 'value', 'low', 'high']
DIFFERENCE_HEADER = ['metric', 'row_a', 'row_b', 'difference', 'low', 'high']


def table_text(header, rows):
    return '\n'.join('\t'.join(map(str, row)) for row in [header, *rows]) + '\n'


def excludes_zero(row):
    low, high = number(row[4]), number(row[5])
    return low is not None and high is not None and (low > 0 or high < 0)


def report_lines(ci_rows, difference_rows, rows, keys, units, resamples, seed, skipped, level=LEVEL):
    """Markdown for report.md: the intervals per metric and row, and the differences whose interval excludes 0."""
    if not ci_rows:
        return []
    ci = {(m, r): (v, lo, hi) for m, r, v, lo, hi in ci_rows}
    shown = [r for r in rows if any((k, r) in ci for k in keys)]
    lines = ['', '## Uncertainty', '',
             '{:g}% bootstrap intervals from {} resamples (seed {}) of {}{}; the interval of a difference between two '
             'methods is paired (both computed on the same resample).'.format(
                 100 * level, resamples, seed, units, ', {} resamples left out (not scorable)'.format(skipped) if skipped else ''), '',
             '| metric | ' + ' | '.join(shown) + ' |', '|---' * (len(shown) + 1) + '|']
    for key in keys:
        cells = []
        for row in shown:
            value, low, high = ci.get((key, row), ('', '', ''))
            cells.append('' if value == '' else '{} ({}–{})'.format(value, low, high) if low != '' else str(value))
        if any(cells):
            lines.append('| {} | '.format(key) + ' | '.join(cells) + ' |')
    clear = [row for row in difference_rows if excludes_zero(row)]
    lines += ['', 'Differences between methods whose interval excludes 0 ({} of {}):'.format(len(clear), len(difference_rows)), '']
    lines += ['- {}: {} − {} = {} ({}–{})'.format(*row) for row in clear] or ['- none']
    return lines
