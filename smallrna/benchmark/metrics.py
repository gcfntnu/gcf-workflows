"""Benchmark metrics on the common table, the same for every method.

Truth of type `counts` (synthetic data) gives reads per mature miRNA and sample; per method:

- `recovered`: share of the true miRNA reads counted on the right miRNA (sum of min(observed, truth) / sum of truth);
- `exact`: share of true miRNAs counted exactly (within EXACT_TOLERANCE, as methods that divide multimapping reads
  give fractional counts whose sums carry floating-point noise);
- `misassigned`: reads counted on miRNAs that are not in the truth;
- `unplaced`: reads on names that could not be matched to a miRBase accession.

Pools of known composition are scored in pool.py.
"""

EXACT_TOLERANCE = 1e-6


def per_mirna(common, truth, methods=None):
    """{(method, sample, mature_id): (truth, observed)} for every truth miRNA and every observed accession."""
    if methods is None:
        methods = sorted({row[0] for row in common})
    observed = {}
    for method, sample, _, mature_id, count in common:
        if mature_id:
            observed[(method, sample, mature_id)] = observed.get((method, sample, mature_id), 0) + count
    table = {}
    for (sample, mature_id), count in truth.items():
        for method in methods:
            table[(method, sample, mature_id)] = (count, observed.get((method, sample, mature_id), 0))
    for key, count in observed.items():
        table.setdefault(key, (0, count))
    return table


def count_metrics(common, truth, methods=None, keys=None):
    """truth: {(sample, mature_id): reads}. Returns {method: {metric: value}}.

    `keys`: the truth pairs that `recovered` and `exact` are computed over instead of all of them (a bootstrap
    resample, with repeats); `misassigned` and `unplaced` are totals over the whole table either way."""
    if methods is None:
        methods = sorted({row[0] for row in common})
    table = per_mirna(common, truth, methods)
    result = {}
    for method in methods:
        rows = [(t, o) for (m, _, _), (t, o) in table.items() if m == method]
        if keys is None:
            true_rows = [(t, o) for t, o in rows if t > 0]
        else:
            true_rows = [table[(method, *key)] for key in keys if truth.get(key, 0) > 0]
        total = sum(t for t, _ in true_rows)
        result[method] = {
            'recovered': round(sum(min(t, o) for t, o in true_rows) / total, 4) if total else float('nan'),
            'exact': round(sum(1 for t, o in true_rows if abs(t - o) <= EXACT_TOLERANCE) / len(true_rows), 4) if true_rows else float('nan'),
            'misassigned': sum(o for t, o in rows if t == 0),
            'unplaced': sum(row[4] for row in common if row[0] == method and not row[3]),
        }
    return result
