"""Benchmark metrics on the common table, the same for every method.

Truth of type `counts` (synthetic data) gives reads per mature miRNA and sample; per method:

- `recovered`: share of the true miRNA reads counted on the right miRNA (sum of min(observed, truth) / sum of truth);
- `exact`: share of true miRNAs counted exactly;
- `misassigned`: reads counted on miRNAs that are not in the truth;
- `unplaced`: reads on names that could not be matched to a miRBase accession.

Truth of type `equimolar` (a pool of miRNAs at equal amounts) gives only the pool's accessions; per method:

- `detected`: share of pool miRNAs with any count;
- `cv`: coefficient of variation of counts across the pool miRNAs (0 when equal, as intended);
- `within_2fold`: share of pool miRNAs within two-fold of the pool median.
"""
import statistics


def per_mirna(common, truth):
    """{(method, sample, mature_id): (truth, observed)} for every truth miRNA and every observed accession."""
    observed = {}
    for method, sample, _, mature_id, count in common:
        if mature_id:
            observed[(method, sample, mature_id)] = observed.get((method, sample, mature_id), 0) + count
    methods = sorted({row[0] for row in common})
    table = {}
    for (sample, mature_id), count in truth.items():
        for method in methods:
            table[(method, sample, mature_id)] = (count, observed.get((method, sample, mature_id), 0))
    for key, count in observed.items():
        table.setdefault(key, (0, count))
    return table


def count_metrics(common, truth):
    """truth: {(sample, mature_id): reads}. Returns {method: {metric: value}}."""
    table = per_mirna(common, truth)
    result = {}
    for method in sorted({row[0] for row in common}):
        rows = [(t, o) for (m, _, _), (t, o) in table.items() if m == method]
        true_rows = [(t, o) for t, o in rows if t > 0]
        total = sum(t for t, _ in true_rows)
        result[method] = {
            'recovered': round(sum(min(t, o) for t, o in true_rows) / total, 4) if total else float('nan'),
            'exact': round(sum(1 for t, o in true_rows if t == o) / len(true_rows), 4) if true_rows else float('nan'),
            'misassigned': sum(o for t, o in rows if t == 0),
            'unplaced': sum(row[4] for row in common if row[0] == method and not row[3]),
        }
    return result


def equimolar_metrics(common, pool_ids):
    """pool_ids: the pool's accessions. Returns {method: {metric: value}} over samples pooled per miRNA."""
    result = {}
    for method in sorted({row[0] for row in common}):
        samples = sorted({row[1] for row in common if row[0] == method})
        values = []
        for sample in samples:
            counts = {}
            for m, s, _, mature_id, count in common:
                if m == method and s == sample and mature_id in pool_ids:
                    counts[mature_id] = counts.get(mature_id, 0) + count
            values.append([counts.get(i, 0) for i in sorted(pool_ids)])
        per_mirna = [statistics.mean(column) for column in zip(*values)] if values else []
        if not per_mirna:
            continue
        mean, median = statistics.mean(per_mirna), statistics.median(per_mirna)
        result[method] = {
            'detected': round(sum(1 for v in per_mirna if v > 0) / len(per_mirna), 4),
            'cv': round(statistics.pstdev(per_mirna) / mean, 4) if mean else float('nan'),
            'within_2fold': round(sum(1 for v in per_mirna if median and median / 2 <= v <= median * 2) / len(per_mirna), 4),
        }
    return result
