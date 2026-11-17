"""Metrics for biological datasets: samples without a known composition, so methods are compared with each other.

Agreement is not accuracy: two methods can agree and both be wrong. What a biological dataset adds to the synthetic
and pool datasets is real sample composition (abundant non-miRNA RNAs, artefacts of the kit and the sequencer) and,
where a library is sequenced in several lanes, the sequencing noise floor.

All functions work on the common table of adapters.py: (method, sample, mature_name, mature_id, count) rows, counts
per mature miRNA accession. Counts per million (CPM) are relative to each method's miRNA reads in the sample.
"""
import collections
import json
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool import spearman  # noqa: E402

MIN_READS = 10          # a miRNA counts as detected in a sample with at least this many reads
TOP = 100               # top-N overlap between methods
LANE_MIN_MEAN = 20      # lane dispersion only for miRNAs with at least this mean count per lane
OVERDISPERSED = 3.0     # dispersion index above which a miRNA is listed


def counts_by(common):
    """{(method, sample): {mature_id: count}} from common-table rows with an accession."""
    out = collections.defaultdict(lambda: collections.defaultdict(float))
    for method, sample, _, mature_id, count in common:
        if mature_id:
            out[(method, sample)][mature_id] += float(count)
    return out


def cpm(counts):
    total = sum(counts.values())
    return {k: 1e6 * v / total for k, v in counts.items()} if total else {}


def method_summary(common, methods, samples):
    """Per method: median miRNA reads, median miRNAs detected per sample, miRNAs detected in every sample, median
    share of the miRNA reads in the 10 most abundant miRNAs."""
    by = counts_by(common)
    out = {}
    for m in methods:
        per = [by.get((m, s), {}) for s in samples]
        detected = [{k for k, v in c.items() if v >= MIN_READS} for c in per]
        top10 = [sum(sorted(c.values(), reverse=True)[:10]) / sum(c.values()) for c in per if sum(c.values())]
        out[m] = {'mirna_reads': round(statistics.median(sum(c.values()) for c in per)),
                  'detected': statistics.median(len(d) for d in detected),
                  'detected_all_samples': len(set.intersection(*detected)) if detected else 0,
                  'top10_share': round(statistics.median(top10), 4) if top10 else ''}
    return out


def agreement(common, methods, samples):
    """Per pair of methods (medians over samples): Spearman of log2(CPM + 1) over miRNAs either detects and over those
    both detect, median |log2 CPM ratio| over miRNAs both detect, overlap of the TOP most abundant, and miRNAs only one
    of them detects (only_a for the first method of the pair, only_b for the second)."""
    by = counts_by(common)
    out = {}
    for i, a in enumerate(methods):
        for b in methods[i + 1:]:
            rho, rho_both, ratio, overlap, only_a, only_b = [], [], [], [], [], []
            for s in samples:
                ca, cb = by.get((a, s), {}), by.get((b, s), {})
                pa, pb = cpm(ca), cpm(cb)
                da = {k for k, v in ca.items() if v >= MIN_READS}
                db = {k for k, v in cb.items() if v >= MIN_READS}
                either, both = sorted(da | db), sorted(da & db)
                if len(either) > 2:
                    rho.append(spearman([math.log2(pa.get(k, 0) + 1) for k in either], [math.log2(pb.get(k, 0) + 1) for k in either]))
                if len(both) > 2:
                    rho_both.append(spearman([math.log2(pa[k] + 1) for k in both], [math.log2(pb[k] + 1) for k in both]))
                if both:
                    ratio.append(statistics.median(abs(math.log2(pa[k] / pb[k])) for k in both))
                ta = set(sorted(ca, key=ca.get, reverse=True)[:TOP])
                tb = set(sorted(cb, key=cb.get, reverse=True)[:TOP])
                if ta and tb:
                    overlap.append(len(ta & tb) / min(len(ta), len(tb)))
                only_a.append(len(da - db))
                only_b.append(len(db - da))
            med = lambda v, n=4: round(statistics.median(v), n) if v else ''
            out[(a, b)] = {'spearman_either': med(rho), 'spearman_both': med(rho_both), 'median_abs_log2_ratio': med(ratio), 'top{}_overlap'.format(TOP): med(overlap),
                           'only_a': med(only_a, 1), 'only_b': med(only_b, 1)}
    return out


def blocked(common, methods, samples, accessions):
    """Median CPM per method of the miRNAs a blocker targets ({name: accession})."""
    by = counts_by(common)
    return {m: {name: round(statistics.median(cpm(by.get((m, s), {})).get(acc, 0) for s in samples), 1)
                for name, acc in accessions.items()} for m in methods}


def filtering(stage_jsons):
    """Read counts through the cutadapt route's three passes, from their JSON reports ({pass: path})."""
    read = {k: json.loads(Path(p).read_text())['read_counts'] for k, p in stage_jsons.items()}
    f = read['filter']['filtered']
    return {'reads_in': read['adapter']['input'], 'with_3p_adapter': read['adapter']['read1_with_adapter'],
            'discarded_5p_adapter': f['discard_trimmed'] or 0, 'too_short': f['too_short'] or 0,
            'too_many_n': f['too_many_n'] or 0, 'too_many_expected_errors': f['too_many_expected_errors'] or 0,
            'discarded_artefacts': read['artefacts']['filtered']['discard_trimmed'] or 0, 'kept': read['artefacts']['output']}


def first_base(fastq, matures):
    """Reads of the given mature miRNAs ({name: sequence}) that match it from the second base on, by first base:
    {name: {base: count}}, and the median quality of canonical and substituted first bases."""
    cores = collections.defaultdict(dict)
    for name, seq in matures.items():
        cores[len(seq) - 1][seq[1:]] = (name, seq[0])
    counts = {name: collections.Counter() for name in matures}
    quality = {'canonical': [], 'substituted': []}
    with open(fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                read = line.strip()
            elif i % 4 == 3:
                for length, table in cores.items():
                    hit = table.get(read[1:1 + length])
                    if hit:
                        name, canonical = hit
                        counts[name][read[0]] += 1
                        q = quality['canonical' if read[0] == canonical else 'substituted']
                        if len(q) < 200000:
                            q.append(ord(line[0]) - 33)
                        break
    return counts, {k: round(float(statistics.median(v)), 1) if v else '' for k, v in quality.items()}


def chi2_expectations(df, threshold=OVERDISPERSED):
    """Median of chi2(df)/df (Wilson-Hilferty) and P(chi2(df)/df > threshold): what pure sampling gives for the
    dispersion index of L = df + 1 lanes. With 4 lanes the median is 0.79, not 1 (the mean is 1)."""
    median = (1 - 2 / (9 * df)) ** 3
    z = ((threshold) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    return median, 0.5 * math.erfc(z / math.sqrt(2))


def lane_replicates(common, methods, lanes):
    """Lane-to-lane variation per method, for samples sequenced in several lanes ({sample: [lane sample, ...]}).

    For each method, sample and miRNA with at least LANE_MIN_MEAN reads per lane on average, the dispersion index
    sum((c_i - T_i p)^2 / (T_i p (1 - p))) / (L - 1), with c_i the lane counts, T_i the method's miRNA reads in lane i
    and p = sum(c) / sum(T): about 1 when lanes differ only by sampling (multinomial within the method's miRNA reads,
    which for a rare miRNA is Poisson; an abundant one varies less). Also the CV of CPM across lanes, binned by mean
    count (powers of two), to find the count from which lanes agree within 10%."""
    by = counts_by(common)
    out = {}
    for m in methods:
        dispersion, cv_bins, high, dfs = [], collections.defaultdict(list), collections.Counter(), []
        for sample, members in lanes.items():
            per = [by.get((m, lane), {}) for lane in members]
            totals = [sum(c.values()) for c in per]
            if len(per) < 2 or not all(totals):
                continue
            for mid in set().union(*per):
                c = [p.get(mid, 0) for p in per]
                mean = sum(c) / len(c)
                if mean < 1:
                    continue
                rates = [ci / ti * 1e6 for ci, ti in zip(c, totals)]
                mu = statistics.mean(rates)
                if mu:
                    cv_bins[int(math.log2(mean))].append(statistics.pstdev(rates) / mu)
                if mean >= LANE_MIN_MEAN:
                    p = sum(c) / sum(totals)
                    d = sum((ci - ti * p) ** 2 / (ti * p * (1 - p)) for ci, ti in zip(c, totals)) / (len(c) - 1)
                    dispersion.append(d)
                    dfs.append(len(c) - 1)
                    if d > OVERDISPERSED:
                        high[mid] += 1
        threshold = next((2 ** b for b in sorted(cv_bins) if statistics.median(cv_bins[b]) < 0.1), '')
        df = statistics.mode(dfs) if dfs else 0
        exp_median, exp_tail = chi2_expectations(df) if df else ('', 0)
        out[m] = {'dispersion_mean': round(statistics.mean(dispersion), 3) if dispersion else '',
                  'dispersion_median': round(statistics.median(dispersion), 3) if dispersion else '',
                  'expected_median': round(exp_median, 3) if df else '',
                  'mirna_sample_pairs': len(dispersion), 'overdispersed': sum(1 for d in dispersion if d > OVERDISPERSED),
                  'expected_overdispersed': round(exp_tail * len(dispersion)) if df else '',
                  'reads_for_cv_below_10pct': threshold, 'overdispersed_mirnas': [k for k, _ in high.most_common(5)]}
    return out
