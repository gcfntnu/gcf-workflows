"""Metrics for biological datasets: samples without a known composition, so methods are compared with each other.

Agreement is not accuracy: two methods can agree and both be wrong. What a biological dataset adds to the synthetic
and pool datasets is real sample composition (abundant non-miRNA RNAs, artefacts of the kit and the sequencer) and,
where a library is sequenced in several lanes, the sequencing noise floor.

All functions work on the common table of adapters.py: (method, sample, mature_name, mature_id, count) rows, counts
per mature miRNA accession. Counts per million (CPM) are relative to each method's miRNA reads in the sample.
"""
import collections
import functools
import importlib.util
import json
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool import spearman  # noqa: E402

MERGE_CUTADAPT = Path(__file__).resolve().parents[1] / 'rules/filter/scripts/merge_cutadapt.py'


@functools.lru_cache(maxsize=None)
def merge_cutadapt():
    """The workflow's merge of the cutadapt pass reports, loaded by path when needed (rules/filter/scripts stays off
    sys.path, and the rest of this module doesn't depend on it)."""
    spec = importlib.util.spec_from_file_location('merge_cutadapt', MERGE_CUTADAPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

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


def collapse(common, id_to_seq):
    """Common-table rows with each accession replaced by the first (sorted) accession of the same mature sequence, and
    the groups of more than one. Methods credit the reads of identical sequences to different members (miRge3 to its
    merge group's first name, UNITAS and mirflowz to one or the other or both: hsa-miR-199a-3p and hsa-miR-199b-3p), so
    comparing methods needs one accession per sequence."""
    by_seq = collections.defaultdict(list)
    for accession, seq in id_to_seq.items():
        by_seq[seq].append(accession)
    alias, groups = {}, []
    for accessions in by_seq.values():
        accessions.sort()
        if len(accessions) > 1:
            groups.append(accessions)
        alias.update((a, accessions[0]) for a in accessions)
    return [(m, s, n, alias.get(i, i), c) for m, s, n, i, c in common], sorted(groups)


def unplaced(common, methods, samples, top=10):
    """Per method: the median share of its miRNA reads on names without a miRBase accession, which every other metric
    leaves out, and the `top` such names by reads over all samples."""
    total, missing = collections.defaultdict(float), collections.defaultdict(float)
    names = {m: collections.Counter() for m in methods}
    for method, sample, name, mature_id, count in common:
        total[(method, sample)] += float(count)
        if not mature_id and method in names:
            missing[(method, sample)] += float(count)
            names[method][name] += float(count)
    share = {m: [missing[(m, s)] / total[(m, s)] for s in samples if total[(m, s)]] for m in methods}
    return {m: {'share': round(statistics.median(share[m]), 4) if share[m] else '',
                'names': sorted(names[m].items(), key=lambda kv: (-kv[1], kv[0]))[:top]} for m in methods}


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


def sample_stats(common, raw, methods, samples):
    """{(method, sample): per-sample values behind method_summary and unplaced}: miRNA reads, detected miRNAs (set),
    share of the 10 most abundant (None without reads), share on names without an accession (None without reads)."""
    by = counts_by(common)
    total, missing = collections.defaultdict(float), collections.defaultdict(float)
    for method, sample, _, mature_id, count in raw:
        total[(method, sample)] += float(count)
        if not mature_id:
            missing[(method, sample)] += float(count)
    stats = {}
    for m in methods:
        for s in set(samples):
            c = by.get((m, s), {})
            reads = sum(c.values())
            stats[(m, s)] = {'reads': reads, 'detected': frozenset(k for k, v in c.items() if v >= MIN_READS),
                             'top10': sum(sorted(c.values(), reverse=True)[:10]) / reads if reads else None,
                             'unplaced': missing[(m, s)] / total[(m, s)] if total[(m, s)] else None}
    return stats


def summary_of(stats, methods, samples, all_samples=True):
    """method_summary with unplaced's share, from sample_stats, for any list of samples (a bootstrap resample may
    repeat one): the same values as method_summary and unplaced over the same samples. `all_samples=False` leaves out
    detected_all_samples (a resample repeating a sample biases it, so the bootstrap does not use it)."""
    out = {}
    for m in methods:
        per = [stats[(m, s)] for s in samples]
        top10 = [p['top10'] for p in per if p['top10'] is not None]
        share = [p['unplaced'] for p in per if p['unplaced'] is not None]
        out[m] = {'mirna_reads': round(statistics.median(p['reads'] for p in per)),
                  'detected': statistics.median(len(p['detected']) for p in per),
                  **({'detected_all_samples': len(frozenset.intersection(*(p['detected'] for p in per))) if per else 0}
                     if all_samples else {}),
                  'top10_share': round(statistics.median(top10), 4) if top10 else '',
                  'unplaced': round(statistics.median(share), 4) if share else ''}
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
                # NaN (a method without reads in the sample, so constant) would make the median depend on sample order
                r = spearman([math.log2(pa.get(k, 0) + 1) for k in either], [math.log2(pb.get(k, 0) + 1) for k in either])
                rho += [] if math.isnan(r) else [r]
                r = spearman([math.log2(pa[k] + 1) for k in both], [math.log2(pb[k] + 1) for k in both])
                rho_both += [] if math.isnan(r) else [r]
                if both:
                    ratio.append(statistics.median(abs(math.log2(pa[k] / pb[k])) for k in both))
                ta = set(sorted((k for k in ca if ca[k] > 0), key=lambda k: (-ca[k], k))[:TOP])
                tb = set(sorted((k for k in cb if cb[k] > 0), key=lambda k: (-cb[k], k))[:TOP])
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


CUTADAPT_PASSES = ('adapter', 'filter', 'artefacts')  # merge_cutadapt.PASSES; filtering() fails if they differ


def filtering(stage_jsons):
    """Read counts through the cutadapt route's three passes ({pass: JSON report path}), as the workflow's merged
    report gives them (rules/filter/scripts/merge_cutadapt.py, which fails unless each pass read what the one before
    wrote): the artefact pass split into 3' artefacts and reagent oligo reads."""
    merger = merge_cutadapt()
    if tuple(stage_jsons) != tuple(merger.PASSES):
        raise ValueError('cutadapt passes {} differ from the workflow\'s {}'.format(tuple(stage_jsons), merger.PASSES))
    paths = {k: str(p) for k, p in stage_jsons.items()}
    reports = {k: json.loads(Path(p).read_text()) for k, p in paths.items()}
    try:
        merged = merger.merge(reports['adapter'], reports['filter'], reports['artefacts'], 'sample', paths)['read_counts']
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError('malformed cutadapt reports ({}: {}): {}'.format(type(error).__name__, error, ', '.join(paths.values())))
    f = merged['filtered']
    return {'reads_in': merged['input'], 'with_3p_adapter': merged['read1_with_adapter'],
            'discarded_5p_adapter': f['discard_trimmed'] or 0, 'too_short': f['too_short'] or 0,
            'too_many_n': f['too_many_n'] or 0, 'too_many_expected_errors': f['too_many_expected_errors'] or 0,
            'discarded_3p_artefacts': f['gcf_three_prime_artefacts'], 'discarded_reagent_oligos': f['gcf_reagent_oligos'],
            'kept': merged['output']}


def first_base(fastq, matures):
    """Reads of the given mature miRNAs ({name: sequence}) that match it from the second base on, by first base:
    {name: {base: count}}, and the median quality of canonical and substituted first bases. Some miRNAs differ only
    at the first base (hsa-miR-17-5p and hsa-miR-106a-5p): a read with either first base is that miRNA's, one with
    another first base is left out, as it cannot be placed."""
    cores = collections.defaultdict(lambda: collections.defaultdict(dict))  # length -> sequence from base 2 -> {first base: name}
    for name, seq in matures.items():
        cores[len(seq) - 1][seq[1:]][seq[0]] = name
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
                        if read[0] not in hit and len(hit) > 1:
                            break
                        name = hit.get(read[0]) or next(iter(hit.values()))
                        counts[name][read[0]] += 1
                        q = quality['canonical' if read[0] in hit else 'substituted']
                        if len(q) < 200000:
                            q.append(ord(line[0]) - 33)
                        break
    return counts, {k: round(float(statistics.median(v)), 1) if v else '' for k, v in quality.items()}


def chi2_cdf(x, df):
    """P(chi2(df)/df <= x), Wilson-Hilferty approximation."""
    z = (x ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    return 0.5 * math.erfc(-z / math.sqrt(2))


def chi2_expectations(dfs, threshold=OVERDISPERSED):
    """What pure sampling gives for dispersion indices with the given degrees of freedom (lanes - 1, one per
    miRNA-sample pair): their median, and how many are expected above threshold. With 4 lanes the median is 0.79, not
    1 (the mean is 1)."""
    lo, hi = 0.0, 50.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if sum(chi2_cdf(mid, df) for df in dfs) < len(dfs) / 2 else (lo, mid)
    return lo, sum(1 - chi2_cdf(threshold, df) for df in dfs)


ABUNDANCE_BINS = (20, 100, 1000, 10000)    # mean reads per lane: overdispersion share per bin


def lane_replicates(common, methods, lanes):
    """Lane-to-lane variation per method, for samples sequenced in several lanes ({sample: [lane sample, ...]}).

    For each method, sample and miRNA with at least LANE_MIN_MEAN reads per lane on average, the dispersion index
    sum((c_i - T_i p)^2 / (T_i p (1 - p))) / (L - 1), with c_i the lane counts, T_i the method's miRNA reads in lane i
    and p = sum(c) / sum(T): about 1 when lanes differ only by sampling (multinomial within the method's miRNA reads,
    which for a rare miRNA is Poisson; an abundant one varies less). Also the CV of CPM across lanes, binned by mean
    count (powers of two), to find the count from which lanes agree within 10% (median sample CV). All methods count
    the same reads, so a lane effect (e.g. a sequencing error rate that differs by lane) raises the index for every
    method, and more for abundant miRNAs, whose sampling noise is small; the share above OVERDISPERSED per bin of mean
    reads (ABUNDANCE_BINS) shows it."""
    by = counts_by(common)
    out = {}
    for m in methods:
        dispersion, cv_bins, high, dfs, bins = [], collections.defaultdict(list), collections.Counter(), [], collections.defaultdict(list)
        for sample, members in lanes.items():
            per = [by.get((m, lane), {}) for lane in members]
            totals = [sum(c.values()) for c in per]
            if len(per) < 2 or not all(totals):
                continue
            for mid in sorted(set().union(*per)):
                c = [p.get(mid, 0) for p in per]
                mean = sum(c) / len(c)
                if mean < 1:
                    continue
                rates = [ci / ti * 1e6 for ci, ti in zip(c, totals)]
                mu = statistics.mean(rates)
                if mu:
                    cv_bins[int(math.log2(mean))].append(statistics.stdev(rates) / mu)
                if mean >= LANE_MIN_MEAN:
                    p = sum(c) / sum(totals)
                    d = sum((ci - ti * p) ** 2 / (ti * p * (1 - p)) for ci, ti in zip(c, totals)) / (len(c) - 1)
                    dispersion.append(d)
                    dfs.append(len(c) - 1)
                    bins[max(b for b in ABUNDANCE_BINS if b <= mean)].append(d > OVERDISPERSED)
                    if d > OVERDISPERSED:
                        high[mid] += 1
        b = next((b for b in sorted(cv_bins) if statistics.median(cv_bins[b]) < 0.1), None)
        exp_median, exp_tail = chi2_expectations(dfs) if dfs else ('', '')
        out[m] = {'dispersion_mean': round(statistics.mean(dispersion), 3) if dispersion else '',
                  'dispersion_median': round(statistics.median(dispersion), 3) if dispersion else '',
                  'expected_median': round(exp_median, 3) if dfs else '',
                  'mirna_sample_pairs': len(dispersion), 'overdispersed': sum(1 for d in dispersion if d > OVERDISPERSED),
                  'expected_overdispersed': round(exp_tail) if dfs else '',
                  'reads_for_cv_below_10pct': '{}-{}'.format(2 ** b, 2 ** (b + 1) - 1) if b is not None else '',
                  'overdispersed_by_mean': {k: (sum(v), len(v)) for k, v in sorted(bins.items())},
                  'overdispersed_mirnas': [k for k, _ in sorted(high.items(), key=lambda kv: (-kv[1], kv[0]))[:5]]}
    return out


KNOCKOUT_MIN_READS = 10  # a miRNA gets a knockout-to-control ratio with at least this many reads on average in the controls
SURVIVING = 0.25         # a miRNA survives a knockout that keeps more than this share of its control level (#452)


def knockout(common, methods, kept, units, control):
    """Knockout cells against their controls (#452): miRNA reads per million filtered reads (`kept`: {sample: the reads
    the filter kept}), not per million miRNA reads, which the knockout changes. `units`: [(sample, group, level)]; each
    knockout sample is compared with the `control` samples of its level (e.g. the parental cells sequenced in the same
    set). Per method:
    - 'control_rpm': the controls' mean miRNA reads per million filtered reads;
    - 'share': {group: the mean over its samples of their miRNA RPM over the mean RPM of their level's controls}
      ('' when those controls have no miRNA reads);
    - 'ratios': {group: {accession: (control reads, knockout reads, ratio)}} for the miRNAs with at least
      KNOCKOUT_MIN_READS reads on average in the controls of every level of the group's samples: the mean reads of
      those controls (levels averaged) and of the group's samples, and the mean over the group's samples of the
      miRNA's RPM over that of its level's controls."""
    by = counts_by(common)
    controls = collections.defaultdict(list)
    for sample, group, level in units:
        if group == control:
            controls[level].append(sample)
    knockouts = collections.defaultdict(list)
    for sample, group, level in units:
        if group != control:
            knockouts[group].append((sample, level))
    mean = lambda values: sum(values) / len(values)
    out = {}
    for m in methods:
        counts = {s: by.get((m, s), {}) for s, _, _ in units}
        rpm = {s: 1e6 * sum(c.values()) / kept[s] for s, c in counts.items()}
        level_rpm = {level: mean([rpm[s] for s in samples]) for level, samples in controls.items()}
        level_reads = {level: collections.Counter() for level in controls}
        level_mirna = {level: collections.Counter() for level in controls}
        for level, samples in controls.items():
            for s in samples:
                for i, v in counts[s].items():
                    level_reads[level][i] += v / len(samples)
                    level_mirna[level][i] += 1e6 * v / kept[s] / len(samples)
        share, ratios = {}, {}
        for group, members in knockouts.items():
            fractions = [rpm[s] / level_rpm[level] for s, level in members if level_rpm[level]]
            share[group] = mean(fractions) if len(fractions) == len(members) else ''
            levels = sorted({level for _, level in members}, key=str)
            eligible = set.intersection(*({i for i, v in level_reads[level].items() if v >= KNOCKOUT_MIN_READS} for level in levels))
            ratios[group] = {i: (mean([level_reads[level][i] for level in levels]),
                                 mean([counts[s].get(i, 0.0) for s, _ in members]),
                                 mean([1e6 * counts[s].get(i, 0.0) / kept[s] / level_mirna[level][i] for s, level in members]))
                             for i in sorted(eligible)}
        out[m] = {'control_rpm': mean([rpm[s] for samples in controls.values() for s in samples]), 'share': share, 'ratios': ratios}
    return out


def survivors(result, methods, group):
    """The miRNAs surviving `group` (ratio above SURVIVING) per method in a knockout() result: {method: {'tested',
    'surviving' (a set), 'only' (those no other method lists), 'only_reads' (their mean reads per knockout sample)}},
    and those every method lists."""
    surviving = {m: {i for i, (_, _, ratio) in result[m]['ratios'].get(group, {}).items() if ratio > SURVIVING} for m in methods}
    out = {}
    for m in methods:
        others = set().union(*(surviving[o] for o in methods if o != m))
        only = surviving[m] - others
        out[m] = {'tested': len(result[m]['ratios'].get(group, {})), 'surviving': surviving[m], 'only': only,
                  'only_reads': sum((result[m]['ratios'][group][i][1] for i in only), 0.0)}
    return out, set.intersection(*surviving.values()) if methods else set()
