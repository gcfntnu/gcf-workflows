"""Synthetic RNA pools of known composition: truth groups, the oracle count and the pool metrics.

Composition tables are those of the Giraldez et al. 2018 study (rspengle/CrossU01_exRNA_Manuscript2017): the equimolar
pool (`Sequence.ID`, `Sequence.type`, `Modification`, `sequence`) and the ratiometric pools A and B (`ratio.seqID`,
`Length`, `Type`, `sequence`, `Modification`, `X.SAMPLE.A.`, `X.SAMPLE.B.` with amounts such as `10x` and `1x`).

**Truth groups.** As in the study, only sequences of 16-25 nt whose only modification is a 5' phosphate, and whose
sequence is unique in the pool, are analysed. They are grouped against the organism's miRBase mature sequences:

- `exact`: the sequence of one or more mature miRNAs (the accessions are kept; identical mature sequences of
  paralogs give several);
- `shifted`: a mature miRNA whose 5' and 3' ends each differ by at most 3 nt (longer or shorter, e.g. an older
  miRBase boundary), which methods should report as an isomiR of that miRNA;
- `foreign`: neither (other species, viral, revised or non-miRNA sequences). Reads from these are expected on no
  organism miRNA unless a method assigns them to a similar one;
- `excluded`: outside the analysed set.

**Oracle.** The workflow's filtered reads counted by exact match to each pool sequence: what is in a sample after the
library preparation's bias. Methods are scored against the oracle on the `exact` group, so a kit's per-sequence bias
does not hide the differences between methods, and on the `shifted` group, by kind of shift (SHIFT_KINDS), for how
they count isomiRs. The same count after each pass of the filter route shows what the route removes.
"""
import collections
import csv
import gzip
import math
import statistics

GROUPS = ('exact', 'shifted', 'foreign', 'excluded')
MAX_SHIFT = 3
POOL_COLUMNS = ('pool_id', 'sequence', 'group', 'mature_ids', 'amount_a', 'amount_b')


def read_composition(path):
    """[{pool_id, sequence, analysed, amount_a, amount_b}] from either composition table."""
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    if not rows:
        raise ValueError('{}: empty composition table'.format(path))
    ratiometric = 'ratio.seqID' in rows[0]
    entries = []
    for row in rows:
        sequence = row['sequence'].strip().upper().replace('U', 'T')
        entries.append({'pool_id': row['ratio.seqID' if ratiometric else 'Sequence.ID'], 'sequence': sequence,
                        'phosphate_only': row['Modification'].strip() == "5'-phosphorylation",
                        'amount_a': parse_amount(row['X.SAMPLE.A.']) if ratiometric else 1.0,
                        'amount_b': parse_amount(row['X.SAMPLE.B.']) if ratiometric else 1.0})
    return entries


def parse_amount(text):
    """'10x' -> 10.0, '1.5x' -> 1.5."""
    text = text.strip().lower()
    if not text.endswith('x'):
        raise ValueError('pool amount {!r} is not of the form <number>x'.format(text))
    return float(text[:-1])


def shift(sequence, mature, min_overlap=15):
    """(start, end) when sequence is mature with its 5' and 3' ends each moved by at most MAX_SHIFT nt (either way), the
    overlap identical and at least min_overlap nt long (e.g. an older miRBase boundary), else None. start is the
    position of sequence's first base on mature (> 0: 5' end shorter), end the difference of the 3' ends (> 0: longer);
    the smallest shift when several fit."""
    if sequence == mature:
        return None
    fits = []
    for start in range(-MAX_SHIFT, MAX_SHIFT + 1):
        end = start + len(sequence) - len(mature)
        if abs(end) > MAX_SHIFT:
            continue
        lo, hi = max(start, 0), min(start + len(sequence), len(mature))
        if hi - lo >= min_overlap and sequence[lo - start:hi - start] == mature[lo:hi]:
            fits.append((abs(start) + abs(end), start, end))
    return min(fits)[1:] if fits else None


def shifted(sequence, mature, min_overlap=15):
    return shift(sequence, mature, min_overlap) is not None


# Kinds of shifted entries: the 3' end shorter or longer with the 5' end in place (the seed unchanged), or the 5' end
# moved (another seed; methods may treat it differently); mixed when the entry's miRNAs give different kinds.
SHIFT_KINDS = ('3p_shorter', '3p_longer', '5p_moved', 'mixed')


def shift_kinds(pool, mature):
    """{pool_id: kind} for shifted entries, against their miRNAs' mature sequences ({accession: sequence})."""
    out = {}
    for e in pool:
        if e['group'] != 'shifted':
            continue
        kinds = set()
        for accession in e['mature_ids']:
            start, end = shift(e['sequence'], mature[accession])
            kinds.add('5p_moved' if start else '3p_shorter' if end < 0 else '3p_longer')
        out[e['pool_id']] = kinds.pop() if len(kinds) == 1 else 'mixed'
    return out


def truth_groups(entries, mature):
    """Rows of POOL_COLUMNS: entries grouped against mature ({name: (accession, sequence)})."""
    by_sequence = {}
    for accession, sequence in mature.values():
        by_sequence.setdefault(sequence, set()).add(accession)
    counts = {}
    for entry in entries:
        counts[entry['sequence']] = counts.get(entry['sequence'], 0) + 1
    rows = []
    for entry in entries:
        sequence = entry['sequence']
        accessions = set()
        if not (16 <= len(sequence) <= 25 and entry['phosphate_only'] and counts[sequence] == 1):
            group = 'excluded'
        elif sequence in by_sequence:
            group, accessions = 'exact', by_sequence[sequence]
        else:
            accessions = {a for s, ids in by_sequence.items() if shifted(sequence, s) for a in ids}
            group = 'shifted' if accessions else 'foreign'
        rows.append((entry['pool_id'], sequence, group, ','.join(sorted(accessions)),
                     '{:g}'.format(entry['amount_a']), '{:g}'.format(entry['amount_b'])))
    return rows


def write_pool(path, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(POOL_COLUMNS) + '\n' + ''.join('\t'.join(row) + '\n' for row in rows))


def read_pool(path):
    with open(path, newline='') as fh:
        return [dict(row, amount_a=float(row['amount_a']), amount_b=float(row['amount_b']),
                     mature_ids=[a for a in row['mature_ids'].split(',') if a])
                for row in csv.DictReader(fh, delimiter='\t')]


def oracle(fastq, sequences, front=0, tail=0):
    """{sequence: reads} for reads in fastq (plain or gzipped) that equal one of sequences exactly, after cutting front
    and tail bases (the kit's random ends, for reads that still carry them)."""
    wanted = set(sequences)
    counts = dict.fromkeys(wanted, 0)
    with open(fastq, 'rb') as probe:
        gzipped = probe.read(2) == b'\x1f\x8b'
    with (gzip.open(fastq, 'rt') if gzipped else open(fastq)) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                read = line.strip()
                if front or tail:
                    read = read[front:len(read) - tail]
                if read in wanted:
                    counts[read] += 1
    return counts


FILTER_MIN_READS = 20   # entries listed for their filter losses need at least this many reads after the first stage
FILTER_LOSS = 0.1       # ... and to lose more than this share of them in the later passes


def filter_losses(stage_counts, pool, stages):
    """What each pass of the filter route removes from the pool sequences. stage_counts: {sample: {stage: {sequence:
    reads}}}, the reads after each pass (stages, in order; the first is the adapter pass, or the input reads of a route
    without one, the last the filtered reads).

    Returns per truth group (exact, shifted, foreign) {'sequences': n, stage: reads summed over samples}, and per entry
    with at least FILTER_MIN_READS reads after the first stage that loses more than FILTER_LOSS of them by the last:
    (entry, {stage: reads}, the pass that removes most), most lost first."""
    totals = {e['pool_id']: {stage: sum(stage_counts[s][stage].get(e['sequence'], 0) for s in stage_counts) for stage in stages}
              for e in pool if e['group'] in GROUPS[:3]}
    groups = {}
    for e in pool:
        if e['pool_id'] in totals:
            g = groups.setdefault(e['group'], dict({'sequences': 0}, **dict.fromkeys(stages, 0)))
            g['sequences'] += 1
            for stage in stages:
                g[stage] += totals[e['pool_id']][stage]
    lossy = []
    for e in pool:
        t = totals.get(e['pool_id'])
        if t and t[stages[0]] >= FILTER_MIN_READS and 1 - t[stages[-1]] / t[stages[0]] > FILTER_LOSS:
            worst = max(stages[1:], key=lambda st: t[stages[stages.index(st) - 1]] - t[st])
            lossy.append((e, t, worst))
    lossy.sort(key=lambda x: x[1][stages[-1]] / x[1][stages[0]])
    return {g: groups[g] for g in GROUPS[:3] if g in groups}, lossy


def shared(pool):
    """{pool_id: [other pool_ids]} for exact and shifted entries that share a miRNA (accession) with another one.

    A method counts a miRNA once, whichever pool sequences its reads came from, so its count cannot be split between
    such entries, while the oracle counts each sequence (in the equimolar pool 9 accessions, 7 of them an exact
    sequence and a shifted variant, e.g. one base shorter)."""
    owners = collections.defaultdict(set)
    for e in pool:
        if e['group'] in ('exact', 'shifted'):
            for a in e['mature_ids']:
                owners[a].add(e['pool_id'])
    out = collections.defaultdict(set)
    for ids in owners.values():
        if len(ids) > 1:
            for i in ids:
                out[i] |= ids - {i}
    return {i: sorted(others) for i, others in out.items()}


def scored(pool, group='exact'):
    """The entries of a group (exact or shifted) scored one by one: those sharing no miRNA with another exact or
    shifted entry (shared)."""
    ambiguous = shared(pool)
    return [e for e in pool if e['group'] == group and e['pool_id'] not in ambiguous]


def method_counts(common, pool):
    """{(method, sample, pool_id): count} for exact and shifted entries (sum over their accessions), and per
    (method, sample) the reads the method reported as miRNAs: {'all': every one, 'outside': on accessions of no exact
    or shifted entry (assigned to miRNAs not in the pool), 'unplaced': on names without a miRBase accession}."""
    entry_of = {}
    for entry in pool:
        if entry['group'] in ('exact', 'shifted'):
            for accession in entry['mature_ids']:
                entry_of.setdefault(accession, []).append(entry['pool_id'])
    counts, reads = {}, collections.defaultdict(lambda: {'all': 0, 'outside': 0, 'unplaced': 0})
    for method, sample, _, mature_id, count in common:
        reads[(method, sample)]['all'] += count
        if not mature_id:
            reads[(method, sample)]['unplaced'] += count
            continue
        if mature_id in entry_of:
            for pool_id in entry_of[mature_id]:
                counts[(method, sample, pool_id)] = counts.get((method, sample, pool_id), 0) + count
        else:
            reads[(method, sample)]['outside'] += count
    return counts, dict(reads)


PER_SEQUENCE = ['sample', 'pool_id', 'group', 'shift', 'shared_with', 'oracle']


def per_sequence(counts, oracle_counts, pool, methods, samples, kinds):
    """Rows of per_sequence.tsv (PER_SEQUENCE + methods) for exact and shifted entries: the kind of shift (kinds, from
    shift_kinds), the oracle count, and each method's count, NA where the entry shares its miRNA with others
    (shared_with), as it cannot be split between them."""
    ambiguous = shared(pool)
    return [[sample, e['pool_id'], e['group'], kinds.get(e['pool_id'], ''), ','.join(ambiguous.get(e['pool_id'], [])),
             oracle_counts[sample].get(e['sequence'], 0)] +
            ['NA' if e['pool_id'] in ambiguous else '{:g}'.format(counts.get((m, sample, e['pool_id']), 0)) for m in methods]
            for sample in samples for e in pool if e['group'] in ('exact', 'shifted')]


def spearman(x, y):
    def ranks(values):
        order = sorted(range(len(values)), key=values.__getitem__)
        result = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            for k in range(i, j + 1):
                result[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return result
    if len(x) < 3:
        return float('nan')
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx, syy = sum((a - mx) ** 2 for a in rx), sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float('nan')


def qcd(values):
    """Quartile coefficient of dispersion (Q3 - Q1) / (Q3 + Q1)."""
    if len(values) < 4:
        return float('nan')
    q1, _, q3 = statistics.quantiles(values, n=4)
    return (q3 - q1) / (q3 + q1) if q3 + q1 else float('nan')


def against_oracle(values, truth):
    """A method's counts against the oracle's, over entries the oracle counts: median log2(method + 1) / (oracle + 1),
    share within 2-fold and beyond 10-fold, Spearman correlation."""
    pairs = [(v, t) for v, t in zip(values, truth) if t > 0]
    ratios = [math.log2((v + 1) / (t + 1)) for v, t in pairs]
    return {'median_log2_vs_oracle': statistics.median(ratios) if ratios else float('nan'),
            'within_2fold': sum(1 for r in ratios if abs(r) <= 1) / len(ratios) if ratios else float('nan'),
            'beyond_10fold': sum(1 for r in ratios if abs(r) > math.log2(10)) / len(ratios) if ratios else float('nan'),
            'spearman_vs_oracle': spearman([v for v, _ in pairs], [t for _, t in pairs])}


def detected(values):
    return sum(1 for v in values if v > 0) / len(values) if values else float('nan')


def agreement(counts, oracle_counts, entries, methods, samples):
    """{row: {metric: value}} over entries, averaged over samples: `detected` for 'oracle' and each method, and
    against_oracle for each method (the shifted group and its kinds of shift in an equimolar pool)."""
    result = {}
    for row in ['oracle'] + list(methods):
        per_sample = []
        for sample in samples:
            truth = [oracle_counts[sample].get(e['sequence'], 0) for e in entries]
            values = truth if row == 'oracle' else [counts.get((row, sample, e['pool_id']), 0) for e in entries]
            per_sample.append(dict({'detected': detected(values)}, **({} if row == 'oracle' else against_oracle(values, truth))))
        result[row] = {k: round(statistics.mean(m[k] for m in per_sample), 4) for k in per_sample[0]} if per_sample else {}
    return result


def equimolar_metrics(counts, reads, oracle_counts, pool, methods, samples):
    """{row: {metric: value}} for 'oracle' and each method, over the scored exact entries, averaged over samples.

    Against the oracle (methods only): median log2(method + 1) / (oracle + 1), share within 2-fold and beyond 10-fold,
    Spearman correlation. For every row: share detected, CV and QCD of counts across the group (spread of an
    equimolar pool, so kit bias for the oracle), and of all the reads a method reported as miRNAs, `outside`: the share
    on miRNAs of no exact or shifted entry (reads from foreign pool sequences assigned to the organism's miRNAs) and
    `unplaced`: the share on names without a miRBase accession.
    """
    exact = scored(pool)
    result = {}
    for row in ['oracle'] + list(methods):
        per_sample = []
        for sample in samples:
            truth = [oracle_counts[sample].get(e['sequence'], 0) for e in exact]
            values = truth if row == 'oracle' else [counts.get((row, sample, e['pool_id']), 0) for e in exact]
            m = {'detected': detected(values),
                 'cv': statistics.pstdev(values) / statistics.mean(values) if values and statistics.mean(values) else float('nan'),
                 'qcd': qcd(values)}
            if row != 'oracle':
                m.update(against_oracle(values, truth))
                r = reads.get((row, sample), {'all': 0})
                m['outside'] = r['outside'] / r['all'] if r['all'] else float('nan')
                m['unplaced'] = r['unplaced'] / r['all'] if r['all'] else float('nan')
            per_sample.append(m)
        result[row] = {k: round(statistics.mean(m[k] for m in per_sample), 4) for k in per_sample[0]} if per_sample else {}
    return result


def ratio_scores(pairs):
    """Observed against designed log2(A / B), over (observed, designed) pairs: median absolute error, Spearman
    correlation, `called` and `false_calls` (ratiometric_metrics)."""
    errors = [abs(o - d) for o, d in pairs]
    apart = [(o, d) for o, d in pairs if abs(d) >= 1]
    close = [o for o, d in pairs if abs(d) < 1]
    return {'median_abs_error': round(statistics.median(errors), 4) if errors else float('nan'),
            'spearman_vs_design': round(spearman([o for o, _ in pairs], [d for _, d in pairs]), 4),
            'called': round(sum(1 for o, d in apart if abs(o) >= 1 and o * d > 0) / len(apart), 4) if apart else float('nan'),
            'false_calls': round(sum(1 for o in close if abs(o) >= 1) / len(close), 4) if close else float('nan')}


def ratiometric_metrics(counts, oracle_counts, pool, methods, sample_groups, pseudo=1.0):
    """{row: {metric: value}} for 'oracle' and each method, over the scored exact entries, and with a `shifted_`
    prefix over the scored shifted entries.

    Counts are scaled to counts per million of the exact group in each sample and averaged per pool (A, B). The observed
    log2(A / B) is centred on the exact sequences designed 1:1 (their median set to 0; the pools' totals differ by
    design) and compared with the designed log2(amount_a / amount_b): median absolute error, Spearman correlation,
    `called` (share of sequences designed at least 2-fold apart that are observed at least 2-fold apart in the right
    direction) and `false_calls` (share of sequences designed less than 2-fold apart observed at least 2-fold apart).
    Shifted entries are scaled and centred as the exact ones, so the exact group's results do not depend on them.
    Also returns {row: {pool_id: observed log2(A / B)}} for both groups.
    """
    exact = scored(pool)
    entries = exact + scored(pool, 'shifted')
    result, observed_all = {}, {}
    for row in ['oracle'] + list(methods):
        means = {}
        for group in ('A', 'B'):
            scaled = []
            for sample in [s for s, g in sample_groups.items() if g == group]:
                values = [oracle_counts[sample].get(e['sequence'], 0) if row == 'oracle'
                          else counts.get((row, sample, e['pool_id']), 0) for e in entries]
                total = sum(values[:len(exact)])
                scaled.append([v * 1e6 / total if total else 0.0 for v in values])
            means[group] = [statistics.mean(column) for column in zip(*scaled)] if scaled else []
        if not means['A'] or not means['B']:
            continue
        observed = [math.log2((a + pseudo) / (b + pseudo)) if a + pseudo > 0 and b + pseudo > 0 else float('nan')
                    for a, b in zip(means['A'], means['B'])]
        designed = [math.log2(e['amount_a'] / e['amount_b']) for e in entries]
        equal = [o for o, d in zip(observed[:len(exact)], designed) if d == 0 and not math.isnan(o)]
        if not equal:
            raise ValueError('ratiometric pool without sequences designed 1:1 in the exact group; cannot centre the ratios')
        offset = statistics.median(equal)
        observed = [o - offset for o in observed]
        pairs = list(zip(observed, designed))
        kept = lambda part: [(o, d) for o, d in part if not math.isnan(o)]
        result[row] = ratio_scores(kept(pairs[:len(exact)]))
        if len(entries) > len(exact):
            result[row].update({'shifted_' + k: v for k, v in ratio_scores(kept(pairs[len(exact):])).items()})
        observed_all[row] = {e['pool_id']: o for e, o in zip(entries, observed)}
    return result, observed_all
