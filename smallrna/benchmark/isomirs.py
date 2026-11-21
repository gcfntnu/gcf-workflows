"""isomiR-level scores (#454): what the routes report per read sequence (adapters.isomir_table) against the synthetic
truth, against the pools' oracle and, as a description only, between routes on biological samples.

A route reports an isomiR when it reports the read's own sequence on the read's miRNA: the read sequence is the key, so
the routes need no shared isomiR naming. variant() describes a sequence relative to its mature miRNA (5' and 3'
offsets, non-templated 3' bases, a substitution); the synthetic scores use the truth's own read classes instead.
Accessions of identical mature sequences count as one (representatives, as biological.collapse).

Rows without a sequence (a mirflowz plain row whose accession has no miRBase mature, adapters.py) cannot be matched by
sequence and are left out here. A route's label is its own category: miRge3's `exact miRNA`, labelled canonical, holds
reads matching its library's sequences, not only the miRBase mature (on the genome-wide synthetic set 82% are the
mature sequence, the rest mostly 3' shorter or longer), so its `isomir_as_canonical` describes that class rather than
a mislabelling.
"""
import collections
import math

import numpy

MAX_OFFSET = 3  # nt from either end of the mature sequence
MAX_ADDED = 3  # non-templated 3' bases
KINDS = ('canonical', 'substitution', '5p', '3p_add', '3p_longer', '3p_templated', '3p_trim', 'other')
MIN_DETECTED = 10  # reads over all samples for a biological isomiR to count as detected
SUBSTITUTED = 'substituted'  # truth reads with a substitution, whatever their class
SYNTHETIC_KEYS = ['isomir_recovered', 'isomir_misassigned', 'isomir_not_in_sample', 'isomir_as_canonical', 'canonical_as_isomir']
POOL_KEYS = ['exact_sequence_recovered', 'shifted_isomir_recovered', 'shifted_isomir_within_2fold', 'artefact_share']
BIOLOGICAL_KEYS = ['isomir_share', 'isomirs_detected']

Variant = collections.namedtuple('Variant', 'start end added substitutions')


def variant(read, mature, hairpin=None):
    """The read as a variant of `mature`: Variant(start, end, added, substitutions), or None beyond MAX_OFFSET,
    MAX_ADDED and one substitution. `start` and `end` are the offsets of the read's templated part from the mature's 5'
    and 3' ends (negative: upstream or shorter), `added` its non-templated 3' bases. The read is placed on `hairpin`
    when given, so an extension that matches the hairpin is templated; without it, 3' bases beyond the mature are
    `added`, templated or not, and a read starting upstream of the mature is None. Of several placements, the one with
    the fewest substitutions, then the fewest added bases, then the smallest offsets."""
    reference = hairpin or mature
    at = reference.find(mature) if mature else -1
    if at < 0:
        return None
    best = None
    for added in range(MAX_ADDED + 1):
        core = read[:len(read) - added]
        for start in range(-MAX_OFFSET, MAX_OFFSET + 1):
            begin, end = at + start, start + len(core) - len(mature)
            if not core or begin < 0 or begin + len(core) > len(reference) or abs(end) > MAX_OFFSET:
                continue
            substitutions = sum(a != b for a, b in zip(core, reference[begin:begin + len(core)]))
            rank = (substitutions, added, abs(start) + abs(end))
            if substitutions <= 1 and (best is None or rank < best[0]):
                best = rank, Variant(start, end, read[len(core):], substitutions)
    return best[1] if best else None


def kind(v, templated_known=True):
    """A Variant's kind (KINDS), the first that applies: other (None), substitution, 5p, 3p_add (non-templated bases;
    3p_longer when no hairpin told templated from added), 3p_templated, 3p_trim, canonical."""
    if v is None:
        return 'other'
    if v.substitutions:
        return 'substitution'
    if v.start:
        return '5p'
    if v.added:
        return '3p_add' if templated_known else '3p_longer'
    if v.end > 0:
        return '3p_templated'
    return '3p_trim' if v.end < 0 else 'canonical'


def representatives(id_to_seq):
    """{accession: the first (sorted) accession of the same mature sequence}, as biological.collapse."""
    by_seq = collections.defaultdict(list)
    for accession, seq in id_to_seq.items():
        by_seq[seq].append(accession)
    return {a: min(accessions) for accessions in by_seq.values() for a in accessions}


class Tally:
    """Named sums per unit ({unit: {name: value}}), so a bootstrap resample of the units (repeats allowed) adds up its
    rows; the `fixed` units count in every resample (e.g. reads on names without an accession)."""

    def __init__(self, rows, names, fixed=()):
        self.names, self.units = list(names), sorted(rows)
        self.index = {unit: i for i, unit in enumerate(self.units)}
        self.matrix = numpy.array([[rows[u].get(n, 0.0) for n in self.names] for u in self.units], dtype=float).reshape(len(self.units), -1)
        self.fixed = [self.index[u] for u in fixed if u in self.index]

    def sum(self, units=None):
        """{name: sum} over `units` (default: every unit once)."""
        if units is None:
            totals = self.matrix.sum(axis=0)
        else:
            rows = [self.index[u] for u in units if u in self.index and self.index[u] not in self.fixed] + self.fixed
            totals = self.matrix[rows].sum(axis=0)
        return dict(zip(self.names, totals.tolist()))


def share(part, whole):
    return round(part / whole, 4) if whole else float('nan')


def group(read_class, errors):
    """The truth class a read is scored under: SUBSTITUTED when it has a substitution."""
    return SUBSTITUTED if errors else read_class


def synthetic_scores(truth, isomirs, samples, methods, alias):
    """({method: {metric: value}}, {class: true reads}) over the mixed samples of the genome-wide synthetic dataset.

    A true read of a miRNA is reported when the method reports its sequence on that miRNA in its sample, up to the true
    count. Per class (canonical, the isomiR classes, SUBSTITUTED): `isomir_recovered_<class>`, the reported share;
    `isomir_recovered`: the same over all isomiR classes, without substitutions. Over the method's reads:
    `isomir_misassigned`, the share on a sequence the sample has under another miRNA or among its random reads;
    `isomir_not_in_sample`, the share on sequences the sample does not have. For methods that label reads (adapters.py),
    over their reads on true sequences of the right miRNA: `isomir_as_canonical`, the share of the reads labelled
    canonical that are not canonical in the truth (an isomiR class or a substitution); `canonical_as_isomir`, the share
    of the reads labelled isomiR that are canonical in the truth."""
    tallies, classes = synthetic_tallies(truth, isomirs, samples, methods, alias)
    return synthetic_values(tallies, classes), classes


def synthetic_tallies(truth, isomirs, samples, methods, alias):
    """({method: (Tally per miRNA, labelled)}, {class: true reads}): synthetic_scores' sums per miRNA (representative
    accession), the truth's reads under their miRNA and a method's under the miRNA it reports them on, so a resample of
    the miRNAs scores by adding them up (synthetic_values). Reads on names without an accession count in every resample."""
    mixed = {row['sample'] for row in samples if row['kind'] == 'mixed'}
    true, in_sample = {}, set()
    for row in truth:
        if row['sample'] not in mixed:
            continue
        in_sample.add((row['sample'], row['sequence']))
        if row['mature_id'] and row['class'] not in ('random', 'carrier'):
            key = (row['sample'], alias.get(row['mature_id'], row['mature_id']), row['sequence'])
            entry = true.setdefault(key, [0.0, group(row['class'], int(row['errors']))])  # a second class of the same read keeps the first
            entry[0] += float(row['count'])
    classes = collections.Counter()
    for count, read_class in true.values():
        classes[read_class] += count
    names = ['true:' + c for c in classes] + ['found:' + c for c in classes] + \
        ['reads', 'elsewhere', 'absent', 'canonical', 'canonical_wrong', 'isomiR', 'isomiR_missed']
    tallies = {}
    for method in methods:
        reported, labels = collections.defaultdict(float), collections.defaultdict(lambda: collections.defaultdict(float))
        for m, sample, _, accession, sequence, label, count in isomirs:
            if m == method and sample in mixed and sequence:
                key = (sample, alias.get(accession, accession), sequence)
                reported[key] += count
                labels[key][label] += count
        rows = collections.defaultdict(collections.Counter)
        for key, (count, read_class) in true.items():
            row = rows[key[1]]
            row['true:' + read_class] += count
            row['found:' + read_class] += min(reported.get(key, 0.0), count)
            by_label = labels.get(key, {})
            row['canonical'] += by_label.get('canonical', 0.0)
            row['isomiR'] += by_label.get('isomiR', 0.0)
            row['canonical_wrong' if read_class != 'canonical' else 'isomiR_missed'] += by_label.get(
                'canonical' if read_class != 'canonical' else 'isomiR', 0.0)
        for key, count in reported.items():
            row = rows[key[1]]
            row['reads'] += count
            if key not in true:
                row['elsewhere' if (key[0], key[2]) in in_sample else 'absent'] += count
        labelled = any(label for by_label in labels.values() for label in by_label)
        tallies[method] = (Tally(rows, names, fixed=('',)), labelled)
    return tallies, dict(classes)


def synthetic_values(tallies, classes, units=None):
    """{method: {metric: value}} of synthetic_scores from synthetic_tallies, over the miRNAs `units` (a bootstrap
    resample; default: all)."""
    out = {}
    for method, (tally, labelled) in tallies.items():
        t = tally.sum(units)
        values = {'isomir_recovered_' + c: share(t['found:' + c], t['true:' + c]) for c in classes}
        isomir_classes = [c for c in classes if c not in ('canonical', SUBSTITUTED)]
        values['isomir_recovered'] = share(sum(t['found:' + c] for c in isomir_classes), sum(t['true:' + c] for c in isomir_classes))
        values.update(isomir_misassigned=share(t['elsewhere'], t['reads']), isomir_not_in_sample=share(t['absent'], t['reads']))
        if labelled:
            values.update(isomir_as_canonical=share(t['canonical_wrong'], t['canonical']), canonical_as_isomir=share(t['isomiR_missed'], t['isomiR']))
        out[method] = values
    return out


def class_order(classes):
    return ['canonical'] + sorted(c for c in classes if c not in ('canonical', SUBSTITUTED)) + [SUBSTITUTED]


def synthetic_keys(classes):
    return SYNTHETIC_KEYS + ['isomir_recovered_' + c for c in class_order(classes) if c in classes]


def text(value):
    return '' if value is None or (isinstance(value, float) and math.isnan(value)) else str(value)


def without_line(without):
    return ['', 'Without isomiR output (features only): {}.'.format(', '.join(without))] if without else []


def synthetic_lines(values, classes, methods, without):
    """report.md's isomiR section of the genome-wide synthetic dataset."""
    lines = ['', '## isomiRs', '',
             'Mixed samples: the true reads of each class a method reports as their own sequence on their miRNA (isomiR table, '
             'isomirs.tsv), up to the true count; `substituted` are reads with a substitution, whatever their class. Over a '
             'method\'s reads: `isomir_misassigned` on a sequence the sample has under another miRNA or as a random read, '
             '`isomir_not_in_sample` on sequences the sample does not have. Its own labels against the truth: '
             '`isomir_as_canonical` (of its reads labelled canonical), `canonical_as_isomir` (of those labelled isomiR); a label is '
             'the method\'s own category, and miRge3\'s `exact miRNA` holds reads matching its library\'s sequences, shorter and '
             'longer ones among them, not only the miRBase mature.', '',
             '| class | true reads | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 2) + '|']
    for c in class_order(classes):
        if c in classes:
            lines.append('| {} | {:.0f} | '.format(c, classes[c]) + ' | '.join(text(values.get(m, {}).get('isomir_recovered_' + c))
                                                                         for m in methods) + ' |')
    lines += ['', '| metric | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 1) + '|']
    lines += ['| {} | '.format(k) + ' | '.join(text(values.get(m, {}).get(k)) for m in methods) + ' |' for k in SYNTHETIC_KEYS]
    return lines + without_line(without)


def pool_scores(isomirs, exact, shifted, entries, oracle, methods, id_to_seq):
    """({method: {metric: value}}, {method: {kind: reads}}) of a pool against its oracle, over all runs.

    `exact_sequence_recovered` and `shifted_isomir_recovered`: of the oracle's reads of each scored exact (shifted)
    sequence, the share the method reports as that sequence on the entry's miRNA, per run up to the oracle's count;
    `shifted_isomir_within_2fold`: the share of shifted sequences whose summed count is within 2-fold of the oracle's
    (one read added to both). `artefact_share`: of the method's reads on the pool's miRNAs, the share on sequences the
    pool does not hold (library or sequencing artefacts, as all methods read the same filtered reads), with their kinds
    (variant() against the mature sequence)."""
    tallies, kinds = pool_tallies(isomirs, exact, shifted, entries, oracle, methods, id_to_seq)
    return pool_values(tallies), kinds


PoolTally = collections.namedtuple('PoolTally', 'columns rows exact shifted reported truth on_pool artefacts')


def pool_tallies(isomirs, exact, shifted, entries, oracle, methods, id_to_seq):
    """({method: PoolTally}, {method: {kind: reads}}): pool_scores' counts per scored entry and run, so a resample of the
    entries and runs scores by adding them up (pool_values). A method's reads on a pool miRNA count under the first
    scored entry (exact, then shifted) that lists that miRNA."""
    sequences = {e['sequence'] for e in entries}
    scored = exact + shifted
    rows = {e['pool_id']: i for i, e in enumerate(scored)}
    columns = {sample: j for j, sample in enumerate(sorted(oracle))}
    owner = {}
    for e in scored:
        for accession in e['mature_ids']:
            owner.setdefault(accession, rows[e['pool_id']])
    shape = (len(scored), len(columns))
    truth = numpy.zeros(shape)
    for e in scored:
        for sample, j in columns.items():
            truth[rows[e['pool_id']], j] = oracle[sample].get(e['sequence'], 0)
    tallies, kinds = {}, {}
    for method in methods:
        reported = collections.defaultdict(float)
        for m, sample, _, accession, sequence, _, count in isomirs:
            if m == method and sequence and sample in columns:
                reported[(sample, accession, sequence)] += count
        counts, on_pool, artefacts = numpy.zeros(shape), numpy.zeros(shape), numpy.zeros(shape)
        for e in scored:
            for sample, j in columns.items():
                counts[rows[e['pool_id']], j] = sum(reported.get((sample, a, e['sequence']), 0.0) for a in e['mature_ids'])
        by_kind = collections.Counter()
        for (sample, accession, sequence), c in reported.items():
            if accession in owner:
                on_pool[owner[accession], columns[sample]] += c
                if sequence not in sequences:
                    artefacts[owner[accession], columns[sample]] += c
                    by_kind[kind(variant(sequence, id_to_seq.get(accession, '')), templated_known=False)] += c
        tallies[method] = PoolTally(columns, rows, [rows[e['pool_id']] for e in exact], [rows[e['pool_id']] for e in shifted],
                                    counts, truth, on_pool, artefacts)
        kinds[method] = dict(by_kind)
    return tallies, kinds


def pool_values(tallies, exact=None, shifted=None, samples=None):
    """{method: {metric: value}} of pool_scores from pool_tallies, over the scored entries `exact` and `shifted` and the
    runs `samples` (a bootstrap resample, repeats allowed; default: all)."""
    out = {}
    for method, t in tallies.items():
        ex = t.exact if exact is None else [t.rows[e['pool_id']] for e in exact]
        sh = t.shifted if shifted is None else [t.rows[e['pool_id']] for e in shifted]
        cols = sorted(t.columns.values()) if samples is None else [t.columns[s] for s in samples]

        def total(matrix, rows):
            return float(matrix[numpy.ix_(rows, cols)].sum()) if rows and cols else 0.0

        def recovered(rows):
            if not rows or not cols:
                return float('nan')
            r, o = t.reported[numpy.ix_(rows, cols)], t.truth[numpy.ix_(rows, cols)]
            return share(float(numpy.minimum(r, o).sum()), float(o.sum()))

        if sh and cols:
            r, o = t.reported[numpy.ix_(sh, cols)].sum(axis=1), t.truth[numpy.ix_(sh, cols)].sum(axis=1)
            within = [0.5 <= (a + 1) / (b + 1) <= 2 for a, b in zip(r.tolist(), o.tolist())]
        else:
            within = []
        both = list(ex) + list(sh)
        out[method] = {'exact_sequence_recovered': recovered(ex), 'shifted_isomir_recovered': recovered(sh),
                       'shifted_isomir_within_2fold': share(sum(within), len(within)),
                       'artefact_share': share(total(t.artefacts, both), total(t.on_pool, both))}
    return out


def kind_table(kinds, methods):
    """Lines of a table of each method's reads per kind, as shares of its reads in the table."""
    shown = [k for k in KINDS if any(kinds.get(m, {}).get(k) for m in methods)]
    lines = ['| kind | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 1) + '|']
    for k in shown:
        lines.append('| {} | '.format(k) + ' | '.join(text(share(kinds.get(m, {}).get(k, 0.0), sum(kinds.get(m, {}).values())))
                                                      for m in methods) + ' |')
    return lines


def pool_lines(values, kinds, methods, without):
    """report.md's isomiR section of a pool."""
    lines = ['', '## isomiRs', '',
             'Against the oracle (exact matches of the pool sequences in the filtered reads), per run: the share of its reads of '
             'the scored `exact` (`exact_sequence_recovered`) and `shifted` sequences (`shifted_isomir_recovered`) that a method '
             'reports as that sequence on the entry\'s miRNA; `shifted_isomir_within_2fold`, the shifted sequences within 2-fold of '
             'the oracle. A synthetic pool holds only its designed sequences, so a method\'s reads on other sequences of its miRNAs '
             'are library or sequencing artefacts (`artefact_share`); as all methods read the same filtered reads, the differences '
             'are the methods\'.', '', '| metric | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 1) + '|']
    lines += ['| {} | '.format(k) + ' | '.join(text(values.get(m, {}).get(k)) for m in methods) + ' |' for k in POOL_KEYS]
    lines += ['', 'The artefact reads by kind, against the mature sequence (no hairpin, so `3p_longer` is templated or not; '
              '`other`: beyond {} nt, {} added bases or one substitution):'.format(MAX_OFFSET, MAX_ADDED), ''] + kind_table(kinds, methods)
    return lines + without_line(without)


def biological_scores(isomirs, methods, id_to_seq, alias):
    """({method: {metric: value}}, {method: {kind: reads}}, {(a, b): Jaccard}) of biological samples, a description:
    `isomir_share`, the share of a method's reads on miRNAs with a known sequence that are not the mature sequence;
    `isomirs_detected`, the sequences with at least MIN_DETECTED reads over all samples (per miRNA, identical mature
    sequences as one); their overlap between methods (Jaccard)."""
    tallies, detected, kinds, overlap = biological_tallies(isomirs, methods, id_to_seq, alias)
    values = biological_values(tallies)
    for method, count in detected.items():
        values[method]['isomirs_detected'] = count
    return values, kinds, overlap


def biological_tallies(isomirs, methods, id_to_seq, alias):
    """({method: Tally per sample of reads and canonical reads}, {method: isomiRs detected}, {method: {kind: reads}},
    {(a, b): Jaccard}) for biological_scores: `isomir_share` from the tallies, so a resample of the samples scores by
    adding them up (biological_values); the detected isomiRs and their overlap are counts over the dataset."""
    tallies, detected_sets, kinds = {}, {}, {}
    for method in methods:
        totals, by_kind = collections.defaultdict(float), collections.Counter()
        per_sample = collections.defaultdict(collections.Counter)
        for m, sample, _, accession, sequence, _, count in isomirs:
            if m == method and sequence and accession in id_to_seq:
                totals[(alias.get(accession, accession), sequence)] += count
                k = kind(variant(sequence, id_to_seq[accession]), templated_known=False)
                by_kind[k] += count
                per_sample[sample]['reads'] += count
                per_sample[sample]['canonical'] += count if k == 'canonical' else 0.0
        detected_sets[method] = {key for key, c in totals.items() if c >= MIN_DETECTED}
        tallies[method] = Tally(per_sample, ['reads', 'canonical'])
        kinds[method] = dict(by_kind)
    overlap = {(a, b): share(len(detected_sets[a] & detected_sets[b]), len(detected_sets[a] | detected_sets[b]))
               for i, a in enumerate(methods) for b in methods[i + 1:]}
    return tallies, {m: len(d) for m, d in detected_sets.items()}, kinds, overlap


def biological_values(tallies, samples=None):
    """{method: {'isomir_share': value}} from biological_tallies over the samples `samples` (a bootstrap resample, repeats
    allowed; default: all)."""
    out = {}
    for method, tally in tallies.items():
        t = tally.sum(samples)
        out[method] = {'isomir_share': share(t['reads'] - t['canonical'], t['reads'])}
    return out


def biological_lines(values, kinds, overlap, methods, without):
    """report.md's isomiR section of a biological dataset."""
    lines = ['', '## isomiRs', '',
             'No known composition: a description. `isomir_share`: the share of a method\'s miRNA reads (isomiR table, isomirs.tsv) '
             'that are not the mature sequence; `isomirs_detected`: sequences with at least {} reads over all samples. Reads by '
             'kind, against the mature sequence (`3p_longer` templated or not):'.format(MIN_DETECTED), ''] + kind_table(kinds, methods)
    if overlap:
        lines += ['', 'Detected isomiRs shared by two methods (Jaccard):', '', '| methods | Jaccard |', '|---|---|']
        lines += ['| {} / {} | {} |'.format(a, b, text(v)) for (a, b), v in overlap.items()]
    return lines + without_line(without)
