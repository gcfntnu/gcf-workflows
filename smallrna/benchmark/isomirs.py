"""isomiR-level scores (#454): what the routes report per read sequence (adapters.isomir_table) against the synthetic
truth, against the pools' oracle and, as a description only, between routes on biological samples.

A route reports an isomiR when it reports the read's own sequence on the read's miRNA: the read sequence is the key, so
the routes need no shared isomiR naming. variant() describes a sequence relative to its mature miRNA (5' and 3'
offsets, non-templated 3' bases, a substitution); the synthetic scores use the truth's own read classes instead.
Accessions of identical mature sequences count as one (representatives, as biological.collapse).

Rows without a sequence (mirflowz's plain rows, adapters.py) cannot be matched by sequence and are left out here.
"""
import collections
import math

MAX_OFFSET = 3  # nt from either end of the mature sequence
MAX_ADDED = 3  # non-templated 3' bases
KINDS = ('canonical', 'substitution', '5p', '3p_add', '3p_longer', '3p_templated', '3p_trim', 'other')
MIN_DETECTED = 10  # reads over all samples for a biological isomiR to count as detected
SUBSTITUTED = 'substituted'  # truth reads with a substitution, whatever their class
SYNTHETIC_KEYS = ['isomir_recovered', 'isomir_misassigned', 'isomir_not_in_sample', 'isomir_as_canonical', 'canonical_as_isomir']
POOL_KEYS = ['canonical_isomir_recovered', 'shifted_isomir_recovered', 'shifted_isomir_within_2fold', 'artefact_share']
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
    out = {}
    for method in methods:
        reported, labels = collections.defaultdict(float), collections.defaultdict(lambda: collections.defaultdict(float))
        for m, sample, _, accession, sequence, label, count in isomirs:
            if m == method and sample in mixed and sequence:
                key = (sample, alias.get(accession, accession), sequence)
                reported[key] += count
                labels[key][label] += count
        found = collections.Counter()
        for key, (count, read_class) in true.items():
            found[read_class] += min(reported.get(key, 0.0), count)
        values = {'isomir_recovered_' + c: share(found[c], n) for c, n in classes.items()}
        isomir_classes = [c for c in classes if c not in ('canonical', SUBSTITUTED)]
        values['isomir_recovered'] = share(sum(found[c] for c in isomir_classes), sum(classes[c] for c in isomir_classes))
        total = sum(reported.values())
        elsewhere = sum(c for key, c in reported.items() if key not in true and (key[0], key[2]) in in_sample)
        absent = sum(c for key, c in reported.items() if (key[0], key[2]) not in in_sample)
        values.update(isomir_misassigned=share(elsewhere, total), isomir_not_in_sample=share(absent, total))
        if any(label for by_label in labels.values() for label in by_label):
            canonical = sum(labels[key]['canonical'] for key in true)
            wrong = sum(labels[key]['canonical'] for key, (_, c) in true.items() if c != 'canonical')
            isomir = sum(labels[key]['isomiR'] for key in true)
            missed = sum(labels[key]['isomiR'] for key, (_, c) in true.items() if c == 'canonical')
            values.update(isomir_as_canonical=share(wrong, canonical), canonical_as_isomir=share(missed, isomir))
        out[method] = values
    return out, dict(classes)


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
             '`isomir_as_canonical` (of its reads labelled canonical), `canonical_as_isomir` (of those labelled isomiR).', '',
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

    `canonical_isomir_recovered` and `shifted_isomir_recovered`: of the oracle's reads of each scored exact (shifted)
    sequence, the share the method reports as that sequence on the entry's miRNA, per run up to the oracle's count;
    `shifted_isomir_within_2fold`: the share of shifted sequences whose summed count is within 2-fold of the oracle's
    (one read added to both). `artefact_share`: of the method's reads on the pool's miRNAs, the share on sequences the
    pool does not hold (library or sequencing artefacts, as all methods read the same filtered reads), with their kinds
    (variant() against the mature sequence)."""
    sequences = {e['sequence'] for e in entries}
    mirnas = {accession for e in exact + shifted for accession in e['mature_ids']}
    out, kinds = {}, {}
    for method in methods:
        reported = collections.defaultdict(float)
        for m, sample, _, accession, sequence, _, count in isomirs:
            if m == method and sequence:
                reported[(sample, accession, sequence)] += count

        def count(sample, entry):
            return sum(reported.get((sample, a, entry['sequence']), 0.0) for a in entry['mature_ids'])

        def recovered(group):
            pairs = [(count(s, e), oracle[s].get(e['sequence'], 0)) for e in group for s in oracle]
            return share(sum(min(r, o) for r, o in pairs), sum(o for _, o in pairs))

        within = [0.5 <= (sum(count(s, e) for s in oracle) + 1) / (sum(oracle[s].get(e['sequence'], 0) for s in oracle) + 1) <= 2
                  for e in shifted]
        on_pool = {key: c for key, c in reported.items() if key[1] in mirnas}
        artefacts = collections.Counter()
        for (_, accession, sequence), c in on_pool.items():
            if sequence not in sequences:
                artefacts[kind(variant(sequence, id_to_seq.get(accession, '')), templated_known=False)] += c
        out[method] = {'canonical_isomir_recovered': recovered(exact), 'shifted_isomir_recovered': recovered(shifted),
                       'shifted_isomir_within_2fold': share(sum(within), len(within)),
                       'artefact_share': share(sum(artefacts.values()), sum(on_pool.values()))}
        kinds[method] = dict(artefacts)
    return out, kinds


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
             'the scored `exact` (`canonical_isomir_recovered`) and `shifted` sequences (`shifted_isomir_recovered`) that a method '
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
    out, kinds, detected = {}, {}, {}
    for method in methods:
        totals, by_kind = collections.defaultdict(float), collections.Counter()
        for m, _, _, accession, sequence, _, count in isomirs:
            if m == method and sequence and accession in id_to_seq:
                totals[(alias.get(accession, accession), sequence)] += count
                by_kind[kind(variant(sequence, id_to_seq[accession]), templated_known=False)] += count
        reads = sum(totals.values())
        detected[method] = {key for key, c in totals.items() if c >= MIN_DETECTED}
        out[method] = {'isomir_share': share(reads - by_kind['canonical'], reads), 'isomirs_detected': len(detected[method])}
        kinds[method] = dict(by_kind)
    overlap = {(a, b): share(len(detected[a] & detected[b]), len(detected[a] | detected[b]))
               for i, a in enumerate(methods) for b in methods[i + 1:]}
    return out, kinds, overlap


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
