"""miRNAs scored by MirGeneDB stratum (#493).

miRBase lists entries that a curated miRNA set leaves out (fragments of rRNA and other RNAs, repeat-derived sequences),
and the routes disagree most on those. With `evaluate --mirgenedb GFF3` (MirGeneDB's annotation, e.g. the one in
miRge3's library, whose `miRNA` rows carry `Alias=MIMAT...` and mark star arms with an ID ending in `*`), each miRBase
accession gets a stratum:

- `mirgenedb`: MirGeneDB lists it, not as a star arm;
- `star`: MirGeneDB lists it as a star arm;
- `mirbase_only`: MirGeneDB does not list it.

Identical sequences count as one accession, as the evaluation collapses them (biological.collapse): the group takes
the first of these strata any member has. The strata split the per-miRNA metrics: for biological datasets the miRNAs
detected, those detected by no other route and the share of the miRNA reads per stratum, and the miRNAs the routes
dispute most; for the genome-wide
synthetic set the recovery per stratum.
"""
import collections
import hashlib
import re
import statistics

STRATA = ('mirgenedb', 'star', 'mirbase_only')
MIN_READS = 10         # a miRNA counts as detected in a sample with at least this many reads (biological.MIN_READS)
DISPUTED_READS = 100   # a miRNA is weighed when some route counts at least this many reads over the samples
DISPUTED_FOLD = 4      # disputed: the routes' totals differ at least this much, (max + 1) / (min + 1)
LISTED = 20            # disputed miRNAs listed in the report
SPECIES = {'homo_sapiens': 'Hsa', 'mus_musculus': 'Mmu', 'rattus_norvegicus': 'Rno'}   # MirGeneDB's ID prefixes


def read_gff(path):
    """({miRBase accession: 'mirgenedb' or 'star'} of a MirGeneDB GFF3's `miRNA` rows, {'release': the header's
    MirGeneDB release, 'species': the ID prefixes}); ValueError without any such row."""
    listed, info = {}, {'release': '', 'species': set()}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                found = re.search(r'MirGeneDB\s+(v?[0-9][0-9.]*)', line)
                info['release'] = info['release'] or (found.group(1) if found else '')
                continue
            if len(fields) < 9 or fields[2] != 'miRNA':
                continue
            prefix = re.search(r'(?:^|;)ID=([A-Za-z]+)-', fields[8])
            if prefix:
                info['species'].add(prefix.group(1))
            star = re.search(r'(?:^|;)ID=[^;]*\*(?:;|$)', fields[8]) is not None
            for accession in re.findall(r'MIMAT\d+', fields[8]):
                if listed.get(accession) != 'mirgenedb':  # an accession listed as an arm and as a star arm is an arm
                    listed[accession] = 'star' if star else 'mirgenedb'
    if not listed:
        raise ValueError('{}: no miRNA row with a miRBase accession (Alias=MIMAT...)'.format(path))
    return listed, info


def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def strata(listed, groups=()):
    """stratum(accession): its stratum, an identical-sequence group's (biological.collapse; `groups`: [accessions], the
    first the group's accession) the first of STRATA any member has."""
    of_group = {}
    for accessions in groups:
        found = [listed.get(a, 'mirbase_only') for a in accessions]
        of_group[accessions[0]] = min(found, key=STRATA.index)
    return lambda accession: of_group.get(accession) or listed.get(accession, 'mirbase_only')


def keys():
    """The biological metrics per stratum: `detected.<stratum>`, `detected_alone.<stratum>` and `reads.<stratum>`."""
    return ['{}.{}'.format(what, s) for what in ('detected', 'detected_alone', 'reads') for s in STRATA]


def biological_values(by, methods, samples, stratum):
    """{method: {detected.<stratum>: median per sample of the miRNAs with at least MIN_READS reads,
    detected_alone.<stratum>: of those, the ones no other method detects in the sample, reads.<stratum>: median share
    of the sample's miRNA reads (on accessions)}} from `by` ({(method, sample): {accession: count}},
    biological.counts_by)."""
    detected = {(m, s): {a for a, v in by.get((m, s), {}).items() if v >= MIN_READS} for m in methods for s in samples}
    out = {}
    for m in methods:
        out[m] = {}
        per = [by.get((m, s), {}) for s in samples]
        alone = [detected[(m, x)] - set().union(*(detected[(o, x)] for o in methods if o != m)) for x in samples]
        for s in STRATA:
            out[m]['detected.' + s] = statistics.median(sum(1 for a, v in c.items() if v >= MIN_READS and stratum(a) == s) for c in per)
            out[m]['detected_alone.' + s] = statistics.median(sum(1 for a in d if stratum(a) == s) for d in alone)
            shares = [sum(v for a, v in c.items() if stratum(a) == s) / sum(c.values()) for c in per if sum(c.values())]
            out[m]['reads.' + s] = round(statistics.median(shares), 4) if shares else ''
    return out


def disputed(by, methods, stratum):
    """({stratum: (weighed, disputed)}, [(fold, accession, {method: total}, stratum)] of the disputed miRNAs, most disputed
    first): the miRNAs some route counts at least DISPUTED_READS reads over the samples, and those whose routes'
    totals differ at least DISPUTED_FOLD-fold."""
    totals = collections.defaultdict(collections.Counter)
    for (m, _), counts in by.items():
        if m in methods:
            for a, v in counts.items():
                totals[a][m] += v
    counts, rows = collections.defaultdict(lambda: [0, 0]), []
    for a, per in totals.items():
        if max(per.values()) < DISPUTED_READS:
            continue
        counts[stratum(a)][0] += 1
        fold = (max(per.get(m, 0) for m in methods) + 1) / (min(per.get(m, 0) for m in methods) + 1)
        if fold >= DISPUTED_FOLD:
            counts[stratum(a)][1] += 1
            rows.append((fold, a, dict(per), stratum(a)))
    rows.sort(key=lambda row: (-row[0], row[1]))
    return {s: tuple(counts[s]) for s in STRATA}, rows


def biological_lines(values, counts, rows, methods, name, source, md_table):
    """report.md's section on the strata of a biological or mixture dataset; `name(accession)` its miRBase name,
    `source` a line on the GFF3 read."""
    lines = ['', '## MirGeneDB strata (#493)', '',
             'Each miRBase accession by MirGeneDB stratum: `mirgenedb` (listed), `star` (listed as a star arm), `mirbase_only` (not '
             'listed); identical sequences count as one, in the first stratum any of them has. Per route: `detected.<stratum>`, the '
             'miRNAs with at least {} reads (median over samples); `detected_alone.<stratum>`, those of them no other route detects '
             'in the sample; and `reads.<stratum>`, the share of the miRNA reads (median over samples). {}'.format(MIN_READS, source), '']
    lines += md_table(['metric'] + methods, [[k] + [values[m][k] for m in methods] for k in keys()])
    lines += ['', '**Disputed miRNAs:** of the miRNAs some route counts at least {} reads over the samples, those whose routes\' '
              'totals differ at least {}-fold, (max + 1) / (min + 1): {}.'.format(
                  DISPUTED_READS, DISPUTED_FOLD, '; '.join('{} {} of {}'.format(s, counts[s][1], counts[s][0]) for s in STRATA))]
    if rows:
        lines += ['', 'The {} most disputed, with each route\'s total reads:'.format(min(LISTED, len(rows))), '']
        lines += md_table(['miRNA', 'stratum', 'fold'] + methods,
                          [[name(a), s, '{:.0f}'.format(fold)] + ['{:.0f}'.format(per.get(m, 0)) for m in methods]
                           for fold, a, per, s in rows[:LISTED]])
    return lines
