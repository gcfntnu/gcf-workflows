"""Decoy reads of other RNA classes for the genome-wide synthetic dataset (#412 PR 2), declared under `synthetic.decoys`.

Decoys are reads of other RNA classes, which no method should count as miRNAs. They are scored, not as classes of
their own, but as reads counted as miRNAs; the truth keeps each read's class and source, so the classes can be scored
later without new data.

- **Units:** a class and the source its molecules come from (UNITS), e.g. `tRNA-gtrnadb`. Each unit reads one filled
  source file (SOURCES; `precursor-mirbase` reads the miRBase hairpins the generator already has).
- **Fragments:** each unit draws a fixed catalogue of `fragments` distinct fragment species once, with lognormal
  weights (`weight_sd`), as degradation favours some fragments; reads are then drawn from the catalogue. Fragment
  types follow the class: tRNA 5' and 3' halves (cut in the anticodon loop), 5'-tRFs, 3'-tRFs ending in the added CCA
  and internal i-tRFs, a share `modification` of them with a mismatch at m1G9 or m1A58; Y RNA and snoRNA mostly 5'- and
  3'-end fragments; piRNAs 24-32 nt, starting with U (T) in a share `pirna_5u`; the other classes uniform fragments of
  `length`. Precursor fragments come from the hairpins outside every mature.
- **No miRBase look-alikes:** a fragment that shares OVERLAP nucleotides with any of the organism's miRBase hairpins
  (precursor fragments: that overlaps a mature), or equals a read a miRNA could give, is drawn again, so these decoys
  measure reads unrelated to miRBase. The miRBase matures that lie in other-class molecules are found instead
  (`embedded`) and probed from both sides: their canonical reads (embedded probes), and fragments of their host
  molecules that cover at least half the mature without being a read any miRNA could give (overlap probes; #412's
  class of fragments overlapping a miRBase hairpin, reported on its own and not as false positives).
- **Mixed samples:** unit u gets `load[u]` decoy reads per miRNA read of the sample, times 2^(`nuisance` x the sample's
  nuisance score), Poisson; reads carry sequencing errors as miRNA reads do.
- **Probe samples** (`probe.reads` per unit and replicate) hold the carrier and one unit's reads, without errors.

Everything is drawn from its own stream (Decoys.rng), so the miRNA reads of a declaration stay the same with or
without decoys.
"""
import collections
import gzip
import re

import numpy as np

SOURCES = ('refseq_rrna', 'gtrnadb', 'ensembl_ncrna', 'ensembl_cdna', 'pirnaclusterdb')
UNITS = {'rRNA-refseq': 'refseq_rrna', 'rRNA-ensembl': 'ensembl_ncrna', 'tRNA-gtrnadb': 'gtrnadb', 'tRNA-ensembl': 'ensembl_ncrna',
         'Y_RNA-ensembl': 'ensembl_ncrna', 'snoRNA-ensembl': 'ensembl_ncrna', 'other_ncRNA-ensembl': 'ensembl_ncrna',
         'mRNA-ensembl': 'ensembl_cdna', 'piRNA-pirnaclusterdb': 'pirnaclusterdb', 'precursor-mirbase': None}
EMBEDDING = ('rRNA', 'tRNA', 'Y_RNA', 'snoRNA', 'other_ncRNA')  # classes searched for miRBase matures inside their molecules
REFSEQ_RRNA = {'NR_003286': '18S', 'NR_003287': '28S', 'NR_003285': '5.8S', 'NR_023363': '5S', 'NR_046235': '45S'}
Y_RNAS = ('RNY1', 'RNY3', 'RNY4', 'RNY5')
OTHER_NCRNA = re.compile(r'RN7SK|RN7SL[12]|VTRNA\d-\d')  # 7SK, 7SL and vault RNAs; snRNAs by biotype
OVERLAP = 18
MAX_LENGTH = 45  # inserts longer than this are not drawn (a 3' tRNA half with CCA is about 40 nt)
TRNA_FRAGMENTS = {'5p_half': 0.3, '3p_half': 0.15, '5p_tRF': 0.25, '3p_tRF': 0.2, 'i_tRF': 0.1}
END_FRAGMENTS = {'Y_RNA': ({'5p_end': 0.6, '3p_end': 0.3, 'internal': 0.1}, (25, 33)),
                 'snoRNA': ({'5p_end': 0.35, '3p_end': 0.35, 'internal': 0.3}, (17, 30))}
TRNA_MODIFIED = (8, 57)  # 0-based positions of m1G9 and m1A58 in the mature tRNA
DEFAULTS = {'fragments': 2000, 'weight_sd': 1.5, 'length': [16, 40], 'modification': 0.2, 'pirna_5u': 0.85, 'nuisance': 0.0}


def unit_class(unit):
    return unit.rsplit('-', 1)[0]


def fasta(path):
    """[(header, DNA sequence)] of a FASTA file, gzipped or not; RNA and lower case become upper-case DNA."""
    out, header, parts = [], None, []
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    out.append((header, ''.join(parts)))
                header, parts = line[1:].strip(), []
            else:
                parts.append(line.strip().upper().replace('U', 'T'))
    if header is not None:
        out.append((header, ''.join(parts)))
    return out


def ensembl_fields(header):
    """{field: value} of an Ensembl FASTA header (`transcript_biotype:snoRNA gene_symbol:SNORD3A ...`), with `id`."""
    tokens = header.split()
    fields = dict(t.split(':', 1) for t in tokens[1:] if ':' in t)
    fields['id'] = tokens[0]
    return fields


def molecules(unit, path, hairpins=None):
    """[(name, sequence, anticodon)] of a unit's molecules (anticodon: tRNAs whose header names it, else None).
    `hairpins` [(name, sequence, mature spans)] gives the precursor unit's."""
    if unit == 'precursor-mirbase':
        return [(name, seq, spans) for name, seq, spans in hairpins]
    records = fasta(path)
    if unit == 'rRNA-refseq':
        by_kind = {}
        for header, seq in records:
            kind = REFSEQ_RRNA.get(header.split()[0].split('.')[0])
            if kind:
                by_kind[kind] = seq
        missing = sorted(set(REFSEQ_RRNA.values()) - set(by_kind))
        if missing:
            raise ValueError('{}: no {} rRNA record'.format(path, ', '.join(missing)))
        out = [(kind, by_kind[kind], None) for kind in ('18S', '5.8S', '28S', '5S')]
        pre, cuts = by_kind['45S'], []
        for kind in ('18S', '5.8S', '28S'):  # the spacers: the 45S outside its matures
            at = pre.find(by_kind[kind])
            if at < 0:
                raise ValueError('{}: the {} rRNA is not in the 45S pre-rRNA'.format(path, kind))
            cuts.append((at, at + len(by_kind[kind])))
        bounds = [0] + [x for cut in sorted(cuts) for x in cut] + [len(pre)]
        out += [(name, pre[bounds[2 * i]:bounds[2 * i + 1]], None) for i, name in enumerate(('5ETS', 'ITS1', 'ITS2', '3ETS'))]
        return [m for m in out if m[1]]
    if unit == 'tRNA-gtrnadb':
        out = []
        for header, seq in records:
            name = header.split()[0].replace('Homo_sapiens_', '')
            codon = re.search(r'\(([ACGTU]{3})\)', header)
            if 'Undet' in name or not codon:
                continue
            out.append((name, seq if seq.endswith('CCA') else seq + 'CCA', codon.group(1).replace('U', 'T')))
        return out
    if unit == 'piRNA-pirnaclusterdb':
        return [(header.split()[0], seq, None) for header, seq in records]
    pick = {'rRNA-ensembl': lambda f: f.get('transcript_biotype') in ('rRNA', 'Mt_rRNA'),
            'tRNA-ensembl': lambda f: f.get('transcript_biotype') == 'Mt_tRNA',
            'Y_RNA-ensembl': lambda f: f.get('gene_symbol') in Y_RNAS,
            'snoRNA-ensembl': lambda f: f.get('transcript_biotype') in ('snoRNA', 'scaRNA'),
            'other_ncRNA-ensembl': lambda f: f.get('transcript_biotype') == 'snRNA' or bool(OTHER_NCRNA.fullmatch(f.get('gene_symbol', ''))),
            'mRNA-ensembl': lambda f: f.get('transcript_biotype') == 'protein_coding'}[unit]
    out = []
    for header, seq in records:
        fields = ensembl_fields(header)
        if pick(fields):
            seq = seq + 'CCA' if unit == 'tRNA-ensembl' and not seq.endswith('CCA') else seq
            out.append((fields.get('gene_symbol') or fields['id'], seq, None))
    return out


def embedded(molecules_by_unit, matures):
    """{class: {key: [molecule names]}}: the matures ({key: sequence}) that lie exactly inside a molecule of a unit of an
    EMBEDDING class ({unit: [(name, sequence, extra)]}); the synthetic generator's embedded matures (#412) and the
    biological evaluation's flag (#415)."""
    by_length = collections.defaultdict(dict)
    for key, seq in matures.items():
        by_length[len(seq)].setdefault(seq, key)
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for unit, mols in molecules_by_unit.items():
        cls = unit_class(unit)
        if cls not in EMBEDDING:
            continue
        for name, seq, _ in mols:
            for length, by_seq in by_length.items():
                for i in range(len(seq) - length + 1):
                    key = by_seq.get(seq[i:i + length])
                    if key and name not in out[cls][key]:
                        out[cls][key].append(name)
    return {cls: dict(keys) for cls, keys in out.items() if keys}


def allowed_differences(length):
    """Differences a mature of this length may have from a molecule and still be flagged as its fragment (#415): none
    under 18 nt, one at 18-19 nt, two from 20 nt; a chance match in the classes' few megabases stays unlikely."""
    return 0 if length < 18 else 1 if length < 20 else 2


GENERIC_LOCUS = re.compile(r'ENST\d+.*|.*(\dP|P\d+)|U\d+(atac)?|Y_RNA|7SK|Metazoa_SRP|snoU\d+\w*')


def functional(name):
    """A named gene that is not a pseudogene or a generic Rfam locus (RNU6-346P, U4, Y_RNA, an unnamed ENST)."""
    return not GENERIC_LOCUS.fullmatch(name)


def near_embedded(molecules_by_unit, matures):
    """{class: {key: [(molecule name, differences)]}}: the matures ({key: sequence}) that lie inside a molecule of an
    EMBEDDING class exactly or with at most allowed_differences: mismatches, or bases beyond either end of the
    molecule (a mature on a tRNA's added CCA end, an entry with an extra base). The biological evaluation's flag (#415):
    miRBase holds entries from such molecules that differ from them by a base or two (sequencing errors, modifications,
    another copy). Differences only against functional genes: pseudogenes and generic loci overlap repeats that give
    miRBase families (miR-548 and RNU6 pseudogenes), so there only exact matches count. Candidates come from three
    disjoint 5-mers of the mature, of which one is exact."""
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    for unit, mols in molecules_by_unit.items():
        cls = unit_class(unit)
        if cls not in EMBEDDING:
            continue
        index = collections.defaultdict(list)
        for m, (_, seq, _) in enumerate(mols):
            for i in range(len(seq) - 4):
                index[seq[i:i + 5]].append((m, i))
        for key, mature in matures.items():
            length, limit = len(mature), allowed_differences(len(mature))
            if length < 15:
                continue
            tried = set()
            for start in (0, 5, 10):
                for m, i in index.get(mature[start:start + 5], ()):
                    offset = i - start
                    if (m, offset) in tried:
                        continue
                    tried.add((m, offset))
                    name, seq, _ = mols[m]
                    if not functional(name):
                        limit_here = 0
                    else:
                        limit_here = limit
                    diff = 0
                    for j, base in enumerate(mature):
                        k = offset + j
                        if not 0 <= k < len(seq) or seq[k] != base:
                            diff += 1
                            if diff > limit_here:
                                break
                    if diff <= limit_here and diff < out[cls][key].get(name, limit + 1):
                        out[cls][key][name] = diff
    return {cls: {key: sorted(names.items()) for key, names in keys.items() if names} for cls, keys in out.items()
            if any(keys.values())}


def contained(molecules_by_unit, sequences):
    """{sequence: [classes]}: the sequences (reads) that lie exactly inside a molecule of a unit of an EMBEDDING class
    (#415). One length at a time: the set of every substring of that length of a unit's molecules, built and dropped."""
    by_length = collections.defaultdict(set)
    for seq in sequences:
        by_length[len(seq)].add(seq)
    out = collections.defaultdict(set)
    for length, wanted in sorted(by_length.items()):
        for unit, mols in molecules_by_unit.items():
            if unit_class(unit) not in EMBEDDING:
                continue
            present = {seq[i:i + length] for _, seq, _ in mols for i in range(len(seq) - length + 1)}
            for seq in wanted & present:
                out[seq].add(unit_class(unit))
    return {seq: sorted(classes) for seq, classes in out.items()}


class Decoys:
    """The declared units' fragment catalogues, the embedded miRNAs, and the decoy rows of mixed and probe samples."""

    def __init__(self, spec, sources, universe, hairpins, reads, seed):
        self.spec = dict(DEFAULTS, **spec)
        self.rng = np.random.default_rng(np.random.SeedSequence(int(seed), spawn_key=(1,)))
        self.reads = reads.fork(self.rng)  # the miRNA reads' error and quality model, on this stream
        self.forbidden, self.min_length = reads.forbidden, reads.min_length
        if self.spec['length'][0] < self.min_length:
            raise ValueError('decoys: length starts at {}, below the filter\'s minimum read length {}; such reads never reach '
                             'a route'.format(self.spec['length'][0], self.min_length))
        self.kmers = {seq[i:i + OVERLAP] for _, seq, _ in hairpins for i in range(len(seq) - OVERLAP + 1)}
        self.units = list(self.spec['load'])
        self.molecules = {u: molecules(u, sources.get(UNITS[u]) if UNITS[u] else None, hairpins) for u in self.units}
        empty = [u for u in self.units if not self.molecules[u]]
        if empty:
            raise ValueError('decoys: no molecules for {}'.format(', '.join(empty)))
        self.catalogue = {u: self.draw_catalogue(u) for u in self.units}
        self.embedded = self.find_embedded(universe)
        self.matures = {rep: universe.mirnas[rep]['sequence'] for reps in self.embedded.values() for rep in reps}

    def lengths(self, low=None, high=None):
        lo, hi = self.spec['length']
        return int(self.rng.integers(low or lo, (high or hi) + 1))

    def fragment(self, unit, name, seq, extra):
        """(fragment, type) of one molecule, or None."""
        cls, rng = unit_class(unit), self.rng
        n = len(seq)
        if cls == 'tRNA':
            kind = rng.choice(list(TRNA_FRAGMENTS), p=list(TRNA_FRAGMENTS.values()))
            near = [i for i in range(28, min(40, n - 3)) if seq[i:i + 3] == extra] if extra else []
            anticodon = min(near, key=lambda i: abs(i - 33)) if near else min(33, n // 2)
            cut = anticodon + int(rng.integers(0, 3))  # angiogenin cuts in the anticodon loop
            start, end = {'5p_half': (0, cut), '3p_half': (cut, n)}.get(kind, (None, None))
            if kind == '5p_tRF':
                start, end = 0, self.lengths(16, 30)
            elif kind == '3p_tRF':
                start, end = n - self.lengths(18, 22), n
            elif kind == 'i_tRF':
                length = self.lengths(16, 30)
                if n - length - 20 < 10:
                    return None
                start = int(rng.integers(10, n - length - 10))
                end = start + length
            fragment = seq[start:end]
            inside = [p - start for p in TRNA_MODIFIED if start <= p < end]
            if inside and rng.random() < float(self.spec['modification']):
                p = int(rng.choice(inside))
                fragment = fragment[:p] + rng.choice([b for b in 'ACGT' if b != fragment[p]]) + fragment[p + 1:]
                kind += ':m1G9' if p + start == TRNA_MODIFIED[0] else ':m1A58'
            return fragment, kind
        if cls in END_FRAGMENTS:
            shares, (lo, hi) = END_FRAGMENTS[cls]
            kind = rng.choice(list(shares), p=list(shares.values()))
            length = min(self.lengths(lo, hi), n)
            start = {'5p_end': 0, '3p_end': n - length}.get(kind)
            if start is None:
                start = int(rng.integers(0, n - length + 1))
            return seq[start:start + length], kind
        if cls == 'piRNA':
            length = self.lengths(24, 32)
            if n < length:
                return None
            if rng.random() < float(self.spec['pirna_5u']):
                ts = [m.start() for m in re.finditer('T', seq[:n - length + 1])]
                start = int(rng.choice(ts)) if ts else 0
            else:
                start = int(rng.integers(0, n - length + 1))
            return seq[start:start + length], '5U' if seq[start] == 'T' else 'other'
        if cls == 'precursor':
            free = []  # stretches of the hairpin outside its matures
            last = 0
            for a, b in sorted(extra) + [(n, n)]:
                if a - last >= 16:
                    free.append((last, a))
                last = max(last, b)
            if not free:
                return None
            a, b = free[int(rng.integers(len(free)))]
            length = min(self.lengths(), b - a)
            start = int(rng.integers(a, b - length + 1))
            return seq[start:start + length], 'loop' if 0 < start and start + length < n else 'flank'
        length = self.lengths()
        if n < length:
            return None
        start = int(rng.integers(0, n - length + 1))
        return seq[start:start + length], 'fragment'

    def acceptable(self, unit, fragment):
        if not self.min_length <= len(fragment) <= MAX_LENGTH or set(fragment) - set('ACGT') or fragment in self.forbidden:
            return False
        if unit == 'precursor-mirbase':
            return True  # outside the matures by construction
        return not any(fragment[i:i + OVERLAP] in self.kmers for i in range(len(fragment) - OVERLAP + 1))

    def draw_catalogue(self, unit):
        """[(fragment, origin)] and their weights: up to `fragments` distinct fragments of a unit."""
        mols, wanted = self.molecules[unit], int(self.spec['fragments'])
        seen, out = set(), []
        for _ in range(50 * wanted):
            if len(out) >= wanted:
                break
            name, seq, extra = mols[int(self.rng.integers(len(mols)))]
            drawn = self.fragment(unit, name, seq, extra)
            if drawn is None or drawn[0] in seen or not self.acceptable(unit, drawn[0]):
                continue
            seen.add(drawn[0])
            out.append((drawn[0], '{}:{}'.format(name, drawn[1])))
        if not out:
            raise ValueError('decoys: no acceptable fragment of {}'.format(unit))
        weights = np.exp(self.rng.normal(0, float(self.spec['weight_sd']), len(out)))
        return out, weights / weights.sum()

    def find_embedded(self, universe):
        """{class: {miRNA representative: [molecule names]}}: the organism's miRBase matures that lie exactly inside a
        molecule of an EMBEDDING class's declared units."""
        return embedded({u: self.molecules[u] for u in self.units},
                        {rep: m['sequence'] for rep, m in universe.mirnas.items() if m['loci']})

    def draw_overlaps(self, excluded=()):
        """{class: ([(fragment, origin, source)], weights)}: up to `fragments` distinct fragments per EMBEDDING class of the
        host molecules of its embedded matures (not those `excluded`, the carrier), each covering at least half of the
        mature and no read a miRNA could give; origin `<molecule>:overlap:<mature representative>`."""
        out = {}
        for cls in sorted(self.embedded):
            hosts = []
            for unit in self.units:
                if unit_class(unit) != cls:
                    continue
                for name, seq, _ in self.molecules[unit]:
                    for rep, names in self.embedded[cls].items():
                        if rep in excluded or name not in names:
                            continue
                        mature = self.matures[rep]
                        at = seq.find(mature)
                        while at >= 0:
                            hosts.append((name, seq, rep, at, at + len(mature), unit.rsplit('-', 1)[1]))
                            at = seq.find(mature, at + 1)
            if not hosts:
                continue
            wanted, seen, found = int(self.spec['fragments']), set(), []
            for _ in range(50 * wanted):
                if len(found) >= wanted:
                    break
                name, seq, rep, a, b, source = hosts[int(self.rng.integers(len(hosts)))]
                length, cover = self.lengths(), (b - a + 1) // 2
                low, high = max(0, a + cover - length), min(len(seq) - length, b - cover)
                if high < low:
                    continue
                start = int(self.rng.integers(low, high + 1))
                fragment = seq[start:start + length]
                if fragment in seen or fragment in self.forbidden or set(fragment) - set('ACGT'):
                    continue
                seen.add(fragment)
                found.append((fragment, '{}:overlap:{}'.format(name, rep), source))
            if found:
                weights = np.exp(self.rng.normal(0, float(self.spec['weight_sd']), len(found)))
                out[cls] = (found, weights / weights.sum())
        return out

    def draw(self, unit, n):
        """{catalogue index: reads} of n reads of a unit."""
        drawn = self.rng.multinomial(n, self.catalogue[unit][1])
        return {i: int(k) for i, k in enumerate(drawn) if k}

    def mixed(self, mirna_reads, nuisance):
        """[(class, source, origin, read, count, error positions)] of a mixed sample with mirna_reads miRNA reads."""
        rows = []
        factor = 2 ** (float(self.spec['nuisance']) * nuisance)
        for unit in self.units:
            n = int(self.rng.poisson(float(self.spec['load'][unit]) * mirna_reads * factor))
            cls, source = unit.rsplit('-', 1)
            for i, k in self.draw(unit, n).items():
                fragment, origin = self.catalogue[unit][0][i]
                rows += [(cls, source, origin, read, c, hits) for read, c, hits in self.reads.with_errors(fragment, k)]
        return rows

    def probe(self, unit, n):
        """[(class, source, origin, read, count, ())] of a probe sample's n reads of one unit, without errors."""
        cls, source = unit.rsplit('-', 1)
        return [(cls, source, self.catalogue[unit][0][i][1], self.catalogue[unit][0][i][0], k, ()) for i, k in self.draw(unit, n).items()]
