"""Genome-wide synthetic human small-RNA reads of known origin (#412), generated from the miRBase files a declaration
fills (`synthetic.mode: genome`); the reduced-reference generator is synthetic.py.

What it draws (all from the declaration's fixed seed):

- **miRNAs.** `mirnas` annotated matures of the organism. Accessions with the same mature sequence count as one miRNA,
  named by their first accession in sorted order, as the evaluation collapses them (biological.collapse). Pairs one
  substitution apart are taken first, kept at a 100:1 ratio either way, or with the partner left out
  (`near_identical`). Each pair is a family of its own: the members' other near-identical relatives are absent and are
  relatives of no other chosen miRNA, so a read on a relative has one source (in large families such as miR-548 most
  pairs are left out). Then matures at several loci, then the rest at random, both without near-identical relatives.
- **Abundances.** Expected reads per million of a mixed sample, by rank: a Zipf curve, then a log-uniform low tail
  (`abundance`); ranks are assigned to the miRNAs at random, a pair's low member taking its partner's abundance / ratio.
- **Counts.** Mixed samples in two groups (`groups`). log2 mu_ij = log2 base_i + u_i1 v_j1 + u_i2 v_j2:
  - factor 1 is the group contrast (v = -1/2, +1/2), with sparse loadings, the true log2 fold changes (`structure.group`);
  - factor 2 is a nuisance (standard normal, correlated with the group), whose loadings miRNAs of one miRBase cluster
    share (`coregulate`).
  The mu are scaled to the sample's reads (lognormal), and counts are negative binomial with dispersion
  asymptotic + extra / expected count (`dispersion`).
- **Reads of a miRNA.** Canonical, 3' trimmed by 1-2 nt, 3' extended by 1-2 hairpin bases, a non-templated A or U, or
  the 5' end moved by one (trimmed, or one hairpin base added). Each miRNA draws its own shares around `isomirs`. A
  class whose read would equal another mature, or a read of any class that another miRNA's loci could give, or that
  the hairpin cannot give, is left out for that miRNA. An added A or U that any of the miRNA's hairpins has next (or
  might have: the mature ends its hairpin) is templated, so that class is left out too.
- **Sequencing errors** (`errors`): substitutions at a per-base rate, raised at the first sequenced base (the insert's
  first base; with UMIs a UMI base, whose errors are not modelled); a substitution that would give another mature or a
  read any miRNA could give is drawn again. Base qualities (`qualities`) are binned NextSeq values, low on substituted
  bases.
- **Random reads** (`random_reads`).
- **Probe samples** (`probes`), each a canonical carrier (miRNAs outside the truth set, without near-identical
  relatives) plus one read class of other miRNAs: each isomiR class, or canonical reads with at least one error. Each
  class has `replicates` samples, and a control has the carrier only. Probe miRNAs are truth miRNAs, so they are also
  present, canonical, in the mixed samples.
- **Decoys** (`decoys`, #412 PR 2; decoys.py): reads of other RNA classes in the mixed samples, a probe sample per class
  and source (`P_decoy_<unit>_<replicate>`, kind `decoy`), probe samples of the miRBase matures that lie inside
  other-class molecules (`P_embedded_<class>_<replicate>`) and of fragments of their host molecules overlapping them
  (`P_overlap_<class>_<replicate>`, kind `overlap`), all drawn from a stream of their own after the miRNA reads, so the
  miRNA reads are the same with or without decoys.

Reads are the insert followed by the kit's adapter (with `umi` random nucleotides at both ends of the insert when the
kit has them). The truth has a row per sample, class and read sequence (with its number of substitutions); the design
tables give each sample's group and factor scores, and each miRNA's abundance, loadings, fold change, dispersion,
cluster and strata.
"""
import collections
import gzip
import math
from pathlib import Path
import re

import numpy as np

import decoys

ISOMIR_CLASSES = ('3p_trim_1', '3p_trim_2', '3p_templated_1', '3p_templated_2', '3p_add_A', '3p_add_U', '5p_shift')
READ_CLASSES = ('canonical',) + ISOMIR_CLASSES
PROBE_CLASSES = ISOMIR_CLASSES + ('errors',)
TRUTH_COLUMNS = ('sample', 'class', 'mature_name', 'mature_id', 'sequence', 'count', 'errors', 'source', 'origin')
SAMPLE_COLUMNS = ('sample', 'kind', 'group', 'probe_class', 'replicate', 'reads', 'factor_group', 'factor_nuisance')
MIRNA_COLUMNS = ('mature_id', 'mature_name', 'accessions', 'sequence', 'base_per_million', 'log2fc', 'de', 'nuisance_loading',
                 'dispersion', 'cluster', 'loci', 'near_identical', 'near_role', 'probe_classes', 'inside')
MIN_LENGTH = 15  # the filter's minimum read length (every kit)


def read_fasta(path, prefix):
    """{name: (accession, DNA sequence)} of a miRBase FASTA (mature.fa, hairpin.fa) for names starting with prefix."""
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                fields = line[1:].split()
                name = fields[0] if fields[0].startswith(prefix) else None
                if name:
                    records[name] = [fields[1], []]
            elif name:
                records[name][1].append(line.strip().upper().replace('U', 'T'))
    return {n: (a, ''.join(s)) for n, (a, s) in records.items()}


def read_gff(path):
    """(mature loci {accession: [(chrom, start, end, strand, hairpin accession)]}, hairpin loci {accession: (chrom,
    start, end, strand)}) from a miRBase GFF3."""
    matures, hairpins = collections.defaultdict(list), {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) != 9:
                continue
            attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            locus = (fields[0], int(fields[3]), int(fields[4]), fields[6])
            if fields[2] == 'miRNA':
                matures[attrs.get('Alias', attrs['ID'])].append(locus + (attrs.get('Derives_from', ''),))
            elif fields[2] == 'miRNA_primary_transcript':
                hairpins[attrs.get('Alias', attrs['ID'])] = locus
    return dict(matures), hairpins


def one_apart(sequences):
    """Pairs of distinct sequences of equal length that differ at one position."""
    pairs = set()
    by_length = collections.defaultdict(list)
    for s in sequences:
        by_length[len(s)].append(s)
    for length, group in by_length.items():
        masked = collections.defaultdict(list)
        for s in group:
            for i in range(length):
                masked[s[:i] + '.' + s[i + 1:]].append(s)
        for members in masked.values():
            members = sorted(set(members))
            pairs.update((a, b) for i, a in enumerate(members) for b in members[i + 1:])
    return sorted(pairs)


def clusters(hairpin_loci, distance):
    """{hairpin accession: cluster number}: hairpins on one chromosome and strand within `distance` of the previous one."""
    out, number = {}, 0
    by_place = collections.defaultdict(list)
    for accession, (chrom, start, end, strand) in hairpin_loci.items():
        by_place[(chrom, strand)].append((start, end, accession))
    for place in sorted(by_place):
        last_end = None
        for start, end, accession in sorted(by_place[place]):
            if last_end is None or start - last_end > distance:
                number += 1
            out[accession] = number
            last_end = max(end, last_end or end)
    return out


class Universe:
    """The organism's miRNAs as the generator and the evaluation see them: one miRNA per mature sequence (its first
    accession sorted), its loci, hairpins (the flanks of the mature in each; the first gives the templated reads),
    cluster and near-identical partners."""

    def __init__(self, mature_fa, hairpin_fa, gff, prefix, cluster_distance=10000):
        mature = read_fasta(mature_fa, prefix)
        hairpins = read_fasta(hairpin_fa, prefix)
        loci, hairpin_loci = read_gff(gff)
        hairpin_seq = {accession: seq for accession, seq in hairpins.values()}
        by_seq = collections.defaultdict(list)
        for name, (accession, seq) in mature.items():
            by_seq[seq].append((accession, name))
        self.mirnas = {}            # representative accession -> dict
        self.alias = {}             # accession -> representative
        for seq, members in by_seq.items():
            members.sort()
            rep = members[0][0]
            accessions = [a for a, _ in members]
            self.alias.update((a, rep) for a in accessions)
            places = [l for a in accessions for l in loci.get(a, [])]
            parents = list(dict.fromkeys(l[4] for l in places if l[4] in hairpin_seq))
            self.mirnas[rep] = {'name': members[0][1], 'accessions': accessions, 'sequence': seq, 'loci': len(places),
                                'hairpin': parents[0] if parents else '', 'hairpins': parents, 'cluster': 0}
        cluster_of = clusters(hairpin_loci, cluster_distance)
        for rep, m in self.mirnas.items():
            m['cluster'] = cluster_of.get(m['hairpin'], 0)
            m['flanks'] = []  # (upstream, downstream) of the mature in each hairpin, None where it is not found
            for parent in m['hairpins']:
                hairpin = hairpin_seq[parent]
                offset = hairpin.find(m['sequence'])
                m['flanks'].append((hairpin[:offset], hairpin[offset + len(m['sequence']):]) if offset >= 0 else None)
            m['upstream'], m['downstream'] = m['flanks'][0] if m['flanks'] and m['flanks'][0] else (None, None)
        rep_of_seq = {m['sequence']: rep for rep, m in self.mirnas.items()}
        self.near = collections.defaultdict(set)
        for a, b in one_apart(rep_of_seq):
            self.near[rep_of_seq[a]].add(rep_of_seq[b])
            self.near[rep_of_seq[b]].add(rep_of_seq[a])
        self.mature_sequences = set(rep_of_seq)

    def annotated(self):
        return sorted(rep for rep, m in self.mirnas.items() if m['loci'])


def variant(m, read_class, shift_sign):
    """The read of a miRNA (Universe entry) of a class, or None when the hairpin cannot give it."""
    seq, up, down = m['sequence'], m['upstream'], m['downstream']
    if read_class == 'canonical':
        return seq
    if read_class.startswith('3p_trim_'):
        return seq[:-int(read_class[-1])]
    if read_class.startswith('3p_templated_'):
        k = int(read_class[-1])
        return seq + down[:k] if down is not None and len(down) >= k else None
    if read_class.startswith('3p_add_'):
        base = read_class[-1].replace('U', 'T')
        # an added base that any hairpin has next is templated: not this class; nor where a hairpin's next base is unknown
        nexts = [f[1][:1] if f else '' for f in m['flanks']]
        return seq + base if nexts and all(b not in ('', base) for b in nexts) else None
    if read_class == '5p_shift':
        if shift_sign > 0:
            return seq[1:]
        return up[-1] + seq if up else None
    raise ValueError('unknown read class {!r}'.format(read_class))


def possible_reads(m):
    """Every read of a miRNA's classes that any of its loci could give: both 5' shifts, either added base (templated or
    not), and each hairpin's templated bases."""
    seq = m['sequence']
    out = {seq, seq[:-1], seq[:-2], seq[1:], seq + 'A', seq + 'T'}
    for up, down in filter(None, m['flanks']):
        out.update(seq + down[:k] for k in (1, 2) if len(down) >= k)
        if up:
            out.add(up[-1] + seq)
    return {s for s in out if len(s) >= MIN_LENGTH}


def abundances(n, spec, rng):
    """Expected reads per million for ranks 1..n: Zipf, then a log-uniform tail."""
    tail = spec.get('tail') or {}
    tail_n = min(int(tail.get('mirnas', 0)), n)
    low, high = tail.get('reads', [1, 50])
    tail_reads = np.sort(np.exp(rng.uniform(math.log(low), math.log(high), tail_n)))[::-1]
    head = np.arange(1, n - tail_n + 1, dtype=float) ** -float(spec['zipf'])
    head_reads = head / head.sum() * (1e6 - tail_reads.sum())
    return np.concatenate([head_reads, tail_reads])


def choose(universe, spec, rng):
    """(truth miRNAs, their abundances {rep: per million}, near roles {rep: (partner, role)}, absent relatives).

    A pair (or a single member with its partner absent) is a family of its own: no member or relative of it is a
    relative of another chosen miRNA, and the relatives other than the partner are absent, so a read on a relative has
    one source. The other miRNAs have no near-identical relatives."""
    n = int(spec['mirnas'])
    near_spec = spec.get('near_identical') or {}
    ratio, alone = float(near_spec.get('ratio', 100)), float(near_spec.get('partner_absent', 0.33))
    annotated = set(universe.annotated())
    pairs = sorted({tuple(sorted((a, b))) for a in annotated for b in universe.near[a] if b in annotated})
    chosen, roles, absent = [], {}, set()
    for a, b in rng.permutation(np.array(pairs, dtype=object)).tolist() if pairs else []:
        if n - len(chosen) < 2:  # a pair is never cut
            break
        if ({a, b} | universe.near[a] | universe.near[b]) & (set(chosen) | absent):
            continue
        high, low = (a, b) if rng.random() < 0.5 else (b, a)
        if rng.random() < alone:
            roles[high] = ('', 'alone')
            members = [high]
        else:
            roles[high], roles[low] = (low, 'high'), (high, 'low')
            members = [high, low]
        chosen += members
        absent |= {p for r in members for p in universe.near[r] if p not in members}
    multi = [r for r in sorted(annotated) if universe.mirnas[r]['loci'] > 1 and not universe.near[r]]
    rest = [r for r in sorted(annotated) if universe.mirnas[r]['loci'] == 1 and not universe.near[r]]
    for r in rng.permutation(np.array(multi, dtype=object)).tolist() + rng.permutation(np.array(rest, dtype=object)).tolist():
        if len(chosen) >= n:
            break
        chosen.append(r)
    family = {r: frozenset({r, roles.get(r, ('', ''))[0]} - {''}) for r in chosen}
    for r in list(chosen) + sorted(absent):
        if len({family[c] for c in ({r} | universe.near[r]) if c in family}) > 1:
            raise RuntimeError('{} is near-identical to miRNAs of two chosen families'.format(universe.mirnas[r]['name']))
    lows = [r for r in chosen if roles.get(r, ('', ''))[1] == 'low']
    units = [r for r in chosen if r not in lows]
    values = abundances(len(units), spec['abundance'], rng)
    base = dict(zip(rng.permutation(np.array(units, dtype=object)).tolist(), values))
    for r in lows:
        base[r] = base[roles[r][0]] / ratio
    return chosen, base, roles, absent


def structure(chosen, base, universe, spec, groups, rng):
    """(sample scores [(group, v_group, v_nuisance)], loadings {rep: (log2 fold change, nuisance loading)})."""
    names = list(groups)
    if len(names) != 2 or min(int(v) for v in groups.values()) < 2:
        raise ValueError('synthetic.groups must name two groups of at least two samples each, got {}'.format(dict(groups)))
    labels = [g for g in names for _ in range(int(groups[g]))]
    v1 = np.array([-0.5 if g == names[0] else 0.5 for g in labels])
    st = spec['structure']
    rho = float(st['nuisance'].get('group_correlation', 0))
    z = (v1 - v1.mean()) / v1.std()
    v2 = rho * z + math.sqrt(1 - rho ** 2) * rng.standard_normal(len(labels))
    v2 = (v2 - v2.mean()) / v2.std()
    group = st['group']
    n_de = int(round(float(group['fraction']) * len(chosen)))
    include = set(group.get('include') or [])
    # The group design draws from a stream of its own (#453): a scenario that changes it keeps the base set's nuisance
    # loadings, which follow on the main stream.
    grng = np.random.default_rng(np.random.SeedSequence(int(spec.get('seed', 1)), spawn_key=(2,)))
    # `abundant` (#453's unbalanced scenario): the most abundant chosen miRNAs change first, so the composition shifts
    abundant = group.get('abundant') or 0
    if isinstance(abundant, bool) or not isinstance(abundant, int) or not 0 <= abundant <= n_de:
        raise ValueError('structure.group.abundant must be a whole number from 0 to the miRNAs that change ({}), got {!r}'.format(n_de, abundant))
    top = sorted(chosen, key=lambda r: (-base[r], r))[:abundant]
    rest = [r for r in chosen if r not in top]
    special = [r for r in rest if ('near_identical' in include and universe.near[r]) or
               ('multi_locus' in include and universe.mirnas[r]['loci'] > 1)]
    others = [r for r in rest if r not in special]
    k = min(len(special), (n_de - len(top)) // 4)
    special = grng.permutation(np.array(special, dtype=object)).tolist()
    others = grng.permutation(np.array(others, dtype=object)).tolist()
    de = top + special[:k] + others[:n_de - len(top) - k]
    de += special[k:k + n_de - len(de)]  # too few others: the rest from the special miRNAs, in their drawn order
    fold = {r: float(grng.choice(group['log2'])) for r in de}
    sd = float(st['nuisance'].get('sd', 0.5))
    shared = float((spec.get('coregulate') or {}).get('shared', 0))
    cluster_effect = {}
    nuisance = {}
    for r in chosen:
        c = universe.mirnas[r]['cluster']
        if c and c not in cluster_effect:
            cluster_effect[c] = rng.standard_normal()
        common = cluster_effect[c] if c else rng.standard_normal()
        nuisance[r] = sd * (math.sqrt(shared) * common + math.sqrt(1 - shared) * rng.standard_normal())
    scores = list(zip(labels, v1.tolist(), v2.tolist()))
    return scores, {r: (fold.get(r, 0.0), nuisance[r]) for r in chosen}


def counts(chosen, base, loadings, scores, spec, rng):
    """{(sample index, rep): count} of the mixed samples, their reads, and each miRNA's dispersion."""
    reads_spec = spec['reads']
    mean, cv = float(reads_spec['mean']), float(reads_spec.get('cv', 0))
    sigma = math.sqrt(math.log(1 + cv ** 2))
    library = mean * np.exp(rng.normal(-sigma ** 2 / 2, sigma, len(scores)))
    b = np.array([base[r] for r in chosen])
    u1 = np.array([loadings[r][0] for r in chosen])
    u2 = np.array([loadings[r][1] for r in chosen])
    v1 = np.array([s[1] for s in scores])
    v2 = np.array([s[2] for s in scores])
    log_mu = np.log2(b)[:, None] + np.outer(u1, v1) + np.outer(u2, v2)
    p = 2 ** log_mu
    p /= p.sum(axis=0)
    expected = p * library[None, :]
    disp = spec['dispersion']
    phi = float(disp['asymptotic']) + float(disp['extra']) / np.maximum(expected.mean(axis=1), 1e-9)
    shape = 1 / phi[:, None]
    lam = rng.gamma(shape, expected / shape)
    drawn = rng.poisson(lam)
    table = {(j, r): int(drawn[i, j]) for i, r in enumerate(chosen) for j in range(len(scores))}
    return table, library, dict(zip(chosen, phi.tolist()))


class Reads:
    """Read classes, errors and qualities of the generator; the sequences no generated read may become."""

    def __init__(self, universe, chosen, spec, rng, umi=0):
        self.universe, self.rng, self.umi = universe, rng, umi
        self.shares, self.reads = {}, {}
        self.sign = {r: (1 if rng.random() < 0.5 else -1) for r in chosen}
        mean = spec['isomirs']
        source = collections.defaultdict(set)  # read -> the miRNAs whose loci could give it in some class
        for rep, m in universe.mirnas.items():
            for s in possible_reads(m):
                source[s].add(rep)
        self.forbidden = set(universe.mature_sequences) | set(source)
        for r in chosen:
            m = universe.mirnas[r]
            reads = {}
            for c in READ_CLASSES:
                s = variant(m, c, self.sign[r])
                if s is None or len(s) < MIN_LENGTH:
                    continue
                # an isomiR equal to another mature, or to a read another miRNA could give, is left out (the canonical
                # read stays)
                if c != 'canonical' and (s in universe.mature_sequences or source[s] - {r}):
                    continue
                reads[c] = s
            self.reads[r] = reads
            weights = np.array([float(mean.get(c, 0)) if c in reads else 0.0 for c in READ_CLASSES])
            alpha = np.maximum(weights * 50, 1e-9)
            share = rng.dirichlet(alpha) * (weights > 0)
            self.shares[r] = share / share.sum()
        self.forbidden |= {s for reads in self.reads.values() for s in reads.values()}
        err = spec['errors']
        self.rate, self.first = float(err['rate']), float(err['first_base'])
        q = spec['qualities']
        self.q_correct = (np.array([int(k) for k in q['correct']]), np.array([float(v) for v in q['correct'].values()]))
        self.q_error = (np.array([int(k) for k in q['error']]), np.array([float(v) for v in q['error'].values()]))
        self.pool = {}

    def split(self, r, n):
        """{class: reads} of a miRNA's n reads."""
        drawn = self.rng.multinomial(n, self.shares[r])
        return {c: int(k) for c, k in zip(READ_CLASSES, drawn) if k}

    def error_rates(self, length):
        """Per-base substitution rates of an insert: the first cycle's raised rate falls on its first base, or with UMIs
        on the UMI."""
        return np.array([self.rate if self.umi else self.first] + [self.rate] * (length - 1))

    def substitute(self, seq, at_least_one=False):
        """(read, positions) with substitutions at the per-base rates, never another mature or miRNA read."""
        rates = self.error_rates(len(seq))
        for _ in range(100):
            hits = np.nonzero(self.rng.random(len(seq)) < rates)[0]
            if at_least_one and not len(hits):
                hits = np.array([int(self.rng.choice(len(seq), p=rates / rates.sum()))])
            read = list(seq)
            for i in hits:
                read[i] = self.rng.choice([b for b in 'ACGT' if b != seq[i]])
            read = ''.join(read)
            if not len(hits) or read not in self.forbidden:
                return read, hits.tolist()
        raise RuntimeError('no substitution of {} in 100 draws avoids the reads of other miRNAs'.format(seq))

    def with_errors(self, seq, n):
        """[(read, count, positions)]: n reads of seq, some with substitutions."""
        p_any = 1 - float(np.prod(1 - self.error_rates(len(seq))))
        k = int(self.rng.binomial(n, p_any))
        out = collections.Counter()
        for _ in range(k):
            read, hits = self.substitute(seq, at_least_one=True)
            out[(read, tuple(hits))] += 1
        rows = [(seq, n - k, ())] if n - k else []
        return rows + [(read, c, hits) for (read, hits), c in out.items()]

    def fork(self, rng):
        """The same error and quality model drawing from another stream (decoys.py), with qualities of its own."""
        other = object.__new__(Reads)
        other.__dict__.update(self.__dict__)
        other.rng, other.pool = rng, {}
        return other

    def quality(self, length, error_positions=()):
        """A quality string: correct bases from a pool of 256 drawn strings per length, substituted bases drawn anew."""
        if length not in self.pool:
            values, probs = self.q_correct
            drawn = self.rng.choice(values, size=(256, length), p=probs)
            self.pool[length] = [''.join(chr(33 + int(v)) for v in row) for row in drawn]
        q = self.pool[length][int(self.rng.integers(256))]
        if error_positions:
            ev, ep = self.q_error
            q = list(q)
            for i, v in zip(error_positions, self.rng.choice(ev, size=len(error_positions), p=ep)):
                q[i] = chr(33 + int(v))
            q = ''.join(q)
        return q


def write_fastq(path, sample, rows, adapter, umi, reads, more=()):
    """Gzipped FASTQ of one sample: rows [(insert, count, error positions)]; [umi]insert[umi]adapter. `more` [(rows,
    reads)] follow, their UMIs and qualities drawn from their own Reads (the decoys' stream)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(path, 'wt', compresslevel=1) as fh:
        chunk = []
        for part_rows, part_reads in [(rows, reads)] + list(more):
            rng = part_reads.rng
            for insert, count, positions in part_rows:
                for _ in range(count):
                    left = ''.join(rng.choice(list('ACGT'), umi)) if umi else ''
                    right = ''.join(rng.choice(list('ACGT'), umi)) if umi else ''
                    seq = left + insert + right + adapter
                    qual = part_reads.quality(len(seq), [p + umi for p in positions])
                    chunk.append('@{}_{}\n{}\n+\n{}\n'.format(sample, n, seq, qual))
                    n += 1
                    if len(chunk) >= 10000:
                        fh.write(''.join(chunk))
                        chunk = []
        fh.write(''.join(chunk))
    return n


def hairpin_records(hairpin_fa, prefix, universe):
    """[(name, sequence, [(start, end) of each mature in it])] of the organism's miRBase hairpins."""
    spans = collections.defaultdict(list)
    for m in universe.mirnas.values():
        for parent, flank in zip(m['hairpins'], m['flanks']):
            if flank:
                spans[parent].append((len(flank[0]), len(flank[0]) + len(m['sequence'])))
    return [(name, seq, sorted(spans.get(accession, []))) for name, (accession, seq) in sorted(read_fasta(hairpin_fa, prefix).items())]


def generate(spec, mature_fa, hairpin_fa, gff, prefix, adapter, umi, fastq_path, decoy_sources=None):
    """Write every sample's FASTQ (fastq_path(sample) -> Path) and return (truth rows, sample rows, miRNA rows).
    `decoy_sources` {source: path} are the filled files of `decoys` (decoys.SOURCES)."""
    rng = np.random.default_rng(int(spec.get('seed', 1)))
    universe = Universe(mature_fa, hairpin_fa, gff, prefix, int((spec.get('coregulate') or {}).get('distance', 10000)))
    chosen, base, roles, absent = choose(universe, spec, rng)
    groups = spec['groups']
    scores, loadings = structure(chosen, base, universe, spec, groups, rng)
    table, library, phi = counts(chosen, base, loadings, scores, spec, rng)
    # the isomiR design and the reads' errors draw from a stream of their own (#453): the counts' draws take a varying
    # number of random numbers, so a scenario would otherwise change every miRNA's isomiRs
    reads = Reads(universe, chosen, spec, np.random.default_rng(np.random.SeedSequence(int(spec.get('seed', 1)), spawn_key=(4,))), umi)
    decoy = None
    if spec.get('decoys'):  # after every miRNA draw, on a stream of its own
        decoy = decoys.Decoys(spec['decoys'], decoy_sources or {}, universe, hairpin_records(hairpin_fa, prefix, universe), reads,
                       int(spec.get('seed', 1)))
    names = {r: universe.mirnas[r]['name'] for r in universe.mirnas}
    truth, samples = [], []
    width = len(str(len(scores)))
    random_spec = spec.get('random_reads') or {}

    def random_rows(total):
        k = int(round(float(random_spec.get('per_million', 0)) * total / 1e6))
        low, high = random_spec.get('length', [18, 26])
        out = collections.Counter()
        for _ in range(k):
            s = ''.join(rng.choice(list('ACGT'), int(rng.integers(low, high + 1))))
            if s not in reads.forbidden:
                out[s] += 1
        return out

    for j, (label, f1, f2) in enumerate(scores):
        sample = 'M{:0{}d}'.format(j + 1, width)
        rows = []
        for r in chosen:
            for c, n in reads.split(r, table[(j, r)]).items():
                for read, k, hits in reads.with_errors(reads.reads[r][c], n):
                    rows.append((read, k, hits))
                    truth.append((sample, c, names[r], r, read, k, len(hits), 'mirbase', ''))
        for s, k in random_rows(library[j]).items():
            rows.append((s, k, ()))
            truth.append((sample, 'random', '', '', s, k, 0, '', ''))
        more = []
        if decoy:
            decoy_rows = decoy.mixed(sum(table[(j, r)] for r in chosen), f2)
            truth += [(sample, cls, '', '', read, k, len(hits), source, origin) for cls, source, origin, read, k, hits in decoy_rows]
            more = [([(read, k, hits) for _, _, _, read, k, hits in decoy_rows], decoy.reads)]
        n = write_fastq(fastq_path(sample), sample, rows, adapter, umi, reads, more)
        samples.append((sample, 'mixed', label, '', '', n, round(f1, 4), round(f2, 4)))

    probes = spec.get('probes') or {}
    probe_of = collections.defaultdict(list)
    if probes:
        # the probe samples draw from a stream of their own (#453), so they are the same whatever the mixed samples' design
        prng = np.random.default_rng(np.random.SeedSequence(int(spec.get('seed', 1)), spawn_key=(3,)))
        preads = reads.fork(prng)
        taken = set(chosen) | absent
        carrier_pool = [r for r in universe.annotated() if r not in taken and not universe.near[r]
                        and len(universe.mirnas[r]['accessions']) == 1]
        carrier = prng.permutation(np.array(carrier_pool, dtype=object)).tolist()[:int(probes['carrier']['mirnas'])]
        carrier_reads = int(probes['carrier']['reads'])
        floor = 20.0  # per million: present, canonical, in every mixed sample (miRge3's project rule)
        eligible = [r for r in chosen if not universe.near[r] and base[r] >= floor]
        probe_reads, probe_n = int(probes['probe']['reads']), int(probes['probe']['mirnas'])
        plan = [('control', 1, [])]
        for c in PROBE_CLASSES:
            pool = [r for r in eligible if (c == 'errors' or c in reads.reads[r])]
            members = prng.permutation(np.array(pool, dtype=object)).tolist()[:probe_n]
            for r in members:
                probe_of[r].append(c)
            plan.append((c, int(probes.get('replicates', 1)), members))
        for c, replicates, members in plan:
            for rep in range(1, replicates + 1):
                sample = 'P_control' if c == 'control' else 'P_{}_{}'.format(c, rep)
                rows = []
                for r in carrier:
                    rows.append((universe.mirnas[r]['sequence'], carrier_reads, ()))
                    truth.append((sample, 'carrier', names[r], r, universe.mirnas[r]['sequence'], carrier_reads, 0, 'mirbase', ''))
                for r in members:
                    if c == 'errors':
                        out = collections.Counter()
                        for _ in range(probe_reads):
                            read, hits = preads.substitute(universe.mirnas[r]['sequence'], at_least_one=True)
                            out[(read, tuple(hits))] += 1
                        for (read, hits), k in out.items():
                            rows.append((read, k, hits))
                            truth.append((sample, 'canonical', names[r], r, read, k, len(hits), 'mirbase', ''))
                    else:
                        rows.append((reads.reads[r][c], probe_reads, ()))
                        truth.append((sample, c, names[r], r, reads.reads[r][c], probe_reads, 0, 'mirbase', ''))
                n = write_fastq(fastq_path(sample), sample, rows, adapter, umi, preads)
                samples.append((sample, 'control' if c == 'control' else 'probe', '', '' if c == 'control' else c,
                                rep if c != 'control' else '', n, '', ''))
        for r in carrier:
            probe_of[r].append('carrier')
        if decoy:
            replicates = int(probes.get('replicates', 1))
            carrier_rows = [(universe.mirnas[r]['sequence'], carrier_reads, ()) for r in carrier]

            def carrier_truth(sample):
                return [(sample, 'carrier', names[r], r, universe.mirnas[r]['sequence'], carrier_reads, 0, 'mirbase', '') for r in carrier]
            decoy_reads = int((spec['decoys'].get('probe') or {}).get('reads', 0))
            for unit in decoy.units:  # (a): one unit's reads and the carrier
                for rep in range(1, replicates + 1):
                    sample = 'P_decoy_{}_{}'.format(unit, rep)
                    rows = decoy.probe(unit, decoy_reads)
                    truth += carrier_truth(sample) + [(sample, cls, '', '', read, k, 0, source, origin) for cls, source, origin, read, k, _ in rows]
                    n = write_fastq(fastq_path(sample), sample, carrier_rows + [(read, k, ()) for _, _, _, read, k, _ in rows],
                                    adapter, umi, decoy.reads)
                    samples.append((sample, 'decoy', '', unit, rep, n, '', ''))
            for cls in sorted(decoy.embedded):  # (b): canonical reads of the miRBase matures inside this class's molecules
                pool = sorted(r for r in decoy.embedded[cls] if r not in carrier)
                members = decoy.rng.permutation(np.array(pool, dtype=object)).tolist()[:probe_n] if pool else []
                if not members:
                    continue
                for r in members:
                    probe_of[r].append('embedded_' + cls)
                for rep in range(1, replicates + 1):
                    sample = 'P_embedded_{}_{}'.format(cls, rep)
                    rows = [(universe.mirnas[r]['sequence'], probe_reads, ()) for r in members]
                    truth += carrier_truth(sample) + [(sample, 'canonical', names[r], r, universe.mirnas[r]['sequence'], probe_reads, 0, 'mirbase',
                                                       ','.join(decoy.embedded[cls][r][:3])) for r in members]
                    n = write_fastq(fastq_path(sample), sample, carrier_rows + rows, adapter, umi, decoy.reads)
                    samples.append((sample, 'probe', '', 'embedded_' + cls, rep, n, '', ''))
            for cls, (catalogue, weights) in sorted(decoy.draw_overlaps(set(carrier)).items()):
                # fragments of the host molecules overlapping an embedded mature (no read a miRNA could give)
                for rep in range(1, replicates + 1):
                    sample = 'P_overlap_{}_{}'.format(cls, rep)
                    drawn = [(catalogue[i], int(k)) for i, k in enumerate(decoy.rng.multinomial(decoy_reads, weights)) if k]
                    truth += carrier_truth(sample) + [(sample, cls, '', '', fragment, k, 0, source, origin)
                                                      for (fragment, origin, source), k in drawn]
                    n = write_fastq(fastq_path(sample), sample, carrier_rows + [(fragment, k, ()) for (fragment, _, _), k in drawn],
                                    adapter, umi, decoy.reads)
                    samples.append((sample, 'overlap', '', cls, rep, n, '', ''))

    mirnas = []
    for r in sorted(set(chosen) | set(probe_of), key=lambda x: names[x]):
        m = universe.mirnas[r]
        fold, nuis = loadings.get(r, (0.0, 0.0))
        partner, role = roles.get(r, ('', ''))
        mirnas.append((r, names[r], ','.join(m['accessions']), m['sequence'], round(float(base.get(r, 0.0)), 4),
                       fold, int(fold != 0), round(nuis, 4), round(phi.get(r, 0.0), 6), m['cluster'], m['loci'],
                       ','.join(sorted(universe.near[r])), role, ','.join(probe_of.get(r, [])),
                       ';'.join('{}:{}'.format(cls, ','.join(mols[:3])) for cls, by in sorted((decoy.embedded if decoy else {}).items())
                                for rr, mols in by.items() if rr == r)))
    return truth, samples, mirnas


def write_table(path, columns, rows):
    Path(path).write_text('\t'.join(columns) + '\n' + ''.join('\t'.join(map(str, row)) + '\n' for row in rows))
