"""Synthetic small-RNA reads of known origin, and their truth table.

A dataset declaration lists miRNAs with read counts (`mirnas`), optionally miRNAs drawn at random from those annotated
on the reference's chromosomes (`random_mirnas`), isomiRs of listed miRNAs (`isomirs`) and random reads
(`random_reads`). Every sample gets the same reads. Reads are the insert followed by the kit's adapter, with `umi`
random nucleotides at both ends of the insert when the kit uses them, so the workflow's own filtering trims them.

The truth table has one row per read sequence: sample, class (`canonical`, `isomir`, `random`), mature miRNA name and
accession (empty for random reads), insert sequence and count.
"""
import gzip
from pathlib import Path
import random
import re

ISOMIR_CHANGES = ('3p_trim_1', '3p_add_A', '5p_trim_1')
TRUTH_COLUMNS = ('sample', 'class', 'mature_name', 'mature_id', 'sequence', 'count')


def read_mature(path, prefix):
    """miRBase mature.fa: {name: (accession, DNA sequence)} for names starting with prefix (e.g. 'mmu-')."""
    mature, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                fields = line[1:].split()
                name = fields[0] if fields[0].startswith(prefix) else None
                if name:
                    mature[name] = [fields[1], '']
            elif name:
                mature[name][1] += line.strip().upper().replace('U', 'T')
    return {k: tuple(v) for k, v in mature.items()}


def gff_mature_names(path):
    """Names of the mature miRNAs annotated in a miRBase GFF3."""
    names = set()
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if not line.startswith('#') and len(fields) == 9 and fields[2] == 'miRNA':
                match = re.search(r'(?:^|;)Name=([^;]+)', fields[8])
                if match:
                    names.add(match.group(1))
    return names


def isomir(sequence, change):
    if change == '3p_trim_1':
        return sequence[:-1]
    if change == '3p_add_A':
        return sequence + 'A'
    if change == '5p_trim_1':
        return sequence[1:]
    raise ValueError('unknown isomiR change {!r}; choose from {}'.format(change, ', '.join(ISOMIR_CHANGES)))


def truth_rows(spec, mature, annotated):
    """Rows (class, name, accession, sequence, count) for one sample; deterministic for a given spec['seed']."""
    rng = random.Random(spec.get('seed', 1))
    rows = []
    for name, count in (spec.get('mirnas') or {}).items():
        if name not in mature:
            raise ValueError('mirnas: {} is not in the mature miRNA FASTA'.format(name))
        rows.append(('canonical', name, mature[name][0], mature[name][1], int(count)))
    drawn = spec.get('random_mirnas')
    if drawn:
        pool = sorted(n for n in annotated if n in mature and n not in (spec.get('mirnas') or {}))
        if drawn['n'] > len(pool):
            raise ValueError('random_mirnas: {} requested, only {} available'.format(drawn['n'], len(pool)))
        low, high = drawn['count']
        for name in rng.sample(pool, drawn['n']):
            rows.append(('canonical', name, mature[name][0], mature[name][1], rng.randint(low, high)))
    for item in spec.get('isomirs') or []:
        name = item['mirna']
        if name not in mature:
            raise ValueError('isomirs: {} is not in the mature miRNA FASTA'.format(name))
        rows.append(('isomir', name, mature[name][0], isomir(mature[name][1], item['change']), int(item['count'])))
    reads = spec.get('random_reads')
    if reads:
        for _ in range(reads['n']):
            rows.append(('random', '', '', ''.join(rng.choice('ACGT') for _ in range(reads.get('length', 22))), 1))
    sequences = [row[3] for row in rows]
    if len(set(sequences)) != len(sequences):
        raise ValueError('synthetic reads repeat a sequence; adjust mirnas, isomirs or the seed')
    return rows


def write_sample(path, sample, rows, adapter, umi, rng):
    """Write one sample's reads as gzipped FASTQ: [umi]insert[umi]adapter."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(path, 'wt') as fh:
        for kind, name, _, sequence, count in rows:
            for _ in range(count):
                left = ''.join(rng.choice('ACGT') for _ in range(umi))
                right = ''.join(rng.choice('ACGT') for _ in range(umi))
                read = left + sequence + right + adapter
                fh.write('@{}_{}_{}\n{}\n+\n{}\n'.format(sample, kind, n, read, 'I' * len(read)))
                n += 1
    return n


def generate(spec, mature_fa, mirna_gff, prefix, samples, adapter, umi, fastq_path):
    """Write every sample's FASTQ (fastq_path(sample) -> Path) and return the truth table rows."""
    mature = read_mature(mature_fa, prefix)
    annotated = gff_mature_names(mirna_gff) if mirna_gff else set(mature)
    rows = truth_rows(spec, mature, annotated)
    truth = []
    for i, sample in enumerate(samples):
        write_sample(fastq_path(sample), sample, rows, adapter, umi, random.Random('{}-{}'.format(spec.get('seed', 1), i)))
        truth += [(sample,) + row for row in rows]
    return truth


def write_truth(path, truth):
    Path(path).write_text('\t'.join(TRUTH_COLUMNS) + '\n' + ''.join('\t'.join(map(str, row)) + '\n' for row in truth))
