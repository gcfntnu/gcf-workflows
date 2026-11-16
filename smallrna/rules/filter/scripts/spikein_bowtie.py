#!/usr/bin/env python
"""Per-spike-in read counts of the bowtie spike-in route (filter.spikein.quantifier: bowtie, #336).

    spikein_bowtie.py fasta SET.fa [SET.fa ...] > spikein.fa
        The configured sets as one bowtie reference: names cut at the first whitespace, sequences in upper case DNA
        (U to T). Fails on a name used twice, within or across sets.
    spikein_bowtie.py counts --fasta spikein.fa --hits SAMPLE.hits > SAMPLE.tsv
        One table per sample: every spike-in of the reference in its order, with the reads matching it (0 if none),
        then `__ambiguous`, reads matching more than one spike-in. HITS is bowtie's default output with -k 2: one line
        per alignment (read name, strand, reference name, ...).
    spikein_bowtie.py matrix SAMPLE.tsv [SAMPLE.tsv ...] > spikein_counts.tsv
        The per-sample tables side by side, one column per sample (named by the file name without .tsv); fails when
        the tables do not list the same spike-ins.
"""
import argparse
import collections
from pathlib import Path
import sys

AMBIGUOUS = '__ambiguous'


def read_fasta(path):
    """[(name, sequence)] of a FASTA file; the name is the header up to the first whitespace."""
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                records.append([line[1:].split()[0], ''])
            elif line:
                records[-1][1] += line
    return [tuple(record) for record in records]


def combine(paths):
    records, seen = [], {}
    for path in paths:
        for name, sequence in read_fasta(path):
            if name in seen:
                raise ValueError('spike-in {!r} is in both {} and {}'.format(name, seen[name], path) if seen[name] != path
                                 else 'spike-in {!r} is listed twice in {}'.format(name, path))
            seen[name] = path
            records.append((name, sequence.upper().replace('U', 'T')))
    return records


def count(names, hits):
    """Reads per spike-in from bowtie's default output; a read matching several spike-ins counts as ambiguous."""
    targets = collections.defaultdict(set)
    for line in hits:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4 or fields[1] != '+':
            raise ValueError('not a forward-strand bowtie alignment line: {!r}'.format(line))
        targets[fields[0]].add(fields[2])
    counts = collections.Counter()
    for read, spikeins in targets.items():
        if len(spikeins) > 1:
            counts[AMBIGUOUS] += 1
            continue
        name = next(iter(spikeins))
        if name not in names:
            raise ValueError('read {} matches {!r}, which is not in the spike-in reference'.format(read, name))
        counts[name] += 1
    return [(name, counts[name]) for name in list(names) + [AMBIGUOUS]]


def read_table(path):
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header != ['spikein', 'count']:
            raise ValueError('{}: not a spike-in count table'.format(path))
        return [(name, int(value)) for name, value in (line.rstrip('\n').split('\t') for line in fh)]


def matrix(paths):
    tables = {Path(p).name[:-len('.tsv')] if p.endswith('.tsv') else Path(p).name: read_table(p) for p in paths}
    rows = [name for name, _ in next(iter(tables.values()))]
    for sample, table in tables.items():
        if [name for name, _ in table] != rows:
            raise ValueError('{}: its spike-ins differ from the other samples\''.format(sample))
    lines = ['\t'.join(['spikein'] + list(tables))]
    for i, name in enumerate(rows):
        lines.append('\t'.join([name] + [str(table[i][1]) for table in tables.values()]))
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('fasta').add_argument('fasta', nargs='+')
    counts_parser = sub.add_parser('counts')
    counts_parser.add_argument('--fasta', required=True)
    counts_parser.add_argument('--hits', required=True)
    sub.add_parser('matrix').add_argument('tables', nargs='+')
    args = parser.parse_args()
    try:
        if args.command == 'fasta':
            sys.stdout.write(''.join('>{}\n{}\n'.format(name, sequence) for name, sequence in combine(args.fasta)))
        elif args.command == 'counts':
            names = [name for name, _ in read_fasta(args.fasta)]
            with open(args.hits) as fh:
                rows = count(names, fh)
            sys.stdout.write('spikein\tcount\n' + ''.join('{}\t{}\n'.format(name, n) for name, n in rows))
        else:
            sys.stdout.write(matrix(args.tables))
    except (OSError, ValueError, IndexError) as error:
        sys.exit('spikein_bowtie.py: {}'.format(error))


if __name__ == '__main__':
    main()
