#!/usr/bin/env python
"""Counts of the miraligner route (smallrna quant.method: miraligner, #426), from miraligner's `.mirna` output, as an
earlier in-house small-RNA pipeline counted them.

count: one sample. Each annotation row (a collapsed sequence annotated to a mature miRNA) is weighted
`freq / ambiguity`, the sequence's reads divided equally between the mature miRNAs it was annotated to. The isomiR
table has every row, keyed `<miRNA>_<sequence>_<mism>;<add>;<t5>;<t3>` as in the pipeline; the mature table sums, per
miRNA, the rows without a substitution (`mism` 0), as the pipeline's reader did. Statistics: reads and unique sequences
in miraligner's input, reads annotated to a mature miRNA, reads with an annotation without a substitution, and the
summed mature counts.

join: the samples' mature tables as one table, miRNA and one column per sample (0 where a sample lacks the miRNA).

Counts are written as integers where they are whole, else with Python's shortest round-trip representation.

Usage:
  miraligner_counts.py count --reads TAB --mirna MIRNA --mature TSV --isomir TSV --stats TSV
  miraligner_counts.py join --samples S1 [S2 ...] --tables MATURE_TSV [MATURE_TSV ...]
"""
import argparse
import math
import sys

COLUMNS = ('seq', 'freq', 'mir', 'mism', 'add', 't5', 't3', 'DB', 'ambiguity')


def fail(message):
    sys.exit('miraligner_counts.py: ' + message)


def number(value):
    return str(int(value)) if float(value).is_integer() else repr(float(value))


def read_input(path):
    """(reads, sequences) of miraligner's tabular input (sequence, count)."""
    reads = sequences = 0
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2 or not fields[1].isdigit():
                fail('{}:{}: expected a sequence and a count, got {!r}'.format(path, i, line[:100]))
            reads += int(fields[1])
            sequences += 1
    return reads, sequences


def count(mirna):
    """({miRNA: count}, {isomiR key: count}, annotated reads, substitution-free reads) of a `.mirna` file."""
    mature, isomir, annotated, perfect = {}, {}, {}, {}
    with open(mirna) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        missing = [c for c in COLUMNS if c not in header]
        if missing:
            fail('{}: not a miraligner -freq output; missing columns {}'.format(mirna, ', '.join(missing)))
        col = {c: header.index(c) for c in COLUMNS}
        for i, line in enumerate(fh, 2):
            fields = line.rstrip('\n').split('\t')
            if len(fields) < len(header):
                fail('{}:{}: {} fields, expected {}'.format(mirna, i, len(fields), len(header)))
            if fields[col['DB']] != 'miRNA':
                continue
            seq, mir, mism = fields[col['seq']], fields[col['mir']], fields[col['mism']]
            try:
                freq, ambiguity = float(fields[col['freq']]), int(fields[col['ambiguity']])
            except ValueError:
                fail('{}:{}: freq {!r} or ambiguity {!r} is not a number'.format(mirna, i, fields[col['freq']], fields[col['ambiguity']]))
            if ambiguity < 1:
                fail('{}:{}: ambiguity {}'.format(mirna, i, ambiguity))
            share = freq / ambiguity
            key = '{}_{}_{}'.format(mir, seq, ';'.join(fields[col[c]] for c in ('mism', 'add', 't5', 't3')))
            isomir.setdefault(key, []).append(share)
            annotated[seq] = freq
            if mism == '0':
                mature.setdefault(mir, []).append(share)
                perfect[seq] = freq
    return ({k: math.fsum(v) for k, v in mature.items()}, {k: math.fsum(v) for k, v in isomir.items()},
            math.fsum(annotated.values()), math.fsum(perfect.values()))


def write_table(path, header, table):
    with open(path, 'w') as fh:
        fh.write(header + '\tcount\n')
        fh.writelines('{}\t{}\n'.format(key, number(table[key])) for key in sorted(table))


def run_count(args):
    reads, sequences = read_input(args.reads)
    mature, isomir, annotated, perfect = count(args.mirna)
    write_table(args.mature, 'miRNA', mature)
    write_table(args.isomir, 'isomiR', isomir)
    with open(args.stats, 'w') as fh:
        fh.write('statistic\tvalue\n')
        for name, value in (('reads', reads), ('sequences', sequences), ('annotated_reads', annotated),
                            ('substitution_free_reads', perfect), ('mature_counts', math.fsum(mature.values()))):
            fh.write('{}\t{}\n'.format(name, number(value)))


def run_join(args):
    if len(args.samples) != len(args.tables) or len(set(args.samples)) != len(args.samples):
        fail('give one table per sample and unique sample names ({} samples, {} tables)'.format(len(args.samples), len(args.tables)))
    counts = {}
    for sample, path in zip(args.samples, args.tables):
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            if header != ['miRNA', 'count']:
                fail('{}: expected the header miRNA, count, got {}'.format(path, header))
            for i, line in enumerate(fh, 2):
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 2:
                    fail('{}:{}: expected a miRNA and a count, got {!r}'.format(path, i, line[:100]))
                name, value = fields
                if sample in counts.get(name, {}):
                    fail('{}:{}: {} twice'.format(path, i, name))
                counts.setdefault(name, {})[sample] = value
    out = sys.stdout
    out.write('\t'.join(['miRNA'] + args.samples) + '\n')
    for name in sorted(counts):
        out.write('\t'.join([name] + [counts[name].get(s, '0') for s in args.samples]) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    one = commands.add_parser('count')
    one.add_argument('--reads', required=True)
    one.add_argument('--mirna', required=True)
    one.add_argument('--mature', required=True)
    one.add_argument('--isomir', required=True)
    one.add_argument('--stats', required=True)
    many = commands.add_parser('join')
    many.add_argument('--samples', nargs='+', required=True)
    many.add_argument('--tables', nargs='+', required=True)
    args = parser.parse_args(argv)
    run_count(args) if args.command == 'count' else run_join(args)


if __name__ == '__main__':
    main()
