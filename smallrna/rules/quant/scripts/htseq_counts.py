#!/usr/bin/env python
"""Counts of the htseq route (smallrna quant.method: htseq, #425), the genome counter of an earlier in-house small-RNA
pipeline, from htseq-count's annotated alignments.

count: one sample, counted as the pipeline counted:
- an alignment record counts when its htseq-count feature (the XF tag, the record's last field) contains the species
  code (e.g. `hsa`) and it passes the pipeline's perfect-alignment test: its MD tag starts with a run of 1 to 99
  matching bases followed by neither a digit nor a base (`MD:Z:[1-9][0-9]?[^0-9ACGTU]`). That is an alignment without
  mismatches, but the test reads only that first run: an insertion leaves no trace in MD, and a deletion after the
  first run (`MD:Z:21^A1`) passes too, as in the pipeline;
- ambiguous feature sets (`__ambiguous[a+b]`) are expanded, and the names that start with the species code are kept;
- each (name, read) pair counts once, with the read's count from its fastx_collapser name (`<rank>-<count>`): a read
  counts in full for every name it overlaps at any of its alignments, never divided;
- isomiRs are keyed `<name>_<sequence>`, the read's sequence as sequenced (minus-strand alignments reverse-complemented
  back).
Writes the mature (miRNA, count) and isomiR (isomiR, count) tables, sorted, and the sample's statistics: reads and
unique sequences in the collapsed reads, reads with at least one alignment, reads with a counted alignment, and the
summed mature counts (more than the counted reads when reads count for several names).

join: the samples' mature tables as one table, miRNA and one column per sample (0 where a sample lacks the miRNA).

Usage:
  htseq_counts.py count --species CODE --collapsed FASTA --annotated SAM[.gz] --mature TSV --isomir TSV --stats TSV
  htseq_counts.py join --samples S1 [S2 ...] --tables MATURE_TSV [MATURE_TSV ...]
"""
import argparse
import gzip
import re
import sys

PERFECT = re.compile(r'MD:Z:[1-9][0-9]?[^0-9ACGTU]')
AMBIGUOUS = re.compile(r'(.*?)__ambiguous\[(.*?)\](.*)')
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def fail(message):
    sys.exit('htseq_counts.py: ' + message)


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def read_count(name):
    """The read count in a fastx_collapser name, `<rank>-<count>`."""
    rank, _, count = name.partition('-')
    if not rank.isdigit() or not count.isdigit():
        fail('read name {!r} is not <rank>-<count> (fastx_collapser)'.format(name))
    return int(count)


def expand(feature):
    """The names of a feature: those of an ambiguous set (`__ambiguous[a+b]`), or the feature itself."""
    match = AMBIGUOUS.match(feature)
    if not match:
        return [feature]
    pre, names, post = match.groups()
    return [pre + name + post for name in names.split('+')]


def collapsed_reads(path):
    """(reads, unique sequences) of a fastx_collapser FASTA."""
    reads = sequences = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                reads += read_count(line[1:].strip())
                sequences += 1
    return reads, sequences


def count(species, annotated):
    """({name: count}, {name_sequence: count}, aligned reads, counted reads) of htseq-count's annotated records."""
    mature, isomir, aligned, counted = set(), set(), {}, {}
    with open_text(annotated) as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            record = line.rstrip('\n')
            fields = record.split('\t')
            if len(fields) < 12 or not fields[-1].startswith('XF:Z:'):
                fail('{}: an alignment record without htseq-count\'s XF tag last: {!r}'.format(annotated, record[:200]))
            read, flag, sequence, feature = fields[0], int(fields[1]), fields[9], fields[-1][5:]
            if not flag & 4:
                aligned[read] = read_count(read)
            if species not in feature or not PERFECT.search(record):
                continue
            if flag & 16:
                sequence = sequence[::-1].translate(COMPLEMENT)
            names = [n for n in expand(feature) if n.startswith(species)]
            if names:
                counted[read] = read_count(read)
            for name in names:
                mature.add((name, read))
                isomir.add(('{}_{}'.format(name, sequence), read))
    return totals(mature), totals(isomir), sum(aligned.values()), sum(counted.values())


def totals(pairs):
    out = {}
    for key, read in pairs:
        out[key] = out.get(key, 0) + read_count(read)
    return out


def write_table(path, header, table):
    with open(path, 'w') as fh:
        fh.write(header + '\tcount\n')
        fh.writelines('{}\t{}\n'.format(key, table[key]) for key in sorted(table))


def run_count(args):
    if not args.species:
        fail('--species is empty: the organism has no miRBase code')
    reads, sequences = collapsed_reads(args.collapsed)
    mature, isomir, aligned, counted = count(args.species, args.annotated)
    write_table(args.mature, 'miRNA', mature)
    write_table(args.isomir, 'isomiR', isomir)
    with open(args.stats, 'w') as fh:
        fh.write('statistic\tvalue\n')
        for name, value in (('reads', reads), ('sequences', sequences), ('aligned_reads', aligned),
                            ('counted_reads', counted), ('mature_counts', sum(mature.values()))):
            fh.write('{}\t{}\n'.format(name, value))


def run_join(args):
    if len(args.samples) != len(args.tables) or len(set(args.samples)) != len(args.samples):
        fail('give one table per sample and unique sample names ({} samples, {} tables)'.format(len(args.samples), len(args.tables)))
    counts = {}
    for sample, path in zip(args.samples, args.tables):
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            if header != ['miRNA', 'count']:
                fail('{}: expected the header miRNA, count, got {}'.format(path, header))
            for line in fh:
                name, value = line.rstrip('\n').split('\t')
                counts.setdefault(name, {})[sample] = value
    out = sys.stdout
    out.write('\t'.join(['miRNA'] + args.samples) + '\n')
    for name in sorted(counts):
        out.write('\t'.join([name] + [counts[name].get(s, '0') for s in args.samples]) + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    one = commands.add_parser('count')
    one.add_argument('--species', required=True)
    one.add_argument('--collapsed', required=True)
    one.add_argument('--annotated', required=True)
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
