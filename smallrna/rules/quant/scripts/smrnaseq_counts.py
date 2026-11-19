#!/usr/bin/env python
"""Project tables from per-sample `mirtop counts --add-extra` output, as nf-core/smrnaseq 2.4.1 makes them (#403).

Reproduces PIVOT_LONGER, PIVOT_WIDER and DATATABLE_MERGE (bin/pivot_longer.R, bin/pivot_wider.R,
bin/collapse_mirtop.r):

- the samples' tables are joined on their 12 annotation columns (UID, Read, miRNA, Variant, iso_5p, iso_3p,
  iso_add3p, iso_snp, iso_5p_nt, iso_3p_nt, iso_add3p_nt, iso_snp_nt), one count column per sample, 0 where a sample
  lacks the row; rows are sorted by UID, ties by the whole line (byte order: `sort -t$'\\t' -k1,1` in the C locale).
  Written as the per-isomiR table (`--joined`);
- of rows with the same UID (one collapsed read annotated to several miRNAs or variants) only the first is kept, so
  such a read counts in full for the miRNA whose line sorts first; the counts are then summed per miRNA (`--mirna`).

Counts are written as integers. The sample columns follow --tables, whose names come from each table's last column.
Each table is written to `<table>.partial` beside it and renamed into place when complete, so an interrupted run never
leaves a partial table at the final path (the benchmark reads mirna.tsv outside Snakemake). The staging name is fixed,
so a rerun after a hard kill overwrites the leftover instead of adding another copy; Snakemake's directory lock keeps
to one writer per table.

Usage: smrnaseq_counts.py --tables S1/mirtop.tsv S2/mirtop.tsv --joined joined_samples_mirtop.tsv --mirna mirna.tsv
"""
import argparse
import csv
import os
import sys

ANNOTATION = ('UID', 'Read', 'miRNA', 'Variant', 'iso_5p', 'iso_3p', 'iso_add3p', 'iso_snp',
              'iso_5p_nt', 'iso_3p_nt', 'iso_add3p_nt', 'iso_snp_nt')


def number(text, path):
    try:
        value = float(text)
    except ValueError:
        raise ValueError('{}: count {!r} is not a number'.format(path, text))
    if value < 0 or value != int(value):
        raise ValueError('{}: count {!r} is not a whole number of reads'.format(path, text))
    return int(value)


def read_table(path):
    """(sample, {annotation tuple: count}) from one sample's mirtop counts table."""
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        if header is None or tuple(header[:len(ANNOTATION)]) != ANNOTATION or len(header) != len(ANNOTATION) + 1:
            raise ValueError('{}: expected the columns {} and one sample, got {}'.format(path, ', '.join(ANNOTATION), header))
        sample, rows = header[-1], {}
        for i, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError('{}:{}: {} fields, expected {}'.format(path, i, len(row), len(header)))
            key = tuple(row[:len(ANNOTATION)])
            if key in rows:
                raise ValueError('{}:{}: row repeated'.format(path, i))
            rows[key] = number(row[-1], path)
    return sample, rows


def join(tables):
    """Header and rows of the joined table, sorted as nf-core sorts it."""
    samples, wide = [], {}
    for sample, rows in tables:
        if sample in samples:
            raise ValueError('sample {} in two tables'.format(sample))
        samples.append(sample)
        for key, count in rows.items():
            wide.setdefault(key, {})[sample] = count
    lines = [list(key) + [str(counts.get(s, 0)) for s in samples] for key, counts in wide.items()]
    lines.sort(key=lambda line: (line[0].encode(), '\t'.join(line).encode()))
    return list(ANNOTATION) + samples, lines


def collapse(header, lines):
    """Counts per miRNA, the first row of each UID only, miRNAs in order of first appearance."""
    samples = header[len(ANNOTATION):]
    seen, totals = set(), {}
    for line in lines:
        if line[0] in seen:
            continue
        seen.add(line[0])
        counts = totals.setdefault(line[2], [0] * len(samples))
        for i, c in enumerate(line[len(ANNOTATION):]):
            counts[i] += int(c)
    return samples, totals


def write_table(path, rows):
    """Rows to path, through `<path>.partial` renamed into place once written."""
    tmp = '{}.partial'.format(path)
    try:
        with open(tmp, 'w', newline='') as fh:
            csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--tables', nargs='+', required=True)
    parser.add_argument('--joined', required=True)
    parser.add_argument('--mirna', required=True)
    args = parser.parse_args(argv)
    try:
        header, lines = join(read_table(p) for p in args.tables)
    except ValueError as e:
        sys.exit('smrnaseq_counts.py: {}'.format(e))
    write_table(args.joined, [header] + lines)
    samples, totals = collapse(header, lines)
    write_table(args.mirna, [['miRNA'] + samples] + [[mirna] + counts for mirna, counts in totals.items()])


if __name__ == '__main__':
    main()
