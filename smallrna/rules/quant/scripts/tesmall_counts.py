#!/usr/bin/env python
"""Join TEsmall's per-sample count_summary.txt tables into one table: ftype, fid, then one count column per sample.

TEsmall writes `fid`, `ftype` and one column named by its sample label. A feature missing from a sample is 0 there.
A table whose count column is not the expected sample, or a feature listed twice in one table, stops with an error.

Usage: tesmall_counts.py --samples S1 S2 ... --tables S1/count_summary.txt S2/count_summary.txt ... > tesmall_counts.tsv
"""
import argparse
import csv
import sys


def read_table(path, sample):
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader, None)
        if header is None or header[:2] != ['fid', 'ftype'] or header[2:] != [sample]:
            sys.exit('tesmall_counts.py: {}: expected columns fid, ftype, {}; got {}'.format(path, sample, header))
        counts = {}
        for row in reader:
            if not row:
                continue
            key = (row[1], row[0])
            if key in counts:
                sys.exit('tesmall_counts.py: {}: {} {} is listed twice'.format(path, *key))
            counts[key] = row[2]
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--samples', nargs='+', required=True)
    parser.add_argument('--tables', nargs='+', required=True)
    args = parser.parse_args()
    if len(args.samples) != len(args.tables):
        sys.exit('tesmall_counts.py: {} samples but {} tables'.format(len(args.samples), len(args.tables)))
    tables = [read_table(path, sample) for sample, path in zip(args.samples, args.tables)]
    keys = sorted(set().union(*tables))
    out = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    out.writerow(['ftype', 'fid'] + args.samples)
    for key in keys:
        out.writerow(list(key) + [table.get(key, '0') for table in tables])


if __name__ == '__main__':
    main()
