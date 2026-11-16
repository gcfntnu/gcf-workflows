#!/usr/bin/env python
"""miRge3 miRNA counts from mapped.csv, with a project-level exact-read rule instead of miRge3's per-sample one (#349).

miRge3 0.1.4 writes miR.Counts.csv in its summary step (mirge/libs/summary.py, `mirge_can`): in each sample it sets a
miRNA's count to 0 when the miRNA has fewer than 2 exact-match reads there, isomiR reads included. The reads stay
assigned in its mapped.csv (columns `exact miRNA` and `isomiR miRNA`, one read sequence per row, read counts per
sample).

This script counts, per sample, each miRNA's exact and isomiR reads, with names folded by the library's merge file
(`<org>_merges_<db>.csv`: group name, then its members; miRNAs with identical sequences and a miRNA's SNP variants) as
miRge3 does, and keeps a miRNA when it has at least MIN_PROJECT_EXACT exact reads in the project (all samples
together); other miRNAs are 0. Requiring exact reads keeps miRge3's protection against isomiR-only assignments, which
on synthetic pools are mostly reads of other sequences; counting it over the project instead of per sample keeps
miRNAs that are present but sparse in a sample. Rows and sample columns are those of miRge3's own table (every miRNA
of the library, unexpressed ones 0), so the output keeps its format, with whole numbers. It stops when a read has both
an exact and an isomiR name, when a miRNA name is not a row of miRge3's table, or when the samples differ.

Usage: mirge3_counts.py --mapped mapped.csv --merges human_merges_miRBase.csv --mirge3-table miR.Counts.mirge3.csv
"""
import argparse
import csv
import sys

MIN_PROJECT_EXACT = 2
# mapped.csv's annotation columns (miRge3 0.1.4 digest.py; `spike-in` only with --spikeIn); the other columns are samples
ANNOTATION_COLUMNS = {'Sequence', 'annotFlag', 'exact miRNA', 'hairpin miRNA', 'mature tRNA', 'primary tRNA', 'snoRNA',
                      'rRNA', 'ncrna others', 'mRNA', 'isomiR miRNA', 'spike-in'}


def read_merges(path):
    """Member name -> group name."""
    groups = {}
    with open(path) as fh:
        for line in fh:
            fields = line.strip().split(',')
            for member in fields[1:]:
                groups[member] = fields[0]
    return groups


def read_mirge3_table(path):
    """miRge3's miR.Counts.csv -> (miRNA names in order, sample names)."""
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh))
    if not rows or rows[0][0] != 'miRNA':
        raise ValueError('{}: not a miRge3 miR.Counts.csv (first column miRNA)'.format(path))
    return [row[0] for row in rows[1:] if row], rows[0][1:]


def count(mapped, merges, names, samples):
    """Per miRNA (miRge3 table names) and sample: exact plus isomiR reads from mapped.csv rows, for miRNAs with at least
    MIN_PROJECT_EXACT exact reads in all samples together; other miRNAs 0."""
    counts = {name: dict.fromkeys(samples, 0) for name in names}
    exact_reads = dict.fromkeys(names, 0)
    reader = csv.DictReader(mapped)
    missing = [s for s in samples if s not in reader.fieldnames]
    extra = [c for c in reader.fieldnames if c not in ANNOTATION_COLUMNS and c not in samples]
    if missing or extra:
        raise ValueError('mapped.csv and miRge3 table samples differ: missing {}, extra {}'.format(missing, extra))
    for row in reader:
        exact, isomir = row['exact miRNA'], row['isomiR miRNA']
        if exact and isomir:
            raise ValueError('read {} has both an exact ({}) and an isomiR ({}) miRNA'.format(row['Sequence'], exact, isomir))
        if not (exact or isomir):
            continue
        name = merges.get(exact or isomir, exact or isomir)
        if name not in counts:
            raise ValueError('miRNA {!r} (read {}) is not a row of the miRge3 table'.format(name, row['Sequence']))
        for sample in samples:
            counts[name][sample] += int(row[sample])
            if exact:
                exact_reads[name] += int(row[sample])
    for name in names:
        if exact_reads[name] < MIN_PROJECT_EXACT:
            counts[name] = dict.fromkeys(samples, 0)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--mapped', required=True)
    parser.add_argument('--merges', required=True)
    parser.add_argument('--mirge3-table', required=True)
    args = parser.parse_args()
    try:
        names, samples = read_mirge3_table(args.mirge3_table)
        with open(args.mapped, newline='') as fh:
            counts = count(fh, read_merges(args.merges), names, samples)
    except (OSError, ValueError, KeyError) as error:
        sys.exit('mirge3_counts.py: {}'.format(error))
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['miRNA'] + samples)
    for name in names:
        writer.writerow([name] + [counts[name][s] for s in samples])


if __name__ == '__main__':
    main()
