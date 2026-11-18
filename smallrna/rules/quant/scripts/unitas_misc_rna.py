#!/usr/bin/env python
"""UNITAS misc_RNA reads by gene family (#379), from the target names in each sample's unitas.hits_per_target.txt.

    unitas_misc_rna.py --summaries S/unitas.annotation_summary.txt [...] --hits-per-target S/unitas.hits_per_target.txt [...]
        > misc_rna_families.tsv

UNITAS reports Ensembl's misc_RNA biotype as one class; in blood-derived libraries it is mostly Y RNA fragments.
A misc_RNA target belongs to the first family with a name prefix that starts its name (the reference's Ensembl gene
name, e.g. `RNY4-ENSG00000252316.1`), ignoring case, so mouse and rat symbols (`Rny1`, `Rn7sk`, `Rn7sl1`) match too.
Names that match none are `other`, including transcript-only names without a gene symbol (`ENST...`).
A read UNITAS splits between targets keeps its fractional counts, so a sample's families must add up to the misc_RNA
class of its annotation summary; a difference, an unpaired file or a table with another header stops the script.
An empty or header-only hits table (a sample without reads) has no targets.

Writes Sample_ID, family, count: every family for every sample, in the order of FAMILIES, samples sorted.
unitas_misc_rna_mqc.py imports FAMILIES and RULES for the MultiQC section.
"""
import argparse
import csv
from pathlib import Path
import sys

# (family, label, name prefixes); the last family takes every other name.
FAMILIES = (('Y_RNA', 'Y RNA', ('RNY', 'Y_RNA')), ('7SK', '7SK', ('RN7SK', '7SK')),
            ('7SL', '7SL (SRP RNA)', ('RN7SL', 'Metazoa_SRP')), ('other', 'other misc_RNA', ()))
NAMES = tuple(name for name, _, _ in FAMILIES)
RULES = '; '.join('{} as {}'.format(' and '.join(p + '...' for p in prefixes), label) for _, label, prefixes in FAMILIES if prefixes) + \
    '; other names, including transcript-only names without a gene symbol, as other misc_RNA'
HEADER = ['TRANSCRIPT_CLASS', 'TRANSCRIPT_NAME', 'READ_COUNT']


def family(target):
    """The misc_RNA family of a UNITAS target name."""
    upper = target.upper()
    return next(name for name, _, prefixes in FAMILIES if not prefixes or upper.startswith(tuple(p.upper() for p in prefixes)))


def misc_rna_families(path):
    """misc_RNA reads of one hits_per_target table by family, every family present."""
    counts = dict.fromkeys(NAMES, 0.0)
    with open(path, newline='') as fh:
        rows = csv.reader(fh, delimiter='\t')
        header = next(rows, None)
        if header is None:
            return counts
        if header != HEADER:
            raise ValueError('{}: not a UNITAS hits_per_target table (expected columns {})'.format(path, ' '.join(HEADER)))
        for row in rows:
            if len(row) != 3:
                raise ValueError('{}: expected 3 columns, got {!r}'.format(path, row))
            if row[0] == 'misc_RNA':
                counts[family(row[1])] += float(row[2])
    return counts


def summary_misc_rna(path):
    """The misc_RNA class count of a UNITAS annotation summary, 0 when it has none."""
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line[0].isspace():
            fields = line.rsplit('\t', 1)
            if len(fields) != 2:
                raise ValueError('{}: expected a class and a count, got {!r}'.format(path, line))
            name, count = fields
            if name == 'misc_RNA':
                return float(count)
    return 0.0


def check_total(counts, misc_rna, path):
    total = sum(counts.values())
    if abs(total - misc_rna) > 0.01 + 1e-6 * misc_rna:
        raise ValueError('{}: misc_RNA targets add up to {:.6g} reads, the annotation summary has {:.6g}'.format(path, total, misc_rna))


def table(summaries, hits):
    """{sample: {family: count}}, samples named by the directories that hold their summary and hits table."""
    by_sample = lambda paths: {Path(p).parent.name: p for p in paths}
    summary, target_hits = by_sample(summaries), by_sample(hits)
    if len(summary) != len(summaries) or len(target_hits) != len(hits) or sorted(summary) != sorted(target_hits):
        raise ValueError('need one summary and one hits table per sample directory, got {} and {}'.format(sorted(summary), sorted(target_hits)))
    result = {}
    for sample in sorted(summary):
        counts = misc_rna_families(target_hits[sample])
        check_total(counts, summary_misc_rna(summary[sample]), target_hits[sample])
        result[sample] = counts
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--summaries', nargs='+', required=True, help='unitas.annotation_summary.txt, one per sample directory')
    parser.add_argument('--hits-per-target', nargs='+', required=True, dest='hits', help='unitas.hits_per_target.txt, one per sample directory')
    args = parser.parse_args()
    try:
        result = table(args.summaries, args.hits)
    except (OSError, ValueError) as error:
        sys.exit('unitas_misc_rna.py: {}'.format(error))
    out = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    out.writerow(['Sample_ID', 'family', 'count'])
    for sample, counts in result.items():
        out.writerows([sample, name, '{:.10g}'.format(counts[name])] for name in NAMES)


if __name__ == '__main__':
    main()
