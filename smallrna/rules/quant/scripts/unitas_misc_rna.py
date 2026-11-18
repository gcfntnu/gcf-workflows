#!/usr/bin/env python
"""UNITAS misc_RNA reads by gene family (#379), from the target names in each sample's unitas.hits_per_target.txt.

    unitas_misc_rna.py --summaries S/unitas.annotation_summary.txt [...] --hits-per-target S/unitas.hits_per_target.txt [...]
        --table misc_rna_families.tsv --mqc unitas_misc_rna_mqc.yaml

--table: Sample_ID, family, count; every family for every sample, in the order of FAMILIES, samples sorted.
--mqc: MultiQC custom content (YAML in JSON notation), a bar per sample with its misc_RNA reads by family, or a short
note when no sample has misc_RNA reads.

UNITAS reports Ensembl's misc_RNA biotype as one class; in blood-derived libraries it is mostly Y RNA fragments.
A misc_RNA target belongs to the first family with a name prefix that starts its name (the reference's Ensembl gene
name, e.g. `RNY4-ENSG00000252316.1`), ignoring case, so mouse and rat symbols (`Rny1`, `Rn7sk`, `Rn7s1`) match too;
`RN7S` follows `RN7SK`, so it takes the 7SL genes (`RN7SL1`, mouse `Rn7s1`) but not 7SK.
Names that match none are `other`, including transcript-only names without a gene symbol (`ENST...`).
A read UNITAS splits between targets keeps its fractional counts, so a sample's families must add up to the misc_RNA
class of its annotation summary, read as unitas_annotations.py reads it (up to the first empty line). A difference,
an unpaired file, a hits table with another header, a malformed row or a count that is not a number stops the script.
An empty or header-only hits table (a sample without reads) has no targets; blank lines in it are skipped.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

# (family, label, name prefixes); the last family takes every other name.
FAMILIES = (('Y_RNA', 'Y RNA', ('RNY', 'Y_RNA')), ('7SK', '7SK', ('RN7SK', '7SK')),
            ('7SL', '7SL (SRP RNA)', ('RN7S', 'Metazoa_SRP')), ('other', 'other misc_RNA', ()))
NAMES = tuple(name for name, _, _ in FAMILIES)
PREFIXES = tuple((name, tuple(p.upper() for p in prefixes)) for name, _, prefixes in FAMILIES if prefixes)
RULES = '; '.join('{} as {}'.format(' and '.join(p + '...' for p in prefixes), label) for _, label, prefixes in FAMILIES if prefixes) + \
    '; other names, including transcript-only names without a gene symbol, as other misc_RNA'
DESCRIPTION = ('Reads of the UNITAS misc_RNA class by gene family, from the names of the Ensembl genes they hit, the first '
               'match counting: {}. Reads split between genes keep their fractional counts.').format(RULES)
HITS_HEADER = ['TRANSCRIPT_CLASS', 'TRANSCRIPT_NAME', 'READ_COUNT']


def family(target):
    """The misc_RNA family of a UNITAS target name."""
    upper = target.upper()
    return next((name for name, prefixes in PREFIXES if upper.startswith(prefixes)), NAMES[-1])


def number(value, path, line):
    try:
        return float(value)
    except ValueError:
        raise ValueError('{}: line {}: expected a number, got {!r}'.format(path, line, value)) from None


def misc_rna_families(path):
    """misc_RNA reads of one hits_per_target table by family, every family present."""
    counts = dict.fromkeys(NAMES, 0.0)
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
        rows = ((reader.line_num, row) for row in reader if row)
        header = next(rows, None)
        if header is None:
            return counts
        if header[1] != HITS_HEADER:
            raise ValueError('{}: line {}: not a UNITAS hits_per_target header (expected {})'.format(path, header[0], ' '.join(HITS_HEADER)))
        for line, row in rows:
            if len(row) != 3:
                raise ValueError('{}: line {}: expected 3 columns, got {!r}'.format(path, line, row))
            if row[0] == 'misc_RNA':
                counts[family(row[1])] += number(row[2], path, line)
    return counts


def summary_misc_rna(path):
    """The misc_RNA class count of a UNITAS annotation summary, 0 when it has none. Read as unitas_annotations.py
    (file2pandas) reads it: up to the first empty line, a row a name and a count, three spaces per sub-level."""
    for line, text in enumerate(Path(path).read_text().split('\n'), start=1):
        if not text:
            break
        fields = text.split('\t')
        if len(fields) != 2:
            raise ValueError('{}: line {}: expected a class and a count, got {!r}'.format(path, line, text))
        if fields[0] == 'misc_RNA':
            return number(fields[1], path, line)
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


def format_table(result):
    lines = ['Sample_ID\tfamily\tcount']
    lines += ['{}\t{}\t{:.10g}'.format(sample, name, counts[name]) for sample, counts in result.items() for name in NAMES]
    return '\n'.join(lines) + '\n'


def mqc(result):
    section = {'id': 'unitas_misc_rna', 'section_name': 'UNITAS: misc_RNA families', 'description': DESCRIPTION}
    if not any(any(counts.values()) for counts in result.values()):
        return {**section, 'plot_type': 'html', 'data': '<p>No sample has misc_RNA reads.</p>'}
    return {**section, 'plot_type': 'bargraph',
            'pconfig': {'id': 'unitas_misc_rna_plot', 'title': 'UNITAS: misc_RNA families', 'ylab': '# Reads'},
            'categories': {name: {'name': label} for name, label, _ in FAMILIES},
            'data': {sample: {name: round(value, 2) for name, value in counts.items()} for sample, counts in result.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--summaries', nargs='+', required=True, help='unitas.annotation_summary.txt, one per sample directory')
    parser.add_argument('--hits-per-target', nargs='+', required=True, dest='hits', help='unitas.hits_per_target.txt, one per sample directory')
    parser.add_argument('--table', required=True, help='output: misc_RNA reads by family per sample (TSV)')
    parser.add_argument('--mqc', required=True, help='output: MultiQC custom content section (YAML)')
    args = parser.parse_args()
    try:
        result = table(args.summaries, args.hits)
    except (OSError, ValueError) as error:
        sys.exit('unitas_misc_rna.py: {}'.format(error))
    Path(args.table).write_text(format_table(result))
    Path(args.mqc).write_text(json.dumps(mqc(result), indent=1) + '\n')


if __name__ == '__main__':
    main()
