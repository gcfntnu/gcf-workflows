#!/usr/bin/env python
"""UNITAS misc_RNA reads by gene family (#379), from the target names in each sample's unitas.hits_per_target.txt.

    unitas_misc_rna.py --annotations annotations.tsv --hits-per-target S/unitas.hits_per_target.txt [...]
        --table misc_rna_families.tsv --mqc unitas_misc_rna_mqc.yaml

--table: Sample_ID, family, count; every family for every sample, in the order of FAMILIES, samples sorted.
--mqc: MultiQC custom content (YAML in JSON notation), a bar per sample with its misc_RNA reads by family, or a short
note when no sample has misc_RNA reads.

UNITAS reports Ensembl's misc_RNA biotype as one class; in blood-derived libraries it is mostly Y RNA fragments.
A misc_RNA target belongs to the first family with a name prefix that starts its name (the reference's Ensembl gene
name, e.g. `RNY4-ENSG00000252316.1`), ignoring case, so mouse and rat symbols (`Rny1`, `Rn7sk`, `Rn7s1`) match too;
`RN7S` follows `RN7SK`, so it takes the 7SL genes (`RN7SL1`, mouse `Rn7s1`) but not 7SK.
Names that match none are `other`, including transcript-only names without a gene symbol (`ENST...`).
A read UNITAS splits between targets keeps its fractional counts, so a sample's families must add up to its misc_RNA
class in annotations.tsv (unitas_annotations.py; Level1 misc_RNA without Level2, 0 when the sample has none).
A difference, a sample missing from either input, a hits table with another header, a data row without class or
target name, or a count that is not a finite number stops the script, naming the file and line.
An empty or header-only hits table (a sample without reads) has no targets; blank lines in it are skipped.
"""
import argparse
import csv
import json
import math
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
ANNOTATIONS_HEADER = ['Level1', 'Level2', 'Level3', 'Count', 'Sample_ID']


def family(target):
    """The misc_RNA family of a UNITAS target name."""
    upper = target.upper()
    return next((name for name, prefixes in PREFIXES if upper.startswith(prefixes)), NAMES[-1])


def number(value, path, line):
    """A finite count, or an error naming the file and line."""
    try:
        count = float(value)
    except ValueError:
        count = math.nan
    if not math.isfinite(count):
        raise ValueError('{}: line {}: expected a finite number, got {!r}'.format(path, line, value))
    return count


def rows(path):
    """(line number, fields) of a TSV without CSV quoting, skipping blank lines (empty or only spaces; a line with a tab
    has fields and is a row)."""
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            if len(row) > 1 or (row and row[0].strip()):
                yield reader.line_num, row


def misc_rna_families(path):
    """misc_RNA reads of one hits_per_target table by family, every family present."""
    counts = dict.fromkeys(NAMES, 0.0)
    lines = rows(path)
    header = next(lines, None)
    if header is None:
        return counts
    if header[1] != HITS_HEADER:
        raise ValueError('{}: line {}: not a UNITAS hits_per_target header (expected {})'.format(path, header[0], ' '.join(HITS_HEADER)))
    for line, row in lines:
        if len(row) != 3 or not row[0].strip() or not row[1].strip():
            raise ValueError('{}: line {}: expected a class, a target name and a count, got {!r}'.format(path, line, row))
        count = number(row[2], path, line)
        if row[0] == 'misc_RNA':
            counts[family(row[1])] += count
    return counts


def annotation_misc_rna(path):
    """{sample: misc_RNA class count} from annotations.tsv, every sample in it (0 without a misc_RNA row)."""
    lines = rows(path)
    header = next(lines, (1, None))
    if header[1] != ANNOTATIONS_HEADER:
        raise ValueError('{}: line {}: not a UNITAS annotations table (expected columns {})'.format(path, header[0], ' '.join(ANNOTATIONS_HEADER)))
    totals = {}
    for line, row in lines:
        if len(row) != 5:
            raise ValueError('{}: line {}: expected 5 columns, got {!r}'.format(path, line, row))
        level1, level2, level3, count, sample = row
        totals.setdefault(sample, 0.0)
        if level1 == 'misc_RNA' and not level2 and not level3:
            totals[sample] = number(count, path, line)
    return totals


def check_total(counts, misc_rna, path):
    total = sum(counts.values())
    if not abs(total - misc_rna) <= 0.01 + 1e-6 * misc_rna:
        raise ValueError('{}: misc_RNA targets add up to {:.6g} reads, annotations.tsv has {:.6g}'.format(path, total, misc_rna))


def table(annotations, hits):
    """{sample: {family: count}}, samples named by the directories that hold their hits table, as in annotations.tsv."""
    target_hits = {Path(p).parent.name: p for p in hits}
    if len(target_hits) != len(hits):
        raise ValueError('need one hits table per sample directory, got {}'.format(hits))
    totals = annotation_misc_rna(annotations)
    missing, extra = sorted(set(target_hits) - set(totals)), sorted(set(totals) - set(target_hits))
    if missing or extra:
        raise ValueError('{}: samples {} have hits tables but are not in it; samples {} are in it without hits tables'.format(annotations, missing, extra))
    result = {}
    for sample in sorted(target_hits):
        counts = misc_rna_families(target_hits[sample])
        check_total(counts, totals[sample], target_hits[sample])
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
    parser.add_argument('--annotations', required=True, help='annotations.tsv from unitas_annotations.py')
    parser.add_argument('--hits-per-target', nargs='+', required=True, dest='hits', help='unitas.hits_per_target.txt, one per sample directory')
    parser.add_argument('--table', required=True, help='output: misc_RNA reads by family per sample (TSV)')
    parser.add_argument('--mqc', required=True, help='output: MultiQC custom content section (YAML)')
    args = parser.parse_args()
    try:
        result = table(args.annotations, args.hits)
    except (OSError, ValueError) as error:
        sys.exit('unitas_misc_rna.py: {}'.format(error))
    Path(args.table).write_text(format_table(result))
    Path(args.mqc).write_text(json.dumps(mqc(result), indent=1) + '\n')


if __name__ == '__main__':
    main()
