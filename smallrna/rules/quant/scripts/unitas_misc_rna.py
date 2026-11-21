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
Y RNA that Ensembl, and so the UNITAS dump, files under another class counts as misc_RNA, family Y RNA (#423,
y_rna_elsewhere): human RNY5, a lncRNA transcript without a gene symbol; a second mouse Rny1 gene, a lncRNA; and rat's
Y RNA genes, which have a class of their own (`Y_RNA`). unitas_annotations.py moves the same reads to misc_RNA in
annotations.tsv, so the three species report Y RNA in one place and each sample's classes still add up.
A read UNITAS splits between targets keeps its fractional counts, so a sample's families must add up to its misc_RNA
class in annotations.tsv (unitas_annotations.py; Level1 misc_RNA without Level2, 0 when the sample has none, and
for a sample without rows there, as from an empty summary): both 0, or both positive and equal up to summation order.
A difference, a second misc_RNA row for a sample, a sample in annotations.tsv without hits table, a hits table with
another header, a data row without class or target name, or a count that is not a non-negative decimal number stops
the script, naming the file and line; a sample without rows in annotations.tsv is checked against 0 like a sample
without a misc_RNA row. Counts are written in their shortest exact form as plain decimals (`3.0`, `0.5`,
`0.00005`), never in scientific notation. The MultiQC section is a note when every plotted value (2 decimals) is 0.
An empty or header-only hits table (a sample without reads) has no targets; blank lines in it are skipped.
"""
import argparse
import csv
import decimal
import json
import math
from pathlib import Path
import re
import sys

# (family, label, name prefixes); the last family takes every other name.
FAMILIES = (('Y_RNA', 'Y RNA', ('RNY', 'Y_RNA')), ('7SK', '7SK', ('RN7SK', '7SK')),
            ('7SL', '7SL (SRP RNA)', ('RN7S', 'Metazoa_SRP')), ('other', 'other misc_RNA', ()))
NAMES = tuple(name for name, _, _ in FAMILIES)
PREFIXES = tuple((name, tuple(p.upper() for p in prefixes)) for name, _, prefixes in FAMILIES if prefixes)
RULES = '; '.join('{} as {}'.format(' and '.join(p + '...' for p in prefixes), label) for _, label, prefixes in FAMILIES if prefixes) + \
    '; other names, including transcript-only names without a gene symbol, as other misc_RNA'
DESCRIPTION = ('Reads of the UNITAS misc_RNA class by gene family, from the names of the Ensembl genes they hit, the first '
               'match counting: {}. Y RNA genes that Ensembl files under another class (human RNY5, a mouse Rny1 gene, the rat '
               'Y RNA genes) count here, as Y RNA; the UNITAS section above shows UNITAS\'s own classes, with those reads under '
               'lncRNA or Y_RNA. Reads split between genes keep their fractional counts.').format(RULES)
# Y RNA outside the misc_RNA class (#423): genes by Ensembl ID, as the dump's target names carry it (`<name>-<gene>.<version>`),
# and whole classes.
Y_RNA_GENES = ('ENSG00000286171', 'ENSMUSG00000138036')  # human RNY5; the mouse Rny1 gene with the biotype lncRNA
Y_RNA_CLASSES = ('Y_RNA',)  # rat: the Ensembl biotype of its Y RNA genes
HITS_HEADER = ['TRANSCRIPT_CLASS', 'TRANSCRIPT_NAME', 'READ_COUNT']
COUNT = re.compile(r'\d+(\.\d+)?([eE][+-]?\d+)?')
ANNOTATIONS_HEADER = ['Level1', 'Level2', 'Level3', 'Count', 'Sample_ID']


def family(target):
    """The misc_RNA family of a UNITAS target name."""
    upper = target.upper()
    return next((name for name, prefixes in PREFIXES if upper.startswith(prefixes)), NAMES[-1])


def y_rna_elsewhere(transcript_class, target):
    """Whether a target outside the misc_RNA class is Y RNA (Y_RNA_GENES, Y_RNA_CLASSES)."""
    if transcript_class == 'misc_RNA':
        return False
    return transcript_class in Y_RNA_CLASSES or any('-{}.'.format(gene) in target or target.endswith('-' + gene) for gene in Y_RNA_GENES)


def number(value, path, line):
    """A count: a non-negative, finite decimal number, or an error naming the file and line."""
    count = float(value) if COUNT.fullmatch(value) else math.nan
    if not math.isfinite(count):
        raise ValueError('{}: line {}: expected a non-negative decimal number, got {!r}'.format(path, line, value))
    return count


def rows(path):
    """(line number, fields) of a TSV without CSV quoting, skipping blank lines (empty or only spaces; a line with a tab
    has fields and is a row)."""
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
        while True:
            try:
                row = next(reader)
            except StopIteration:
                return
            except csv.Error as error:
                raise ValueError('{}: line {}: {}'.format(path, reader.line_num, error)) from None
            if len(row) > 1 or (row and row[0].strip()):
                yield reader.line_num, row


def hits(path):
    """(class, target name, count) of a hits_per_target table's rows; nothing for an empty file."""
    lines = rows(path)
    header = next(lines, None)
    if header is None:
        return
    if header[1] != HITS_HEADER:
        raise ValueError('{}: line {}: not a UNITAS hits_per_target header (expected {})'.format(path, header[0], ' '.join(HITS_HEADER)))
    for line, row in lines:
        if len(row) != 3 or not row[0].strip() or not row[1].strip():
            raise ValueError('{}: line {}: expected a class, a target name and a count, got {!r}'.format(path, line, row))
        yield row[0], row[1], number(row[2], path, line)


def misc_rna_families(path):
    """misc_RNA reads of one hits_per_target table by family, every family present, with the Y RNA of other classes
    (y_rna_elsewhere); summed with math.fsum."""
    counts = {name: [] for name in NAMES}
    for transcript_class, target, count in hits(path):
        if transcript_class == 'misc_RNA':
            counts[family(target)].append(count)
        elif y_rna_elsewhere(transcript_class, target):
            counts[NAMES[0]].append(count)
    return {name: math.fsum(values) for name, values in counts.items()}


def y_rna_by_class(path):
    """{class: reads} of one hits_per_target table's Y RNA outside the misc_RNA class: what unitas_annotations.py moves
    to misc_RNA."""
    moved = {}
    for transcript_class, target, count in hits(path):
        if y_rna_elsewhere(transcript_class, target):
            moved.setdefault(transcript_class, []).append(count)
    return {transcript_class: math.fsum(values) for transcript_class, values in moved.items()}


def annotation_misc_rna(path):
    """{sample: misc_RNA class count} from annotations.tsv, every sample in it (0 without a misc_RNA row, an error with two)."""
    lines = rows(path)
    header = next(lines, (1, None))
    if header[1] != ANNOTATIONS_HEADER:
        raise ValueError('{}: line {}: not a UNITAS annotations table (expected columns {})'.format(path, header[0], ' '.join(ANNOTATIONS_HEADER)))
    totals = {}
    for line, row in lines:
        if len(row) != 5:
            raise ValueError('{}: line {}: expected 5 columns, got {!r}'.format(path, line, row))
        level1, level2, level3, count, sample = row
        totals.setdefault(sample, None)
        if level1 == 'misc_RNA' and not level2 and not level3:
            if totals[sample] is not None:
                raise ValueError('{}: line {}: a second misc_RNA row for sample {}'.format(path, line, sample))
            totals[sample] = number(count, path, line)
    return {sample: 0.0 if total is None else total for sample, total in totals.items()}


# UNITAS writes the summary at full precision (15 significant digits), and a read split between many targets gives each
# a small share (1/6400 seen), so only summation order may differ: 4e-14 relative at most in development runs.
TOTAL_REL_TOL = 1e-9


def check_total(counts, misc_rna, path, source='annotations.tsv has'):
    """The families must add up to the sample's misc_RNA class: both 0, or both positive and equal within TOTAL_REL_TOL."""
    total = math.fsum(counts.values())
    if (total > 0) != (misc_rna > 0) or not math.isclose(total, misc_rna, rel_tol=TOTAL_REL_TOL):
        raise ValueError('{}: misc_RNA targets add up to {!r} reads, {} {!r}'.format(path, total, source, misc_rna))


def table(annotations, hits):
    """{sample: {family: count}}, samples named by the directories that hold their hits table, as in annotations.tsv."""
    target_hits = {Path(p).parent.name: p for p in hits}
    if len(target_hits) != len(hits):
        raise ValueError('need one hits table per sample directory, got {}'.format(hits))
    totals = annotation_misc_rna(annotations)
    extra = sorted(set(totals) - set(target_hits))
    if extra:
        raise ValueError('{}: samples {} are in it without hits tables'.format(annotations, extra))
    result = {}
    for sample in sorted(target_hits):
        counts = misc_rna_families(target_hits[sample])
        # A sample without reads may have an empty summary, so no rows in annotations.tsv: its misc_RNA total is 0, as for
        # a sample there without a misc_RNA row.
        source = 'annotations.tsv has' if sample in totals else 'annotations.tsv has no rows for the sample, so'
        check_total(counts, totals.get(sample, 0.0), target_hits[sample], source)
        result[sample] = counts
    return result


def decimal_count(value):
    """The shortest exact form of a count as a plain decimal, never in scientific notation (3.0, 0.5, 0.00005)."""
    text = format(decimal.Decimal(repr(value)), 'f')
    return text if '.' in text else text + '.0'  # 1e16 and above: repr has no decimal point


def format_table(result):
    lines = ['Sample_ID\tfamily\tcount']
    lines += ['{}\t{}\t{}'.format(sample, name, decimal_count(counts[name])) for sample, counts in result.items() for name in NAMES]
    return '\n'.join(lines) + '\n'


def mqc(result):
    section = {'id': 'unitas_misc_rna', 'section_name': 'UNITAS: misc_RNA families', 'description': DESCRIPTION}
    data = {sample: {name: round(value, 2) for name, value in counts.items()} for sample, counts in result.items()}
    if not any(any(counts.values()) for counts in data.values()):
        return {**section, 'plot_type': 'html', 'data': '<p>No sample has misc_RNA reads to plot (all below 0.005).</p>'}
    return {**section, 'plot_type': 'bargraph',
            'pconfig': {'id': 'unitas_misc_rna_plot', 'title': 'UNITAS: misc_RNA families', 'ylab': '# Reads'},
            'categories': {name: {'name': label} for name, label, _ in FAMILIES}, 'data': data}


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
