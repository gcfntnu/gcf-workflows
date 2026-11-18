#!/usr/bin/env python
"""UNITAS misc_RNA reads by gene family for the MultiQC report (#379), from misc_rna_families.tsv.

    unitas_misc_rna_mqc.py misc_rna_families.tsv > unitas_misc_rna_mqc.yaml

Writes MultiQC custom content as YAML in JSON notation: a bar per sample with its misc_RNA reads by family, or a
short note when no sample has misc_RNA reads. Families, labels and rules come from quant/scripts/unitas_misc_rna.py.
The UNITAS annotation plot keeps misc_RNA as one class; this section shows what it holds.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'quant' / 'scripts'))
from unitas_misc_rna import FAMILIES, NAMES, RULES  # noqa: E402

DESCRIPTION = ('Reads of the UNITAS misc_RNA class by gene family, from the names of the Ensembl genes they hit: {}. '
               'Reads split between genes keep their fractional counts.').format(RULES)


def families(path):
    """{sample: {family: count}} from misc_rna_families.tsv."""
    data = {}
    with open(path, newline='') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        if rows.fieldnames != ['Sample_ID', 'family', 'count']:
            raise ValueError('{}: not a misc_RNA family table (columns Sample_ID, family, count)'.format(path))
        for row in rows:
            if row['family'] not in NAMES:
                raise ValueError('{}: unknown misc_RNA family {!r}'.format(path, row['family']))
            data.setdefault(row['Sample_ID'], dict.fromkeys(NAMES, 0.0))[row['family']] = float(row['count'])
    return data


def mqc(data):
    section = {'id': 'unitas_misc_rna', 'section_name': 'UNITAS: misc_RNA families', 'description': DESCRIPTION}
    if not any(any(counts.values()) for counts in data.values()):
        return {**section, 'plot_type': 'html', 'data': '<p>No sample has misc_RNA reads.</p>'}
    return {**section, 'plot_type': 'bargraph',
            'pconfig': {'id': 'unitas_misc_rna_plot', 'title': 'UNITAS: misc_RNA families', 'ylab': '# Reads'},
            'categories': {name: {'name': label} for name, label, _ in FAMILIES},
            'data': {sample: {name: round(value, 2) for name, value in counts.items()} for sample, counts in sorted(data.items())}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('families', help='misc_rna_families.tsv from unitas_misc_rna.py')
    args = parser.parse_args()
    try:
        data = families(args.families)
    except (OSError, ValueError) as error:
        sys.exit('unitas_misc_rna_mqc.py: {}'.format(error))
    json.dump(mqc(data), sys.stdout, indent=1)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
