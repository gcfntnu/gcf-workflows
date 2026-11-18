#!/usr/bin/env python
"""UNITAS misc_RNA reads by gene family for the MultiQC report (#379), from the misc_RNA rows of annotations.tsv.

    unitas_misc_rna_mqc.py annotations.tsv > unitas_misc_rna_mqc.yaml

Writes MultiQC custom content (plot_type: bargraph) as YAML in JSON notation: one bar per sample with its misc_RNA
reads split into Y RNA, 7SK, 7SL and other misc_RNA (family rules in quant/scripts/unitas_misc_rna.py).
The UNITAS annotation plot keeps misc_RNA as one class; this section shows what it holds.
"""
import argparse
import csv
import json
import sys

CATEGORIES = {'Y_RNA': 'Y RNA', '7SK': '7SK', '7SL': '7SL (SRP RNA)', 'other': 'other misc_RNA'}
DESCRIPTION = ('Reads of the UNITAS misc_RNA class by gene family, from the names of the Ensembl genes they hit: '
               'RNY and Y_RNA genes as Y RNA, RN7SK and 7SK as 7SK, RN7SL and Metazoa_SRP as 7SL; other names as other '
               'misc_RNA. Reads split between genes keep their fractional counts.')


def families(path):
    """misc_RNA family counts per sample; every sample in the table, families it lacks as 0."""
    data = {}
    with open(path, newline='') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        missing = {'Level1', 'Level2', 'Level3', 'Count', 'Sample_ID'} - set(rows.fieldnames or [])
        if missing:
            raise ValueError('{}: not a UNITAS annotations table (no {})'.format(path, ', '.join(sorted(missing))))
        for row in rows:
            sample = data.setdefault(row['Sample_ID'], dict.fromkeys(CATEGORIES, 0.0))
            if row['Level1'] != 'misc_RNA' or not row['Level2']:
                continue
            if row['Level2'] not in CATEGORIES:
                raise ValueError('{}: unknown misc_RNA family {!r}'.format(path, row['Level2']))
            sample[row['Level2']] = float(row['Count'])
    return data


def mqc(data):
    return {'id': 'unitas_misc_rna', 'section_name': 'UNITAS: misc_RNA families', 'description': DESCRIPTION,
            'plot_type': 'bargraph',
            'pconfig': {'id': 'unitas_misc_rna_plot', 'title': 'UNITAS: misc_RNA families', 'ylab': '# Reads'},
            'categories': {key: {'name': name} for key, name in CATEGORIES.items()},
            'data': {sample: {key: round(value, 2) for key, value in counts.items()} for sample, counts in sorted(data.items())}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('annotations', help='annotations.tsv from unitas_annotations.py')
    args = parser.parse_args()
    try:
        data = families(args.annotations)
    except (OSError, ValueError) as error:
        sys.exit('unitas_misc_rna_mqc.py: {}'.format(error))
    json.dump(mqc(data), sys.stdout, indent=1)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
