#!/usr/bin/env python
"""Record the selected miRge3 library and check that it holds the requested miRNA database.

The library must hold the files miRge3 0.1.4 opens for `quant.mirge3.db` (annotation, FASTA and bowtie indexes of
that database); a missing file stops with an error naming it. The miRNA release is read from the library's
`gcf_mirna_db.tsv`, other sources from `gcf_manifest.tsv`, and the archive check from the `mirge3_library.py`
record. Libraries without these files (the SourceForge library) are recorded as `NA`; nothing is filled in from
current settings. miRge3 counts are not joined to another annotation, so no other release is checked.

Writes a three-column TSV: field, value and the file the value comes from.

Usage: mirge3_reference.py --library DIR --mirna-db mirbase|mirgenedb --release LABEL [--verification TSV] --output TSV
"""
import argparse
import csv
import os
from pathlib import Path
import sys

DATABASES = {'mirbase': 'miRBase', 'mirgenedb': 'MirGeneDB'}
MANIFEST = ('ensembl_release', 'genome_assembly', 'mirbase_gff_assembly', 'mirgenedb_gff_assembly', 'bowtie', 'builder_image')


def fail(message):
    sys.exit('mirge3_reference.py: ' + message)


def database_files(org, db):
    return ['annotation.Libs/{}_{}.gff3'.format(org, db), 'annotation.Libs/{}_merges_{}.csv'.format(org, db),
            'annotation.Libs/{}_miRNAs_in_repetitive_element_{}.csv'.format(org, db),
            'fasta.Libs/{}_mature_{}.fa'.format(org, db), 'fasta.Libs/{}_mirna_SNP_pseudo_{}.fa'.format(org, db),
            'index.Libs/{}_hairpin_{}.1.ebwt'.format(org, db), 'index.Libs/{}_mirna_{}.1.ebwt'.format(org, db)]


def read_table(path, key, value):
    if not path.is_file():
        return None
    with open(path) as fh:
        return {row[key]: row[value] for row in csv.DictReader(fh, delimiter='\t')}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--library', required=True, type=Path, help='the miRge3 library, <libraries-path>/<org>')
    parser.add_argument('--mirna-db', required=True, choices=sorted(DATABASES))
    parser.add_argument('--release', required=True, help='db.mirge3.release')
    parser.add_argument('--verification', type=Path, help='record of mirge3_library.py; absent for the SourceForge library')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()

    if not args.library.is_dir():
        fail('miRge3 library not found: {}'.format(args.library))
    org, db = args.library.name, DATABASES[args.mirna_db]
    missing = [name for name in database_files(org, db) if not (args.library / name).is_file()]
    if missing:
        fail('quant.mirge3.db is {}, but the miRge3 library {} (db.mirge3.release {}) lacks:\n  {}'.format(
            args.mirna_db, args.library, args.release, '\n  '.join(missing)))

    rows = [('release', args.release, 'db.mirge3.release'), ('library', str(args.library), 'db.mirge3.release'),
            ('mirna_database', args.mirna_db, 'quant.mirge3.db')]
    mirna = read_table(args.library / 'gcf_mirna_db.tsv', 'database', 'release')
    if mirna is None:
        rows.append(('mirna_release', 'NA', 'NA'))
    elif args.mirna_db not in mirna:
        fail('{}/gcf_mirna_db.tsv lists {}, not quant.mirge3.db {}'.format(args.library, ', '.join(mirna) or 'nothing', args.mirna_db))
    else:
        rows.append(('mirna_release', mirna[args.mirna_db], 'gcf_mirna_db.tsv'))
    manifest = read_table(args.library / 'gcf_manifest.tsv', 'key', 'value') or {}
    rows += [(key, manifest[key], 'gcf_manifest.tsv') if key in manifest else (key, 'NA', 'NA') for key in MANIFEST]
    sources = args.library / 'gcf_sources.tsv'
    if sources.is_file():
        with open(sources) as fh:
            rows.append(('source_checksums', '{} files with SHA-256'.format(sum(1 for _ in fh) - 1), 'gcf_sources.tsv'))
    else:
        rows.append(('source_checksums', 'NA', 'NA'))
    if args.verification:
        record = read_table(args.verification, 'field', 'value')
        if record is None:
            fail('verification record not found: {}'.format(args.verification))
        rows += [('archive_status', record.get('status', 'NA'), args.verification.name),
                 ('archive_published_sha256', record.get('published_sha256', 'NA'), args.verification.name)]
    else:
        rows += [('archive_status', 'not verified: library without published checksums', 'NA'),
                 ('archive_published_sha256', 'NA', 'NA')]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name('.{}.{}.part'.format(args.output.name, os.getpid()))
    partial.write_text('field\tvalue\tsource\n' + ''.join('{}\t{}\t{}\n'.format(*row) for row in rows))
    os.replace(partial, args.output)


if __name__ == '__main__':
    main()
