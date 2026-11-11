#!/usr/bin/env python
"""Record the selected UNITAS reference dump and check its miRNA annotation is compatible.

Reads what the dump records about itself: `gcf_mirna_db.tsv` (dumps from 05-11-2026 on with MirGeneDB, written
by unitas_build_refdump.sh), UNITAS' own `db_versions.info` (every dump) and `gcf_sources.tsv` (dumps from
04-11-2026 on), plus the archive verification record of `unitas_verify_dump.py`. Nothing missing is filled in
from current defaults: an absent record is written as `NA`.

The miRNA feature annotation (miRBase GFF for `db.mirbase.release`, or the MirGeneDB GFF in the dump) is joined
to UNITAS' miRNA counts by name, so it must come from the same miRNA database and release as the dump. The script
stops with an error naming both when they differ, when the dump is not from `db.unitas.mirna_db`, or when the
dump's miRNA release cannot be read. Other sources (Ensembl, GtRNAdb, ...) are recorded, not checked: UNITAS'
reference is independent of the genome alignment reference.

Writes a three-column TSV: field, value and the file the value comes from.

Usage: unitas_reference.py --dump DIR --verification TSV --mirna-db mirbase|mirgenedb --annotation GFF
                           --release LABEL [--mirbase-releases 20,21,22,23.1] --output TSV
"""
import argparse
import csv
from pathlib import Path
import re
import sys

NAMES = {'mirbase': 'miRBase', 'mirgenedb': 'MirGeneDB'}
# db_versions.info: "<name> version/date: ....... <value>", e.g. "Ensembl version/date: ...... Release 116"
VERSION_LINE = re.compile(r'^(?P<name>.+?)\s+version/date:\s*\.+\s*(?P<value>.*?)\s*$')
FIELDS = {'Ensembl': 'ensembl_release', 'EnsemblGenomes': 'ensemblgenomes_release', 'Genomic tRNA database': 'gtrnadb',
          'piRNA cluster database': 'pirna_clusters', 'tRF-1 sequence data': 'trf1', 'tRNA-leader sequence data': 'trna_leader',
          'SILVA rRNA (SSU) database': 'silva_ssu', 'SILVA rRNA (LSU) database': 'silva_lsu', 'SeqMap': 'seqmap'}


def fail(message):
    sys.exit('unitas_reference.py: ' + message)


def mirna_from_versions(value):
    """(database, release) from the db_versions.info miRBase line: 'Release 23.1', '22.1' or 'MirGeneDB 3.0'."""
    match = re.fullmatch(r'MirGeneDB\s+(\S+)', value)
    if match:
        return 'mirgenedb', match.group(1)
    match = re.fullmatch(r'(?:Release\s+)?(\d+(?:\.\d+)?)', value)
    return ('mirbase', match.group(1)) if match else None


def mirna_from_gff(path):
    with open(path) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            match = re.match(r'^# microRNAs:\s+(miRBase|MirGeneDB) v(\S+)\s*$', line)
            if match:
                return match.group(1).lower(), match.group(2)
    fail('{} has no "# microRNAs: miRBase v<release>" or "MirGeneDB v<release>" header line'.format(path))


def describe(mirna):
    return '{} {}'.format(NAMES[mirna[0]], mirna[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--dump', required=True, type=Path)
    parser.add_argument('--verification', required=True, type=Path)
    parser.add_argument('--mirna-db', required=True, choices=sorted(NAMES))
    parser.add_argument('--annotation', required=True, type=Path)
    parser.add_argument('--release', required=True, help='selected dump release (db.unitas.release or mirgenedb_release)')
    parser.add_argument('--mirbase-releases', default='', help='comma-separated miRBase releases gcfdb can annotate from')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()

    rows = [('dump', str(args.dump), 'configuration'), ('dump_release', args.release, 'configuration')]
    with open(args.verification) as fh:
        verification = dict(row for row in csv.reader(fh, delimiter='\t'))
    for field in ('status', 'published_sha256', 'checked'):
        rows.append(('archive_' + field, verification.get(field, 'NA'), args.verification.name))

    versions, info = {}, args.dump / 'db_versions.info'
    if info.exists():
        for line in info.read_text().splitlines():
            match = VERSION_LINE.match(line.strip())
            if match:
                versions[match.group('name')] = match.group('value')
    from_versions = mirna_from_versions(versions['miRBase database']) if 'miRBase database' in versions else None

    recorded, mirna_tsv = None, args.dump / 'gcf_mirna_db.tsv'
    if mirna_tsv.exists():
        table = list(csv.DictReader(mirna_tsv.open(), delimiter='\t'))
        if len(table) != 1 or table[0].get('database') not in NAMES or not table[0].get('release'):
            fail('{} must have one row with database (mirbase or mirgenedb) and release'.format(mirna_tsv))
        recorded = (table[0]['database'], table[0]['release'])
    if recorded and from_versions and recorded != from_versions:
        fail('{} records {}, but db_versions.info records {}'.format(mirna_tsv, describe(recorded), describe(from_versions)))
    dump_mirna = recorded or from_versions
    if dump_mirna is None:
        fail('cannot read the miRNA database release of {} (no gcf_mirna_db.tsv, and db_versions.info has no readable '
             'miRBase line), so its compatibility with the miRNA annotation cannot be checked'.format(args.dump))
    source = mirna_tsv.name if recorded else info.name
    rows += [('mirna_database', dump_mirna[0], source), ('mirna_release', dump_mirna[1], source)]
    for name, field in FIELDS.items():
        value = versions.get(name, 'NA')
        value = re.sub(r'^Release\s+', '', value)
        date = re.fullmatch(r'(\d\d)\.(\d\d)\.(\d{4}) \(dd\.mm\.yyyy\)', value)
        if date:
            value = '{}-{}-{}'.format(date.group(3), date.group(2), date.group(1))
        rows.append((field, value, info.name if name in versions else 'NA'))
    sources = args.dump / 'gcf_sources.tsv'
    if sources.exists():
        count = sum(1 for _ in csv.DictReader(sources.open(), delimiter='\t'))
        rows.append(('source_checksums', '{} files with SHA-256'.format(count), sources.name))
    else:
        rows.append(('source_checksums', 'NA', 'NA'))

    if dump_mirna[0] != args.mirna_db:
        fail('db.unitas.mirna_db is {}, but the selected dump {} was built from {}'.format(args.mirna_db, args.dump, describe(dump_mirna)))
    annotation = mirna_from_gff(args.annotation)
    rows += [('annotation', str(args.annotation), 'configuration'),
             ('annotation_database', annotation[0], args.annotation.name), ('annotation_release', annotation[1], args.annotation.name)]
    if annotation != dump_mirna:
        served = [r for r in args.mirbase_releases.split(',') if r]
        if args.mirna_db == 'mirgenedb':
            hint = 'the dump and its gcf_mirgenedb.gff3 disagree; the dump is inconsistent'
        elif dump_mirna[1] in served:
            hint = 'select a dump built from {} (db.unitas.release) or set db.mirbase.release to {}'.format(describe(annotation), dump_mirna[1])
        else:
            hint = 'select a dump built from {} (db.unitas.release); no annotation for {} is available (db.mirbase.release: {})'.format(
                describe(annotation), describe(dump_mirna), ', '.join(served) or 'NA')
        fail('the UNITAS dump {} (release {}) was built from {}, but the miRNA feature annotation {} is {}. '
             'UNITAS counts are joined to the annotation by miRNA name, which differs between releases: {}'.format(
                 args.dump, args.release, describe(dump_mirna), args.annotation, describe(annotation), hint))

    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['field', 'value', 'source'])
        writer.writerows(rows)


if __name__ == '__main__':
    main()
