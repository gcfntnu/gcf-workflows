#!/usr/bin/env python
"""Download, verify and unpack a versioned miRge3 library.

The release directory of the store (`<winecellar>/mirge3/<release>/`) publishes `miRge3_lib_<species>.tar.gz`
archives and `SHA256SUMS`, written by gcfntnu/gcf-docker `mirge3-build`. Each archive holds one top-level `<org>/`
directory, miRge3's `--organism-name`.

The archive's SHA-256 must equal the published one before anything is unpacked. The unpacked library must hold
every file miRge3 0.1.4 opens for the options gcf-workflows uses (both miRNA databases; tRNA annotation for human)
and the provenance files `gcf_mirna_db.tsv`, `gcf_manifest.tsv` and `gcf_sources.tsv`, which stay in the library.
The library's own spike-in index (`index.Libs/<org>_spike-in.*`) is removed: spike-ins are a library prep choice,
and smallrna builds that index from `filter.spikein.ref` before each miRge3 run.

The library is moved into `<libraries-path>/<org>/` only when every check passed. An existing `<org>/` without
this record is not replaced: the script stops and names it. A missing `SHA256SUMS`, an unlisted archive or any
failed check stops with an error.

Writes a two-column TSV (field, value), replaced atomically.

Usage: mirge3_library.py --release-url URL --release LABEL --species homo_sapiens --organism-name human
                         --libraries-path DIR --output TSV
"""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request

CHUNK = 1 << 20
DATABASES = ('miRBase', 'MirGeneDB')
INDEXES = ('genome', 'mature_trna', 'pre_trna', 'mrna', 'ncrna_others', 'rrna', 'snorna')
TRNA_ANNOTATION = ('trna.str', 'trna_aminoacid_anticodon.csv', 'trna_deduplicated_list.csv', 'tRF_infor.csv', 'tRF_merges.csv')
PROVENANCE = ('gcf_mirna_db.tsv', 'gcf_manifest.tsv', 'gcf_sources.tsv')


def fail(message):
    sys.exit('mirge3_library.py: ' + message)


def library_files(org):
    """Files of <org>/ that miRge3 0.1.4 opens with --gff-out --isoform-entropy --AtoI (--tRNA-frag for human)."""
    files = []
    for db in DATABASES:
        files += ['annotation.Libs/{}_{}.gff3'.format(org, db), 'annotation.Libs/{}_merges_{}.csv'.format(org, db),
                  'annotation.Libs/{}_miRNAs_in_repetitive_element_{}.csv'.format(org, db),
                  'fasta.Libs/{}_mature_{}.fa'.format(org, db), 'fasta.Libs/{}_mirna_SNP_pseudo_{}.fa'.format(org, db),
                  'index.Libs/{}_hairpin_{}.1.ebwt'.format(org, db), 'index.Libs/{}_mirna_{}.1.ebwt'.format(org, db)]
    files += ['index.Libs/{}_{}.1.ebwt'.format(org, index) for index in INDEXES]
    if org == 'human':
        files += ['annotation.Libs/{}_{}'.format(org, name) for name in TRNA_ANNOTATION]
    return files + list(PROVENANCE)


def fetch(url, destination=None):
    """Return the body of url (or write it to destination and return its SHA-256); None when it does not exist."""
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            if destination is None:
                return response.read().decode()
            digest = hashlib.sha256()
            with open(destination, 'wb') as out:
                for block in iter(lambda: response.read(CHUNK), b''):
                    digest.update(block)
                    out.write(block)
            return digest.hexdigest()
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        fail('cannot download {}: HTTP {}'.format(url, error.code))
    except urllib.error.URLError as error:
        if url.startswith('file:') and isinstance(error.reason, FileNotFoundError):
            return None
        fail('cannot download {}: {}'.format(url, error.reason))
    except ValueError as error:
        fail('cannot download {} (is winecellar.url set?): {}'.format(url, error))


def published_checksum(sums, archive):
    for line in sums.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1].lstrip('*') == archive:
            return fields[0].lower()
    return None


def extract(archive, directory, org):
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            parts = Path(member.name).parts
            if member.name.startswith('/') or '..' in parts or not (member.isfile() or member.isdir()):
                fail('{} has an unexpected member {!r}'.format(archive.name, member.name))
            if parts and parts[0] != org:
                fail('{} holds {!r}, not only the {}/ directory of --organism-name {}'.format(archive.name, member.name, org, org))
        if hasattr(tarfile, 'data_filter'):
            tar.extractall(directory, filter='data')
        else:
            tar.extractall(directory)


def read_table(path, key, value):
    with open(path) as fh:
        return {row[key]: row[value] for row in csv.DictReader(fh, delimiter='\t')}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--release-url', required=True, help='store directory of the release')
    parser.add_argument('--release', required=True, help='release label, e.g. 09-11-2026')
    parser.add_argument('--species', required=True, help='e.g. homo_sapiens')
    parser.add_argument('--organism-name', required=True, help="miRge3's --organism-name, e.g. human")
    parser.add_argument('--libraries-path', required=True, type=Path, help="miRge3's --libraries-path")
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()

    org = args.organism_name
    target = args.libraries_path / org
    if target.exists():
        fail('{} exists but has no verification record {}; it was not installed by this release check. '
             'Remove it to download release {} again'.format(target, args.output, args.release))
    archive = 'miRge3_lib_{}.tar.gz'.format(args.species)
    base = args.release_url.rstrip('/')
    sums_url, archive_url = base + '/SHA256SUMS', base + '/' + archive

    sums = fetch(sums_url)
    if sums is None:
        fail('release {} publishes no SHA256SUMS ({})'.format(args.release, sums_url))
    expected = published_checksum(sums, archive)
    if expected is None:
        fail('{} does not list {}'.format(sums_url, archive))

    args.libraries_path.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='.mirge3-', dir=args.libraries_path))
    try:
        got = fetch(archive_url, work / archive)
        if got is None:
            fail('{} is listed in {} but cannot be downloaded'.format(archive, sums_url))
        if got != expected:
            fail('{} has SHA-256 {}, but {} publishes {}'.format(archive_url, got, sums_url, expected))
        extract(work / archive, work / 'extracted', org)
        library = work / 'extracted' / org
        missing = [name for name in library_files(org) if not (library / name).is_file()]
        if missing:
            fail('{} lacks {} file(s) miRge3 or the provenance record needs:\n  {}'.format(
                archive_url, len(missing), '\n  '.join(missing[:20])))
        removed = sorted(path.name for path in (library / 'index.Libs').glob('{}_spike-in.*'.format(org)))
        for name in removed:
            (library / 'index.Libs' / name).unlink()
        mirna = read_table(library / 'gcf_mirna_db.tsv', 'database', 'release')
        manifest = read_table(library / 'gcf_manifest.tsv', 'key', 'value')
        os.replace(library, target)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    record = {'release': args.release, 'species': args.species, 'mirge3_organism': org, 'library': str(target),
              'archive_url': archive_url, 'checksums_url': sums_url, 'published_sha256': expected, 'status': 'verified',
              'files_checked': str(len(library_files(org))), 'mirbase_release': mirna.get('mirbase', 'NA'),
              'mirgenedb_release': mirna.get('mirgenedb', 'NA'), 'ensembl_release': manifest.get('ensembl_release', 'NA'),
              'removed_spikein_index': ','.join(removed) or 'none',
              'checked': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name('.{}.{}.part'.format(args.output.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, args.output)


if __name__ == '__main__':
    main()
