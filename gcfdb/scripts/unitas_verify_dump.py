#!/usr/bin/env python
"""Verify a cached UNITAS reference dump against its published archive checksum.

The release directory of the store (`<winecellar>/unitas/<release>/`) may publish `SHA256SUMS` for its
`UNITAS_refdump_<organism>.tar.gz` archives. When the archive is listed there, it is downloaded, its SHA-256
must equal the published one, and every file in the cached dump must be identical to the archive's copy: same
file names, same contents (apart from Snakemake's `.snakemake_timestamp` marker). A mismatch stops with an error naming the files. The cache itself is never changed.

When the release publishes no `SHA256SUMS` (releases before 04-11-2026), or does not list the archive, nothing
can be verified: the record says so (`status`), and the dump is still usable. Other download errors stop.

Writes a two-column TSV (field, value), replaced atomically.

Usage: unitas_verify_dump.py --dump DIR --release-url URL --output TSV
"""
import argparse
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
# Snakemake marks every directory() output with this file; it is not part of the dump.
SNAKEMAKE_MARKER = '.snakemake_timestamp'


def fail(message):
    sys.exit('unitas_verify_dump.py: ' + message)


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


def file_digests(root):
    digests = {}
    for path in sorted(root.rglob('*')):
        if path.is_file() or path.is_symlink():
            digest = hashlib.sha256()
            with open(path, 'rb') as fh:
                for block in iter(lambda: fh.read(CHUNK), b''):
                    digest.update(block)
            digests[str(path.relative_to(root))] = digest.hexdigest()
    return digests


def extract(archive, directory):
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.name.startswith('/') or '..' in Path(member.name).parts or not (member.isfile() or member.isdir()):
                fail('{} has an unexpected member {!r}'.format(archive, member.name))
        if hasattr(tarfile, 'data_filter'):
            tar.extractall(directory, filter='data')
        else:
            tar.extractall(directory)


def compare(dump, extracted):
    cached, published = file_digests(dump), file_digests(extracted)
    cached.pop(SNAKEMAKE_MARKER, None)
    problems = (['missing from the cache: ' + name for name in sorted(set(published) - set(cached))]
                + ['not in the archive: ' + name for name in sorted(set(cached) - set(published))]
                + ['differs from the archive: ' + name for name in sorted(set(cached) & set(published))
                   if cached[name] != published[name]])
    return len(published), problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--dump', required=True, type=Path, help='cached UNITAS_refdump_<organism> directory')
    parser.add_argument('--release-url', required=True, help='store directory of the release')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()

    if not args.dump.is_dir():
        fail('reference dump not found: {}'.format(args.dump))
    archive = args.dump.name + '.tar.gz'
    base = args.release_url.rstrip('/')
    sums_url, archive_url = base + '/SHA256SUMS', base + '/' + archive
    record = {'dump': str(args.dump), 'archive_url': archive_url, 'checksums_url': sums_url,
              'published_sha256': 'NA', 'files_compared': 'NA'}

    sums = fetch(sums_url)
    expected = published_checksum(sums, archive) if sums is not None else None
    if sums is None:
        record['status'] = 'not verified: the release publishes no SHA256SUMS'
    elif expected is None:
        record['status'] = 'not verified: SHA256SUMS does not list {}'.format(archive)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix='.verify-', dir=args.output.parent))
        try:
            got = fetch(archive_url, work / archive)
            if got is None:
                fail('{} is listed in {} but cannot be downloaded'.format(archive, sums_url))
            if got != expected:
                fail('{} has SHA-256 {}, but {} publishes {}'.format(archive_url, got, sums_url, expected))
            extract(work / archive, work / 'extracted')
            extracted = work / 'extracted' / args.dump.name
            if not extracted.is_dir():
                fail('{} does not hold {}/'.format(archive_url, args.dump.name))
            count, problems = compare(args.dump, extracted)
            if problems:
                fail('the cached dump {} is not the published archive {} ({} differences):\n  {}'.format(
                    args.dump, archive_url, len(problems), '\n  '.join(problems[:20])))
        finally:
            shutil.rmtree(work, ignore_errors=True)
        record.update(status='verified', published_sha256=expected, files_compared=str(count))
    record['checked'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_name('.{}.{}.part'.format(args.output.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, args.output)


if __name__ == '__main__':
    main()
