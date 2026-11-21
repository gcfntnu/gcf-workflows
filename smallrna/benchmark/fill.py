"""Fill a public dataset into EXT_DIR: download its runs (and a pool's composition table) once, checked against checksums.

A pool declaration names a data folder relative to EXT_DIR (`data`), the composition table (`composition`: file, url,
sha256) and the runs (`runs`: run, url, md5, ...); a biological declaration has no composition and lists its runs
per sample (`samples`: sample, group, runs); a genome-wide synthetic declaration (`synthetic.mode: genome`, #412) has
no runs, only the miRBase files its reads are generated from. Files go to `<EXT_DIR>/<data>/`: runs as
`fastq/<run>.fastq.gz`, the table as `composition/<file>`, the miRBase mature sequences (`mirbase`: release, url,
sha256) as `mirbase/<release>/mature.fa`, and a synthetic declaration's hairpins and GFF (`mirbase.hairpin`,
`mirbase.gff`: url, sha256) as `mirbase/<release>/hairpin.fa` and `genes.gff3`. Each download is written to its own
temporary file in the same folder and renamed only after its checksum matches, so a file in place is always complete,
also when fills of declarations sharing a data folder run at the same time; files already in place are skipped. A
transfer the server cuts before the announced length is resumed from where it stopped, on a new connection (one at a
time, with a pause before each). This is
the only step of the benchmark that writes to EXT_DIR, and it is run deliberately when a dataset is added; `prepare` and
`evaluate` only read.
"""
import hashlib
import http.client
import os
from pathlib import Path
import shutil
import tempfile
import time
import urllib.error
import urllib.request

CHUNK = 1 << 20
CONNECTIONS = 20   # per file: the first, and one for each time the transfer is cut and resumed
PAUSE = 30         # seconds before a resuming connection, times its number, at most 300
UMASK = os.umask(0)
os.umask(UMASK)


def targets(spec, ext_dir):
    """[(url, algorithm, digest, final path)] for every file of a pool, biological or genome-wide synthetic declaration."""
    data = Path(ext_dir) / spec['data']
    mirbase = spec['mirbase']
    files = []
    if spec['kind'] == 'pool':
        comp = spec['composition']
        files.append((comp['url'], 'sha256', comp['sha256'], data / 'composition' / comp['file']))
    files.append((mirbase['url'], 'sha256', mirbase['sha256'], data / 'mirbase' / str(mirbase['release']) / 'mature.fa'))
    if spec['kind'] == 'synthetic':
        for key, name in (('hairpin', 'hairpin.fa'), ('gff', 'genes.gff3')):
            files.append((mirbase[key]['url'], 'sha256', mirbase[key]['sha256'], data / 'mirbase' / str(mirbase['release']) / name))
        sources = (spec['synthetic'].get('decoys') or {}).get('sources') or {}  # #412 PR 2
        files += [(value['url'], 'sha256', value['sha256'], decoy_path(spec, ext_dir, key, value['url'])) for key, value in sorted(sources.items())]
        return files
    runs = spec['runs'] if spec['kind'] == 'pool' else [run for sample in spec['samples'] for run in sample['runs']]
    files += [(run['url'], 'md5', run['md5'], data / 'fastq' / '{}.fastq.gz'.format(run['run'])) for run in runs]
    return files


def transfer(url, out, checksum, log=print):
    """Write url's bytes to `out` and into `checksum`. A server can close a long transfer early, which reads like the
    end of the file: when fewer bytes arrived than the response announced (Content-Length), or the connection broke,
    the rest is asked for on a new connection (a Range request from the bytes already written), up to CONNECTIONS in
    all, one at a time. An HTTP error status is not retried. Without an announced length there is nothing to compare
    with, and the caller's checksum decides."""
    received, total = 0, None
    for connection in range(CONNECTIONS):
        if connection:
            log('{}: cut at {} of {} bytes; resuming'.format(url, received, total if total is not None else 'unknown'))
            time.sleep(min(PAUSE * connection, 300))
        request = urllib.request.Request(url, headers={'Range': 'bytes={}-'.format(received)}) if received else url
        try:
            response = urllib.request.urlopen(request, timeout=600)
        except urllib.error.HTTPError:
            raise
        except (OSError, http.client.HTTPException) as error:  # no connection: retried only for a transfer under way
            if not received:
                raise
            log('{}: {}'.format(url, error))
            continue
        with response:
            if received and getattr(response, 'status', None) != 206:
                raise ValueError('{}: no partial content for a resumed download (HTTP {})'.format(url, getattr(response, 'status', None)))
            headers = getattr(response, 'headers', None)
            length = headers.get('Content-Length') if headers is not None else None
            if length is not None:
                total = received + int(length)
            while True:
                try:
                    block = response.read(CHUNK)
                except (OSError, http.client.HTTPException) as error:  # a reset or timed-out connection, an incomplete read
                    if total is None:
                        raise  # no announced length to resume against
                    log('{}: {}'.format(url, error))
                    break
                if not block:
                    break
                out.write(block)  # a local error (e.g. a full disk) is not a cut transfer: it ends the download
                checksum.update(block)
                received += len(block)
        if total is None or received >= total:
            return
    raise ValueError('{}: incomplete after {} connections ({} of {} bytes)'.format(url, CONNECTIONS, received, total))


def decoy_path(spec, ext_dir, key, url):
    """Where a genome-wide synthetic declaration's decoy source `key` is filled: <data>/decoys/<key>.fa, or .fa.gz for a
    gzipped download."""
    return Path(ext_dir) / spec['data'] / 'decoys' / '{}.fa{}'.format(key, '.gz' if url.split('?')[0].endswith('.gz') else '')


def download(url, algorithm, digest, path, log=print):
    """Download url to path through a temporary file of this call only; fail (and remove that file) when the checksum
    differs. A concurrent download of the same file writes its own temporary file, so neither can change the other's."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.partial-{}.'.format(path.name), dir=path.parent)
    partial = Path(name)
    checksum = hashlib.new(algorithm)
    try:
        with os.fdopen(fd, 'wb') as out:
            transfer(url, out, checksum, log)
        if checksum.hexdigest() != digest:
            raise ValueError('{}: {} {} does not match the declared {}'.format(url, algorithm, checksum.hexdigest(), digest))
        os.chmod(partial, 0o666 & ~UMASK)  # mkstemp creates it readable by the owner only
        os.replace(partial, path)
    finally:
        if partial.exists():
            partial.unlink()
    log('{}: {} checked'.format(path, algorithm))


def fill(spec, ext_dir, log=print, fetch=download):
    """Download the declaration's missing files into EXT_DIR. Returns the number of files downloaded."""
    genome = spec.get('kind') == 'synthetic' and (spec.get('synthetic') or {}).get('mode') == 'genome'
    if spec.get('kind') not in ('pool', 'biological') and not genome:
        raise ValueError('{}: only pool, biological and genome-wide synthetic datasets are filled; a reduced synthetic dataset is '
                         'generated by prepare from its reference'.format(spec['name']))
    done = 0
    for url, algorithm, digest, path in targets(spec, ext_dir):
        if path.exists():
            log('{}: in place'.format(path))
            continue
        fetch(url, algorithm, digest, path, log)
        done += 1
    free = shutil.disk_usage(Path(ext_dir) / spec['data']).free
    log('{}: {} file(s) downloaded, {:.0f} GB free'.format(spec['name'], done, free / 1e9))
    return done
