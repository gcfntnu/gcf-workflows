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
time, with a pause before each).

Runs can also come from NCBI's copy of the Sequence Read Archive on AWS (`fill(source='sra')`), for runs that declare
`sra: {md5, bytes}` (the SRA file's, from NCBI's locator, https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=<run>).
ENA allows one slow connection; this copy takes many range requests at once. The SRA file is downloaded in parts,
checked against its MD5, converted with sra-tools' fasterq-dump and written as `fastq/<run>.fastq.gz`, with
`fastq/<run>.sra.json` recording where it came from. Its reads are the ones in ENA's file, in the same order, with
other read names: the file is checked by its number of reads and, when the run declares `content_md5`, by the MD5 over
its sequence and quality lines (`content`), which is the same for a file from either archive. Single-end runs only.

This is
the only step of the benchmark that writes to EXT_DIR, and it is run deliberately when a dataset is added; `prepare` and
`evaluate` only read.
"""
import concurrent.futures
import gzip
import hashlib
import http.client
import itertools
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import time
import urllib.error
import urllib.request

CHUNK = 1 << 20
CONNECTIONS = 20   # per file: the first, and one for each time the transfer is cut and resumed
PAUSE = 30         # seconds before a resuming connection, times its number, at most 300
SRA_URL = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}'  # NCBI's open-data copy of public runs, on AWS
PART = 32 << 20    # bytes per range request of an SRA file
PART_ATTEMPTS = 5  # requests per part
UMASK = os.umask(0)
os.umask(UMASK)
USER_AGENT = 'gcf-workflows-benchmark'  # the client name downloads send


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
    files += [(run['url'], 'md5', run['md5'], data / 'fastq' / '{}.fastq.gz'.format(run['run'])) for run in runs(spec)]
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
        # a client name of its own on every request: some servers (GtRNAdb) refuse Python's default
        headers = dict({'User-Agent': USER_AGENT}, **({'Range': 'bytes={}-'.format(received)} if received else {}))
        request = urllib.request.Request(url, headers=headers)
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


def content(lines, out=None):
    """(reads, MD5 over the sequence and quality lines) of FASTQ `lines` (an iterable of byte blocks, cut anywhere). The
    MD5 does not depend on read names or compression, so it is the same for a run's file from ENA and from SRA. With
    `out`, the records are also written there, with a bare `+` as third line."""
    digest, reads, rest = hashlib.md5(), 0, b''
    for block in itertools.chain(lines, [b'\n']):  # a last line without its newline still ends its record
        rows = (rest + block).split(b'\n')
        rest = rows.pop()
        whole = len(rows) // 4 * 4
        rest = b'\n'.join(rows[whole:] + [rest])
        rows = rows[:whole]
        if not rows:
            continue
        if any(not name.startswith(b'@') for name in rows[0::4]) or any(not plus.startswith(b'+') for plus in rows[2::4]):
            raise ValueError('not FASTQ: a record without `@` or `+` at read {}'.format(reads + 1))
        if any(not seq or len(seq) != len(quality) for seq, quality in zip(rows[1::4], rows[3::4])):
            raise ValueError('not FASTQ: a record without bases, or with other than one quality per base, after read {}'.format(reads))
        pairs = [None] * (whole // 2)
        pairs[0::2], pairs[1::2] = rows[1::4], rows[3::4]
        digest.update(b'\n'.join(pairs) + b'\n')
        reads += whole // 4
        if out is not None:
            rows[2::4] = [b'+'] * (whole // 4)
            out.write(b'\n'.join(rows) + b'\n')
    if rest.strip():
        raise ValueError('not FASTQ: an incomplete record after read {}'.format(reads))
    return reads, digest.hexdigest()


def file_content(path):
    """content() of a gzipped FASTQ file."""
    with gzip.open(path, 'rb') as source:
        return content(iter(lambda: source.read(CHUNK * 16), b''))


def fetch_parts(url, size, path, connections, log=print):
    """Download url (`size` bytes) to path with range requests of PART bytes, `connections` at a time; each part is
    asked for up to PART_ATTEMPTS times. An HTTP error status is not retried."""
    spans = [(start, min(start + PART, size) - 1) for start in range(0, size, PART)]
    with open(path, 'wb') as out:
        out.truncate(size)
    fd = os.open(path, os.O_WRONLY)

    def part(span):
        start, end = span
        for attempt in range(PART_ATTEMPTS):
            if attempt:
                time.sleep(min(PAUSE * attempt, 300))
            try:
                request = urllib.request.Request(url, headers={'User-Agent': USER_AGENT, 'Range': 'bytes={}-{}'.format(start, end)})
                with urllib.request.urlopen(request, timeout=300) as response:
                    if getattr(response, 'status', None) != 206:
                        raise ValueError('{}: no partial content for bytes {}-{} (HTTP {})'.format(url, start, end, getattr(response, 'status', None)))
                    data = response.read()
            except urllib.error.HTTPError:
                raise
            except (OSError, http.client.HTTPException) as error:
                log('{}: bytes {}-{}: {}'.format(url, start, end, error))
                continue
            if len(data) == end - start + 1:
                os.pwrite(fd, data, start)
                return
            log('{}: bytes {}-{}: {} bytes received'.format(url, start, end, len(data)))
        raise ValueError('{}: bytes {}-{} not received in {} requests'.format(url, start, end, PART_ATTEMPTS))

    try:
        with concurrent.futures.ThreadPoolExecutor(connections) as pool:
            list(pool.map(part, spans))
    finally:
        os.close(fd)


def file_md5(path):
    digest = hashlib.md5()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(CHUNK * 16), b''):
            digest.update(block)
    return digest.hexdigest()


def sra_run(run, path, connections=8, fasterq_dump='fasterq-dump', threads=4, log=print):
    """Fill one run from NCBI's copy of SRA on AWS: download `run`'s SRA file (declared `sra`: md5, bytes) in parts,
    check its MD5, convert it with fasterq-dump (`fasterq_dump`: the command, e.g. through a container) and write the
    reads to `path`, gzipped, after their number matches the declared `reads` and, when declared, their content() MD5
    matches `content_md5`. `<run>.sra.json` beside it records the source. Everything is done in a temporary folder
    next to `path`, removed afterwards. Returns the content MD5."""
    name, declared = run['run'], run.get('sra') or {}
    if not declared.get('md5') or not declared.get('bytes') or 'reads' not in run:
        raise ValueError('{}: a run filled from SRA declares `reads` and `sra: {{md5, bytes}}`; fill it from ENA, or declare them'.format(name))
    url = SRA_URL.format(run=name)
    path.parent.mkdir(parents=True, exist_ok=True)
    # absolute, as fasterq-dump runs inside it, but not resolved: a store reached through a link keeps the path the
    # caller named, which is the one a container has bound
    work = Path(os.path.abspath(tempfile.mkdtemp(prefix='.partial-{}.sra.'.format(name), dir=path.parent)))
    try:
        sra = work / '{}.sra'.format(name)  # fasterq-dump names the reads after the file
        fetch_parts(url, int(declared['bytes']), sra, connections, log)
        found = file_md5(sra)
        if found != declared['md5']:
            raise ValueError('{}: md5 {} does not match the declared {}'.format(url, found, declared['md5']))
        command = [*shlex.split(fasterq_dump), str(sra), '--stdout', '--split-spot', '--skip-technical', '--threads', str(threads),
                   '--temp', str(work)]
        packed = work / path.name
        messages = work / 'fasterq-dump.log'
        with open(messages, 'wb') as errors, subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, cwd=work) as process, \
                gzip.open(packed, 'wb', compresslevel=4) as out:
            try:
                reads, digest = content(iter(lambda: process.stdout.read(CHUNK * 16), b''), out)
            except ValueError as error:  # a stream that ends early: fasterq-dump's own message says why
                reads, digest, broken = 0, None, error
            else:
                broken = None
        errors = messages.read_text(errors='replace').strip()
        if process.returncode:
            raise ValueError('{}: fasterq-dump failed ({}): {}'.format(name, process.returncode, errors[-500:]))
        if broken:
            raise ValueError('{}: {}'.format(name, broken))
        if reads != int(run['reads']):
            raise ValueError('{}: {} reads from SRA, {} declared (a paired run gives two per spot; only single-end runs are '
                             'filled from SRA)'.format(name, reads, run['reads']))
        if run.get('content_md5') and digest != run['content_md5']:
            raise ValueError('{}: content md5 {} does not match the declared {}'.format(name, digest, run['content_md5']))
        version = subprocess.run([*shlex.split(fasterq_dump), '--version'], capture_output=True, text=True).stdout.split()
        record = {'run': name, 'source': url, 'sra_md5': found, 'sra_bytes': int(declared['bytes']), 'reads': reads,
                  'content_md5': digest, 'content_md5_declared': bool(run.get('content_md5')),
                  'fasterq_dump': version[-1] if version else None, 'filled': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
        note = work / 'record.json'
        note.write_text(json.dumps(record, indent=1) + '\n')
        for item in (packed, note):
            os.chmod(item, 0o666 & ~UMASK)
        os.replace(packed, path)  # the reads first: a record never stands beside a file it does not describe
        os.replace(note, path.parent / '{}.sra.json'.format(name))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    log('{}: from SRA, {} reads, content md5 {}{}'.format(path, reads, digest, '' if run.get('content_md5') else ' (not declared)'))
    return digest


def runs(spec):
    return spec['runs'] if spec['kind'] == 'pool' else [run for sample in spec['samples'] for run in sample['runs']]


def fill(spec, ext_dir, log=print, fetch=download, source='ena', connections=8, jobs=4, fasterq_dump='fasterq-dump'):
    """Download the declaration's missing files into EXT_DIR. Returns the number of files downloaded. `source`: where
    the runs come from, `ena` (their declared url and md5) or `sra` (sra_run: NCBI's copy on AWS, `connections` range
    requests per run and `jobs` runs at a time); every other file comes from its declared url."""
    genome = spec.get('kind') == 'synthetic' and (spec.get('synthetic') or {}).get('mode') == 'genome'
    if spec.get('kind') not in ('pool', 'biological') and not genome:
        raise ValueError('{}: only pool, biological and genome-wide synthetic datasets are filled; a reduced synthetic dataset is '
                         'generated by prepare from its reference'.format(spec['name']))
    if source not in ('ena', 'sra'):
        raise ValueError('source must be ena or sra, got {!r}'.format(source))
    from_sra = {}
    if source == 'sra' and not genome:
        fastq = Path(ext_dir) / spec['data'] / 'fastq'
        from_sra = {fastq / '{}.fastq.gz'.format(run['run']): run for run in runs(spec)}
    done, waiting = 0, []
    for url, algorithm, digest, path in targets(spec, ext_dir):
        if path.exists():
            log('{}: in place'.format(path))
        elif path in from_sra:
            waiting.append(path)
        else:
            fetch(url, algorithm, digest, path, log)
            done += 1
    if waiting:
        with concurrent.futures.ThreadPoolExecutor(jobs) as pool:  # the work is in downloads, fasterq-dump and zlib
            results = [pool.submit(sra_run, from_sra[path], path, connections, fasterq_dump, 4, log) for path in waiting]
            failed = [str(result.exception()) for result in results if result.exception()]
        done += len(waiting) - len(failed)
        if failed:
            raise ValueError('{} of {} runs not filled from SRA: {}'.format(len(failed), len(waiting), '; '.join(failed)))
    free = shutil.disk_usage(Path(ext_dir) / spec['data']).free
    log('{}: {} file(s) downloaded, {:.0f} GB free'.format(spec['name'], done, free / 1e9))
    return done
