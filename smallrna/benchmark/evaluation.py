"""Publish and read a completed evaluation's metrics and report together."""
import csv
from contextlib import contextmanager
import fcntl
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory


REPORT_MARKER = '<!-- benchmark evaluation: see evaluation.json -->'


@contextmanager
def evaluation_lock(bench, *, writer):
    # NFS requires write access on the descriptor for exclusive locks. Keep
    # this inode stable across publications. Legacy readers create no files.
    try:
        lock = (bench / '.evaluation.lock').open('ab' if writer else 'rb')
    except FileNotFoundError:
        if writer:
            raise
        yield False
        return
    with lock:
        fcntl.flock(lock, fcntl.LOCK_EX if writer else fcntl.LOCK_SH)
        yield True


def publish(bench, header, rows, lines):
    table = '\n'.join('\t'.join(map(str, row)) for row in [header, *rows]) + '\n'
    files = {'metrics.tsv': table.encode(), 'report.md': ('\n'.join([*lines, '', REPORT_MARKER]) + '\n').encode()}
    digests = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    files['evaluation.json'] = (json.dumps(digests, sort_keys=True) + '\n').encode()
    with TemporaryDirectory(dir=bench) as directory:
        stage = Path(directory)
        for name, data in files.items():
            (stage / name).write_bytes(data)
        pending = stage / 'pending.json'
        pending.write_text('{"status":"publishing"}\n')
        with evaluation_lock(bench, writer=True):
            pending.replace(bench / 'evaluation.json')
            for name in files:  # The completion record is installed last.
                (stage / name).replace(bench / name)


def evaluation_files(bench):
    snapshot = {}
    for name in ('metrics.tsv', 'report.md', 'evaluation.json'):
        try:
            snapshot[name] = (bench / name).read_bytes()
        except FileNotFoundError:
            pass
    return snapshot


def read_metrics(bench):
    message = f'{bench}: incomplete or inconsistent evaluation; run `benchmark.py evaluate --project {bench.parent}`'
    with evaluation_lock(bench, writer=False) as locked:
        files = evaluation_files(bench)
    # A first publisher may have installed the lock during an unlocked legacy
    # read. Discard that snapshot and wait for its completed publication.
    if not locked and (bench / '.evaluation.lock').exists():
        with evaluation_lock(bench, writer=False):
            files = evaluation_files(bench)
    if 'metrics.tsv' not in files or 'report.md' not in files:
        raise ValueError(message)
    data, report = files['metrics.tsv'], files['report.md']
    table = list(csv.reader(io.StringIO(data.decode()), delimiter='\t'))
    if not table or table[0][:1] != ['metric'] or len(table) < 2 or any(len(row) != len(table[0]) for row in table):
        raise ValueError(message)
    if 'evaluation.json' in files:
        if json.loads(files['evaluation.json']) != {'metrics.tsv': hashlib.sha256(data).hexdigest(),
                                                    'report.md': hashlib.sha256(report).hexdigest()}:
            raise ValueError(message)
    else:
        # Historical writers put the primary metrics first. A later table can
        # have the same heading but different metrics (e.g. lane agreement).
        if REPORT_MARKER.encode() in report:
            raise ValueError(message)
        lines = ['| ' + ' | '.join(row) + ' |' for row in table]
        lines.insert(1, '|---' * len(table[0]) + '|')
        report_lines = report.decode().splitlines()
        first = next((i for i, line in enumerate(report_lines) if line.startswith('| metric |')), None)
        if first is None or report_lines[first:first + len(lines)] != lines:
            raise ValueError(message)
        end = first + len(lines)
        if end < len(report_lines) and report_lines[end].startswith('|'):
            raise ValueError(message)
    return table
