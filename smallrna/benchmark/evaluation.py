"""Publish and read a completed evaluation's metrics and report together."""
import csv
from contextlib import contextmanager
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory


@contextmanager
def directory_lock(directory, mode):
    fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(fd, mode)
        yield
    finally:
        os.close(fd)


def publish(bench, header, rows, lines):
    table = '\n'.join('\t'.join(map(str, row)) for row in [header, *rows]) + '\n'
    files = {'metrics.tsv': table.encode(), 'report.md': ('\n'.join(lines) + '\n').encode()}
    digests = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    files['evaluation.json'] = (json.dumps(digests, sort_keys=True) + '\n').encode()
    with TemporaryDirectory(dir=bench) as directory:
        stage = Path(directory)
        for name, data in files.items():
            (stage / name).write_bytes(data)
        with directory_lock(bench, fcntl.LOCK_EX):
            for name in files:  # The completion record is installed last.
                (stage / name).replace(bench / name)


def read_metrics(bench):
    message = f'{bench}: incomplete or inconsistent evaluation; run `benchmark.py evaluate --project {bench.parent}`'
    with directory_lock(bench, fcntl.LOCK_SH):
        try:
            data = (bench / 'metrics.tsv').read_bytes()
            report = (bench / 'report.md').read_bytes()
        except FileNotFoundError as error:
            raise ValueError(message) from error
        marker = bench / 'evaluation.json'
        table = list(csv.reader(io.StringIO(data.decode()), delimiter='\t'))
        if not table or table[0][:1] != ['metric'] or len(table) < 2 or any(len(row) != len(table[0]) for row in table):
            raise ValueError(message)
        if marker.exists():
            if json.loads(marker.read_text()) != {'metrics.tsv': hashlib.sha256(data).hexdigest(),
                                                  'report.md': hashlib.sha256(report).hexdigest()}:
                raise ValueError(message)
        else:
            # Existing releases embed the entire metrics table in the report. Check
            # that exact table before accepting a pair without a completion record.
            lines = ['| ' + ' | '.join(row) + ' |' for row in table]
            lines.insert(1, '|---' * len(table[0]) + '|')
            report_lines = report.decode().splitlines()
            matching = [i for i, line in enumerate(report_lines) if line == lines[0]]
            if not any(report_lines[i:i + len(lines)] == lines
                       and (i + len(lines) == len(report_lines) or not report_lines[i + len(lines)].startswith('|'))
                       for i in matching):
                raise ValueError(message)
    return table
