"""Publish and read a completed evaluation's metrics and report together."""
import csv
from contextlib import contextmanager
import fcntl
import hashlib
import io
import json
import os
import stat
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml


REPORT_MARKER = '<!-- benchmark evaluation: see evaluation.json -->'


@contextmanager
def evaluation_lock(bench, *, writer):
    # NFS requires write access on the descriptor for exclusive locks. Keep
    # this inode stable across publications. Legacy readers create no files.
    path = bench / '.evaluation.lock'
    mode = os.O_WRONLY | os.O_CREAT if writer else os.O_RDONLY
    try:
        fd = os.open(path, mode | os.O_NOFOLLOW | os.O_NONBLOCK, 0o666)
    except FileNotFoundError:
        if writer:
            raise
        yield False
        return
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):
            raise ValueError(f'expected a regular evaluation lock file: {path}')
        fcntl.flock(fd, fcntl.LOCK_EX if writer else fcntl.LOCK_SH)
        yield True
    finally:
        os.close(fd)


def reads_processed(project):
    """Reads the project's run filtered, summed over its configured samples: the input of each sample's first cutadapt
    pass, from the JSON report the run wrote (`data/tmp/smallrna/filter/cutadapt_stages/<sample>.adapter.json`, kept
    with `--notemp`). None when a sample has no report (another filter route, or the reports were removed). A report
    without a whole read count raises ValueError."""
    samples = yaml.safe_load((project / 'config.yaml').read_text())['samples']
    total = 0
    for sample in samples:
        path = project / 'data/tmp/smallrna/filter/cutadapt_stages' / f'{sample}.adapter.json'
        if not path.is_file():
            return None
        try:
            count = json.loads(path.read_text())['read_counts']['input']
        except (ValueError, KeyError, TypeError) as error:
            raise ValueError(f'{path}: no read_counts.input in the cutadapt report ({type(error).__name__}: {error})') from None
        if type(count) is not int or count < 0:
            raise ValueError(f'{path}: read_counts.input is {count!r}, not a read count')
        total += count
    return total


# Tables published with some evaluations (bootstrap intervals, #413; runtime and memory, #386); listed with their
# SHA-256 in evaluation.json.
EXTRA_FILES = ('metrics_ci.tsv', 'differences.tsv', 'routes.tsv', 'resources.tsv')


def publish(bench, header, rows, lines, reads_processed=None, extra=None, bootstrap=None):
    """Publish metrics.tsv and report.md with evaluation.json, the completion record: their SHA-256 and the reads the
    evaluated run processed (`reads_processed`, as read when evaluating; null when unknown). `extra`: {name: text} of
    EXTRA_FILES published with them (an earlier evaluation's EXTRA_FILES not in it are removed), and `bootstrap`
    ({'resamples': n, 'seed': s}) recorded with them."""
    table = '\n'.join('\t'.join(map(str, row)) for row in [header, *rows]) + '\n'
    files = {'metrics.tsv': table.encode(), 'report.md': ('\n'.join([*lines, '', REPORT_MARKER]) + '\n').encode()}
    extra = extra or {}
    unknown = sorted(set(extra) - set(EXTRA_FILES))
    if unknown:
        raise ValueError('not an evaluation file: {}'.format(', '.join(unknown)))
    files.update({name: text.encode() for name, text in extra.items()})
    record = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    record['reads_processed'] = reads_processed
    if bootstrap is not None:
        record['bootstrap'] = bootstrap
    files['evaluation.json'] = (json.dumps(record, sort_keys=True) + '\n').encode()
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
            for name in set(EXTRA_FILES) - set(files):
                (bench / name).unlink(missing_ok=True)


def evaluation_files(bench):
    snapshot = {}
    for name in ('metrics.tsv', 'report.md', 'evaluation.json', *EXTRA_FILES):
        try:
            snapshot[name] = (bench / name).read_bytes()
        except FileNotFoundError:
            pass
    return snapshot


def read_metrics(bench):
    return read_evaluation(bench)[0]


def read_evaluation(bench):
    """(metrics table, completion record, intervals) from one locked snapshot; the record is {} for an older
    evaluation without evaluation.json, and lacks `reads_processed` for one published before it was recorded.
    intervals: the rows of metrics_ci.tsv as dicts when the record lists it (checksum checked), else []."""
    message = f'{bench}: incomplete or inconsistent evaluation; run `benchmark.py evaluate --project {bench.parent}`'
    with evaluation_lock(bench, writer=False) as locked:
        files = evaluation_files(bench)
    # A first publisher may have installed the lock during an unlocked legacy
    # read. Discard that snapshot and wait for its completed publication.
    if not locked and (bench / '.evaluation.lock').exists():
        with evaluation_lock(bench, writer=False):
            files = evaluation_files(bench)
    if 'metrics.tsv' not in files and 'evaluation.json' not in files:
        raise ValueError(f'{bench.parent}: no benchmark/metrics.tsv; run `benchmark.py evaluate --project {bench.parent}`')
    if 'metrics.tsv' not in files or 'report.md' not in files:
        raise ValueError(message)
    data, report = files['metrics.tsv'], files['report.md']
    table = list(csv.reader(io.StringIO(data.decode()), delimiter='\t'))
    if not table or table[0][:1] != ['metric'] or len(table) < 2 or any(len(row) != len(table[0]) for row in table):
        raise ValueError(message)
    record = {}
    if 'evaluation.json' in files:
        try:
            record = json.loads(files['evaluation.json'])
        except ValueError:
            raise ValueError(message) from None
        digests = {'metrics.tsv': hashlib.sha256(data).hexdigest(), 'report.md': hashlib.sha256(report).hexdigest()}
        if not isinstance(record, dict) or set(record) - {'metrics.tsv', 'report.md', 'reads_processed', 'bootstrap', *EXTRA_FILES} or \
                {k: record.get(k) for k in digests} != digests:
            raise ValueError(message)
        for name in EXTRA_FILES:
            if name in record and (name not in files or hashlib.sha256(files[name]).hexdigest() != record[name]):
                raise ValueError(message)
        count = record.get('reads_processed')
        if count is not None and (type(count) is not int or count < 0):
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
    intervals = []
    headers = {'metrics_ci.tsv': ['metric', 'row', 'value', 'low', 'high'],
               'differences.tsv': ['metric', 'row_a', 'row_b', 'difference', 'low', 'high']}
    for name, header in headers.items():
        if name in record:
            reader = csv.DictReader(io.StringIO(files[name].decode()), delimiter='\t')
            rows = list(reader)
            # a long row has its extra fields under None, a short one None values
            if reader.fieldnames != header or any(None in row or None in row.values() for row in rows):
                raise ValueError(message)
            if name == 'metrics_ci.tsv':
                intervals = rows
    return table, record, intervals
