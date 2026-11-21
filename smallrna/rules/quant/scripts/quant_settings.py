#!/usr/bin/env python
"""A quantification route's counting settings (#417): quant/<route>/settings.tsv with field, value and source, read by
the benchmark's evaluation. It is written after the route's count table and carries that table's SHA-256
(`table_sha256`), so a record is accepted only for the very table it describes, as quant_reference's is.

The settings come as a JSON list of [field, value, source] (each route's *_settings() in rules/quant/<route>.smk).
`--recorded FILE:KEY,...` adds keys of a configuration file the route itself wrote, in the `key: <JSON value>` lines of
run_mirflowz.py, with the file's name as their source, so they are recorded as the tool ran with them.

Usage: quant_settings.py --table TABLE --settings JSON [--recorded FILE:KEY,KEY...] --output TSV
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def fail(message):
    sys.exit('quant_settings.py: ' + message)


def recorded(spec):
    """[(field, value, source)] of the keys `spec` (FILE:KEY,KEY...) names in a run's `key: <JSON value>` file."""
    path, _, keys = spec.rpartition(':')
    values = {}
    for line in Path(path).read_text().splitlines():
        key, sep, value = line.partition(': ')
        if sep:
            values[key] = json.loads(value)
    missing = [k for k in keys.split(',') if k not in values]
    if missing:
        fail('{} has no {}'.format(path, ', '.join(missing)))
    return [(k, json.dumps(values[k]) if not isinstance(values[k], str) else values[k], Path(path).name) for k in keys.split(',')]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--table', required=True, type=Path, help="the route's count table")
    parser.add_argument('--settings', required=True, help='JSON list of [field, value, source]')
    parser.add_argument('--recorded', help='FILE:KEY,KEY... of a configuration the route wrote')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    rows = [tuple(str(v) for v in row) for row in json.loads(args.settings)]
    if any(len(row) != 3 for row in rows):
        fail('each setting is [field, value, source], got {}'.format(args.settings))
    if args.recorded:
        rows += recorded(args.recorded)
    fields = [row[0] for row in rows]
    if len(set(fields)) != len(fields) or 'table_sha256' in fields:
        fail('fields must be unique and not table_sha256, got {}'.format(', '.join(fields)))
    if any(c in v for row in rows for v in row for c in '\t\n'):
        fail('a setting holds a tab or a line break')
    sha = hashlib.sha256(args.table.read_bytes()).hexdigest()
    text = 'field\tvalue\tsource\n' + ''.join('\t'.join(row) + '\n' for row in rows) + 'table_sha256\t{}\t{}\n'.format(sha, args.table)
    args.output.write_text(text)


if __name__ == '__main__':
    main()
