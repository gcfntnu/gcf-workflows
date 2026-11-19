#!/usr/bin/env python
"""Move one miRge3 run's outputs into the project's miRge3 folder (smallrna mirge3_quant, #404).

miRge3 writes each run into <run-root>/miRge.<date>_<time>/. Its own miR.Counts.csv goes straight to --mirtable
(miR.Counts.mirge3.csv): miR.Counts.csv in --output is mirge3_counts' table (#349), so miRge3's table never lands
there, not even when a later step fails. Every other entry replaces the entry of the same name in --output, file or
directory, so an earlier run's entries (e.g. the directory tRFs.samples.tmp, which no rule declares) cannot block a
rerun; the run root is removed afterwards.

Checked before anything is moved: the run root holds exactly one run, the run has miR.Counts.csv, and no entry has a
name another rule owns in --output (--reserved).

Usage: mirge3_collect.py --run-root DIR --output DIR --mirtable CSV [--reserved NAME ...]
"""
import argparse
import os
from pathlib import Path
import shutil
import sys

MIRGE3_TABLE = 'miR.Counts.csv'


def fail(message):
    sys.exit('mirge3_collect.py: ' + message)


def remove(path):
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def collect(run_root, output, mirtable, reserved=()):
    """Move the single run in run_root into output (miRge3's table to mirtable); returns the moved names."""
    runs = sorted(run_root.glob('miRge.*'))
    if len(runs) != 1:
        fail('expected one miRge3 run in {}, found {}'.format(run_root, [p.name for p in runs] or 'none'))
    run = runs[0]
    entries = sorted(run.iterdir())
    if not (run / MIRGE3_TABLE).is_file():
        fail('miRge3 wrote no {} in {}'.format(MIRGE3_TABLE, run))
    clash = sorted(p.name for p in entries if p.name in set(reserved) | {mirtable.name})
    if clash:
        fail('miRge3 wrote {} in {}, owned by other rules in {}'.format(', '.join(clash), run, output))
    os.replace(run / MIRGE3_TABLE, mirtable)
    moved = []
    for entry in entries:
        if entry.name == MIRGE3_TABLE:
            continue
        target = output / entry.name
        remove(target)
        os.replace(entry, target)
        moved.append(entry.name)
    shutil.rmtree(run_root)
    return moved


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--run-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--mirtable', required=True, type=Path)
    parser.add_argument('--reserved', nargs='*', default=[])
    args = parser.parse_args()
    collect(args.run_root, args.output, args.mirtable, args.reserved)


if __name__ == '__main__':
    main()
