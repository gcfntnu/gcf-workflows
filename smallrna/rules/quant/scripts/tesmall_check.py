#!/usr/bin/env python
"""Fail a TEsmall run that reported success but lost a step (#419).

TEsmall 2.1.0 does not check its own steps: when `sort` cannot write its temporary files, or a Bowtie step cannot open
an index, it goes on and ends with "TEsmall run completed", and the count table silently lacks every annotated class
(only the 3' tRF counts remain). This check, run in the sample's work directory after TEsmall, stops the job when:

- TEsmall's log reports a failed step: `write failed`, `Disk quota exceeded`, `No space left`, a Python traceback, or an
  htslib error (`[E::`);
- one of its Bowtie logs (`<sample>.{rRNA,genome,tRNA}.log`) is missing, reports an error or a missing index, or lacks
  Bowtie's `# reads processed` summary (a normal log also says how many reads "failed to align", which is not an error);
- reads aligned (`<sample>.aligned.rinfo` has rows) but the annotation (`<sample>.anno`) has none.

Usage: tesmall_check.py --work DIR --sample NAME --log TESMALL_LOG
"""
import argparse
from pathlib import Path
import re
import sys

LOG_FAILURES = re.compile(r'write failed|Disk quota exceeded|No space left|Traceback \(most recent call last\)|\[E::')
BOWTIE_FAILURES = re.compile(r'Could not locate a Bowtie index|Error|error:', re.I)
BOWTIE_STEPS = ('rRNA', 'genome', 'tRNA')


def rows(path):
    """Data rows of a file with a header line (0 for a missing or header-only file)."""
    if not path.is_file():
        return 0
    with open(path) as fh:
        return max(sum(1 for line in fh if line.strip()) - 1, 0)


def problems(work, sample, log):
    found = []
    text = log.read_text(errors='replace') if log.is_file() else ''
    if not text:
        found.append('{}: missing or empty'.format(log))
    for line in text.splitlines():
        if LOG_FAILURES.search(line):
            found.append('{}: {}'.format(log.name, line.strip()))
    for step in BOWTIE_STEPS:
        path = work / '{}.{}.log'.format(sample, step)
        if not path.is_file():
            found.append('{}: missing'.format(path.name))
            continue
        content = path.read_text(errors='replace')
        bad = [line.strip() for line in content.splitlines() if BOWTIE_FAILURES.search(line) and 'failed to align' not in line]
        if bad:
            found.append('{}: {}'.format(path.name, bad[0]))
        elif '# reads processed' not in content:
            found.append('{}: no "# reads processed" summary'.format(path.name))
    aligned, annotated = rows(work / '{}.aligned.rinfo'.format(sample)), rows(work / '{}.anno'.format(sample))
    if aligned and not annotated:
        found.append('{} reads aligned but {}.anno has no annotation'.format(aligned, sample))
    return found


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--work', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--log', required=True)
    args = parser.parse_args(argv)
    found = problems(Path(args.work), args.sample, Path(args.log))
    if found:
        sys.exit('tesmall_check.py: TEsmall run of {} is incomplete:\n  {}'.format(args.sample, '\n  '.join(found)))


if __name__ == '__main__':
    main()
