#!/usr/bin/env python
"""Fail a TEsmall run that reported success but lost a step (#419).

TEsmall 2.1.0 does not check its own steps: when `sort` cannot write its temporary files, or a Bowtie step cannot open
an index, it goes on and ends with "TEsmall run completed", and the count table silently lacks every annotated class
(only the 3' tRF counts remain). This check, run in the sample's work directory after TEsmall, stops the job when:

- TEsmall's log reports anything failing: any `sort:` diagnostic (a `write failed` or `read failed` can leave a
  partial annotation), an error, a failure, a full disk or quota, a killed process, a Python traceback or an htslib
  error (`[E::`). A successful TEsmall log contains none of these words; TEsmall does not pass its subprocesses' exit
  status on, so its log is the only record of a failed step;
- one of its Bowtie logs (`<sample>.{rRNA,genome,tRNA}.log`) is missing, has a fatal diagnostic (a line starting with
  `Error`, a missing index, an exception or abort), or lacks Bowtie's `# reads processed` summary (a normal log also
  says how many reads "failed to align", which is not an error);
- the genome Bowtie log reports aligned reads, but `<sample>.aligned.rinfo` or the annotation (`<sample>.anno`) is
  missing or has no rows.

Usage: tesmall_check.py --work DIR --sample NAME --log TESMALL_LOG
"""
import argparse
from pathlib import Path
import re
import sys

LOG_FAILURES = re.compile(r'^sort: |error|fail|quota|no space left|killed|segmentation fault|traceback|\[E::', re.I)
BOWTIE_FAILURES = re.compile(r'^\s*(Error|Could not locate a Bowtie index|terminate called|Exception|Aborted)')
BOWTIE_STEPS = ('rRNA', 'genome', 'tRNA')
GENOME_ALIGNED = re.compile(r'^# reads with at least one (?:reported )?alignment: (\d+)', re.M)


def rows(path):
    """Data rows of a file with a header line (0 for a header-only file, None for a missing one)."""
    if not path.is_file():
        return None
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
        bad = [line.strip() for line in content.splitlines() if BOWTIE_FAILURES.search(line)]
        if bad:
            found.append('{}: {}'.format(path.name, bad[0]))
        elif '# reads processed' not in content:
            found.append('{}: no "# reads processed" summary'.format(path.name))
        elif step == 'genome':
            match = GENOME_ALIGNED.search(content)
            aligned = int(match.group(1)) if match else None
            if aligned is None:
                found.append('{}: no aligned-read count'.format(path.name))
            elif aligned:
                for name in ('{}.aligned.rinfo'.format(sample), '{}.anno'.format(sample)):
                    n = rows(work / name)
                    if not n:
                        found.append('{} reads aligned to the genome but {} is {}'.format(aligned, name, 'missing' if n is None else 'empty'))
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
