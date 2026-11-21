#!/usr/bin/env python
"""The miRNAs expected to survive a DROSHA knockout: the list a knockout declaration names (`knockout.expected`, #478).

Canonical miRNAs need DROSHA; a few are made without it. Two sources:
- Kim et al. 2016 (PNAS 113:E1881, doi:10.1073/pnas.1602532113), Fig. 3B: the miRNAs among the 200 most abundant in
  HCT116 cells that keep more than 1% of their reads in the DROSHA knockout (5' capped pre-miRNAs and a mirtron); fixed
  here by miRBase accession (PAPER).
- mirtronDB (version 2, doi:10.6084/m9.figshare.29344775): mirtrons are spliced out of introns and need no DROSHA.
  Every human miRBase mature that lies exactly inside a mirtronDB human precursor is listed, with the PubMed IDs of the
  papers behind those precursors.

The list has `#` lines naming the sources and inputs (with SHA-256), then the columns mature_id, mirna, reason (`5'
capped` or `mirtron`) and source; one row per mature, sorted by name.

mirtronDB's precursor table has no fixed URL; its download page exports it from a form:
  curl -X POST -d 'btn_prec_tabular=19&specie_hidden=19' http://mirtrondb.cp.utfpr.edu.br/functions.php

Usage: knockout_expected.py --precursors TSV --exported YYYY-MM-DD --mature FA --output TSV
"""
import argparse
import csv
import hashlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic  # noqa: E402

COLUMNS = ['mature_id', 'mirna', 'reason', 'source']
PAPER_SOURCE = 'Kim et al. 2016 Fig. 3B'
# (accession, the paper's name, reason); the paper's miR-484-3p is miRBase's hsa-miR-484-5p, whose hairpin the paper finds
# mis-annotated (Fig. S3A). Its miR-1254-5p, from an Alu element, is not in miRBase 23.1 (neither mature nor hairpin).
PAPER = [('MIMAT0030021', 'miR-7706-3p', "5' capped"), ('MIMAT0017994', 'miR-3615-3p', "5' capped"),
         ('MIMAT0000510', 'miR-320a-3p', "5' capped"), ('MIMAT0002174', 'miR-484-3p', "5' capped"),
         ('MIMAT0005792', 'miR-320b-3p', "5' capped"), ('MIMAT0004949', 'miR-877-5p', 'mirtron')]
PAPER_ABSENT = 'hsa-miR-1254'


def fail(message):
    sys.exit('knockout_expected.py: ' + message)


def precursors(path):
    """[(name, PubMed ID, DNA sequence)] of a mirtronDB tabular export of human precursors."""
    rows = [[field.strip() for field in row] for row in csv.reader(open(path, newline=''), delimiter='\t')
            if row and not row[0].startswith('#')]
    if not rows or rows[0][:4] != ['specie', 'id', 'name', 'paper'] or rows[0][-1] != 'sequence':
        fail('{}: not a mirtronDB tabular export (specie, id, name, paper, ..., sequence)'.format(path))
    body = rows[1:]
    if not body or any(len(row) != len(rows[0]) or row[0] != 'H. sapiens' or row[4] != 'precursor' for row in body):
        fail('{}: expected the human precursors only (specie H. sapiens, type precursor)'.format(path))
    return [(row[2], row[3], row[-1].upper().replace('U', 'T')) for row in body]


def expected_rows(mirtrons, mature):
    """The rows of the list: PAPER's miRNAs, and every mature (miRBase {name: (accession, sequence)}) exactly inside a
    mirtron precursor ([(name, PubMed ID, sequence)]); one per accession, sorted by name."""
    by_id = {accession: name for name, (accession, _) in mature.items()}
    missing = [a for a, _, _ in PAPER if a not in by_id]
    if missing:
        fail('{} of Kim et al. 2016 not in the mature sequences; check the release'.format(', '.join(missing)))
    if any(name == PAPER_ABSENT or name.startswith(PAPER_ABSENT + '-') for name in mature):
        fail('{} is in the mature sequences: add Kim et al. 2016\'s miR-1254-5p to PAPER'.format(PAPER_ABSENT))
    rows = {a: [a, by_id[a], reason, [PAPER_SOURCE]] for a, _, reason in PAPER}
    for name, (accession, seq) in mature.items():
        papers = sorted({pmid for _, pmid, precursor in mirtrons if seq in precursor}, key=int)
        if papers:
            row = rows.setdefault(accession, [accession, name, 'mirtron', []])
            if row[2] != 'mirtron':
                fail('{} is {} in Kim et al. 2016, but inside a mirtron precursor'.format(name, row[2]))
            row[3].append('mirtronDB v2 PMID ' + ','.join(papers))
    return [[a, name, reason, '; '.join(sources)] for a, name, reason, sources in sorted(rows.values(), key=lambda r: r[1])]


def read_list(path):
    """{accession: (mirna, reason, source)} of a list written by this script; exits on another layout."""
    lines = [line.rstrip('\n').split('\t') for line in open(path) if line.strip() and not line.startswith('#')]
    if not lines or lines[0] != COLUMNS or any(len(row) != len(COLUMNS) or not row[0] or not row[2] for row in lines[1:]):
        fail('{}: expected the columns {} and a mature_id and reason on every row'.format(path, ', '.join(COLUMNS)))
    if len({row[0] for row in lines[1:]}) != len(lines) - 1:
        fail('{}: a mature_id is listed twice'.format(path))
    return {row[0]: tuple(row[1:]) for row in lines[1:]}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--precursors', required=True, help="mirtronDB's tabular export of the human precursors")
    parser.add_argument('--exported', required=True, help='the date of that export, YYYY-MM-DD')
    parser.add_argument('--mature', required=True, help="miRBase's mature.fa (the release the benchmark uses)")
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    mirtrons = precursors(args.precursors)
    rows = expected_rows(mirtrons, synthetic.read_mature(args.mature, 'hsa-'))
    counts = {reason: sum(row[2] == reason for row in rows) for reason in ("5' capped", 'mirtron')}
    header = [
        '# miRNAs expected to survive a DROSHA knockout (smallrna benchmark, knockout.expected); written by knockout_expected.py.',
        '# Kim et al. 2016 (PNAS 113:E1881, doi:10.1073/pnas.1602532113), Fig. 3B: the miRNAs HCT116 cells keep without DROSHA;',
        "#   the paper's miR-484-3p is miRBase's hsa-miR-484-5p (Fig. S3A); its miR-1254-5p is not in this miRBase release.",
        '# mirtronDB v2 (doi:10.6084/m9.figshare.29344775): {} human precursors, exported {}, sha256'.format(len(mirtrons), args.exported),
        '#   {}; a miRBase mature exactly inside a precursor'.format(sha256(args.precursors)),
        '#   is a mirtron; PMID: the papers behind those precursors.',
        '# miRBase mature.fa sha256 {}; {} rows: {} 5\' capped, {} mirtrons.'.format(
            sha256(args.mature), len(rows), counts["5' capped"], counts['mirtron'])]
    with open(args.output, 'w', newline='') as fh:
        fh.write('\n'.join(header) + '\n')
        csv.writer(fh, delimiter='\t', lineterminator='\n').writerows([COLUMNS] + rows)
    print('{}: {} rows ({} 5\' capped, {} mirtrons)'.format(args.output, len(rows), counts["5' capped"], counts['mirtron']))


if __name__ == '__main__':
    main()
