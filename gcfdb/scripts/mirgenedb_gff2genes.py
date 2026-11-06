#!/usr/bin/env python
"""Feature info for UNITAS miRNA tables built from a MirGeneDB reference dump.

Reads a MirGeneDB GFF3 (gcf_mirgenedb.gff3 in the dump) and writes one row per mature or star
miRNA (feature `miRNA`), indexed by `mirna_id`, the name UNITAS reports in a MirGeneDB dump: the gene
without its species code and the arm, joined by "-" (Hsa-Let-7-P1b_3p* -> Let-7-P1b-3p). UNITAS gives
reads on a hairpin arm that same name, so the star "*" is not part of it; see `star`.

Columns: mirna_id, seqid, source, feature, start, end, score, strand, phase, id (MirGeneDB ID),
alias (miRBase accession; empty when MirGeneDB lists none), precursor (MirGeneDB precursor ID),
arm (5p or 3p) and star (True for the star arm).

Usage: mirgenedb_gff2genes.py genes.gff3 > genes.tsv
"""
import csv
import re
import sys

ID = re.compile(r'^(?P<code>[A-Za-z]+)-(?P<name>(?P<gene>.+)_(?P<arm>[35]p)(?P<star>\*?))$')
COLUMNS = ['mirna_id', 'seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase',
           'id', 'alias', 'precursor', 'arm', 'star']


def rows(lines):
    seen = set()
    for number, line in enumerate(lines, 1):
        line = line.rstrip('\r\n')
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != 9:
            sys.exit('mirgenedb_gff2genes.py: line {} has {} fields, expected 9'.format(number, len(fields)))
        if fields[2] != 'miRNA':
            continue
        attrs = dict(a.split('=', 1) for a in fields[8].split(';') if a)
        match = ID.match(attrs.get('ID', ''))
        if not match:
            sys.exit('mirgenedb_gff2genes.py: line {}: not a MirGeneDB miRNA ID: {!r}'.format(number, attrs.get('ID')))
        name = '{}-{}'.format(match.group('gene'), match.group('arm'))
        if name in seen:
            sys.exit('mirgenedb_gff2genes.py: duplicate miRNA {}'.format(name))
        seen.add(name)
        yield [name] + fields[:8] + [attrs['ID'], attrs.get('Alias', ''),
                                     '{}-{}_pre'.format(match.group('code'), match.group('gene')),
                                     match.group('arm'), str(bool(match.group('star')))]


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__.strip().splitlines()[-1])
    with open(sys.argv[1]) as fh:
        lines = fh.readlines()
    if not lines or lines[0].strip() != '##gff-version 3':
        sys.exit('mirgenedb_gff2genes.py: {} is not GFF3'.format(sys.argv[1]))
    if not any(re.match(r'^# microRNAs:\s+MirGeneDB v', line) for line in lines):
        sys.exit('mirgenedb_gff2genes.py: {} is not a MirGeneDB GFF'.format(sys.argv[1]))
    out = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    out.writerow(COLUMNS)
    count = 0
    for row in rows(lines):
        out.writerow(row)
        count += 1
    if not count:
        sys.exit('mirgenedb_gff2genes.py: no miRNA features in {}'.format(sys.argv[1]))
