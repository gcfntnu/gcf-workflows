#!/usr/bin/env python
"""Make a reduced reference for the small-RNA benchmark: a few chromosomes of a genome, with matching annotations.

Writes to `--output`: `genome.fa.gz` (the listed chromosomes), `genes.gtf.gz` and `mirna.gff3` (their features only;
GFF header lines kept, e.g. `genome-build-id`), `mature.fa` (the miRBase mature sequences of `--species`) and
`sources.tsv` (each input with its SHA-256). Inputs are paths or URLs, plain or gzipped. GFF chromosome names may be
UCSC style (`chr19`) while the genome uses Ensembl names (`19`); both forms of each listed chromosome are kept.

mirflowz uses the reduced reference through `db.mirflowz.custom`; UNITAS and miRge3 match reads by sequence and use
their normal references. A reduced reference makes mirflowz's mapping fit a workstation; it cannot show behaviour that
depends on the rest of the genome (e.g. multimapping to other chromosomes).

Usage: reduce_reference.py --genome FA --gtf GTF --mirna-gff GFF3 --mature FA --species mmu --chromosomes 19 [...]
                           --output DIR
"""
import argparse
import gzip
import hashlib
import io
import os
from pathlib import Path
import shutil
import sys
import tempfile
import urllib.request


def fail(message):
    sys.exit('reduce_reference.py: ' + message)


def fetch(source, work):
    """Local copy of a path or URL; returns (path, sha256)."""
    if '://' in source:
        target = work / Path(source.split('?')[0]).name
        with urllib.request.urlopen(source, timeout=300) as response, open(target, 'wb') as out:
            shutil.copyfileobj(response, out, 1 << 20)
    else:
        target = Path(source)
        if not target.is_file():
            fail('input not found: {}'.format(source))
    digest = hashlib.sha256()
    with open(target, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return target, digest.hexdigest()


def open_text(path):
    with open(path, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    return io.TextIOWrapper(gzip.open(path)) if gzipped else open(path)


def keep_names(chromosomes):
    names = set()
    for name in chromosomes:
        names |= {name, 'chr' + name, name[3:] if name.startswith('chr') else name}
    names |= {'MT', 'chrM'} if {'MT', 'chrM', 'M'} & set(chromosomes) else set()
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('--genome', '--gtf', '--mirna-gff', '--mature', '--species', '--output'):
        parser.add_argument(name, required=True)
    parser.add_argument('--chromosomes', required=True, nargs='+')
    args = parser.parse_args()
    keep = keep_names(args.chromosomes)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix='.reduce-', dir=output.parent) as work:
        work = Path(work)
        staged = work / 'reference'
        staged.mkdir()
        sources = {}
        genome, sources['genome'] = fetch(args.genome, work)
        kept, writing = [], False
        with open_text(genome) as src, gzip.open(staged / 'genome.fa.gz', 'wt', compresslevel=6) as out:
            for line in src:
                if line.startswith('>'):
                    writing = line[1:].split()[0] in keep
                    if writing:
                        kept.append(line[1:].split()[0])
                if writing:
                    out.write(line)
        missing = [c for c in args.chromosomes if not ({c, 'chr' + c, c[3:]} & set(kept))]
        if missing:
            fail('chromosome(s) not in the genome: {}'.format(', '.join(missing)))
        gtf, sources['gtf'] = fetch(args.gtf, work)
        with open_text(gtf) as src, gzip.open(staged / 'genes.gtf.gz', 'wt') as out:
            for line in src:
                if line.startswith('#') or line.split('\t', 1)[0] in keep:
                    out.write(line)
        gff, sources['mirna_gff'] = fetch(args.mirna_gff, work)
        features = 0
        with open_text(gff) as src, open(staged / 'mirna.gff3', 'w') as out:
            for line in src:
                if line.startswith('#'):
                    out.write(line)
                elif line.split('\t', 1)[0] in keep:
                    out.write(line)
                    features += 1
        if not features:
            fail('no miRNA features on {} in {}'.format(', '.join(args.chromosomes), args.mirna_gff))
        mature, sources['mature'] = fetch(args.mature, work)
        with open_text(mature) as src, open(staged / 'mature.fa', 'w') as out:
            keep_record = False
            for line in src:
                if line.startswith('>'):
                    keep_record = line[1:].startswith(args.species + '-')
                if keep_record:
                    out.write(line)
        (staged / 'sources.tsv').write_text('file\tsource\tsha256\n' + ''.join(
            '{}\t{}\t{}\n'.format(key, getattr(args, key), digest) for key, digest in sources.items()) +
            'chromosomes\t{}\t\n'.format(','.join(kept)))
        if output.exists():
            shutil.rmtree(output)
        os.replace(staged, output)


if __name__ == '__main__':
    main()
