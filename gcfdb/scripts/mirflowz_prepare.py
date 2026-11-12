#!/usr/bin/env python
"""Prepare a mirflowz reference with mirflowz's own preparation rules, run unmodified inside the mirflowz image.

mirflowz's top-level Snakefile deletes its intermediates directory, where the prepared reference lives, after every
successful run. This script therefore writes its own small Snakefile that includes only
`$MIRFLOWZ_DIR/workflow/rules/prepare.smk`, runs its `finish_prepare` target in a work directory next to `--output`,
checks that every file mirflowz's mapping and quantification read is present, and moves them into `--output`.

mirflowz reads the genome and GTF with `zcat` and the GFF as plain text: plain inputs are compressed in the work
directory. The GFF is always rewritten there, in the miRBase layout mirflowz expects: the GFFs miRBase now serves (e.g.
23.1) repeat a mature miRNA's `ID` at every locus, where older files gave later occurrences a replica suffix; mirflowz's
annotation step (gffutils) stops on repeated IDs. The first occurrence keeps its `ID`, later ones get `_1`, `_2`, ...,
in file order; `Alias` and all other fields are unchanged. The work directory is removed afterwards; the inputs are
not changed.

Writes a two-column TSV (field, value): releases, assembly, inputs with SHA-256, and the mirflowz version.

Usage: mirflowz_prepare.py --genome FA[.gz] --gtf GTF[.gz] --mirna-gff GFF3[.gz] --chr-map TSV --mirbase-release R
                           --ensembl-release R --assembly A --threads N --output DIR --record TSV
"""
import argparse
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

EXTENSION = 6
# Files of the prepared reference that mirflowz 0.10.0's mapping and quantification rules read.
PREPARED = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
            'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv',
            'extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION), 'extended_primir_annotation_{}_nt.gff3'.format(EXTENSION))


def fail(message):
    sys.exit('mirflowz_prepare.py: ' + message)


def is_gzip(path):
    with open(path, 'rb') as fh:
        return fh.read(2) == b'\x1f\x8b'


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_mirflowz_gff(source, target):
    """Copy the GFF with repeated IDs suffixed _1, _2, ... after their first occurrence; return how many were renamed."""
    seen, renamed = {}, 0
    with (gzip.open(source, 'rt') if is_gzip(source) else open(source)) as src, open(target, 'w') as out:
        for line in src:
            fields = line.rstrip('\n').split('\t')
            if not line.startswith('#') and len(fields) == 9:
                attributes = fields[8].split(';')
                for i, attribute in enumerate(attributes):
                    if attribute.startswith('ID='):
                        identifier = attribute[3:]
                        count = seen.get(identifier, 0)
                        seen[identifier] = count + 1
                        if count:
                            attributes[i] = 'ID={}_{}'.format(identifier, count)
                            renamed += 1
                fields[8] = ';'.join(attributes)
                line = '\t'.join(fields) + '\n'
            out.write(line)
    return renamed


def mirflowz_config(work, genome, gtf, gff, chr_map, intermediates, samples, mirflowz_dir):
    """mirflowz's configuration; every key of its template must be present."""
    return {'samples': str(samples), 'genome_file': str(genome), 'gtf_file': str(gtf), 'mirna_file': str(gff),
            'bed_file': str(work / 'empty.bed'), 'map_chr_file': str(chr_map),
            'output_dir': str(work / 'results'), 'pileups_dir': str(work / 'results/pileups'),
            'intermediates_dir': str(intermediates), 'local_log': str(work / 'logs/local'),
            'cluster_log': str(work / 'logs/cluster'), 'scripts_dir': str(mirflowz_dir / 'scripts') + '/',
            'extension': EXTENSION, 'q_value': 10, 'p_value': 50, 'error_rate': 0.1, 'minimum_length': 15,
            'overlap': 3, 'max_n': 0, 'max_length_reads': 30, 'nh': 100, 'mir_list': ['isomir', 'mirna', 'pri-mir'],
            'lib_dict': {}}


def write_yaml(path, data):
    """Write the flat configuration as YAML (JSON-compatible values)."""
    import json
    path.write_text(''.join('{}: {}\n'.format(key, json.dumps(value)) for key, value in data.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    for name in ('--genome', '--gtf', '--mirna-gff', '--chr-map', '--output', '--record'):
        parser.add_argument(name, required=True, type=Path)
    for name in ('--mirbase-release', '--ensembl-release', '--assembly'):
        parser.add_argument(name, required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    mirflowz_dir = Path(os.environ.get('MIRFLOWZ_DIR', ''))
    if not (mirflowz_dir / 'workflow/rules/prepare.smk').is_file():
        fail('MIRFLOWZ_DIR does not hold mirflowz ({}); run in the mirflowz image'.format(mirflowz_dir))
    for path in (args.genome, args.gtf, args.mirna_gff, args.chr_map):
        if not path.is_file():
            fail('input not found: {}'.format(path))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.mirflowz-prepare-', dir=args.output.parent) as work:
        work = Path(work).resolve()
        inputs = {}
        for key, path in (('genome', args.genome), ('gtf', args.gtf)):
            if is_gzip(path):
                inputs[key] = path.resolve()
            else:
                inputs[key] = work / (key + '.gz')
                with open(path, 'rb') as src, gzip.open(inputs[key], 'wb', compresslevel=1) as dst:
                    shutil.copyfileobj(src, dst, 1 << 20)
        inputs['mirna_gff'] = work / 'mirna.gff3'
        renamed = write_mirflowz_gff(args.mirna_gff, inputs['mirna_gff'])
        (work / 'empty.bed').touch()
        (work / 'samples.tsv').write_text('sample\tsample_file\tadapter\tformat\n')
        write_yaml(work / 'config.yaml', mirflowz_config(work, inputs['genome'], inputs['gtf'], inputs['mirna_gff'],
                                                         args.chr_map.resolve(), work / 'prepared', work / 'samples.tsv', mirflowz_dir))
        (work / 'Snakefile').write_text('include: {!r}\n'.format(str(mirflowz_dir / 'workflow/rules/prepare.smk')))
        env = dict(os.environ, XDG_CACHE_HOME=str(work / '.cache'))
        result = subprocess.run(['snakemake', '-s', str(work / 'Snakefile'), '--configfile', str(work / 'config.yaml'),
                                 '--directory', str(work), '--cores', str(args.threads), '--nolock',
                                 # The work directory is removed afterwards, so no metadata is kept; Snakemake then also
                                 # skips reading the conda files that mirflowz's rules name (conda is not used).
                                 '--drop-metadata', 'finish_prepare'],
                                capture_output=True, text=True, env=env)
        if result.returncode != 0:
            logs = sorted((work / 'logs/local').glob('*.log'), key=os.path.getmtime)
            tail = ''.join('\n--- {}\n{}'.format(p.name, p.read_text()[-1500:]) for p in logs[-3:])
            fail('mirflowz preparation failed:\n{}{}'.format(result.stderr[-3000:], tail))
        missing = [name for name in PREPARED if not (work / 'prepared' / name).is_file()]
        if missing:
            fail('mirflowz preparation wrote no {}'.format(', '.join(missing)))
        if args.output.exists():
            shutil.rmtree(args.output)
        os.replace(work / 'prepared', args.output)

    record = {'mirbase_release': args.mirbase_release, 'ensembl_release': args.ensembl_release, 'assembly': args.assembly,
              'genome': str(args.genome.resolve()), 'genome_sha256': sha256(args.genome),
              'gtf': str(args.gtf.resolve()), 'gtf_sha256': sha256(args.gtf),
              'mirna_gff': str(args.mirna_gff.resolve()), 'mirna_gff_sha256': sha256(args.mirna_gff),
              'gff_ids_suffixed': str(renamed), 'chr_map_sha256': sha256(args.chr_map), 'mirflowz_dir': str(mirflowz_dir),
              'prepared_files': ','.join(PREPARED), 'extension': str(EXTENSION)}
    args.record.parent.mkdir(parents=True, exist_ok=True)
    partial = args.record.with_name('.{}.{}.part'.format(args.record.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, args.record)


if __name__ == '__main__':
    main()
