#!/usr/bin/env python
"""Prepare a mirflowz reference with mirflowz's own preparation rules, run unmodified inside the mirflowz image.

Two parts, one subcommand each (#440):
- `genome`: from the genome and GTF, the files mirflowz maps with (processed genome, transcriptome, their segemehl
  indexes, exons, the genome's SAM header and chromosome sizes), shared by every miRNA annotation on the reference;
- `annotation`: from the miRNA GFF and the chromosome map, with the genome part's chromosome sizes, mirflowz's
  extended miRNA and pri-miRNA annotations, one per corrected miRBase GFF3 (gcfdb/mirflowz.db).

mirflowz's top-level Snakefile deletes its intermediates directory, where the prepared reference lives, after every
successful run. This script therefore writes its own small Snakefile that includes only
`$MIRFLOWZ_DIR/workflow/rules/prepare.smk`, runs the part's files as targets in a work directory next to `--output`,
checks that every file mirflowz's mapping and quantification read is present, and moves them into `--output`. The
annotation part runs with `--allowed-rules map_chr_names extend_mirs_annotations` and the chromosome sizes linked into
the work directory, so nothing of the genome part is built again.

mirflowz reads the genome and GTF with `zcat` and the GFF as plain text: plain genome and GTF files are compressed in
the work directory. The GFF is always rewritten there, in the miRBase layout mirflowz expects: the GFFs miRBase now
serves (e.g. 23.1) repeat a mature miRNA's `ID` at every locus, where older files gave later occurrences a replica
suffix; mirflowz's annotation step (gffutils) stops on repeated IDs. The first occurrence keeps its `ID`, later ones
get `_1`, `_2`, ..., in file order; `Alias` and all other fields are unchanged. The work directory is removed
afterwards; the inputs are not changed.

mirflowz declares threads and memory (GB) as custom resources, which `--cores` does not limit: each rule's values are
capped at `--threads` and `--mem-mb` (the outer job's allocation, in whole GB), and `--resources` limits their totals.

Writes a two-column TSV (field, value) per part: releases, assembly, inputs with SHA-256, and the mirflowz version.

Usage: mirflowz_prepare.py genome --genome FA[.gz] --gtf GTF[.gz] --mirbase-release R --ensembl-release R --assembly A
                                  --threads N --mem-mb N --output DIR --record TSV
       mirflowz_prepare.py annotation --mirna-gff GFF3[.gz] --chr-map TSV --chr-size TSV --mirbase-release R
                                      --output DIR --record TSV
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
# Resources (threads, mem in GB) that mirflowz 0.10.0's preparation rules declare.
RESOURCES = {'generate_segemehl_index_genome': (8, 60), 'generate_segemehl_index_transcriptome': (8, 10)}
# Files of the genome part that mirflowz 0.10.0's mapping and quantification rules read.
GENOME = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
          'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv')
# Files of the annotation part they read, and the preparation rules that write them from the GFF.
ANNOTATION = ('extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION), 'extended_primir_annotation_{}_nt.gff3'.format(EXTENSION))
ANNOTATION_RULES = ('map_chr_names', 'extend_mirs_annotations')


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


def mem_budget_gb(mem_mb):
    """The memory budget in whole GB, the unit of mirflowz's own declarations; at least 1 GB."""
    if mem_mb < 1024:
        fail('memory budget {} MB is below 1 GB; mirflowz declares memory in GB'.format(mem_mb))
    return mem_mb // 1024


def resource_args(resources, threads, mem_gb):
    """Snakemake options that cap each rule's declared threads and memory and limit their totals."""
    set_threads, set_resources = [], []
    for rule, (rule_threads, rule_mem) in resources.items():
        capped = min(rule_threads, threads)
        set_threads.append('{}={}'.format(rule, capped))
        set_resources.append('{}:threads={}'.format(rule, capped))
        if rule_mem is not None:
            set_resources.append('{}:mem={}'.format(rule, min(rule_mem, mem_gb)))
    return ['--set-threads', *set_threads, '--set-resources', *set_resources,
            '--resources', 'threads={}'.format(threads), 'mem={}'.format(mem_gb)]


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


def prepare(work, config, targets, threads, options, output, mirflowz_dir):
    """Run prepare.smk for `targets` (names in the intermediates directory) in `work`, check them and move the
    intermediates directory to `output`."""
    write_yaml(work / 'config.yaml', config)
    (work / 'Snakefile').write_text('include: {!r}\n'.format(str(mirflowz_dir / 'workflow/rules/prepare.smk')))
    intermediates = Path(config['intermediates_dir'])
    env = dict(os.environ, XDG_CACHE_HOME=str(work / '.cache'))
    result = subprocess.run(['snakemake', '-s', str(work / 'Snakefile'), '--configfile', str(work / 'config.yaml'),
                             '--directory', str(work), '--cores', str(threads), '--nolock', '--workflow-profile', 'none',
                             *[str(intermediates / name) for name in targets],
                             # The work directory is removed afterwards, so no metadata is kept; Snakemake then also
                             # skips reading the conda files that mirflowz's rules name (conda is not used).
                             '--drop-metadata',
                             # Last: these options take several values each, so a target after them would be read as one.
                             *options],
                            capture_output=True, text=True, env=env)
    if result.returncode != 0:
        logs = sorted((work / 'logs/local').glob('*.log'), key=os.path.getmtime)
        tail = ''.join('\n--- {}\n{}'.format(p.name, p.read_text()[-1500:]) for p in logs[-3:])
        fail('mirflowz preparation failed:\n{}{}'.format(result.stderr[-3000:], tail))
    missing = [name for name in targets if not (intermediates / name).is_file()]
    if missing:
        fail('mirflowz preparation wrote no {}'.format(', '.join(missing)))
    for path in intermediates.iterdir():  # the inputs linked in (the annotation part's chromosome sizes)
        if path.is_symlink():
            path.unlink()
    if output.exists():
        shutil.rmtree(output)
    os.replace(intermediates, output)


def write_record(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name('.{}.{}.part'.format(path.name, os.getpid()))
    partial.write_text('field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))
    os.replace(partial, path)


def genome_part(args, mirflowz_dir):
    """The genome part: prepare.smk's genome and transcriptome files from the genome and GTF."""
    for path in (args.genome, args.gtf):
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
        (work / 'empty.bed').touch()
        (work / 'samples.tsv').write_text('sample\tsample_file\tadapter\tformat\n')
        config = mirflowz_config(work, inputs['genome'], inputs['gtf'], 'NA', 'NA', work / 'prepared', work / 'samples.tsv', mirflowz_dir)
        prepare(work, config, GENOME, args.threads, resource_args(RESOURCES, args.threads, mem_budget_gb(args.mem_mb)), args.output, mirflowz_dir)
    write_record(args.record, {
        'mirbase_release': args.mirbase_release, 'ensembl_release': args.ensembl_release, 'assembly': args.assembly,
        'genome': str(args.genome.resolve()), 'genome_sha256': sha256(args.genome), 'gtf': str(args.gtf.resolve()),
        'gtf_sha256': sha256(args.gtf), 'mirflowz_dir': str(mirflowz_dir), 'prepared_files': ','.join(GENOME)})


def annotation_part(args, mirflowz_dir):
    """The annotation part: prepare.smk's extended annotations from the miRNA GFF, the chromosome map and the genome
    part's chromosome sizes, with only the rules that read the GFF allowed."""
    for path in (args.mirna_gff, args.chr_map, args.chr_size):
        if not path.is_file():
            fail('input not found: {}'.format(path))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.mirflowz-annotation-', dir=args.output.parent) as work:
        work = Path(work).resolve()
        (work / 'annotation').mkdir()
        os.symlink(args.chr_size.resolve(), work / 'annotation' / 'chr_size.tsv')
        gff = work / 'mirna.gff3'
        renamed = write_mirflowz_gff(args.mirna_gff, gff)
        (work / 'empty.bed').touch()
        (work / 'samples.tsv').write_text('sample\tsample_file\tadapter\tformat\n')
        config = mirflowz_config(work, 'NA', 'NA', gff, args.chr_map.resolve(), work / 'annotation', work / 'samples.tsv', mirflowz_dir)
        prepare(work, config, ANNOTATION, 1, ['--allowed-rules', *ANNOTATION_RULES], args.output, mirflowz_dir)
    write_record(args.record, {
        'mirbase_release': args.mirbase_release, 'mirna_gff': str(args.mirna_gff.resolve()), 'mirna_gff_sha256': sha256(args.mirna_gff),
        'gff_ids_suffixed': str(renamed), 'chr_map_sha256': sha256(args.chr_map), 'chr_size_sha256': sha256(args.chr_size),
        'mirflowz_dir': str(mirflowz_dir), 'annotation_files': ','.join(ANNOTATION), 'extension': str(EXTENSION)})


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parts = parser.add_subparsers(dest='part', required=True)
    genome = parts.add_parser('genome', help='the genome and transcriptome files, from the genome and GTF')
    for name in ('--genome', '--gtf'):
        genome.add_argument(name, required=True, type=Path)
    for name in ('--ensembl-release', '--assembly'):
        genome.add_argument(name, required=True)
    genome.add_argument('--threads', type=int, default=1)
    genome.add_argument('--mem-mb', type=int, required=True, help="the outer job's memory allocation (MB)")
    annotation = parts.add_parser('annotation', help='the extended miRNA annotations, from the miRNA GFF')
    for name in ('--mirna-gff', '--chr-map', '--chr-size'):
        annotation.add_argument(name, required=True, type=Path)
    for part in (genome, annotation):
        part.add_argument('--mirbase-release', required=True)
        for name in ('--output', '--record'):
            part.add_argument(name, required=True, type=Path)
    args = parser.parse_args()

    mirflowz_dir = Path(os.environ.get('MIRFLOWZ_DIR', ''))
    if not (mirflowz_dir / 'workflow/rules/prepare.smk').is_file():
        fail('MIRFLOWZ_DIR does not hold mirflowz ({}); run in the mirflowz image'.format(mirflowz_dir))
    (genome_part if args.part == 'genome' else annotation_part)(args, mirflowz_dir)


if __name__ == '__main__':
    main()
