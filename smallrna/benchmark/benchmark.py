#!/usr/bin/env python
"""Benchmark small-RNA quantification methods on datasets with a known truth.

The benchmark runs the workflow's own `quant_all` routes; it never calls the tools itself. Three steps:

1. `prepare` builds a smallrna project from a dataset declaration (`datasets/<name>.yaml`): configuration, sample
   table and reads (generated for synthetic datasets), and the truth table, kept in `<project>/benchmark/`.
2. Run `snakemake quant_all` in the project as for any project on that machine (containers, cores, resources).
3. `evaluate` turns each method's output into one common table (adapters.py), computes the same metrics for every
   method (metrics.py) and writes `<project>/benchmark/report.md` with the tables behind it.

References and public datasets live under `EXT_DIR` (`ext_dir` in the project configuration); a reduced reference for
mirflowz is made with reduce_reference.py and named in the declaration.

Usage: benchmark.py prepare --dataset YAML --project DIR --ext-dir DIR
       benchmark.py evaluate --project DIR [--mirflowz-gff GFF3]
"""
import argparse
import csv
from pathlib import Path
import shutil
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapters  # noqa: E402
import metrics  # noqa: E402
import synthetic  # noqa: E402

PREFIX = {'homo_sapiens': 'hsa-', 'mus_musculus': 'mmu-', 'rattus_norvegicus': 'rno-'}
EXTENSION = 6


def fail(message):
    sys.exit('benchmark.py: ' + message)


def load_dataset(path):
    spec = yaml.safe_load(Path(path).read_text())
    for key in ('name', 'kind', 'organism', 'libprepkit', 'methods'):
        if key not in spec:
            fail('{}: missing {!r}'.format(path, key))
    if spec['organism'] not in PREFIX:
        fail('{}: organism must be one of {}'.format(path, ', '.join(PREFIX)))
    if spec['kind'] != 'synthetic':
        fail('{}: only synthetic datasets are supported so far, got kind {!r}'.format(path, spec['kind']))
    return spec


def reference_dir(spec, ext_dir):
    return Path(ext_dir) / spec['reference']['dir'] if spec.get('reference') else None


def prepare(args):
    spec = load_dataset(args.dataset)
    project, ext_dir = Path(args.project), Path(args.ext_dir).resolve()
    if project.exists():
        fail('{} exists; choose a new project directory'.format(project))
    ref = reference_dir(spec, ext_dir)
    if ref is None or not (ref / 'mature.fa').is_file():
        fail('reference {} lacks mature.fa; make it with reduce_reference.py'.format(ref))
    samples = ['S{:02d}'.format(i) for i in range(1, spec['synthetic'].get('samples', 2) + 1)]
    fastq = lambda sample: project / 'data/raw/fastq' / spec['name'] / '{}_R1.fastq.gz'.format(sample)
    truth = synthetic.generate(spec['synthetic'], ref / 'mature.fa', ref / 'mirna.gff3', PREFIX[spec['organism']], samples,
                               spec['adapter'], int(spec.get('umi', 0)), fastq)
    read_length = max(len(row[4]) for row in truth) + 2 * int(spec.get('umi', 0)) + len(spec['adapter'])

    config = {'workflow': 'smallrna', 'project_id': ['BENCH-0000-001'], 'organism': spec['organism'],
              'libprepkit': spec['libprepkit'], 'read_geometry': [read_length], 'machine': '', 'ext_dir': str(ext_dir),
              'fastq_dir': 'data/raw/fastq',
              'quant': {'method': ','.join(spec['methods']), 'batch': {'method': 'skip'}},
              'samples': {s: {'R1': '{}/{}_R1.fastq.gz'.format(spec['name'], s), 'R2': '', 'paired_end': 0, 'Sample_ID': s}
                          for s in samples}}
    if 'mirflowz' in spec['methods']:
        config['db'] = {'mirflowz': {'custom': {'name': spec['reference']['name'], 'genome': str(ref / 'genome.fa.gz'),
                                                'gtf': str(ref / 'genes.gtf.gz'), 'mirna_gff': str(ref / 'mirna.gff3')}}}
    config.update(spec.get('config') or {})
    (project / 'pep').mkdir(parents=True, exist_ok=True)
    (project / 'config.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
    (project / 'pep/pep_config.yaml').write_text('pep_version: 2.0.0\nsample_table: sample_table.csv\n')
    (project / 'pep/sample_table.csv').write_text('sample_name,R1,paired_end\n' + ''.join(
        '{},{},0\n'.format(s, config['samples'][s]['R1']) for s in samples))
    workflow = Path(__file__).resolve().parents[1] / 'smallrna.smk'
    (project / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\ninclude: {!r}\n'.format(str(workflow)))
    bench = project / 'benchmark'
    bench.mkdir()
    shutil.copy(args.dataset, bench / 'dataset.yaml')
    synthetic.write_truth(bench / 'truth.tsv', truth)
    print('{}: {} samples, {} reads each; run `snakemake quant_all` there, then `benchmark.py evaluate --project {}`'.format(
        project, len(samples), sum(row[5] for row in truth if row[0] == samples[0]), project))


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, header, rows):
    Path(path).write_text('\t'.join(header) + '\n' + ''.join('\t'.join(str(v) for v in row) + '\n' for row in rows))


def evaluate(args):
    project = Path(args.project)
    bench = project / 'benchmark'
    spec = load_dataset(bench / 'dataset.yaml')
    config = yaml.safe_load((project / 'config.yaml').read_text())
    ext_dir = Path(config['ext_dir'])
    prefix = PREFIX[spec['organism']]
    mature = synthetic.read_mature(reference_dir(spec, ext_dir) / 'mature.fa', prefix)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    gff = args.mirflowz_gff
    if gff is None and 'mirflowz' in spec['methods']:
        gff = ext_dir / 'mirflowz' / 'custom-{}'.format(spec['reference']['name']) / spec['organism'] / 'prepared' / \
            'extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION)
    common = adapters.common_table(spec['methods'], project / 'data/tmp/smallrna/quant', name_to_id, prefix, gff)
    write_tsv(bench / 'common.tsv', adapters.COMMON_COLUMNS, common)

    truth = {}
    for row in read_tsv(bench / 'truth.tsv'):
        if row['mature_id']:
            key = (row['sample'], row['mature_id'])
            truth[key] = truth.get(key, 0) + int(row['count'])
    summary = metrics.count_metrics(common, truth)
    table = metrics.per_mirna(common, truth)
    methods = list(spec['methods'])
    names = {}
    for _, _, name, mature_id, _ in common:
        names.setdefault(mature_id, name)
    names.update({mature_id: name for name, mature_id in name_to_id.items() if any(k[1] == mature_id for k in truth)})
    samples = sorted({k[0] for k in truth})
    per = []
    for mature_id in sorted({k[2] for k in table if k[2]}, key=lambda i: names.get(i, i)):
        t = truth.get((samples[0], mature_id), 0)
        per.append([names.get(mature_id, ''), mature_id, t] + [table.get((m, samples[0], mature_id), (t, 0))[1] for m in methods])
    write_tsv(bench / 'per_mirna.tsv', ['mature_name', 'mature_id', 'truth'] + methods, per)
    keys = ['recovered', 'exact', 'misassigned', 'unplaced']
    write_tsv(bench / 'metrics.tsv', ['metric'] + methods, [[k] + [summary.get(m, {}).get(k, '') for m in methods] for k in keys])

    lines = ['# {}'.format(spec['name']), '', '{} ({}), {} samples; truth per sample.'.format(
        spec['organism'], spec['libprepkit'], len(samples)), '', '## Metrics', '',
        '| metric | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 1) + '|']
    lines += ['| {} | '.format(k) + ' | '.join(str(summary.get(m, {}).get(k, '')) for m in methods) + ' |' for k in keys]
    lines += ['', '## Counts per miRNA (sample {})'.format(samples[0]), '',
              '| miRNA | truth | ' + ' | '.join(methods) + ' |', '|---' * (len(methods) + 2) + '|']
    lines += ['| {} | {} | '.format(row[0], row[2]) + ' | '.join('{:g}'.format(v) for v in row[3:]) + ' |'
              for row in per if row[2] or any(row[3:])]
    unplaced = sorted({(row[0], row[2]) for row in common if not row[3]})
    if unplaced:
        lines += ['', '## Names without a miRBase accession', ''] + ['- {}: `{}`'.format(m, n) for m, n in unplaced]
    (bench / 'report.md').write_text('\n'.join(lines) + '\n')
    print(bench / 'report.md')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--dataset', required=True)
    p.add_argument('--project', required=True)
    p.add_argument('--ext-dir', required=True, help='EXT_DIR holding the references named in the declaration')
    e = sub.add_parser('evaluate')
    e.add_argument('--project', required=True)
    e.add_argument('--mirflowz-gff', type=Path, help="mirflowz's extended miRNA annotation; default from the declaration")
    args = parser.parse_args()
    prepare(args) if args.command == 'prepare' else evaluate(args)


if __name__ == '__main__':
    main()
