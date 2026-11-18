#!/usr/bin/env python
"""Benchmark small-RNA quantification methods on datasets with a known truth.

The benchmark runs the workflow's own `quant_all` routes; it never calls the tools itself. Steps:

0. `fill` (public datasets only) downloads a dataset's runs, composition table and miRBase sequences into
   `EXT_DIR/<data>/`, checked against the declared checksums (fill.py). It is the only step that writes to EXT_DIR.
1. `prepare` builds a smallrna project from a dataset declaration (`datasets/<name>.yaml`): configuration, sample
   table, reads (generated for synthetic datasets; linked, or the first `--reads` copied, for pool datasets) and the
   truth (`truth.tsv` per read for synthetic data; `pool.tsv` with truth groups for pools, pool.py), kept in
   `<project>/benchmark/`.
2. Run `snakemake quant_all` in the project as for any project on that machine (containers, cores, resources); for a
   pool, with `--notemp`, since the oracle count reads the workflow's filtered reads.
3. `evaluate` turns each method's output into one common table (adapters.py), computes the same metrics for every
   method (metrics.py, pool.py) and writes `<project>/benchmark/report.md` with the tables behind it.

References and public datasets live under `EXT_DIR` (`ext_dir` in the project configuration); a reduced reference for
mirflowz is made with reduce_reference.py and named in the declaration.

Usage: benchmark.py fill --dataset YAML --ext-dir DIR
       benchmark.py prepare --dataset YAML --project DIR --ext-dir DIR [--methods LIST] [--reads N]
       benchmark.py evaluate --project DIR [--mirflowz-gff GFF3]
"""
import argparse
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapters  # noqa: E402
import fill as filler  # noqa: E402
import metrics  # noqa: E402
import pool  # noqa: E402
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
    if spec['kind'] not in ('synthetic', 'pool'):
        fail('{}: kind must be synthetic or pool, got {!r}'.format(path, spec['kind']))
    required = {'synthetic': ('adapter', 'synthetic', 'reference'), 'pool': ('pool', 'data', 'composition', 'mirbase', 'runs')}[spec['kind']]
    for key in required:
        if key not in spec:
            fail('{}: a {} dataset needs {!r}'.format(path, spec['kind'], key))
    if spec['kind'] == 'pool' and spec['pool'] not in ('equimolar', 'ratiometric'):
        fail('{}: pool must be equimolar or ratiometric, got {!r}'.format(path, spec['pool']))
    return spec


def reference_dir(spec, ext_dir):
    return Path(ext_dir) / spec['reference']['dir'] if spec.get('reference') else None


def write_project(project, spec, ext_dir, samples, read_length, methods, mirflowz_reference=None):
    """smallrna project: config.yaml, PEP sample table and a Snakefile including the worktree's workflow."""
    config = {'workflow': 'smallrna', 'project_id': ['BENCH-0000-001'], 'organism': spec['organism'],
              'libprepkit': spec['libprepkit'], 'read_geometry': [read_length], 'machine': '', 'ext_dir': str(ext_dir),
              'fastq_dir': 'data/raw/fastq',
              'quant': {'method': ','.join(methods), 'batch': {'method': 'skip'}},
              'samples': {s: {'R1': '{}/{}_R1.fastq.gz'.format(spec['name'], s), 'R2': '', 'paired_end': 0, 'Sample_ID': s}
                          for s in samples}}
    if 'mirflowz' in methods and mirflowz_reference:
        ref = mirflowz_reference
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
    (bench / 'dataset.yaml').write_text(yaml.safe_dump(dict(spec, methods=list(methods)), sort_keys=False))
    return bench


def project_fastq(project, spec, sample):
    return project / 'data/raw/fastq' / spec['name'] / '{}_R1.fastq.gz'.format(sample)


def prepare(args):
    spec = load_dataset(args.dataset)
    project, ext_dir = Path(args.project), Path(args.ext_dir).resolve()
    if project.exists():
        fail('{} exists; choose a new project directory'.format(project))
    methods = args.methods.split(',') if args.methods else list(spec['methods'])
    unknown = [m for m in methods if m not in adapters.ADAPTERS]
    if unknown:
        fail('--methods: unknown {}; choose from {}'.format(', '.join(unknown), ', '.join(adapters.ADAPTERS)))
    (prepare_pool if spec['kind'] == 'pool' else prepare_synthetic)(spec, project, ext_dir, methods, args)


def prepare_synthetic(spec, project, ext_dir, methods, args):
    ref = reference_dir(spec, ext_dir)
    if ref is None or not (ref / 'mature.fa').is_file():
        fail('reference {} lacks mature.fa; make it with reduce_reference.py'.format(ref))
    samples = ['S{:02d}'.format(i) for i in range(1, spec['synthetic'].get('samples', 2) + 1)]
    truth = synthetic.generate(spec['synthetic'], ref / 'mature.fa', ref / 'mirna.gff3', PREFIX[spec['organism']], samples,
                               spec['adapter'], int(spec.get('umi', 0)), lambda sample: project_fastq(project, spec, sample))
    read_length = max(len(row[4]) for row in truth) + 2 * int(spec.get('umi', 0)) + len(spec['adapter'])
    bench = write_project(project, spec, ext_dir, samples, read_length, methods, ref)
    synthetic.write_truth(bench / 'truth.tsv', truth)
    print('{}: {} samples, {} reads each; run `snakemake quant_all` there, then `benchmark.py evaluate --project {}`'.format(
        project, len(samples), sum(row[5] for row in truth if row[0] == samples[0]), project))


def pool_paths(spec, ext_dir):
    data = Path(ext_dir) / spec['data']
    return (data / 'composition' / spec['composition']['file'], data / 'mirbase' / str(spec['mirbase']['release']) / 'mature.fa',
            {run['run']: data / 'fastq' / '{}.fastq.gz'.format(run['run']) for run in spec['runs']})


def head_fastq(source, target, reads):
    """Copy the first `reads` reads of a gzipped FASTQ."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, 'rt') as src, gzip.open(target, 'wt', compresslevel=3) as out:
        for i, line in enumerate(src):
            if i >= 4 * reads:
                break
            out.write(line)


def prepare_pool(spec, project, ext_dir, methods, args):
    composition, mature_fa, fastqs = pool_paths(spec, ext_dir)
    missing = [str(p) for p in [composition, mature_fa] + list(fastqs.values()) if not p.is_file()]
    if missing:
        fail('not filled: {}; run `benchmark.py fill --dataset ...` first'.format(', '.join(missing)))
    samples = list(fastqs)
    for sample, source in fastqs.items():
        target = project_fastq(project, spec, sample)
        if args.reads:
            head_fastq(source, target, args.reads)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(source)
    with gzip.open(fastqs[samples[0]], 'rt') as fh:
        fh.readline()
        read_length = len(fh.readline().strip())
    bench = write_project(project, spec, ext_dir, samples, read_length, methods)
    rows = pool.truth_groups(pool.read_composition(composition), synthetic.read_mature(mature_fa, PREFIX[spec['organism']]))
    pool.write_pool(bench / 'pool.tsv', rows)
    write_tsv(bench / 'samples.tsv', ['sample', 'group', 'lab', 'reads'],
              [[r['run'], r['group'], r.get('lab', ''), args.reads or r.get('reads', '')] for r in spec['runs']])
    groups = {g: sum(1 for r in rows if r[2] == g) for g in pool.GROUPS}
    print('{}: {} runs{}; pool {}; run `snakemake quant_all --notemp` there (the oracle reads the filtered reads), '
          'then `benchmark.py evaluate --project {}`'.format(project, len(samples), ', {} reads each'.format(args.reads) if args.reads else '',
                                                             ', '.join('{} {}'.format(n, g) for g, n in groups.items()), project))


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
    if spec['kind'] == 'pool':
        return evaluate_pool(args, project, bench, spec, ext_dir)
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


def mirflowz_annotation(ext_dir, spec):
    """mirflowz's extended annotation for the organism's full reference; fails unless exactly one is prepared."""
    found = sorted(Path(ext_dir).glob('mirflowz/mirbase-{}_ensembl-*/{}/prepared/extended_mirna_annotation_{}_nt.gff3'.format(
        spec['mirbase']['release'], spec['organism'], EXTENSION)))
    if len(found) != 1:
        fail('{} prepared mirflowz references for {} under {}; name one with --mirflowz-gff'.format(
            len(found), spec['organism'], Path(ext_dir) / 'mirflowz'))
    return found[0]


# The reads the oracle counts (rules/filter/cliptrim.smk): with Bowtie spike-ins (filter.spikein.ref lists sets, a
# project setting the workflow defaults to none, and filter.spikein.quantifier is bowtie) the prefilter's reads,
# which every method reads; otherwise the trimming route's (cutadapt, the kit default; fastp; fastq with
# filter.trim.quantifier: skip). With native spike-ins UNITAS and miRge3 read the trimmed reads, and the prefilter
# runs only for the other methods, so the trimmed reads are the ones every selection leaves. The prefilter reads the
# trimming route's reads, so with --notemp both remain.
SPIKEIN_STAGE = 'spikein/bowtie'
TRIM_STAGES = ('cutadapt', 'fastp', 'fastq')
WORKFLOW_CONFIG = Path(__file__).resolve().parents[1] / 'smallrna.config'


def bowtie_spikeins(project):
    """Whether the project removes spike-ins with the Bowtie prefilter: filter.spikein.ref lists sets (as
    rules/quant/spikein.smk reads it) and filter.spikein.quantifier, from the project or the workflow default, is bowtie."""
    config = yaml.safe_load((project / 'config.yaml').read_text()) or {}
    spikein = (config.get('filter') or {}).get('spikein') or {}
    default = yaml.safe_load(WORKFLOW_CONFIG.read_text())['filter']['spikein']['quantifier']
    return any(name.strip() for name in str(spikein.get('ref') or '').split(',')) and spikein.get('quantifier', default) == 'bowtie'


def filtered_fastq(project, sample):
    """The reads of a sample the oracle counts (temporary in the workflow: run quant_all with --notemp): the spike-in
    prefilter's with Bowtie spike-ins, otherwise those of the one trimming route that left reads; reads of two trimming
    routes mean an earlier run with other settings, so it fails."""
    filtered = project / 'data/tmp/smallrna/filter'
    name = '{}_R1.fastq'.format(sample)
    stages = [SPIKEIN_STAGE] if bowtie_spikeins(project) else TRIM_STAGES
    found = [path for path in (filtered / stage / name for stage in stages) if path.is_file()]
    if not found:
        fail('no filtered reads for {} under {}; run `snakemake quant_all --notemp`'.format(sample, filtered / '|'.join(stages)))
    if len(found) > 1:
        fail('filtered reads for {} from more than one trimming route: {}; remove those of the earlier run'.format(
            sample, ', '.join(str(path) for path in found)))
    return found[0]


def evaluate_pool(args, project, bench, spec, ext_dir):
    _, mature_fa, _ = pool_paths(spec, ext_dir)
    prefix = PREFIX[spec['organism']]
    name_to_id = {name: accession for name, (accession, _) in synthetic.read_mature(mature_fa, prefix).items()}
    methods = list(spec['methods'])
    gff = args.mirflowz_gff or (mirflowz_annotation(ext_dir, spec) if 'mirflowz' in methods else None)
    common = adapters.common_table(methods, project / 'data/tmp/smallrna/quant', name_to_id, prefix, gff)
    write_tsv(bench / 'common.tsv', adapters.COMMON_COLUMNS, common)
    entries = pool.read_pool(bench / 'pool.tsv')
    sample_groups = {row['sample']: row['group'] for row in read_tsv(bench / 'samples.tsv')}
    sequences = [e['sequence'] for e in entries if e['group'] != 'excluded']
    oracle = {sample: pool.oracle(filtered_fastq(project, sample), sequences) for sample in sample_groups}
    counts, reads = pool.method_counts(common, entries)
    write_tsv(bench / 'per_sequence.tsv', pool.PER_SEQUENCE + methods, pool.per_sequence(counts, oracle, entries, methods, sample_groups))
    rows = ['oracle'] + methods
    if spec['pool'] == 'equimolar':
        summary = pool.equimolar_metrics(counts, reads, oracle, entries, methods, list(sample_groups))
        keys = ['detected', 'cv', 'qcd', 'median_log2_vs_oracle', 'within_2fold', 'beyond_10fold', 'spearman_vs_oracle', 'outside', 'unplaced']
        observed = {}
    else:
        summary, observed = pool.ratiometric_metrics(counts, oracle, entries, methods, sample_groups)
        keys = ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls']
    write_tsv(bench / 'metrics.tsv', ['metric'] + rows, [[k] + [summary.get(r, {}).get(k, '') for r in rows] for k in keys])

    groups = {g: sum(1 for e in entries if e['group'] == g) for g in pool.GROUPS}
    exact = pool.scored(entries)
    lines = ['# {}'.format(spec['name']), '',
             '{} pool, {} ({}), {} runs. Pool sequences against miRBase {}: {}.'.format(
                 spec['pool'], spec['organism'], spec['libprepkit'], len(sample_groups), spec['mirbase']['release'],
                 ', '.join('{} {}'.format(n, g) for g, n in groups.items())),
             'Metrics over {} `exact` sequences ({} sharing their miRNA with a `shifted` one are left out); `oracle` is the '
             'filtered reads counted by exact match to the pool sequences.'.format(len(exact), groups['exact'] - len(exact)), '',
             '| metric | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 1) + '|']
    lines += ['| {} | '.format(k) + ' | '.join(str(summary.get(r, {}).get(k, '')) for r in rows) + ' |' for k in keys]
    if observed:
        designs = sorted({(e['amount_a'], e['amount_b']) for e in exact}, key=lambda ab: ab[0] / ab[1])
        lines += ['', '## Median observed log2(A/B) per designed ratio', '',
                  '| A:B | designed | n | ' + ' | '.join(r for r in rows if r in observed) + ' |', '|---' * (3 + len(observed)) + '|']
        for a, b in designs:
            ids = [e['pool_id'] for e in exact if (e['amount_a'], e['amount_b']) == (a, b)]
            lines.append('| {:g}:{:g} | {:.2f} | {} | '.format(a, b, math.log2(a / b), len(ids)) + ' | '.join(
                '{:.2f}'.format(statistics.median(observed[r][i] for i in ids)) for r in rows if r in observed) + ' |')
    (bench / 'report.md').write_text('\n'.join(lines) + '\n')
    print(bench / 'report.md')


def fill(args):
    spec = load_dataset(args.dataset)
    try:
        filler.fill(spec, Path(args.ext_dir))
    except (OSError, ValueError) as error:
        fail(str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--dataset', required=True)
    p.add_argument('--project', required=True)
    p.add_argument('--ext-dir', required=True, help='EXT_DIR holding the references and data named in the declaration')
    p.add_argument('--methods', help="comma-separated subset of the declaration's methods, e.g. unitas,mirge3")
    p.add_argument('--reads', type=int, help='pool datasets: copy the first N reads of each run instead of linking the runs')
    f = sub.add_parser('fill', help='download a public dataset into EXT_DIR (the only step that writes there)')
    f.add_argument('--dataset', required=True)
    f.add_argument('--ext-dir', required=True)
    e = sub.add_parser('evaluate')
    e.add_argument('--project', required=True)
    e.add_argument('--mirflowz-gff', type=Path, help="mirflowz's extended miRNA annotation; default from the declaration")
    args = parser.parse_args()
    {'prepare': prepare, 'fill': fill, 'evaluate': evaluate}[args.command](args)


if __name__ == '__main__':
    main()
