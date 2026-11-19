#!/usr/bin/env python
"""Benchmark small-RNA quantification methods on datasets with a known truth, and compare them on biological samples.

The benchmark runs the workflow's own `quant_all` routes; it never calls the tools itself. Steps:

0. `fill` (public datasets only) downloads a dataset's runs, composition table and miRBase sequences into
   `EXT_DIR/<data>/`, checked against the declared checksums (fill.py). It is the only step that writes to EXT_DIR.
1. `prepare` builds a smallrna project from a dataset declaration (`datasets/<name>.yaml`): configuration, sample
   table, reads (generated for synthetic datasets; linked, or the first `--reads` copied, for pool and biological
   datasets) and the truth (`truth.tsv` per read for synthetic data; `pool.tsv` with truth groups for pools, pool.py;
   none for biological samples), kept in `<project>/benchmark/`. A biological sample's runs (e.g. the lanes of one
   library) are merged by the workflow; with `--lanes` each run is its own sample, for the lane-replicate metric.
2. Run `snakemake quant_all` in the project as for any project on that machine (containers, cores, resources); for a
   pool or biological dataset, with `--notemp`, since the oracle count and the first-base check read the workflow's
   filtered reads (and the filtering table its cutadapt reports).
3. `evaluate` turns each method's output into one common table (adapters.py), computes the same metrics for every
   method (metrics.py, pool.py; biological.py compares methods with each other) and writes
   `<project>/benchmark/report.md` with the tables behind it.
4. `summary` reads several evaluated projects and writes one overview, a section per dataset kind with each project's
   key metrics per method (summary.py), as Markdown or a self-contained HTML file; it writes nothing into the projects.

References and public datasets live under `EXT_DIR` (`ext_dir` in the project configuration); a reduced reference for
mirflowz is made with reduce_reference.py and named in the declaration.

Usage: benchmark.py fill --dataset YAML --ext-dir DIR
       benchmark.py prepare --dataset YAML --project DIR --ext-dir DIR [--methods LIST] [--reads N] [--lanes]
       benchmark.py evaluate --project DIR [--mirflowz-gff GFF3] [--unitas-ref DIR]
       benchmark.py summary --projects DIR [DIR ...] --output FILE.md|FILE.html
"""
import argparse
import collections
import csv
import gzip
import json
import math
from pathlib import Path
import shlex
import statistics
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapters  # noqa: E402
import biological  # noqa: E402
import evaluation  # noqa: E402
import fill as filler  # noqa: E402
import metrics  # noqa: E402
import pool  # noqa: E402
import summary as summarizer  # noqa: E402
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
    if spec['kind'] not in ('synthetic', 'pool', 'biological'):
        fail('{}: kind must be synthetic, pool or biological, got {!r}'.format(path, spec['kind']))
    required = {'synthetic': ('adapter', 'synthetic', 'reference'), 'pool': ('pool', 'data', 'composition', 'mirbase', 'runs'),
                'biological': ('data', 'mirbase', 'samples')}[spec['kind']]
    for key in required:
        if key not in spec:
            fail('{}: a {} dataset needs {!r}'.format(path, spec['kind'], key))
    if spec['kind'] == 'pool' and spec['pool'] not in ('equimolar', 'ratiometric'):
        fail('{}: pool must be equimolar or ratiometric, got {!r}'.format(path, spec['pool']))
    if spec['kind'] == 'biological':
        names = [s['sample'] for s in spec['samples']]
        if len(set(names)) != len(names) or any(not s.get('runs') for s in spec['samples']):
            fail('{}: samples need unique names and at least one run each'.format(path))
        runs = [run['run'] for s in spec['samples'] for run in s['runs']]
        repeated = sorted({r for r in runs if runs.count(r) > 1})
        if repeated:
            fail('{}: runs listed more than once: {}'.format(path, ', '.join(repeated)))
        lanes = [run['lane'] for s in spec['samples'] for run in s['runs'] if 'lane' in run]
        if any(type(lane) is not int for lane in lanes):
            fail('{}: a run\'s lane must be a whole number, e.g. `lane: 1`'.format(path))
    return spec


def reference_dir(spec, ext_dir):
    return Path(ext_dir) / spec['reference']['dir'] if spec.get('reference') else None


def write_project(project, spec, ext_dir, samples, read_length, methods, mirflowz_reference=None, r1=None):
    """smallrna project: config.yaml, PEP sample table and a Snakefile including the worktree's workflow. `r1` gives a
    sample's read files (comma-separated, merged by the workflow); by default `<name>/<sample>_R1.fastq.gz`."""
    r1 = r1 or {s: '{}/{}_R1.fastq.gz'.format(spec['name'], s) for s in samples}
    config = {'workflow': 'smallrna', 'project_id': ['BENCH-0000-001'], 'organism': spec['organism'],
              'libprepkit': spec['libprepkit'], 'read_geometry': [read_length], 'machine': '', 'ext_dir': str(ext_dir),
              'fastq_dir': 'data/raw/fastq',
              'quant': {'method': ','.join(methods), 'batch': {'method': 'skip'}},
              'samples': {s: {'R1': r1[s], 'R2': '', 'paired_end': 0, 'Sample_ID': s} for s in samples}}
    if {'mirflowz', 'tesmall'} & set(methods) and mirflowz_reference:  # TEsmall's database is built on the mirflowz reference
        ref = mirflowz_reference
        config['db'] = {'mirflowz': {'custom': {'name': spec['reference']['name'], 'genome': str(ref / 'genome.fa.gz'),
                                                'gtf': str(ref / 'genes.gtf.gz'), 'mirna_gff': str(ref / 'mirna.gff3')}}}
    config.update(spec.get('config') or {})
    (project / 'pep').mkdir(parents=True, exist_ok=True)
    (project / 'config.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
    (project / 'pep/pep_config.yaml').write_text('pep_version: 2.0.0\nsample_table: sample_table.csv\n')
    (project / 'pep/sample_table.csv').write_text('sample_name,R1,paired_end\n' + ''.join(
        '{},{},0\n'.format(s, '"{}"'.format(r1[s]) if ',' in r1[s] else r1[s]) for s in samples))
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
    if args.reads is not None and args.reads <= 0:
        fail('--reads must be a positive integer')
    if getattr(args, 'lanes', False) and spec['kind'] != 'biological':
        fail('--lanes is for biological datasets')
    {'pool': prepare_pool, 'biological': prepare_biological, 'synthetic': prepare_synthetic}[spec['kind']](spec, project, ext_dir, methods, args)


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
    """Copy up to the limit and return the number of complete four-line records."""
    target.parent.mkdir(parents=True, exist_ok=True)
    copied = 0
    with gzip.open(source, 'rb') as src, gzip.open(target, 'wb', compresslevel=3) as out:
        while copied < reads:
            header = src.readline()
            if not header:
                break
            record = [header, src.readline(), src.readline(), src.readline()]
            if not all(record):
                fail(f'{source}: incomplete FASTQ record {copied + 1}')
            out.writelines(record)
            copied += 1
    return copied


def prepare_runs(project, spec, fastqs, runs, limit):
    counts = {}
    for run in runs:
        name = run['run']
        source, target = fastqs[name], project_fastq(project, spec, name)
        if limit is None:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(source)
            counts[name] = run.get('reads')
        else:
            counts[name] = head_fastq(source, target, limit)
    return counts


def prepare_pool(spec, project, ext_dir, methods, args):
    composition, mature_fa, fastqs = pool_paths(spec, ext_dir)
    missing = [str(p) for p in [composition, mature_fa] + list(fastqs.values()) if not p.is_file()]
    if missing:
        fail('not filled: {}; run `benchmark.py fill --dataset ...` first'.format(', '.join(missing)))
    samples = list(fastqs)
    counts = prepare_runs(project, spec, fastqs, spec['runs'], args.reads)
    with gzip.open(fastqs[samples[0]], 'rt') as fh:
        fh.readline()
        read_length = len(fh.readline().strip())
    bench = write_project(project, spec, ext_dir, samples, read_length, methods)
    rows = pool.truth_groups(pool.read_composition(composition), synthetic.read_mature(mature_fa, PREFIX[spec['organism']]))
    pool.write_pool(bench / 'pool.tsv', rows)
    write_tsv(bench / 'samples.tsv', ['sample', 'group', 'lab', 'reads', 'reads_source'],
              [[r['run'], r['group'], r.get('lab', ''), counts[r['run']] if counts[r['run']] is not None else '',
                'copied' if args.reads else 'declaration'] for r in spec['runs']])
    groups = {g: sum(1 for r in rows if r[2] == g) for g in pool.GROUPS}
    print('{}: {} runs{}; pool {}; run `snakemake quant_all --notemp` there (the oracle reads the filtered reads), '
          'then `benchmark.py evaluate --project {}`'.format(project, len(samples), ', up to {} reads each'.format(args.reads) if args.reads else '',
                                                             ', '.join('{} {}'.format(n, g) for g, n in groups.items()), project))


def biological_paths(spec, ext_dir):
    data = Path(ext_dir) / spec['data']
    return (data / 'mirbase' / str(spec['mirbase']['release']) / 'mature.fa',
            {run['run']: data / 'fastq' / '{}.fastq.gz'.format(run['run']) for s in spec['samples'] for run in s['runs']})


def fastq_lanes(path, reads=1000):
    """The lanes in the first read names of a gzipped FASTQ, from the Illumina read name
    (`instrument:run:flowcell:lane:tile:x:y`, first or second field, as SRA writes it), or None when a read name has none."""
    lanes = set()
    with gzip.open(path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i >= 4 * reads:
                break
            if i % 4 == 0:
                parts = next((f.split(':') for f in line[1:].split() if f.count(':') >= 6), None)
                if not parts or not parts[3].isdigit():
                    return None
                lanes.add(int(parts[3]))
    return lanes


def prepare_biological(spec, project, ext_dir, methods, args):
    mature_fa, fastqs = biological_paths(spec, ext_dir)
    missing = [str(p) for p in [mature_fa] + list(fastqs.values()) if not p.is_file()]
    if missing:
        fail('not filled: {}; run `benchmark.py fill --dataset ...` first'.format(', '.join(missing)))
    r1, rows = {}, []
    for s in spec['samples']:
        runs = sorted(s['runs'], key=lambda r: r.get('lane', 0))
        for run in runs:
            if 'lane' not in run:
                continue
            found = fastq_lanes(fastqs[run['run']])
            if found is None:
                print('{}: no lane in the read names; declared lane {} not checked'.format(run['run'], run['lane']), file=sys.stderr)
            elif found != {run['lane']}:
                fail('{}: declared lane {}, but its reads are from lane {}'.format(run['run'], run['lane'], ', '.join(map(str, sorted(found)))))
        files = ['{}/{}_R1.fastq.gz'.format(spec['name'], run['run']) for run in runs]
        if args.lanes:
            for run, f in zip(runs, files):
                name = '{}_L{}'.format(s['sample'], run.get('lane', run['run']))
                if name in r1:
                    fail('{}: two runs of {} in lane {} (e.g. two flow cells); --lanes needs one run per lane'.format(
                        args.dataset, s['sample'], run['lane']))
                r1[name] = f
                rows.append([name, s['sample'], s.get('group', ''), run.get('lane', ''), [run['run']]])
        else:
            r1[s['sample']] = ','.join(files)
            rows.append([s['sample'], s['sample'], s.get('group', ''), 'merged', [r['run'] for r in runs]])
    with gzip.open(fastqs[spec['samples'][0]['runs'][0]['run']], 'rt') as fh:
        fh.readline()
        read_length = len(fh.readline().strip())
    counts = prepare_runs(project, spec, fastqs, (r for sample in spec['samples'] for r in sample['runs']), args.reads)
    for row in rows:
        run_ids = row.pop()
        values = [counts[name] for name in run_ids]
        count = sum(int(n) for n in values) if all(n not in (None, '') for n in values) else ''
        row.extend([','.join(run_ids), count, 'copied' if args.reads else 'declaration'])
    bench = write_project(project, spec, ext_dir, list(r1), read_length, methods, r1=r1)
    write_tsv(bench / 'samples.tsv', ['sample', 'library', 'group', 'lane', 'runs', 'reads', 'reads_source'], rows)
    print('{}: {} {}{}; run `snakemake quant_all --notemp` there (the first-base check reads the filtered reads), then '
          '`benchmark.py evaluate --project {}`'.format(project, len(r1), 'lane samples' if args.lanes else 'samples',
                                                        ', first {} reads of each run'.format(args.reads) if args.reads else '', project))


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
    if spec['kind'] == 'biological':
        return evaluate_biological(args, project, bench, spec, ext_dir)
    prefix = PREFIX[spec['organism']]
    mature = synthetic.read_mature(reference_dir(spec, ext_dir) / 'mature.fa', prefix)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    gff = args.mirflowz_gff
    if gff is None and 'mirflowz' in spec['methods']:
        gff = ext_dir / 'mirflowz' / 'custom-{}'.format(spec['reference']['name']) / spec['organism'] / 'prepared' / \
            'extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION)
    common = common_table(args, project, bench, spec, ext_dir, spec['methods'], mature, gff)

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
    metric_rows = [[k] + [summary.get(m, {}).get(k, '') for m in methods] for k in keys]

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
    evaluation.publish(bench, ['metric'] + methods, metric_rows, lines)
    print(bench / 'report.md')


MAIN_CONFIG = Path(__file__).resolve().parents[2] / 'main.config'


def unitas_refdump(project, ext_dir, organism):
    """The UNITAS reference dump the project's UNITAS run read: db.unitas from the project, else main.config."""
    unitas = dict(yaml.safe_load(MAIN_CONFIG.read_text())['db']['unitas'])
    unitas.update(((yaml.safe_load((project / 'config.yaml').read_text()) or {}).get('db') or {}).get('unitas') or {})
    if unitas.get('mirna_db', 'mirbase') != 'mirbase':
        fail('db.unitas.mirna_db is {}: the benchmark matches miRBase names'.format(unitas['mirna_db']))
    refdump = Path(ext_dir) / 'unitas' / str(unitas['release']) / 'UNITAS_refdump_{}'.format(organism)
    if not refdump.is_dir():
        fail('no UNITAS reference dump at {}; name one with --unitas-ref'.format(refdump))
    return refdump


def common_table(args, project, bench, spec, ext_dir, methods, mature, gff):
    """Every method's counts in one table (adapters.py), written to common.tsv; with UNITAS, also how each UNITAS name
    was placed (unitas_names.tsv: by name, by hairpin and arm, or not)."""
    prefix = PREFIX[spec['organism']]
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    hairpins = adapters.unitas_hairpins(args.unitas_ref or unitas_refdump(project, ext_dir, spec['organism'])) if 'unitas' in methods else None
    common = adapters.common_table(methods, project / 'data/tmp/smallrna/quant', name_to_id, prefix, gff, hairpins, mature)
    write_tsv(bench / 'common.tsv', adapters.COMMON_COLUMNS, common)
    if hairpins is not None:
        reads = collections.Counter()
        for method, _, name, _, count in common:
            if method == 'unitas':
                reads[name] += float(count)
        placed = {name: adapters.unitas_place(name, name_to_id, prefix, hairpins, mature) for name in reads}
        write_tsv(bench / 'unitas_names.tsv', ['name', 'step', 'mature_id', 'reads'],
                  [[n, placed[n][1] or 'unplaced', placed[n][0], '{:g}'.format(reads[n])] for n in sorted(reads)])
        steps = collections.defaultdict(lambda: [0, 0.0])
        for n, (_, step) in placed.items():
            steps[step or 'unplaced'][0] += 1
            steps[step or 'unplaced'][1] += reads[n]
        labels = {'name': 'by name', 'hairpin': 'by hairpin and arm', 'unplaced': 'unplaced'}
        print('UNITAS names: ' + ', '.join('{} {} ({:.0f} reads)'.format(steps[k][0], v, steps[k][1]) for k, v in labels.items() if k in steps))
    return common


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


def trim_route(project, sample):
    """The one trimming route that left reads for a sample; reads of two routes mean an earlier run with other
    settings, so it fails."""
    filtered = project / 'data/tmp/smallrna/filter'
    name = '{}_R1.fastq'.format(sample)
    found = [stage for stage in TRIM_STAGES if (filtered / stage / name).is_file()]
    if not found:
        fail('no filtered reads for {} under {}; run `snakemake quant_all --notemp`'.format(sample, filtered / '|'.join(TRIM_STAGES)))
    if len(found) > 1:
        fail('filtered reads for {} from more than one trimming route: {}; remove those of the earlier run'.format(
            sample, ', '.join(str(filtered / stage / name) for stage in found)))
    return found[0]


def filtered_fastq(project, sample):
    """The reads of a sample the oracle counts (temporary in the workflow: run quant_all with --notemp): the spike-in
    prefilter's with Bowtie spike-ins, otherwise those of the trimming route (trim_route)."""
    filtered = project / 'data/tmp/smallrna/filter'
    name = '{}_R1.fastq'.format(sample)
    if bowtie_spikeins(project):
        if not (filtered / SPIKEIN_STAGE / name).is_file():
            fail('no filtered reads for {} under {}; run `snakemake quant_all --notemp`'.format(sample, filtered / SPIKEIN_STAGE))
        return filtered / SPIKEIN_STAGE / name
    return filtered / trim_route(project, sample) / name


# The passes of each trimming route (rules/filter/cutadapt.smk, cliptrim.smk), by the folder of the reads after them.
# The reads after the adapter pass still carry the kit's random ends, which the next pass cuts (cutadapt -u, fastp
# -f/-t, as its JSON report records); with filter.trim.quantifier: skip the workflow trims nothing (`input`).
ROUTE_PASSES = {'cutadapt': [('adapter', 'cutadapt_trimmed'), ('filter', 'cutadapt_filtered'), ('artefacts', 'cutadapt')],
                'fastp': [('adapter', 'fastp_trimmed'), ('filter', 'fastp')], 'fastq': [('input', 'fastq')]}
END_OPTIONS = {'cutadapt': {'-u': 'cut', '--cut': 'cut'},
               'fastp': {'-f': 'front', '--trim_front1': 'front', '-t': 'tail', '--trim_tail1': 'tail'}}


def end_cuts(project, route, sample):
    """(front, tail): the bases the pass after the adapter pass cuts from each end of a read, from its JSON report."""
    report = project / 'data/tmp/smallrna/filter' / '{}_stages'.format(route) / '{}.filter.json'.format(sample)
    if not report.is_file():
        fail('no {} filter report for {}: {}; run `snakemake quant_all --notemp`'.format(route, sample, report))
    data = json.loads(report.read_text())
    args = data['command_line_arguments'] if route == 'cutadapt' else shlex.split(data['command'])
    options, cuts = END_OPTIONS[route], {'front': 0, 'tail': 0}
    for i, arg in enumerate(args):
        name, _, value = arg.partition('=')
        if name in options:
            value = int(value or args[i + 1])
            kind = options[name] if options[name] != 'cut' else 'front' if value >= 0 else 'tail'
            cuts[kind] += abs(value)
    return cuts['front'], cuts['tail']


def filter_stages(project, sample):
    """[(stage, reads, (front, tail))]: a sample's reads after each pass of the filter route, the last being
    filtered_fastq, with the bases to cut from each end before matching (the adapter pass's reads only). All are
    temporary in the workflow (--notemp)."""
    filtered = project / 'data/tmp/smallrna/filter'
    name = '{}_R1.fastq'.format(sample)
    route = trim_route(project, sample)
    passes = ROUTE_PASSES[route]
    stages = [(stage, filtered / folder / name) for stage, folder in passes]
    final = filtered_fastq(project, sample)
    if not stages or stages[-1][1] != final:
        stages.append(('spikein' if SPIKEIN_STAGE in str(final) else 'filtered', final))
    missing = [str(path) for _, path in stages if not path.is_file()]
    if missing:
        fail('reads of the filter passes missing: {}; run `snakemake quant_all --notemp`'.format(', '.join(missing)))
    return [(stage, path, end_cuts(project, route, sample) if route in END_OPTIONS and i == 0 else (0, 0)) for i, (stage, path) in enumerate(stages)]


def evaluate_pool(args, project, bench, spec, ext_dir):
    _, mature_fa, _ = pool_paths(spec, ext_dir)
    prefix = PREFIX[spec['organism']]
    mature = synthetic.read_mature(mature_fa, prefix)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    methods = list(spec['methods'])
    gff = args.mirflowz_gff or (mirflowz_annotation(ext_dir, spec) if 'mirflowz' in methods else None)
    common = common_table(args, project, bench, spec, ext_dir, methods, mature, gff)
    entries = pool.read_pool(bench / 'pool.tsv')
    kinds = pool.shift_kinds(entries, dict(mature.values()))
    sample_groups = {row['sample']: row['group'] for row in read_tsv(bench / 'samples.tsv')}
    sequences = [e['sequence'] for e in entries if e['group'] != 'excluded']
    stages = {sample: filter_stages(project, sample) for sample in sample_groups}
    stage_counts = {sample: {stage: pool.oracle(path, sequences, *cuts) for stage, path, cuts in stages[sample]} for sample in sample_groups}
    stage_names = [stage for stage, _, _ in next(iter(stages.values()))]
    oracle = {sample: stage_counts[sample][stage_names[-1]] for sample in sample_groups}
    write_tsv(bench / 'filter_stages.tsv', ['sample', 'pool_id', 'group'] + stage_names,
              [[sample, e['pool_id'], e['group']] + [stage_counts[sample][st].get(e['sequence'], 0) for st in stage_names]
               for sample in sample_groups for e in entries if e['group'] != 'excluded'])
    counts, reads = pool.method_counts(common, entries)
    write_tsv(bench / 'per_sequence.tsv', pool.PER_SEQUENCE + methods, pool.per_sequence(counts, oracle, entries, methods, sample_groups, kinds))
    rows = ['oracle'] + methods
    shifted = pool.scored(entries, 'shifted')
    by_kind = {kind: [e for e in shifted if kinds[e['pool_id']] == kind] for kind in pool.SHIFT_KINDS}
    by_kind = {kind: group for kind, group in by_kind.items() if group}
    if spec['pool'] == 'equimolar':
        summary = pool.equimolar_metrics(counts, reads, oracle, entries, methods, list(sample_groups))
        keys = ['detected', 'cv', 'qcd', 'median_log2_vs_oracle', 'within_2fold', 'beyond_10fold', 'spearman_vs_oracle', 'outside', 'unplaced']
        kind_keys = ['detected', 'median_log2_vs_oracle', 'within_2fold']
        if shifted:
            for row, values in pool.agreement(counts, oracle, shifted, methods, list(sample_groups)).items():
                summary[row].update({'shifted_' + k: v for k, v in values.items()})
            keys += ['shifted_' + k for k in ['detected', 'median_log2_vs_oracle', 'within_2fold', 'beyond_10fold', 'spearman_vs_oracle']]
        per_kind = {kind: pool.agreement(counts, oracle, group, methods, list(sample_groups)) for kind, group in by_kind.items()}
        observed = {}
    else:
        summary, observed = pool.ratiometric_metrics(counts, oracle, entries, methods, sample_groups)
        keys = ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls']
        keys += ['shifted_' + k for k in keys] if shifted else []
        kind_keys = ['median_abs_error', 'called']
        per_kind = {kind: {r: pool.ratio_scores([(observed[r][e['pool_id']], math.log2(e['amount_a'] / e['amount_b']))
                                                 for e in group if not math.isnan(observed[r][e['pool_id']])])
                           for r in observed} for kind, group in by_kind.items()}
    metric_rows = [[k] + [summary.get(r, {}).get(k, '') for r in rows] for k in keys]

    groups = {g: sum(1 for e in entries if e['group'] == g) for g in pool.GROUPS}
    exact = pool.scored(entries)
    lines = ['# {}'.format(spec['name']), '',
             '{} pool, {} ({}), {} runs. Pool sequences against miRBase {}: {}.'.format(
                 spec['pool'], spec['organism'], spec['libprepkit'], len(sample_groups), spec['mirbase']['release'],
                 ', '.join('{} {}'.format(n, g) for g, n in groups.items())),
             'Metrics over {} `exact` sequences ({} sharing their miRNA with a `shifted` one are left out), and with a `shifted_` '
             'prefix over {} `shifted` sequences ({} left out); `oracle` is the filtered reads counted by exact match to the pool '
             'sequences.'.format(len(exact), groups['exact'] - len(exact), len(shifted), groups['shifted'] - len(shifted)), '',
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
    lines += shifted_section(per_kind, {kind: len(group) for kind, group in by_kind.items()}, kind_keys,
                             [r for r in rows if any(r in v for v in per_kind.values())])
    lines += filter_section(*pool.filter_losses(stage_counts, entries, stage_names), stage_names, trim_route(project, next(iter(sample_groups))),
                            next(iter(stages.values()))[0][2])
    evaluation.publish(bench, ['metric'] + rows, metric_rows, lines)
    print(bench / 'report.md')


SHIFT_TEXT = {'3p_shorter': "3' end shorter", '3p_longer': "3' end longer", '5p_moved': "5' end moved", 'mixed': 'mixed'}


def shifted_section(per_kind, sizes, keys, rows):
    """Report lines: the shifted group's metrics per kind of shift (pool.SHIFT_KINDS), with the number of sequences."""
    if not per_kind:
        return []
    lines = ['', '## Shifted sequences by kind of shift', '',
             'A method should count an isomiR on its miRNA; its count is set against the oracle (equimolar) or the design '
             '(ratiometric). The 3\' end shorter or longer keeps the seed; a moved 5\' end changes it.', '',
             '| kind | n | metric | ' + ' | '.join(rows) + ' |', '|---' * (len(rows) + 3) + '|']
    for kind, values in per_kind.items():
        for i, k in enumerate(keys):
            lines.append('| {} | {} | {} | '.format(SHIFT_TEXT[kind] if i == 0 else '', sizes[kind] if i == 0 else '', k) +
                         ' | '.join(str(values.get(r, {}).get(k, '')) for r in rows) + ' |')
    return lines


def filter_section(groups, lossy, stages, route, cuts):
    """Report lines: what each pass of the filter route removes from the pool sequences (pool.filter_losses)."""
    if len(stages) < 2:
        return ['', '## Filter route', '', 'The `{}` route has no passes to compare.'.format(route)]
    passes = stages[1:]
    cut = ' after cutting {} and {} bases, the kit\'s random ends'.format(*cuts) if any(cuts) else ''
    lines = ['', '## Filter route', '',
             'Pool reads after the `{}` route\'s `{}` step (summed over runs; matched exactly{}), and the share of them each later '
             'pass removes.'.format(route, stages[0], cut), '',
             '| group | sequences | after {} | '.format(stages[0]) + ' | '.join('removed: {}'.format(p) for p in passes) + ' | kept |',
             '|---' * (len(passes) + 4) + '|']
    share = lambda part, whole: '{:.2%}'.format(part / whole) if whole else ''
    for group, g in groups.items():
        lines.append('| {} | {} | {} | '.format(group, g['sequences'], g[stages[0]]) + ' | '.join(
            share(g[stages[i]] - g[stages[i + 1]], g[stages[0]]) for i in range(len(passes))) + ' | {} |'.format(share(g[stages[-1]], g[stages[0]])))
    lines += ['', '{} sequences with at least {} reads after the `{}` step lose more than {:.0%} of them{}'.format(
        len(lossy), pool.FILTER_MIN_READS, stages[0], pool.FILTER_LOSS, ':' if lossy else '.')]
    if lossy:
        lines += ['', '| pool_id | group | length | GC | after {} | kept | removed most by |'.format(stages[0]), '|---' * 7 + '|']
        for e, t, worst in lossy[:LOSSY_SHOWN]:
            seq = e['sequence']
            lines.append('| {} | {} | {} | {:.0%} | {} | {} | {} |'.format(e['pool_id'], e['group'], len(seq), (seq.count('G') + seq.count('C')) / len(seq),
                                                                       t[stages[0]], share(t[stages[-1]], t[stages[0]]), worst))
        if len(lossy) > LOSSY_SHOWN:
            lines.append('\nThe {} most affected are shown; all are in `filter_stages.tsv`.'.format(LOSSY_SHOWN))
    return lines


LOSSY_SHOWN = 20


# Human miRNAs the NEXTflex blood miRNA blockers target (#361), by accession (names change between miRBase releases:
# miR-451a is hsa-miR-451a-5p in 23.1); reported to show whether blockers were used.
BLOCKED = {'miR-451a': 'MIMAT0001631', 'miR-486-5p': 'MIMAT0002177', 'miR-92a-3p': 'MIMAT0000092'}


def md_table(header, rows):
    return (['| ' + ' | '.join(map(str, header)) + ' |', '|---' * len(header) + '|'] +
            ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def evaluate_biological(args, project, bench, spec, ext_dir):
    mature_fa, _ = biological_paths(spec, ext_dir)
    prefix = PREFIX[spec['organism']]
    mature = synthetic.read_mature(mature_fa, prefix)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    id_to_name = {accession: name for name, (accession, _) in mature.items()}
    id_to_seq = {accession: seq for _, (accession, seq) in mature.items()}
    methods = list(spec['methods'])
    table = read_tsv(bench / 'samples.tsv')
    samples = [row['sample'] for row in table]
    lane_mode = any(row['lane'] != 'merged' for row in table)
    # check the temporary inputs (run with --notemp) before writing anything
    reads = {s: filtered_fastq(project, s) for s in samples}
    filt_dir = project / 'data/tmp/smallrna/filter'
    cutadapt = [s for s in samples if (filt_dir / 'cutadapt' / '{}_R1.fastq'.format(s)).is_file()]
    if cutadapt and len(cutadapt) != len(samples):
        fail('cutadapt reads for some samples only ({}); remove those of an earlier run'.format(', '.join(cutadapt)))
    stage_jsons = {s: {k: filt_dir / 'cutadapt_stages' / '{}.{}.json'.format(s, k) for k in biological.CUTADAPT_PASSES} for s in cutadapt}
    missing = [str(p) for paths in stage_jsons.values() for p in paths.values() if not p.is_file()]
    if missing:
        fail('cutadapt reports missing: {}; run `snakemake quant_all --notemp`'.format(', '.join(missing)))
    try:
        filt = {s: biological.filtering(paths) for s, paths in stage_jsons.items()}
    except ValueError as e:
        fail(str(e))

    gff = args.mirflowz_gff or (mirflowz_annotation(ext_dir, spec) if 'mirflowz' in methods else None)
    raw = common_table(args, project, bench, spec, ext_dir, methods, mature, gff)
    common, groups = biological.collapse(raw, id_to_seq)
    unplaced = biological.unplaced(raw, methods, samples)

    summary = biological.method_summary(common, methods, samples)
    for m in methods:
        summary[m]['unplaced'] = unplaced[m]['share']
    keys = ['mirna_reads', 'unplaced', 'detected', 'detected_all_samples', 'top10_share']
    metric_rows = [[k] + [summary[m][k] for m in methods] for k in keys]
    pairs = biological.agreement(common, methods, samples)
    pair_keys = list(next(iter(pairs.values()), {}))
    write_tsv(bench / 'agreement.tsv', ['method_a', 'method_b'] + pair_keys,
              [[a, b] + [v.get(k, '') for k in pair_keys] for (a, b), v in pairs.items()])
    seen = {i for _, _, _, i, _ in raw}
    shown = [g for g in groups if seen & set(g)]
    lines = ['# {}'.format(spec['name']), '',
             '{} ({}), {} {}. No known composition: methods are compared with each other; agreement is not accuracy.'.format(
                 spec['organism'], spec['libprepkit'], len(samples), 'lane samples' if lane_mode else 'samples'), '',
             '## Methods', '',
             '`mirna_reads` are the reads on miRNAs with a miRBase accession; `unplaced` is the median share of a method\'s',
             'miRNA reads on names without one, which every other metric leaves out. miRNAs with the same mature sequence',
             'count as one ({} such groups among the reported miRNAs, e.g. {}).'.format(
                 len(shown), ' = '.join(id_to_name[i] for i in shown[0]) if shown else 'none'), ''] + md_table(
        ['metric'] + methods, metric_rows)
    named = [(m, n, c) for m in methods for n, c in unplaced[m]['names']]
    if named:
        lines += ['', 'Most reads without an accession: ' + '; '.join(
            '{}: {}'.format(m, ', '.join('`{}` ({:.0f})'.format(n, c) for mm, n, c in named if mm == m)) for m in methods if unplaced[m]['names'])]
    lines += ['', '## Agreement between methods (medians over samples)', ''] + md_table(
        ['methods'] + pair_keys, [['{} / {}'.format(a, b)] + [v.get(k, '') for k in pair_keys] for (a, b), v in pairs.items()])

    targets = {n: acc for n, acc in BLOCKED.items() if acc in id_to_name} if spec['organism'] == 'homo_sapiens' else {}
    if targets:
        blocked = biological.blocked(common, methods, samples, targets)
        lines += ['', '## miRNAs targeted by blood miRNA blockers (median CPM)', ''] + md_table(
            ['miRNA'] + methods, [[n] + [blocked[m][n] for m in methods] for n in targets])

    if filt:
        fk = list(next(iter(filt.values())))
        lines += ['', '## Filtering (cutadapt route)', ''] + md_table(['sample'] + fk, [[s] + [filt[s][k] for k in fk] for s in filt])

    abundance = collections.Counter()
    for counts in biological.counts_by(common).values():
        abundance.update(biological.cpm(counts))
    top = {id_to_name[i]: id_to_seq[i] for i, _ in sorted(abundance.items(), key=lambda kv: (-kv[1], kv[0]))[:20] if i in id_to_seq}
    shares = []
    for s in samples:
        counts, quality = biological.first_base(reads[s], top)
        total = sum(sum(c.values()) for c in counts.values())
        canonical = sum(c[top[n][0]] for n, c in counts.items())
        shares.append([s, total, '{:.2%}'.format(1 - canonical / total) if total else '', quality['canonical'], quality['substituted']])
    lines += ['', '## First-base substitutions in the 20 most abundant miRNAs', '',
              'Reads that match a mature miRNA from its second base on but have another first base, and the median base',
              'quality of canonical and substituted first bases (a low quality points to a sequencing error).', ''] + md_table(
        ['sample', 'reads', 'substituted', 'Q canonical', 'Q substituted'], shares)

    if lane_mode:
        lanes = collections.defaultdict(list)
        for row in table:
            lanes[row['library']].append(row['sample'])
        rep = biological.lane_replicates(common, methods, lanes)
        rk = ['dispersion_mean', 'dispersion_median', 'expected_median', 'mirna_sample_pairs', 'overdispersed', 'expected_overdispersed',
              'reads_for_cv_below_10pct']
        lines += ['', '## Lane replicates (one library in several lanes: sequencing noise only)', '',
                  'Sampling alone gives a dispersion index with mean 1, the expected median (chi-square) and the expected number',
                  'above {} by chance. All methods count the same reads, so a lane effect (e.g. a sequencing error rate that'.format(
                      biological.OVERDISPERSED),
                  'differs by lane) raises every method\'s values, most for abundant miRNAs; compare the methods with each other.',
                  '`reads_for_cv_below_10pct`: the mean reads per lane (bin) from which the median CV of CPM across lanes is below 10%.',
                  ''] + md_table(['metric'] + methods, [[k] + [rep[m][k] for m in methods] for k in rk])
        edges = biological.ABUNDANCE_BINS
        labels = {lo: '{}-{}'.format(lo, hi - 1) for lo, hi in zip(edges, edges[1:])}
        labels[edges[-1]] = '>= {}'.format(edges[-1])
        share = lambda m, lo: '{:.1%} of {}'.format(rep[m]['overdispersed_by_mean'][lo][0] / rep[m]['overdispersed_by_mean'][lo][1],
                                                    rep[m]['overdispersed_by_mean'][lo][1]) if lo in rep[m]['overdispersed_by_mean'] else ''
        lines += ['', 'Share of miRNA-sample pairs above {} by mean reads per lane:'.format(biological.OVERDISPERSED), ''] + md_table(
            ['mean reads per lane'] + methods, [[labels[lo]] + [share(m, lo) for m in methods] for lo in edges])
        listed = ['{}: {}'.format(m, id_to_name.get(i, i)) for m in methods for i in rep[m]['overdispersed_mirnas']]
        if listed:
            lines += ['', 'Most often overdispersed: ' + '; '.join(listed)]
    evaluation.publish(bench, ['metric'] + methods, metric_rows, lines)
    print(bench / 'report.md')


def fill(args):
    spec = load_dataset(args.dataset)
    try:
        filler.fill(spec, Path(args.ext_dir))
    except (OSError, ValueError) as error:
        fail(str(error))


def summary(args):
    try:
        summarizer.write(args.projects, args.output)
    except (OSError, ValueError) as error:
        fail(str(error))
    print(args.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--dataset', required=True)
    p.add_argument('--project', required=True)
    p.add_argument('--ext-dir', required=True, help='EXT_DIR holding the references and data named in the declaration')
    p.add_argument('--methods', help="comma-separated subset of the declaration's methods, e.g. unitas,mirge3")
    p.add_argument('--reads', type=int, help='pool and biological datasets: copy the first N reads of each run instead of linking the runs')
    p.add_argument('--lanes', action='store_true', help='biological datasets: each run (lane) its own sample, for the lane-replicate metric')
    f = sub.add_parser('fill', help='download a public dataset into EXT_DIR (the only step that writes there)')
    f.add_argument('--dataset', required=True)
    f.add_argument('--ext-dir', required=True)
    e = sub.add_parser('evaluate')
    e.add_argument('--project', required=True)
    e.add_argument('--mirflowz-gff', type=Path, help="mirflowz's extended miRNA annotation; default from the declaration")
    e.add_argument('--unitas-ref', type=Path, help="UNITAS reference dump (its miRBase hairpins place locus names); default from "
                   "the project's db.unitas under EXT_DIR")
    s = sub.add_parser('summary', help='one overview of several evaluated projects (Markdown or HTML by the extension)')
    s.add_argument('--projects', nargs='+', required=True)
    s.add_argument('--output', required=True, help='.md or .html; written outside the projects')
    args = parser.parse_args()
    {'prepare': prepare, 'fill': fill, 'evaluate': evaluate, 'summary': summary}[args.command](args)


if __name__ == '__main__':
    main()
