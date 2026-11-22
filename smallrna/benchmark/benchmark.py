#!/usr/bin/env python
"""Benchmark small-RNA quantification methods on datasets with a known truth, and compare them on biological samples.

The benchmark runs the workflow's own `quant_all` routes, and de.smk for differential expression; it never calls the
tools itself. Steps:

0. `fill` (public pool and biological datasets, and the genome-wide synthetic dataset's miRBase files) downloads a
   dataset's runs, composition table and miRBase files into `EXT_DIR/<data>/`, checked against the declared checksums
   (fill.py). It is the only step that writes to EXT_DIR.
1. `prepare` builds a smallrna project from a dataset declaration (`datasets/<name>.yaml`): configuration, sample
   table, reads (generated for synthetic datasets; linked, or the first `--reads` copied, for pool and biological
   datasets) and the truth (`truth.tsv` per read for synthetic data; `pool.tsv` with truth groups for pools, pool.py;
   none for biological samples and tissue mixtures), kept in `<project>/benchmark/`. A biological sample's runs (e.g. the lanes of one
   library) are merged by the workflow; with `--lanes` each run is its own sample, for the lane-replicate metric.
2. Run `snakemake quant_all` in the project as for any project on that machine (containers, cores, resources); for a
   pool or biological dataset, with `--notemp`, since the oracle count and the first-base check read the workflow's
   filtered reads (and the filtering table its cutadapt reports).
3. `evaluate` turns each method's output into one common table (adapters.py), computes the same metrics for every
   method (metrics.py, pool.py; biological.py compares methods with each other; mixture.py scores each method on
   tissue mixtures against its own counts of the pure samples, #451) and writes
   `<project>/benchmark/report.md` with the tables behind it. For a dataset that declares `de:` (a ratiometric pool,
   #414; the genome-wide synthetic dataset or a biological dataset, #447) it also writes the inputs of differential
   expression to `<project>/benchmark/de/input/`, and those of its other DE sets to folders beside it (with `within`
   each lab's to `benchmark/de/lab/<lab>/input/`, #444); `snakemake -s <repo>/smallrna/benchmark/de.smk -d <project>
   --use-apptainer` tests them with edgeR's exact and GLM tests and limma-voom (de.R), and `evaluate --de` adds the
   scores of their calls (de.py).
4. `summary` reads several evaluated projects and writes one overview, a section per dataset kind with each project's
   key metrics per method (summary.py), as Markdown or a self-contained HTML file; it writes nothing into the projects.

References and public datasets live under `EXT_DIR` (`ext_dir` in the project configuration); a reduced reference for
mirflowz is made with reduce_reference.py and named in the declaration.

Usage: benchmark.py fill --dataset YAML --ext-dir DIR [--mirbase-release R]
       benchmark.py prepare --dataset YAML --project DIR --ext-dir DIR [--mirbase-release R] [--methods LIST] [--reads N] [--lanes]
                            [--set quant.ROUTE.SETTING=VALUE ...]
       benchmark.py evaluate --project DIR [--de] [--mirflowz-gff GFF3] [--unitas-ref DIR] [--bootstrap N] [--seed S]
       benchmark.py summary --projects DIR [DIR ...] --output FILE.md|FILE.html
       benchmark.py depth --project DIR [--output DIR] [--depths 500k,1M,...] [--seed S] [--plot]
"""
import argparse
import collections
import copy
import csv
import errno
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import shlex
import shutil
import statistics
import sys
from tempfile import TemporaryDirectory

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapters  # noqa: E402
import biological  # noqa: E402
import decoys  # noqa: E402
import depth as depth_series  # noqa: E402
import de  # noqa: E402
import evaluation  # noqa: E402
import fill as filler  # noqa: E402
import genome_metrics  # noqa: E402
import imbalance as imbalance_series  # noqa: E402
import isomirs  # noqa: E402
import knockout_expected  # noqa: E402
import metrics  # noqa: E402
import mixture  # noqa: E402
import pool  # noqa: E402
import resources  # noqa: E402
import summary as summarizer  # noqa: E402
import uncertainty  # noqa: E402
import synthetic  # noqa: E402
import synthetic_genome  # noqa: E402

PREFIX = {'homo_sapiens': 'hsa-', 'mus_musculus': 'mmu-', 'rattus_norvegicus': 'rno-'}
EXTENSION = 6


def fail(message):
    sys.exit('benchmark.py: ' + message)


def merged(base, change):
    """`base` with `change` applied: mappings merged key by key, anything else (lists included) replaced."""
    out = dict(base)
    for key, value in change.items():
        out[key] = merged(out[key], value) if isinstance(value, dict) and isinstance(out.get(key), dict) else value
    return out


def load_dataset(path, mirbase_release=None):
    """The declaration, with `mirbase` the release a run uses: the declared one, or with `mirbase_release` one of those
    listed under `mirbase_releases` (#431)."""
    spec = yaml.safe_load(Path(path).read_text())
    if 'extends' in spec:  # a scenario (#453): another declaration of the same folder with some values changed
        if Path(str(spec['extends'])).name != str(spec['extends']):
            fail('{}: extends {}: name a declaration in the same folder by its file name'.format(path, spec['extends']))
        fixed = sorted(set(spec) & {'kind', 'organism', 'libprepkit', 'data', 'mirbase', 'mirbase_releases'} |
                       {'synthetic.' + k for k in (spec.get('synthetic') or {}) if k in ('mode', 'decoys')})
        if fixed:  # the data and references are the base's: a scenario changes the design only
            fail('{}: a declaration that extends another keeps its dataset and references; it changes {}'.format(path, ', '.join(fixed)))
        base_path = Path(path).parent / str(spec['extends'])
        if not base_path.is_file():
            fail('{}: extends {}, which is not a declaration beside it'.format(path, spec['extends']))
        base = yaml.safe_load(base_path.read_text())
        if 'extends' in base:
            fail('{}: extends {}, which extends another; one level only'.format(path, spec['extends']))
        if 'name' not in spec or spec['name'] == base.get('name'):
            fail('{}: a declaration that extends another needs a name of its own'.format(path))
        spec = merged(base, {k: v for k, v in spec.items() if k != 'extends'})
    for key in ('name', 'kind', 'organism', 'libprepkit', 'methods'):
        if key not in spec:
            fail('{}: missing {!r}'.format(path, key))
    if spec['organism'] not in PREFIX:
        fail('{}: organism must be one of {}'.format(path, ', '.join(PREFIX)))
    if spec['kind'] not in ('synthetic', 'pool', 'biological', 'mixture'):
        fail('{}: kind must be synthetic, pool, biological or mixture, got {!r}'.format(path, spec['kind']))
    genome = spec['kind'] == 'synthetic' and synthetic_mode(spec) == 'genome'
    required = {'synthetic': ('adapter', 'synthetic') + (('data', 'mirbase') if genome else ('reference',)),
                'pool': ('pool', 'data', 'composition', 'mirbase', 'runs'), 'biological': ('data', 'mirbase', 'samples'),
                'mixture': ('data', 'mirbase', 'samples', 'mixture')}[spec['kind']]
    for key in required:
        if key not in spec:
            fail('{}: a {} dataset needs {!r}'.format(path, spec['kind'], key))
    if mirbase_release is not None:
        spec = select_mirbase(spec, str(mirbase_release), path)
    if genome and not all(isinstance(spec['mirbase'].get(k), dict) and {'url', 'sha256'} <= set(spec['mirbase'][k]) for k in ('hairpin', 'gff')):
        fail('{}: a genome-wide synthetic dataset needs mirbase.hairpin and mirbase.gff, each with url and sha256'.format(path))
    if genome and spec['synthetic'].get('decoys') is not None:
        check_decoys(spec['synthetic'], path)
    if spec['kind'] == 'biological' and 'other_rna' in spec:
        check_other_rna(spec['other_rna'], path)
    if spec['kind'] == 'pool' and spec['pool'] not in ('equimolar', 'ratiometric'):
        fail('{}: pool must be equimolar or ratiometric, got {!r}'.format(path, spec['pool']))
    if spec['kind'] == 'mixture':  # #451: two tissues of each subject mixed at known fractions
        components = (spec['mixture'] or {}).get('components')
        if not isinstance(components, list) or len(components) != 2 or len({str(c) for c in components}) != 2:
            fail('{}: mixture.components must name the two tissues, the one at fraction 0 first'.format(path))
        try:
            mixture.Design(spec['samples'])
        except ValueError as error:
            fail('{}: {}'.format(path, error))
    if 'de' in spec:
        check_de(spec, path)
    if 'knockout' in spec:
        check_knockout(spec, path)
    if spec['kind'] in ('biological', 'mixture'):
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


DE_KINDS = 'ratiometric pools, the genome-wide synthetic dataset, mixtures and biological datasets'
DeSet = collections.namedtuple('DeSet', 'suffix folder design model samples arm', defaults=(None,))  # arm: de.ARMS (#453), else None


def de_kind(spec):
    """How a dataset can declare `de:` (#414, #447): 'pool' (a ratiometric pool), 'synthetic' (the genome-wide synthetic
    dataset), 'biological' or 'mixture'; None for the others."""
    if spec['kind'] == 'pool':
        return 'pool' if spec['pool'] == 'ratiometric' else None
    if spec['kind'] == 'synthetic':
        return 'synthetic' if synthetic_mode(spec) == 'genome' else None
    return spec['kind']  # biological, mixture


def de_contrast(spec, units=None):
    """(numerator, denominator) of a dataset's DE contrast, so that a positive log2 fold change is more in the
    numerator: pool A against pool B; a biological dataset's declared `de.contrast`; for the synthetic dataset, from its
    `units` (de_units), the group with the positive score of the group factor against the other, as its designed fold
    changes are."""
    kind = de_kind(spec)
    if kind == 'pool':
        return 'A', 'B'
    if kind == 'biological':
        contrast = (spec['de'] or {}).get('contrast')  # None with null comparisons only (#449)
        return tuple(str(group) for group in contrast) if contrast else None  # as text, as samples.tsv has the groups
    signs = {group: {float(fields['factor_group']) > 0 for _, g, fields in units if g == group} for _, group, _ in units}
    if len(signs) != 2 or sorted(signs.values(), key=sorted) != [{False}, {True}]:
        fail('benchmark/samples.tsv: the mixed samples need two groups, one with a positive and one with a negative factor_group')
    return tuple(sorted(signs, key=lambda group: signs[group] == {False}))


def check_knockout(spec, path):
    """A biological dataset's `knockout:` (#452): the `control` group (e.g. the parental cells) that the samples of
    every other group are compared with, the knockout group whose surviving miRNAs are listed (`survival`), and
    optionally `within`, a field of the samples whose levels pair each knockout sample with the controls of its own
    level (e.g. the set the libraries were made in), and `expected` (#478), a list of the miRNAs made without DROSHA (a file
    next to the declaration, written by knockout_expected.py), which splits each method's surviving miRNAs."""
    if spec['kind'] != 'biological':
        fail('{}: `knockout` (#452) is for biological datasets'.format(path))
    declared = spec['knockout']
    if not isinstance(declared, dict) or set(declared) - {'control', 'survival', 'within', 'expected'}:
        fail('{}: `knockout` takes only `control`, `survival`, `within` and `expected`'.format(path))
    groups = [str(s.get('group', '')) for s in spec.get('samples') or []]
    control, survival = (str(declared.get(k, '')) for k in ('control', 'survival'))
    if control not in groups or survival not in groups or control == survival:
        fail('{}: knockout.control and knockout.survival must name two different groups of the samples, got {!r} and {!r}'.format(
            path, declared.get('control'), declared.get('survival')))
    within = declared.get('within')
    if within is not None:
        if any(within not in s for s in spec['samples']):
            fail('{}: knockout.within is {!r}, but not every sample has it'.format(path, within))
        levels = {str(s[within]) for s in spec['samples'] if str(s['group']) == control}
        alone = [s['sample'] for s in spec['samples'] if str(s['group']) != control and str(s[within]) not in levels]
        if alone:
            fail('{}: no {} sample with the same {} as {}'.format(path, control, within, ', '.join(alone)))
    expected = declared.get('expected')
    if expected is not None:
        if not isinstance(expected, str) or not expected or Path(expected).name != expected:
            fail('{}: knockout.expected names a file next to the declaration, got {!r}'.format(path, expected))
        if not (Path(path).parent / expected).is_file():
            fail('{}: knockout.expected: {} not found; prepare copies the list next to a project\'s dataset.yaml (for a project '
                 'prepared before #478, copy it from smallrna/benchmark/datasets/)'.format(path, Path(path).parent / expected))
        knockout_expected.read_list(Path(path).parent / expected)


def knockout_units(spec, table):
    """[(sample, group, level)] of the samples in benchmark/samples.tsv for knockout(): the group from samples.tsv, the
    level from the declaration's `knockout.within` field (one level without it)."""
    within = spec['knockout'].get('within')
    level = {s['sample']: str(s[within]) if within else '' for s in spec['samples']}
    return [(row['sample'], row['group'], level[row['library']]) for row in table]


def check_de(spec, path):
    """A dataset's `de:` (#414, #447). A ratiometric pool and a biological dataset take `covariates`, columns of their
    runs or samples for the design of edgeR's GLM test and voom, and `within` (#444), a column whose levels are also
    tested on their own, each with the other covariates; a biological dataset names its `contrast`, [numerator,
    denominator] of its samples' groups, and `null` (#449), {splits, seed} of null comparisons within its groups (a
    contrast, null comparisons or both). The genome-wide synthetic dataset and a mixture take no option: their
    contrasts and sets come from their design (de_contrast and de_sets; mixture_de_sets). The whole dataset's design,
    each level's and each null split's pass check_design."""
    kind = de_kind(spec)
    if kind is None:
        fail('{}: `de` (differential expression, #414) is for {}'.format(path, DE_KINDS))
    declared = spec['de'] or {}
    options = {'pool': ['covariates', 'within'], 'synthetic': [], 'mixture': [],
               'biological': ['contrast', 'covariates', 'within', 'null']}[kind]
    if isinstance(declared, dict) and None in declared and 'null' in options:  # YAML reads a bare `null:` key as null (None)
        fail("{}: `de` has a key that YAML reads as null; write the null comparisons' key quoted, `'null':` (#449)".format(path))
    if not isinstance(declared, dict) or set(declared) - set(options):
        names = ['`{}`'.format(o) for o in options]
        fail('{}: `de` takes {}'.format(path, 'only ' + ' and '.join([', '.join(names[:-1])] * (len(names) > 1) + names[-1:]) if names else
                                        'no options for this dataset (`de: {}`): its design gives the contrasts'))
    if kind == 'mixture':
        mixture_de_sets(spec)  # every paired contrast needs two subjects with both fractions: checked here, before a run
    if kind in ('synthetic', 'mixture'):
        return
    unit = 'run' if kind == 'pool' else 'sample'
    # Groups as text: samples.tsv, which the evaluation reads, has them as text, and YAML may give numbers (100, 0).
    units = [(None if u.get('group') is None else str(u['group']), u) for u in (spec['runs'] if kind == 'pool' else spec['samples'])]
    if kind == 'biological' and 'contrast' not in declared and 'null' not in declared:
        fail('{}: a biological dataset\'s `de` needs a `contrast`, `null` comparisons (#449) or both'.format(path))
    covariates = declared.get('covariates') or []
    if not isinstance(covariates, list) or len(set(covariates)) != len(covariates) or \
            not all(isinstance(c, str) and de.NAME.fullmatch(c) and c not in ('sample', 'group', 'row') for c in covariates):
        fail('{}: de.covariates must be distinct column names (letters, digits and `_`, starting with a letter; not sample, group '
             'or row), got {!r}'.format(path, covariates))
    within = declared.get('within')
    if within is not None and (not isinstance(within, str) or not de.NAME.fullmatch(within)):
        fail('{}: de.within must be a column name (letters, digits and `_`, starting with a letter), got {!r}'.format(path, within))
    if kind == 'biological' and 'null' in declared:
        check_null(declared['null'], units, covariates, covariates + ([within] if within and within not in covariates else []), path)
    if kind == 'biological' and 'contrast' not in declared:
        return
    if kind == 'biological':
        contrast = declared.get('contrast')
        present = {group for group, _ in units}
        labels = [str(group) for group in contrast] if isinstance(contrast, list) and None not in contrast else []
        if len(labels) != 2 or labels[0] == labels[1] or not set(labels) <= present or \
                any(not label.strip() or set(label) & set('\t\n') for label in labels):
            fail('{}: de.contrast must name two different groups of the samples, [numerator, denominator], got {!r}'.format(path, contrast))
    contrast = de_contrast(spec)
    units = [(group, fields) for group, fields in units if group in contrast]
    check_design(units, contrast, covariates, path, unit)
    if within is None:
        return
    values = [fields.get(within) for _, fields in units]
    if not all(isinstance(v, str) and de.LEVEL_NAME.fullmatch(v) for v in values) or len(set(values)) < 2:
        fail('{}: de.within {}: every {} needs a value that starts with a letter or digit and has only letters, digits, `.`, '
             '`_` or `-` (it names a folder and metrics), with at least two different ones'.format(path, within, unit))
    if len({value.lower() for value in values}) != len(set(values)):
        fail('{}: de.within {}: values that differ only in case would share a folder: {}'.format(path, within, ', '.join(sorted(set(values)))))
    for level in sorted(set(values)):
        check_design([(group, fields) for group, fields in units if fields[within] == level], contrast, [c for c in covariates if c != within],
                     '{}: de.within {} {}'.format(path, within, level), unit)


def null_covariates(members, covariates):
    """The covariates of a group's null comparisons (#449): those that vary within the group; one the group's samples
    share carries nothing about its halves, and would leave the design without full rank."""
    return [c for c in covariates if len({str(fields.get(c)) for _, _, fields in members}) > 1]


def check_null(declared, units, covariates, strata, path):
    """A biological dataset's `de.null` (#449): {splits: 1 to 500, seed: a whole number}, at least one group with
    de.NULL_MIN samples, and every split (balanced on the `strata` columns) passing check_design with the covariates
    that vary in its group."""
    whole = lambda v: isinstance(v, int) and not isinstance(v, bool)
    if not isinstance(declared, dict) or set(declared) != {'splits', 'seed'} or not whole(declared['seed']) or \
            not whole(declared['splits']) or not 1 <= declared['splits'] <= 500:
        fail('{}: de.null takes `splits` (a whole number from 1 to 500) and `seed` (a whole number), got {!r}'.format(path, declared))
    members = [(fields['sample'], group, fields) for group, fields in units if group is not None]
    splits = de.null_splits(members, declared['splits'], declared['seed'], strata)
    if not splits:
        fail('{}: de.null needs a group with at least {} samples, two in each half'.format(path, de.NULL_MIN))
    for group, halves in splits.items():
        in_group = [m for m in members if m[1] == group]
        names = null_covariates(in_group, covariates)
        fields = {sample: f for sample, _, f in in_group}
        for x, y in halves:
            check_design([('x', fields[s]) for s in x] + [('y', fields[s]) for s in y], ('x', 'y'), names,
                         '{}: de.null, a split of {}'.format(path, group), 'sample')


def check_design(units, contrast, covariates, where, unit='run'):
    """A DE design's runs or samples (#414), `units`: [(group, {column: value})] of the two groups of `contrast`: at
    least two in each group, each covariate with a value for every one and at least two different ones, and residual
    degrees of freedom for ~ covariates + group. `where` begins the messages."""
    for covariate in covariates:
        values = [fields.get(covariate) for _, fields in units]
        if not all(v is not None and str(v).strip() and not set(str(v)) & set('\t\n') for v in values) or len(set(values)) < 2:
            fail('{}: de.covariates {}: every {} needs a value, with at least two different ones'.format(where, covariate, unit))
    sizes = collections.Counter(group for group, _ in units)
    if min(sizes[group] for group in contrast) < 2:
        fail('{}: differential expression needs at least two {}s in each group, got {} {} and {} {}'.format(
            where, unit, sizes[contrast[0]], contrast[0], sizes[contrast[1]], contrast[1]))
    coefficients = 2 + sum(len({fields[c] for _, fields in units}) - 1 for c in covariates)  # intercept, covariate levels, group
    if len(units) - coefficients < 1:
        fail('{}: the design ~ {} of the GLM test and voom has {} coefficients for {} {}s, so no residual degrees of freedom'.format(
            where, ' + '.join(covariates + ['group']), coefficients, len(units), unit))


def select_mirbase(spec, release, path):
    """The declaration with `mirbase` set to `release`: its own, or one listed under `mirbase_releases` with url and
    sha256 (and, for a genome-wide synthetic dataset, hairpin and gff, checked by load_dataset)."""
    if 'mirbase' not in spec:
        fail('{}: --mirbase-release needs a declaration with miRBase files (mirbase)'.format(path))
    declared = str(spec['mirbase']['release'])
    if release == declared:
        return spec
    if spec['kind'] == 'synthetic' and synthetic_mode(spec) == 'reduced':
        fail('{}: a reduced synthetic dataset\'s reference fixes its miRBase release; --mirbase-release is for pools, '
             'biological and genome-wide synthetic datasets'.format(path))
    others = {str(r): files for r, files in (spec.get('mirbase_releases') or {}).items()}
    files = others.get(release)
    if not isinstance(files, dict) or not {'url', 'sha256'} <= set(files):
        fail('{}: miRBase {} is not declared with url and sha256; declared: {}'.format(path, release, ', '.join([declared] + sorted(others))))
    return dict(spec, mirbase=dict(files, release=release))


def synthetic_mode(spec):
    """`reduced` (reads on a reduced reference, the quick check) or `genome` (genome-wide, from filled miRBase sources)."""
    mode = (spec.get('synthetic') or {}).get('mode', 'reduced')
    if mode not in ('reduced', 'genome'):
        fail('{}: synthetic.mode must be reduced or genome, got {!r}'.format(spec.get('name'), mode))
    return mode


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
    if {'mirflowz', 'tesmall', 'tinyrna', 'biocore_htseq'} & set(methods) and mirflowz_reference:  # the routes on the mirflowz reference
        ref = mirflowz_reference
        config['db'] = {'mirflowz': {'custom': {'name': spec['reference']['name'], 'genome': str(ref / 'genome.fa.gz'),
                                                'gtf': str(ref / 'genes.gtf.gz'), 'mirna_gff': str(ref / 'mirna.gff3')}}}
    config.update(copy.deepcopy(spec.get('config') or {}))  # a variant below must not reach the saved declaration
    for key, value in (spec.get('variant') or {}).items():  # a settings variant (#417); the workflow fills in the other settings
        _, route, setting = key.split('.')
        config['quant'].setdefault(route, {})[setting] = value
    if 'mirbase' in spec:  # every method uses the declaration's release (#407), also when another is selected (#431)
        release = str(spec['mirbase']['release'])
        db = config.setdefault('db', {})
        configured = str((db.get('mirbase') or {}).get('release', release))
        if configured != release:
            fail('the declaration\'s config sets db.mirbase.release {}, but the run uses miRBase {}'.format(configured, release))
        db.setdefault('mirbase', {})['release'] = release
    (project / 'pep').mkdir(parents=True, exist_ok=True)
    (project / 'config.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
    (project / 'pep/pep_config.yaml').write_text('pep_version: 2.0.0\nsample_table: sample_table.csv\n')
    (project / 'pep/sample_table.csv').write_text('sample_name,R1,paired_end\n' + ''.join(
        '{},{},0\n'.format(s, '"{}"'.format(r1[s]) if ',' in r1[s] else r1[s]) for s in samples))
    workflow = Path(__file__).resolve().parents[1] / 'smallrna.smk'
    (project / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\ninclude: {!r}\n'.format(str(workflow)))
    bench = project / 'benchmark'
    bench.mkdir()
    saved = {k: v for k, v in spec.items() if k != 'mirbase_releases'}  # the project uses one release
    saved['methods'] = list(methods)
    left_out = [m for m in spec['methods'] if m not in methods]  # declared, not selected (#448): `not run` in routes.tsv
    if left_out:
        saved['not_selected'] = left_out
    (bench / 'dataset.yaml').write_text(yaml.safe_dump(saved, sort_keys=False))
    return bench


def check_variant(spec):
    """A settings variant (#417) is scored only when its route's settings record, bound to the count table, has the
    setting at the variant's value: the record's row whose source is the setting's key."""
    for key, value in (spec.get('variant') or {}).items():
        route = key.split('.')[1]
        if route not in spec['route_settings']:
            fail('settings variant {}: {} is not scored'.format(key, route))
        rows = spec['route_settings'][route]
        if rows is None:
            fail('settings variant {}: {} has no settings record (quant/{}/settings.tsv), so its setting cannot be confirmed'.format(key, route, route))
        recorded = [v for _, v, source in rows if source == key]
        if [str(v).strip().lower() for v in recorded] != [str(value).strip().lower()]:
            fail('settings variant {}: {} counted with {}, not {}'.format(key, route, ', '.join(recorded) or 'no such setting recorded', value))


def settings_variant(assignments, methods):
    """{quant.ROUTE.SETTING: value} of `prepare --set` (#417): each a setting smallrna.config has for a route the project
    runs, set once, to a value (read as YAML) other than its default."""
    defaults = (yaml.safe_load((Path(__file__).resolve().parents[1] / 'smallrna.config').read_text()) or {}).get('quant') or {}
    variant = {}
    for assignment in assignments:
        key, equals, text = assignment.partition('=')
        parts = key.split('.')
        if not equals or len(parts) != 3 or parts[0] != 'quant':
            fail('--set {}: expected quant.ROUTE.SETTING=VALUE'.format(assignment))
        if parts[1] not in methods:
            fail('--set {}: {} is not one of the project\'s routes ({})'.format(assignment, parts[1], ', '.join(methods)))
        if parts[2] not in (defaults.get(parts[1]) or {}):
            fail('--set {}: smallrna.config has no setting {} for {}'.format(assignment, parts[2], parts[1]))
        value = yaml.safe_load(text)
        if key in variant or value == defaults[parts[1]][parts[2]]:
            fail('--set {}: {}'.format(assignment, 'set twice' if key in variant else 'that is the default'))
        variant[key] = value
    return variant


def select_per_group(spec, count):
    """The declaration with the first `count` samples of each group, in declared order (#448): a quick run of a
    biological dataset, e.g. `--per-group 2` for GSE317995's two lowest-numbered samples per group."""
    if spec['kind'] != 'biological':
        fail('--per-group is for biological datasets')
    if count <= 0:
        fail('--per-group must be a positive integer')
    taken = collections.Counter()
    samples = []
    for sample in spec['samples']:
        if taken[sample.get('group', '')] < count:
            taken[sample.get('group', '')] += 1
            samples.append(sample)
    selected = dict(spec, samples=samples)
    declared = spec.get('de') or {}
    if 'null' in declared and max(taken.values()) < de.NULL_MIN:  # the null comparisons need larger groups (#449)
        print('--per-group {}: no group has {} samples, so the null comparisons (de.null) are left out'.format(count, de.NULL_MIN), file=sys.stderr)
        rest = {k: v for k, v in declared.items() if k != 'null'}
        selected = dict(selected, de=rest) if rest else {k: v for k, v in selected.items() if k != 'de'}
    return selected


def project_fastq(project, spec, sample):
    return project / 'data/raw/fastq' / spec['name'] / '{}_R1.fastq.gz'.format(sample)


def prepare(args):
    spec = load_dataset(args.dataset, args.mirbase_release)
    project, ext_dir = Path(args.project), Path(args.ext_dir).resolve()
    if project.exists() or project.is_symlink():
        fail('{} exists; choose a new project directory'.format(project))
    methods = args.methods.split(',') if args.methods else list(spec['methods'])
    unknown = [m for m in methods if m not in adapters.ADAPTERS]
    if unknown:
        fail('--methods: unknown {}; choose from {}'.format(', '.join(unknown), ', '.join(adapters.ADAPTERS)))
    if args.set:
        spec = dict(spec, variant=settings_variant(args.set, methods))
    if args.reads is not None and args.reads <= 0:
        fail('--reads must be a positive integer')
    if getattr(args, 'lanes', False) and spec['kind'] != 'biological':
        fail('--lanes is for biological datasets')
    if getattr(args, 'per_group', None) is not None:
        spec = select_per_group(spec, args.per_group)
    if getattr(args, 'synthetic_seed', None) is not None:  # another draw of a synthetic set (#455); its dataset.yaml keeps it
        if spec['kind'] != 'synthetic':
            fail('--synthetic-seed is for synthetic datasets')
        spec = dict(spec, synthetic=dict(spec['synthetic'], seed=args.synthetic_seed))
    prepare_dataset = {'pool': prepare_pool, 'biological': prepare_biological, 'mixture': prepare_biological,
                       'synthetic': prepare_synthetic}[spec['kind']]
    project.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=project.parent) as temporary:
        staged = Path(temporary) / 'project'
        message = prepare_dataset(spec, staged, ext_dir, methods, args)
        expected = (spec.get('knockout') or {}).get('expected')
        if expected:  # the project keeps the list it is evaluated with, next to its dataset.yaml
            shutil.copyfile(Path(args.dataset).parent / expected, staged / 'benchmark' / expected)
        if project.exists() or project.is_symlink():
            fail('{} exists; choose a new project directory'.format(project))
        # A completed project is nonempty. A concurrent prepare can win the
        # rename, but POSIX rename cannot replace that winner's directory.
        try:
            staged.rename(project)
        except OSError as error:
            if error.errno in (errno.EEXIST, errno.ENOTEMPTY):
                fail('{} exists; choose a new project directory'.format(project))
            raise
    print(message)


def check_decoys(syn, path):
    """`synthetic.decoys` of a genome-wide synthetic declaration (#412 PR 2; decoys.py): `load` {unit: decoy reads per miRNA
    read} of units in decoys.UNITS, `sources` {source: {url, sha256}} for every source a declared unit reads, `probe.reads`
    (a probe sample per unit needs `synthetic.probes`), and the optional fragment settings."""
    d = syn['decoys']
    number = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool)
    where = '{}: synthetic.decoys'.format(path)
    known = {'load', 'sources', 'probe'} | set(decoys.DEFAULTS)
    if not isinstance(d, dict) or set(d) - known:
        fail('{} takes {}, got {!r}'.format(where, ', '.join(sorted(known)), sorted(d) if isinstance(d, dict) else d))
    load = d.get('load')
    if not isinstance(load, dict) or not load or set(load) - set(decoys.UNITS) or not all(number(v) and v >= 0 for v in load.values()):
        fail('{}.load must give decoy reads per miRNA read (a number of at least 0) for units of {}, got {!r}'.format(
            where, ', '.join(decoys.UNITS), load))
    sources = d.get('sources') or {}
    needed = sorted({decoys.UNITS[u] for u in load if decoys.UNITS[u]})
    bad = [s for s in needed if not (isinstance(sources.get(s), dict) and {'url', 'sha256'} <= set(sources[s]))]
    if bad or set(sources) - set(decoys.SOURCES):
        fail('{}.sources needs url and sha256 for {} (sources: {}), got {!r}'.format(where, ', '.join(needed), ', '.join(decoys.SOURCES),
                                                                                 sorted(sources)))
    reads = (d.get('probe') or {}).get('reads', 0)
    if type(reads) is not int or reads < 0 or (reads and not syn.get('probes')):
        fail('{}.probe.reads must be a whole number of at least 0, and needs synthetic.probes (their carrier and control), '
             'got {!r}'.format(where, reads))
    length = d.get('length', decoys.DEFAULTS['length'])
    if not (isinstance(length, list) and len(length) == 2 and all(type(v) is int for v in length) and 15 <= length[0] <= length[1] <= decoys.MAX_LENGTH):
        fail('{}.length must be [shortest, longest] between 15 and {}, got {!r}'.format(where, decoys.MAX_LENGTH, length))
    for key in ('modification', 'pirna_5u'):
        if not (number(d.get(key, 0)) and 0 <= d.get(key, 0) <= 1):
            fail('{}.{} must be a share between 0 and 1, got {!r}'.format(where, key, d[key]))
    if not (type(d.get('fragments', 1)) is int and d.get('fragments', 1) > 0) or not number(d.get('weight_sd', 0)) or \
            d.get('weight_sd', 0) < 0 or not number(d.get('nuisance', 0)):
        fail('{}: fragments must be a whole number above 0, weight_sd a number of at least 0 and nuisance a number'.format(where))


OTHER_RNA_SOURCES = ('refseq_rrna', 'gtrnadb', 'ensembl_ncrna')  # the sources of the embedding classes' molecules (#415)


def check_other_rna(declared, path):
    """A biological declaration's `other_rna` (#415): `sources` {source: {url, sha256}} of OTHER_RNA_SOURCES, at least
    one, whose molecules flag the miRBase entries that lie inside them."""
    sources = declared.get('sources') if isinstance(declared, dict) else None
    if not isinstance(declared, dict) or set(declared) != {'sources'} or not isinstance(sources, dict) or not sources or \
            set(sources) - set(OTHER_RNA_SOURCES) or not all(isinstance(v, dict) and {'url', 'sha256'} <= set(v) for v in sources.values()):
        fail('{}: other_rna takes `sources`, each of {} with url and sha256, got {!r}'.format(path, ', '.join(OTHER_RNA_SOURCES), declared))


def other_rna_files(spec, ext_dir):
    """{source: filled path} of a biological declaration's other_rna sources (fill.decoy_path, folder other_rna)."""
    sources = (spec.get('other_rna') or {}).get('sources') or {}
    return {key: filler.decoy_path(spec, ext_dir, key, value['url'], 'other_rna') for key, value in sources.items()}


def other_rna_molecules(spec, ext_dir):
    """{unit: molecules} of the embedding classes whose source a biological declaration's other_rna fills (#415)."""
    files = other_rna_files(spec, ext_dir)
    missing = sorted(str(p) for p in files.values() if not p.is_file())
    if missing:
        fail('not filled: {}; run `benchmark.py fill --dataset ...` first'.format(', '.join(missing)))
    return {unit: decoys.molecules(unit, files[source]) for unit, source in decoys.UNITS.items()
            if source in files and decoys.unit_class(unit) in decoys.EMBEDDING}


def other_rna_section(flag, iso_rows, methods, no_isomirs, summary, common, molecules, id_to_name):
    """#415 on a biological dataset: the miRBase entries inside molecules of other RNA classes (`flag`, {accession:
    [(class, molecules)]}), each method's read_in_other_class (from its isomiR rows), and the report lines."""
    within = decoys.contained(molecules, {row[4] for row in iso_rows if row[4]})
    reads, inside = collections.Counter(), collections.Counter()
    for method, _, _, _, sequence, _, count in iso_rows:
        if sequence:
            reads[method] += float(count)
            if sequence in within:
                inside[method] += float(count)
    for m in methods:
        summary[m]['read_in_other_class'] = round(inside[m] / reads[m], 4) if m not in no_isomirs and reads[m] else ''
    by_class = collections.defaultdict(list)
    for accession, places in flag.items():
        for cls, names in places:
            by_class[cls].append('{} ({})'.format(id_to_name.get(accession, accession), ', '.join(names[:3]) + (', ...' if len(names) > 3 else '')))
    on = collections.defaultdict(collections.Counter)
    for m, _, _, accession, count in common:
        if accession in flag:
            on[m][accession] += float(count)
    lines = ['', '## miRNAs inside other RNA classes (#415)', '',
             'miRBase entries whose mature sequence lies exactly inside a molecule of another RNA class (the declaration\'s '
             '`other_rna` sources: RefSeq and Ensembl rRNA, GtRNAdb and mitochondrial tRNAs, Ensembl Y RNA, snoRNA and other '
             'ncRNA) are flagged ({} entries, in `other_rna.tsv`); fragments overlapping a mature only partly are not. Per '
             'method: `flagged_reads`, the median share of its miRNA reads on flagged entries; `flagged_detected` and '
             '`detected_unflagged`, the median detected miRNAs among and without them; `read_in_other_class`, the share of its '
             'miRNA reads (isomiR table, #454) whose read sequence lies exactly inside another class\'s molecule (empty for '
             'methods without isomiR output: {}). Both are lower bounds: five classes (no mRNA, introns or repeats) and exact '
             'matches only, so a read with a non-templated end or an error does not match.'.format(len(flag), ', '.join(no_isomirs) or 'none'), '']
    lines += md_table(['metric'] + methods, [[k] + [summary[m].get(k, '') for m in methods]
                                             for k in ('detected', 'flagged_detected', 'detected_unflagged', 'flagged_reads', 'read_in_other_class')])
    lines += [''] + ['- {}: {}.'.format(cls, '; '.join(sorted(entries))) for cls, entries in sorted(by_class.items())]
    top = ['- {}: '.format(m) + ', '.join('{} ({:.0f})'.format(id_to_name.get(a, a), v) for a, v in on[m].most_common(5)) for m in methods if on[m]]
    if top:
        lines += ['', 'The flagged miRNAs with the most reads, per method (all samples):', ''] + top
    return lines


def decoy_files(spec, ext_dir):
    """{source: filled path} of a genome-wide synthetic declaration's decoy sources (fill.decoy_path)."""
    sources = ((spec['synthetic'].get('decoys') or {}).get('sources') or {})
    return {key: filler.decoy_path(spec, ext_dir, key, value['url']) for key, value in sources.items()}


def synthetic_files(spec, ext_dir):
    """The filled miRBase files of a genome-wide synthetic declaration: mature.fa, hairpin.fa, genes.gff3."""
    folder = Path(ext_dir) / spec['data'] / 'mirbase' / str(spec['mirbase']['release'])
    return [folder / name for name in ('mature.fa', 'hairpin.fa', 'genes.gff3')]


def prepare_synthetic_genome(spec, project, ext_dir, methods, args):
    files = synthetic_files(spec, ext_dir)
    sources = decoy_files(spec, ext_dir)
    missing = [str(p) for p in files + list(sources.values()) if not p.is_file()]
    if missing:
        fail('not filled: {}; run `benchmark.py fill --dataset ...` first'.format(', '.join(missing)))
    syn = dict(spec['synthetic'])
    if args.reads:  # --reads sets the mean reads of a mixed sample, for a quick run
        syn['reads'] = dict(syn['reads'], mean=args.reads)
    umi = int(spec.get('umi', 0))
    truth, samples, mirnas = synthetic_genome.generate(syn, *files, PREFIX[spec['organism']], spec['adapter'], umi,
                                                       lambda sample: project_fastq(project, spec, sample), sources)
    read_length = max(len(row[4]) for row in truth) + 2 * umi + len(spec['adapter'])
    bench = write_project(project, spec, ext_dir, [s[0] for s in samples], read_length, methods)
    synthetic_genome.write_table(bench / 'truth.tsv', synthetic_genome.TRUTH_COLUMNS, truth)
    synthetic_genome.write_table(bench / 'samples.tsv', synthetic_genome.SAMPLE_COLUMNS, samples)
    synthetic_genome.write_table(bench / 'mirnas.tsv', synthetic_genome.MIRNA_COLUMNS, mirnas)
    mixed = [s for s in samples if s[1] == 'mixed']
    return ('{}: {} mixed samples ({:,} reads on average), {} probe samples, {} truth miRNAs; run `snakemake quant_all` '
            'there, then `benchmark.py evaluate --project {}`'.format(args.project, len(mixed), round(sum(s[5] for s in mixed) / len(mixed)),
                                                                     len(samples) - len(mixed), sum(1 for m in mirnas if float(m[4]) > 0), args.project))


def prepare_synthetic(spec, project, ext_dir, methods, args):
    if synthetic_mode(spec) == 'genome':
        return prepare_synthetic_genome(spec, project, ext_dir, methods, args)
    ref = reference_dir(spec, ext_dir)
    if ref is None or not (ref / 'mature.fa').is_file():
        fail('reference {} lacks mature.fa; make it with reduce_reference.py'.format(ref))
    samples = ['S{:02d}'.format(i) for i in range(1, spec['synthetic'].get('samples', 2) + 1)]
    truth = synthetic.generate(spec['synthetic'], ref / 'mature.fa', ref / 'mirna.gff3', PREFIX[spec['organism']], samples,
                               spec['adapter'], int(spec.get('umi', 0)), lambda sample: project_fastq(project, spec, sample))
    read_length = max(len(row[4]) for row in truth) + 2 * int(spec.get('umi', 0)) + len(spec['adapter'])
    bench = write_project(project, spec, ext_dir, samples, read_length, methods, ref)
    synthetic.write_truth(bench / 'truth.tsv', truth)
    return '{}: {} samples, {} reads each; run `snakemake quant_all` there, then `benchmark.py evaluate --project {}`'.format(
        args.project, len(samples), sum(row[5] for row in truth if row[0] == samples[0]), args.project)


def pool_paths(spec, ext_dir):
    data = Path(ext_dir) / spec['data']
    return (data / 'composition' / spec['composition']['file'], data / 'mirbase' / str(spec['mirbase']['release']) / 'mature.fa',
            {run['run']: data / 'fastq' / '{}.fastq.gz'.format(run['run']) for run in spec['runs']})


def head_fastq(source, target, reads):
    """Copy up to the limit and return the number of complete four-line records; a damaged gzip source stops."""
    target.parent.mkdir(parents=True, exist_ok=True)
    copied = 0
    try:
        with gzip.open(source, 'rb') as src, gzip.open(target, 'wb', compresslevel=3) as out:
            while copied < reads:
                header = src.readline()
                if not header:
                    break
                record = [header, src.readline(), src.readline(), src.readline()]
                if not all(record):
                    fail(f'{source}: incomplete FASTQ record {copied + 1}')
                name, sequence, plus, quality = (line.rstrip(b'\r\n') for line in record)
                if not (name.startswith(b'@') and len(name) > 1 and plus.startswith(b'+') and sequence and len(sequence) == len(quality)):
                    fail(f'{source}: malformed FASTQ record {copied + 1}')
                out.writelines(record)
                copied += 1
            if copied == reads:
                src.readline()  # at the limit, read on so that a source ending here still has its gzip trailer checked
    except (OSError, EOFError) as error:
        fail(f'{source}: damaged gzip file ({error})')
    return copied


def declared_reads(run):
    count = run.get('reads')
    return int(count) if str(count).isdecimal() else None


def prepare_runs(project, spec, fastqs, runs, limit):
    """Link or copy (first `limit` reads) each run into the project; {run: reads copied, or the declared count of a
    linked run, None when not declared}. The summary reports the reads a run processed from its own reports instead."""
    counts = {}
    for run in runs:
        name = run['run']
        source, target = fastqs[name], project_fastq(project, spec, name)
        if limit is None:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.symlink_to(source)
            counts[name] = declared_reads(run)
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
    write_tsv(bench / 'samples.tsv', ['sample', 'group', 'lab', 'reads'],
              [[r['run'], r['group'], r.get('lab', ''), '' if counts[r['run']] is None else counts[r['run']]] for r in spec['runs']])
    groups = {g: sum(1 for r in rows if r[2] == g) for g in pool.GROUPS}
    return ('{}: {} runs{}; pool {}; run `snakemake quant_all --notemp` there (the oracle reads the filtered reads), '
          'then `benchmark.py evaluate --project {}`'.format(args.project, len(samples), ', up to {} reads each'.format(args.reads) if args.reads else '',
                                                             ', '.join('{} {}'.format(n, g) for g, n in groups.items()), args.project))


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
        count = sum(values) if None not in values else None
        row.extend([','.join(run_ids), '' if count is None else count])
    bench = write_project(project, spec, ext_dir, list(r1), read_length, methods, r1=r1)
    write_tsv(bench / 'samples.tsv', ['sample', 'library', 'group', 'lane', 'runs', 'reads'], rows)
    return ('{}: {} {}{}; run `snakemake quant_all --notemp` there (the first-base check reads the filtered reads), then '
          '`benchmark.py evaluate --project {}`'.format(args.project, len(r1), 'lane samples' if args.lanes else 'samples',
                                                        ', up to {} reads of each run'.format(args.reads) if args.reads else '', args.project))


def bootstrap(args, point, compute, draw, rows, keys, units, shown=None):
    """Bootstrap intervals for an evaluation (#413): (extra files for evaluation.publish, the record of resamples and
    seed, report lines). Nothing when --bootstrap is 0. The files have every metric in `keys`, the report's section the
    metrics in `shown` (default: all); a lab's DE metrics have their own report section (#444)."""
    n, seed = args.bootstrap, args.seed
    if not n:
        return {}, None, []
    ci, differences, skipped = uncertainty.intervals(point, compute, draw, rows, keys, n, seed)
    extra = {'metrics_ci.tsv': uncertainty.table_text(uncertainty.CI_HEADER, ci),
             'differences.tsv': uncertainty.table_text(uncertainty.DIFFERENCE_HEADER, differences)}
    record = {'resamples': n, 'seed': seed, 'skipped': skipped}
    shown = set(keys if shown is None else shown)
    return extra, record, uncertainty.report_lines([row for row in ci if row[0] in shown], [row for row in differences if row[0] in shown],
                                                   rows, [key for key in keys if key in shown], units, n, seed, skipped)


def run_resources(project, bench, spec, reads):
    """(resources.tsv text, report lines) of the methods' runtime and memory (resources.py, #386); a malformed benchmark
    file stops the evaluation."""
    try:
        return resources.summarize(project / 'data/tmp/smallrna/quant', spec['methods'], sorted(scored_samples(spec, bench)), reads)
    except ValueError as error:
        fail(str(error))


def route_settings(project, route):
    """[(field, value, source)] of a route's counting settings, recorded with its count table by quant_settings
    (quant/<route>/settings.tsv, #417), or None for a project quantified before they were recorded."""
    path = project / 'data/tmp/smallrna/quant' / route / 'settings.tsv'
    if not path.is_file():
        return None
    return [(row['field'], row['value'], row['source']) for row in read_tsv(path) if row['field'] != 'table_sha256']


def settings_lines(settings, variant=None):
    """report.md's lines on the evaluated routes' counting settings (#417), and settings.tsv's text; `variant`, the
    settings `prepare --set` changed."""
    lines = ['', 'Counting settings (#417), as each route recorded them with its count table (`quant/<route>/settings.tsv`; all in '
             '`settings.tsv`); `fixed` ones are set in the rule or script, the others by the configuration or the kit:', '']
    if variant:
        lines[1:1] = ['', '**A settings variant** (`prepare --set`): {}. `summary` compares it with a project of the same dataset run '
                      'with the defaults.'.format('; '.join('{} `{}`'.format(k, v) for k, v in sorted(variant.items())))]
    for route, rows in settings.items():
        lines.append('- {}: {}'.format(route, 'not recorded (quantified before #417)' if rows is None else '; '.join(
            '{}{} `{}`'.format(field, ' (fixed)' if source == 'fixed' else '', value) for field, value, source in rows)))
    text = 'route\tfield\tvalue\tsource\n' + ''.join('{}\t{}\t{}\t{}\n'.format(route, *row) for route, rows in settings.items() for row in rows or [])
    return lines, text


def publish(project, bench, spec, header, rows, lines, extra, record):
    """evaluation.publish with the run's reads processed, the declared routes' states (#448) and the methods' runtime
    and memory as the last report sections, with routes.tsv and resources.tsv."""
    reads = run_reads(project)
    text, resource_lines = run_resources(project, bench, spec, reads)
    states = spec['route_states']
    route_lines = ['', '## Routes', '', 'Every declared route (#448): `evaluated`; `not run`, not selected at `prepare` or recorded in '
                   '`benchmark/not_run.tsv`; `stopped`, a check before scoring failed, and the other routes are scored without it.',
                   ''] + md_table(['route', 'state', 'reason'], [list(s) for s in states])
    routes = '\n'.join('\t'.join(row) for row in [('route', 'state', 'reason'), *states]) + '\n'
    setting_lines, settings = settings_lines(spec.get('route_settings') or {}, spec.get('variant'))
    evaluation.publish(bench, header, rows, lines + route_lines + setting_lines + resource_lines, reads,
                       {**extra, 'routes.tsv': routes, 'settings.tsv': settings, 'resources.tsv': text}, record,
                       spec['synthetic'].get('seed') if spec['kind'] == 'synthetic' else None)  # the draw of a synthetic set (#455)


def nonnegative(text):
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError('must be 0 or more, got {}'.format(value))
    return value


def run_reads(project):
    """The reads the evaluated run processed, recorded with the evaluation (evaluation.reads_processed)."""
    try:
        return evaluation.reads_processed(Path(project))
    except ValueError as error:
        fail(str(error))


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, header, rows):
    Path(path).write_text('\t'.join(header) + '\n' + ''.join('\t'.join(str(v) for v in row) + '\n' for row in rows))


def scored_samples(spec, bench):
    """The samples scoring uses: the truth's (synthetic) or benchmark/samples.tsv's (pools, biological)."""
    if spec['kind'] == 'synthetic':
        return {row['sample'] for row in read_tsv(bench / 'truth.tsv')}
    return {row['sample'] for row in read_tsv(bench / 'samples.tsv')}


def check_methods(project, spec, config):
    """The routes that cannot be scored, {method: reason} (#448): a count table missing or without exactly the samples
    scoring uses (column order free; a sample missing would be scored as zero), or, for datasets declaring a miRBase
    release (pools, biological), a route that annotated with another miRNA database or release than the one the truth
    and names come from (#407). The evaluation scores the others without them. Stops when the project's configured
    samples differ from the benchmark's, which concerns every route."""
    renamed = [m for m in spec['methods'] if m in adapters.RENAMED]
    if renamed:
        fail('the project names {} in benchmark/dataset.yaml, since renamed {} (#434): rerun those routes under the new names, '
             'rename them in benchmark/dataset.yaml and evaluate again'.format(', '.join(renamed), ', '.join(adapters.RENAMED[m] for m in renamed)))
    quant = project / 'data/tmp/smallrna/quant'
    samples = scored_samples(spec, project / 'benchmark')
    configured = set(config.get('samples') or {})
    if configured != samples:
        fail('the project\'s configured samples differ from the benchmark\'s{}{}; prepare the project again'.format(
            '; not configured ' + ', '.join(sorted(samples - configured)) if samples - configured else '',
            '; not in the benchmark ' + ', '.join(sorted(configured - samples)) if configured - samples else ''))
    stopped = {}
    for method in spec['methods']:
        if method not in adapters.TABLES:
            fail('no sample and release check for method {}: add its count table to adapters.TABLES'.format(method))
        problems = []
        try:
            found = adapters.table_samples(method, quant)
        except FileNotFoundError as e:
            stopped[method] = 'count table missing ({}); run `snakemake quant_all` there'.format(Path(e.filename).relative_to(quant))
            continue
        if len(set(found)) != len(found) or set(found) != samples:
            missing, extra = sorted(samples - set(found)), sorted(set(found) - samples)
            problems.append('table samples differ from the benchmark\'s{}{}{}'.format(
                '; missing ' + ', '.join(missing) if missing else '', '; extra ' + ', '.join(extra) if extra else '',
                '; repeated columns' if len(set(found)) != len(found) else ''))
        if 'mirbase' in spec:
            try:
                database, release = adapters.method_mirbase(method, quant)
                recorded = (adapters.method_record(method, quant, ('table_sha256',))['table_sha256']
                            if method in adapters.TABLE_BOUND_RECORDS else None)
            except FileNotFoundError as e:
                stopped[method] = 'reference record missing ({}); run `snakemake quant_all` there'.format(Path(e.filename).relative_to(quant))
                continue
            except ValueError as e:
                stopped[method] = '{}; run `snakemake quant_all` there'.format(e)
                continue
            if recorded is not None and recorded != adapters.table_sha256(method, quant):
                problems.append('reference record was written for another count table; run `snakemake quant_all` there')
            expected = str(spec['mirbase']['release'])
            if database != 'mirbase' or release != expected:
                problems.append('annotated with {} {}, but the dataset declares miRBase {}'.format(database, release, expected))
        settings = quant / method / 'settings.tsv'  # the route's counting settings (#417), bound to its table like the record
        if settings.is_file():
            try:
                if adapters.method_record(method, quant, ('table_sha256',), 'settings.tsv')['table_sha256'] != adapters.table_sha256(method, quant):
                    problems.append('settings record was written for another count table; run `snakemake quant_all` there')
            except ValueError as e:
                problems.append(str(e))
        if problems:
            stopped[method] = '; '.join(problems)
    return stopped


ROUTE_STATES = ('evaluated', 'not run', 'stopped')


def route_states(bench, spec, stopped):
    """[(route, state, reason)] for every declared route (#448): `evaluated`; `not run`, not selected at `prepare`
    (the project's `not_selected`) or recorded by hand in benchmark/not_run.tsv (route, reason: e.g. a run that could
    not be done for resources); `stopped`, a check of check_methods failed. A route in not_run.tsv that the evaluation
    scores counts as evaluated."""
    declared = list(spec['methods']) + list(spec.get('not_selected') or [])
    by_hand = {}
    path = bench / 'not_run.tsv'
    if path.is_file():
        rows = read_tsv(path)
        if not rows or set(rows[0]) != {'route', 'reason'}:
            fail('{}: expected the columns route and reason, one row per route not run'.format(path))
        by_hand = {row['route']: row['reason'] for row in rows}
        unknown = sorted(set(by_hand) - set(declared))
        if unknown:
            fail('{}: {} not declared for this project'.format(path, ', '.join(unknown)))
    states = []
    for route in declared:
        if route in spec['methods'] and route not in stopped:
            states.append((route, 'evaluated', ''))
        elif route in by_hand:
            states.append((route, 'not run', by_hand[route]))
        elif route in stopped:
            states.append((route, 'stopped', stopped[route]))
        else:
            states.append((route, 'not run', 'not selected at prepare'))
    return states


def evaluate(args):
    project = Path(args.project)
    bench = project / 'benchmark'
    spec = load_dataset(bench / 'dataset.yaml')
    config = yaml.safe_load((project / 'config.yaml').read_text())
    ext_dir = Path(config['ext_dir'])
    stopped = check_methods(project, spec, config)
    states = route_states(bench, spec, stopped)
    scored = [route for route, state, _ in states if state == 'evaluated']
    if not scored:
        fail('no route can be scored: {}'.format('; '.join('{}: {}'.format(r, why) for r, _, why in states)))
    spec = dict(spec, methods=scored, route_states=states,  # the stopped routes are left out of every table
                route_settings={route: route_settings(project, route) for route in scored})
    check_variant(spec)
    for route, state, reason in states:  # visible to a script that reads only exit codes and the last lines
        if state == 'stopped':
            print('benchmark.py: {} stopped, scored without it: {}'.format(route, reason), file=sys.stderr)
    if args.de:
        check_de_results(project, bench, spec)
    run_resources(project, bench, spec, None)  # benchmark files checked before any scoring table is written
    if spec['kind'] == 'pool':
        return evaluate_pool(args, project, bench, spec, ext_dir)
    if spec['kind'] == 'biological':
        return evaluate_biological(args, project, bench, spec, ext_dir)
    if spec['kind'] == 'mixture':
        return evaluate_mixture(args, project, bench, spec, ext_dir)
    if synthetic_mode(spec) == 'genome':
        return evaluate_synthetic_genome(args, project, bench, spec, ext_dir)
    prefix = PREFIX[spec['organism']]
    mature = synthetic.read_mature(reference_dir(spec, ext_dir) / 'mature.fa', prefix)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    gff = args.mirflowz_gff or (mirflowz_annotation(project) if 'mirflowz' in spec['methods'] else None)
    common = common_table(args, project, bench, spec, ext_dir, spec['methods'], mature, gff)

    truth = {}
    for row in read_tsv(bench / 'truth.tsv'):
        if row['mature_id']:
            key = (row['sample'], row['mature_id'])
            truth[key] = truth.get(key, 0) + int(row['count'])
    methods = list(spec['methods'])
    summary = metrics.count_metrics(common, truth, methods)
    table = metrics.per_mirna(common, truth, methods)
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
    # Every sample gets the same reads (synthetic.py), so the unit is the miRNA with all its samples' truth pairs.
    by_mirna = {}
    for key in sorted(k for k, v in truth.items() if v > 0):
        by_mirna.setdefault(key[1], []).append(key)
    mirnas = sorted(by_mirna)
    extra, record, ci_lines = bootstrap(
        args, summary, lambda keys: metrics.truth_scores(table, truth, methods, keys),
        lambda rng: ([key for m in uncertainty.resample(rng, mirnas) for key in by_mirna[m]],), methods, ['recovered', 'exact'],
        'the truth miRNAs, each with all its samples (`misassigned` and `unplaced` are totals, without an interval)')
    publish(project, bench, spec, ['metric'] + methods, metric_rows, lines + ci_lines, extra, record)
    print(bench / 'report.md')


def check_de_results(project, bench, spec):
    """`evaluate --de` stops before writing anything unless the project's dataset can and does declare `de:` and de.smk
    has written the results of each DE set (de_sets)."""
    if de_kind(spec) is None:
        fail('--de (differential expression, #414) is for {}, not this {} dataset'.format(DE_KINDS, spec.get('pool', spec['kind'])))
    if 'de' not in spec:
        fail("{} declares no `de:`; add the dataset's `de:` (datasets/{}.yaml) to it, run `benchmark.py evaluate` and de.smk, "
             "then evaluate --de".format(bench / 'dataset.yaml', spec['name']))
    sample_rows = read_tsv(bench / 'samples.tsv')
    if de_kind(spec) == 'biological' and any(row['lane'] != 'merged' for row in sample_rows):
        fail('--de needs the samples\' merged runs: {} was prepared with --lanes, whose samples are lanes of one library'.format(project))
    sets = [s for s, *_ in mixture_de_sets(spec)] if de_kind(spec) == 'mixture' else de_sets(spec, de_units(spec, sample_rows))
    missing = [str(Path(s.folder) / name) for s in sets
               for name in ('results.tsv', 'record.tsv') if not (bench / 'de' / s.folder / name).is_file()]
    if missing:
        fail('no {} in {}: run `benchmark.py evaluate --project {}`, then `snakemake -s {} -d {} --use-apptainer`, then evaluate '
             '--de'.format(' or '.join(missing), bench / 'de', project, Path(__file__).resolve().with_name('de.smk'), project))


def unitas_refdump(project):
    """The UNITAS reference dump the project's UNITAS run read: the `dump` row of its reference record (written by the
    workflow), not today's configuration or defaults, which may name another dump since."""
    record = project / 'data/tmp/smallrna/quant/unitas/reference.tsv'
    if not record.is_file():
        fail('no UNITAS reference record {}: run quant_all again'.format(record))
    rows = {row['field']: row['value'] for row in read_tsv(record)}
    if rows.get('mirna_database', 'mirbase') != 'mirbase':
        fail('the UNITAS dump was built from {}: the benchmark matches miRBase names'.format(rows['mirna_database']))
    if not rows.get('dump'):
        fail('{} has no dump row: run quant_all again'.format(record))
    refdump = Path(rows['dump']) if Path(rows['dump']).is_absolute() else project / rows['dump']
    if not refdump.is_dir():
        fail('no UNITAS reference dump at {} (named by {}); name one with --unitas-ref'.format(refdump, record))
    return refdump


def common_table(args, project, bench, spec, ext_dir, methods, mature, gff, pending=None):
    """Every method's counts in one table (adapters.py), written to common.tsv; with UNITAS, also how each UNITAS name
    was placed (unitas_names.tsv: by name, by hairpin and arm, or not). With `pending`, the tables are added to it
    (path, header, rows) for the caller to write later."""
    write = write_tsv if pending is None else lambda *table: pending.append(table)
    prefix = PREFIX[spec['organism']]
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    hairpins = adapters.unitas_hairpins(args.unitas_ref or unitas_refdump(project)) if 'unitas' in methods else None
    common = adapters.common_table(methods, project / 'data/tmp/smallrna/quant', name_to_id, prefix, gff, hairpins, mature)
    write(bench / 'common.tsv', adapters.COMMON_COLUMNS, common)
    if hairpins is not None:
        reads = collections.Counter()
        for method, _, name, _, count in common:
            if method == 'unitas':
                reads[name] += float(count)
        placed = {name: adapters.unitas_place(name, name_to_id, prefix, hairpins, mature) for name in reads}
        write(bench / 'unitas_names.tsv', ['name', 'step', 'mature_id', 'reads'],
              [[n, placed[n][1] or 'unplaced', placed[n][0], '{:g}'.format(reads[n])] for n in sorted(reads)])
        steps = collections.defaultdict(lambda: [0, 0.0])
        for n, (_, step) in placed.items():
            steps[step or 'unplaced'][0] += 1
            steps[step or 'unplaced'][1] += reads[n]
        labels = {'name': 'by name', 'hairpin': 'by hairpin and arm', 'unplaced': 'unplaced'}
        print('UNITAS names: ' + ', '.join('{} {} ({:.0f} reads)'.format(steps[k][0], v, steps[k][1]) for k, v in labels.items() if k in steps))
    return common


def isomir_table(args, project, bench, spec, methods, mature, gff, common, pending=None):
    """The isomiR table of the methods that report read sequences (adapters.isomir_table, #454), checked against their
    counts in the common table and written to isomirs.tsv (with `pending`, added to it): (rows, methods without isomiR
    output)."""
    write = write_tsv if pending is None else lambda *table: pending.append(table)
    name_to_id = {name: accession for name, (accession, _) in mature.items()}
    hairpins = adapters.unitas_hairpins(args.unitas_ref or unitas_refdump(project)) if 'unitas' in methods else None
    try:
        rows, without = adapters.isomir_table(methods, project / 'data/tmp/smallrna/quant', name_to_id, PREFIX[spec['organism']], gff,
                                              hairpins, mature)
        adapters.check_isomirs(rows, common, [m for m in methods if m not in without])
    except (OSError, ValueError) as error:
        fail(str(error))
    write(bench / 'isomirs.tsv', adapters.ISOMIR_COLUMNS, sorted(rows))
    return rows, without


def mirflowz_annotation(project):
    """mirflowz's extended annotation that the project's mirflowz run counted with, from its quant/mirflowz/reference.tsv:
    `<reference>/gff-<mirna_gff_id>/annotation/` (#440); a record written before #440 has no mirna_gff_id, and its run
    read the annotation from the reference's `prepared/`."""
    quant = project / 'data/tmp/smallrna/quant'
    try:
        reference = Path(adapters.method_record('mirflowz', quant, ('reference',))['reference'])
    except (OSError, ValueError) as e:
        fail('{}; name the annotation with --mirflowz-gff'.format(e))
    ids = [row['value'] for row in read_tsv(quant / 'mirflowz/reference.tsv') if row['field'] == 'mirna_gff_id']
    name = 'extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION)
    path = reference / 'gff-{}'.format(ids[0]) / 'annotation' / name if ids else reference / 'prepared' / name
    if not path.is_file():
        fail('no mirflowz annotation at {}, where the project\'s mirflowz run read it; name one with --mirflowz-gff'.format(path))
    return path


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
    gff = args.mirflowz_gff or (mirflowz_annotation(project) if 'mirflowz' in methods else None)
    pending = []  # tables written once the differential-expression checks (evaluate --de) have passed
    common = common_table(args, project, bench, spec, ext_dir, methods, mature, gff, pending)
    iso_rows, no_isomirs = isomir_table(args, project, bench, spec, methods, mature, gff, common, pending)
    entries = pool.read_pool(bench / 'pool.tsv')
    kinds = pool.shift_kinds(entries, dict(mature.values()))
    sample_rows = read_tsv(bench / 'samples.tsv')
    sample_groups = {row['sample']: row['group'] for row in sample_rows}
    sequences = [e['sequence'] for e in entries if e['group'] != 'excluded']
    stages = {sample: filter_stages(project, sample) for sample in sample_groups}
    stage_counts = {sample: {stage: pool.oracle(path, sequences, *cuts) for stage, path, cuts in stages[sample]} for sample in sample_groups}
    stage_names = [stage for stage, _, _ in next(iter(stages.values()))]
    oracle = {sample: stage_counts[sample][stage_names[-1]] for sample in sample_groups}
    pending.append((bench / 'filter_stages.tsv', ['sample', 'pool_id', 'group'] + stage_names,
                    [[sample, e['pool_id'], e['group']] + [stage_counts[sample][st].get(e['sequence'], 0) for st in stage_names]
                     for sample in sample_groups for e in entries if e['group'] != 'excluded']))
    counts, reads = pool.method_counts(common, entries)
    pending.append((bench / 'per_sequence.tsv', pool.PER_SEQUENCE + methods, pool.per_sequence(counts, oracle, entries, methods, sample_groups, kinds)))
    rows = ['oracle'] + methods
    shifted = pool.scored(entries, 'shifted')
    by_kind = {kind: [e for e in shifted if kinds[e['pool_id']] == kind] for kind in pool.SHIFT_KINDS}
    by_kind = {kind: group for kind, group in by_kind.items() if group}
    de_texts = de_results = de_reports = table_keys = None
    de_lines = []
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
        summary, observed = pool.ratiometric_metrics(counts, oracle, entries, methods, list(sample_groups.items()))
        keys = ['median_abs_error', 'spearman_vs_design', 'called', 'false_calls']
        keys += ['shifted_' + k for k in keys] if shifted else []
        kind_keys = ['median_abs_error', 'called']
        per_kind = {kind: {r: pool.ratio_scores([(observed[r][e['pool_id']], math.log2(e['amount_a'] / e['amount_b']))
                                                 for e in group if not math.isnan(observed[r][e['pool_id']])])
                           for r in observed} for kind, group in by_kind.items()}
        if 'de' in spec:
            sets = de_sets(spec, de_units(spec, sample_rows))
            sizes = de.library_sizes(reads, oracle, entries, rows, list(sample_groups))
            count = de.pool_count(counts, oracle, entries)
            de_features = de.pool_features(pool.scored(entries))
            de_ids = [feature for feature, _ in de_features]
            de_way = de.direction(de_features)  # the full design's, for every bootstrap draw (#453)
            set_sizes = filtered_sizes(project, sets, rows, sizes)
            de_texts = {s.folder: de.input_tables(count, rows, de_ids, s.design, set_sizes[s.folder], s.model) for s in sets}
            de_results, de_reports = differential_expression(args, bench, sets, de_texts, rows, de_ids)
        if de_results:
            for suffix, results in de_results.items():
                for row, values in de.scores(results, de_features, rows, suffix).items():
                    summary.setdefault(row, {}).update(values)
            keys += de.metric_keys()
            table_keys = list(keys)  # the report's metrics table: each level's DE metrics have a section of their own
            keys += de.metric_keys([suffix for suffix in de_results if suffix])
            de_lines = de.report_lines(de_results[''], de_features, rows, *de_reports[''][:2], de.POOL_INTRO.format(len(de_features), de.FDR))
            de_lines += ['', rounding_line(count, methods, de_ids, list(sample_groups), 'sequences')]
            de_lines += de.level_lines((spec['de'] or {}).get('within'), [(s.suffix[1:], de_reports[s.suffix][2], de_reports[s.suffix][1],
                                                                          de_reports[s.suffix][0]) for s in sets if s.suffix and s.arm is None],
                                       summary, rows)
            de_lines += de.arm_lines([(s.suffix[1:], s.folder, de_reports[s.suffix][2], de_reports[s.suffix][1], de_reports[s.suffix][0])
                                      for s in sets if s.suffix == '' or s.arm], summary, rows)
    with_isomirs = [m for m in methods if m not in no_isomirs]
    iso_tallies, iso_kinds = isomirs.pool_tallies(iso_rows, pool.scored(entries), shifted, entries, oracle, with_isomirs, dict(mature.values()))
    iso = isomirs.pool_values(iso_tallies)
    for row, values in iso.items():
        summary.setdefault(row, {}).update(values)
    keys += isomirs.POOL_KEYS
    if table_keys is not None:
        table_keys += isomirs.POOL_KEYS
    for path, header, table in pending:
        write_tsv(path, header, table)
    if de_texts:
        write_de_inputs(bench, de_texts)
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
    lines += ['| {} | '.format(k) + ' | '.join(str(summary.get(r, {}).get(k, '')) for r in rows) + ' |' for k in table_keys or keys]
    if observed:
        designs = sorted({(e['amount_a'], e['amount_b']) for e in exact}, key=lambda ab: ab[0] / ab[1])
        lines += ['', '## Median observed log2(A/B) per designed ratio', '',
                  '| A:B | designed | n | ' + ' | '.join(r for r in rows if r in observed) + ' |', '|---' * (3 + len(observed)) + '|']
        for a, b in designs:
            ids = [e['pool_id'] for e in exact if (e['amount_a'], e['amount_b']) == (a, b)]
            lines.append('| {:g}:{:g} | {:.2f} | {} | '.format(a, b, math.log2(a / b), len(ids)) + ' | '.join(
                '{:.2f}'.format(statistics.median(observed[r][i] for i in ids)) for r in rows if r in observed) + ' |')
    lines += de_lines
    lines += shifted_section(per_kind, {kind: len(group) for kind, group in by_kind.items()}, kind_keys,
                             [r for r in rows if any(r in v for v in per_kind.values())])
    lines += isomirs.pool_lines(iso, iso_kinds, with_isomirs, no_isomirs)
    lines += filter_section(*pool.filter_losses(stage_counts, entries, stage_names), stage_names, trim_route(project, next(iter(sample_groups))),
                            next(iter(stages.values()))[0][2])
    if spec['pool'] == 'equimolar':
        def compute(ex, sh, samples):
            result = pool.equimolar_metrics(counts, reads, oracle, entries, methods, samples, exact=ex)
            if sh:
                for row, values in pool.agreement(counts, oracle, sh, methods, samples).items():
                    result[row].update({'shifted_' + k: v for k, v in values.items()})
            for row, values in isomirs.pool_values(iso_tallies, ex, sh, samples).items():  # the isomiR table's tallies (#459)
                result.setdefault(row, {}).update(values)
            return result
        draw = lambda rng: (uncertainty.resample(rng, exact), uncertainty.resample(rng, shifted), uncertainty.resample(rng, list(sample_groups)))
        units = 'the scored pool sequences (exact and shifted) and the runs'
    else:
        def compute(ex, sh, pairs):
            try:
                result = pool.ratiometric_metrics(counts, oracle, entries, methods, pairs, exact=ex, shifted=sh)[0]
            except ValueError:  # a draw without a 1:1 sequence: no ratio metrics, but the DE scores need none
                if not de_results:
                    raise
                result = {}
            for suffix, results in (de_results or {}).items():  # the DE fits stay fixed: their scores resample the sequences only
                for row, values in de.scores(results, de.pool_features(ex), rows, suffix, way=de_way).items():
                    result.setdefault(row, {}).update(values)
            for row, values in isomirs.pool_values(iso_tallies, ex, sh, [sample for sample, _ in pairs]).items():  # (#459)
                result.setdefault(row, {}).update(values)
            return result
        draw = lambda rng: (uncertainty.resample(rng, exact), uncertainty.resample(rng, shifted),
                            uncertainty.resample_groups(rng, list(sample_groups.items())))
        units = 'the scored pool sequences (exact and shifted) and the runs within each pool (A, B)'
        if de_results:
            units += '; the differential-expression metrics resample the sequences only, on the fixed fits'
    extra, record, ci_lines = bootstrap(args, summary, compute, draw, rows, keys, units, table_keys)
    publish(project, bench, spec, ['metric'] + rows, metric_rows, lines + ci_lines, extra, record)
    print(bench / 'report.md')


def de_units(spec, sample_rows):
    """[(sample, group, {column: value})] of a dataset that declares `de:`, in the order of benchmark/samples.tsv: a
    pool's runs and a biological dataset's samples with the fields of their declaration, the synthetic dataset's mixed
    samples with the true scores of its two factors (`nuisance`, a covariate of its adjusted set; `factor_group`, which
    gives the contrast)."""
    kind = de_kind(spec)
    if kind == 'synthetic':
        return [(row['sample'], row['group'], {'nuisance': row['factor_nuisance'], 'factor_group': row['factor_group']})
                for row in sample_rows if row['kind'] == 'mixed']
    if kind == 'pool':
        fields = {run['run']: run for run in spec['runs']}
        return [(row['sample'], row['group'], fields[row['sample']]) for row in sample_rows]
    fields = {s['sample']: s for s in spec['samples']}
    return [(row['sample'], row['group'], fields[row['library']]) for row in sample_rows]


def de_sets(spec, units):
    """The DE sets (DeSet) of a dataset that declares `de:` (#414, #444, #447), over the `units` (de_units) of the
    contrast's two groups: first the whole dataset's (suffix '', folder '' of benchmark/de), then the others, each with
    a metric suffix and a folder `<column>/<level>`:
    - with `within`, each level's (suffix `.<level>`, folder `<within>/<level>`), tested with the declared covariates
      other than `within`;
    - for the synthetic dataset, the one with the true nuisance scores as a numeric covariate (suffix `.adjusted`,
      folder `adjusted/nuisance`);
    - for pools and the synthetic dataset, the whole design again in each normalisation arm (#453, de.ARMS: suffix
      `.<arm>`, folder `arm/<arm>`; `fixed` gets the samples' filtered reads as library sizes from the caller, `tmm`
      a model with `normalisation: tmm`).
    design: [(sample, group, {covariate: value})]; model: de.model_text, with the least mean count a row tests (none for
    a pool's designed sequences; the synthetic and biological datasets' detection thresholds); samples: the set's
    samples per group, for the report."""
    kind = de_kind(spec)
    declared = spec['de'] or {}
    min_mean = {'pool': 0, 'synthetic': genome_metrics.DE_MIN, 'biological': biological.MIN_READS}[kind]
    covariates = list(declared.get('covariates') or [])
    within = declared.get('within')
    nulls = null_sets(declared['null'], units, covariates, covariates + ([within] if within and within not in covariates else []), min_mean) \
        if kind == 'biological' and 'null' in declared else []
    contrast = de_contrast(spec, units)
    if contrast is None:
        return nulls
    numerator, denominator = contrast
    units = [unit for unit in units if unit[1] in (numerator, denominator)]

    def de_set(suffix, folder, members, names, numeric=(), arm=None):
        sizes = collections.Counter(group for _, group, _ in members)
        if min(sizes[numerator], sizes[denominator]) < 2:  # e.g. a project prepared from another declaration
            fail('benchmark/samples.tsv has {} {} and {} {} sample(s) for the DE set {}; differential expression needs at least '
                 'two in each group: prepare the project again from its declaration'.format(
                     sizes[numerator], numerator, sizes[denominator], denominator, folder or 'of the whole dataset'))
        return DeSet(suffix, folder, [(sample, group, {c: fields[c] for c in names}) for sample, group, fields in members],
                     de.model_text(numerator, denominator, numeric, min_mean, 'tmm' if arm == 'tmm' else None),
                     '{} {} and {} {} {}'.format(sizes[numerator], numerator, sizes[denominator], denominator,
                                                 'runs' if kind == 'pool' else 'samples'), arm)

    whole = [] if kind == 'synthetic' else covariates
    arms = [de_set('.' + arm, 'arm/' + arm, units, whole, arm=arm) for arm in de.ARMS] if kind in ('pool', 'synthetic') else []
    if kind == 'synthetic':
        return [de_set('', '', units, []), de_set('.adjusted', 'adjusted/nuisance', units, ['nuisance'], ['nuisance'])] + arms
    sets = [de_set('', '', units, covariates)]
    if within:
        others = [c for c in covariates if c != within]
        sets += [de_set('.' + level, '{}/{}'.format(within, level), [unit for unit in units if unit[2][within] == level], others)
                 for level in sorted({fields[within] for _, _, fields in units})]
    return sets + arms + nulls


def null_sets(declared, units, covariates, strata, min_mean):
    """The DE sets of a biological dataset's null comparisons (#449): each split of de.null_splits, within a group of
    `units`, balanced on the `strata` columns (the covariates and `within`), as pseudo-groups x and y with the
    covariates that vary in the group; suffix `.null.<group>-<k>`, folder `null/<group>-<k>` (de.folder_name), k from 1,
    zero-padded."""
    splits = de.null_splits(units, declared['splits'], declared['seed'], strata)
    names = {group: de.folder_name(group) for group in splits}
    if len(set(names.values())) != len(names):
        fail('benchmark/samples.tsv: groups {} share a folder name for their null comparisons'.format(', '.join(sorted(names))))
    sets = []
    for group, halves in splits.items():
        in_group = [unit for unit in units if unit[1] == group]
        fields = {sample: f for sample, _, f in in_group}
        names_used = null_covariates(in_group, covariates)
        width = len(str(len(halves)))
        for k, (x, y) in enumerate(halves, 1):
            design = [(s, 'x' if s in x else 'y', {c: fields[s][c] for c in names_used}) for s, _, _ in in_group]
            level = '{}-{}'.format(names[group], str(k).zfill(width))
            sets.append(DeSet('.null.' + level, '{}/{}'.format(de.NULL, level), design, de.model_text('x', 'y', (), min_mean),
                              '{} x and {} y samples of {}'.format(len(x), len(y), group)))
    return sets


def differential_expression(args, bench, sets, texts, rows, features, root=None, command='evaluate'):
    """With --de, the results of each DE set (de_sets): ({suffix: results}, {suffix: (record, results SHA-256, the set's
    samples)}), after checking for every set that de.smk's results come from its inputs `texts[folder]` and the current
    de.R, and that the results table is the one de.R recorded, with a line for every row, test and feature; it verifies,
    parses and reports one read of each. Without --de: (None, None). The caller writes each set's inputs of de.smk to
    benchmark/de/<folder>/input/, and benchmark/de/levels.tsv, with its other tables (write_de_inputs). `root`: the
    sets' DE folder instead of benchmark/de, `command` the benchmark.py command that writes their inputs (`depth`,
    #414)."""
    if not args.de:
        return None, None
    root = bench / 'de' if root is None else root
    results, reports = {}, {}
    for de_set in sets:
        set_features = features[de_set.folder] if isinstance(features, dict) else features
        path, record_path = root / de_set.folder / 'results.tsv', root / de_set.folder / 'record.tsv'
        record = de.read_record(record_path)
        changed = de.stale(record, texts[de_set.folder])
        if changed:
            fail('{} comes from other inputs or another de.R ({}); nothing was written: run `benchmark.py {}` without --de '
                 '(it writes the current inputs), then de.smk, then {} --de'.format(path, ', '.join(changed), command, command))
        data = path.read_bytes()
        if de.md5(data) != record.get('md5_results'):
            fail('{} is not the table de.R recorded in {}; nothing was written: run de.smk again, then {} --de'.format(
                path, record_path, command))
        try:
            results[de_set.suffix] = de.read_results(data.decode(), rows, set_features, path)
        except ValueError as error:
            fail(str(error))
        reports[de_set.suffix] = (record, hashlib.sha256(data).hexdigest(), de_set.samples)
    return results, reports


def filtered_sizes(project, sets, rows, sizes):
    """{folder: library sizes} of each DE set: the given `sizes` ({(row, sample): size}), except for the `fixed` arm
    (#453), every row's library size the sample's filtered reads, from the cutadapt route's reports."""
    fixed = [s for s in sets if s.arm == 'fixed']
    if not fixed:
        return {s.folder: sizes for s in sets}
    samples = [sample for sample, _, _ in fixed[0].design]
    stages = project / 'data/tmp/smallrna/filter/cutadapt_stages'
    paths = {s: {k: stages / '{}.{}.json'.format(s, k) for k in biological.CUTADAPT_PASSES} for s in samples}
    missing = [str(p) for ps in paths.values() for p in ps.values() if not p.is_file()]
    if missing:
        fail('the `fixed` normalisation arm (#453) takes each sample\'s filtered reads from the cutadapt route\'s reports: {} '
             'missing; run `snakemake quant_all --notemp`'.format(', '.join(missing[:3]) + (' ...' if len(missing) > 3 else '')))
    try:
        kept = {s: float(biological.filtering(p)['kept']) for s, p in paths.items()}
    except ValueError as error:
        fail(str(error))
    by_row = {(row, s): kept[s] for row in rows for s in samples}
    return {s.folder: by_row if s.arm == 'fixed' else sizes for s in sets}


def biological_de_features(sets, methods, common):
    """{folder: features} of a biological dataset's DE sets (de_sets): the miRNAs any method reports in the set's
    samples in `common`; each null comparison its own split's, the contrast's sets (the whole and each level of
    `within`) the whole contrast's."""
    by = biological.counts_by(common)
    reported = lambda design: sorted({i for m in methods for sample, _, _ in design for i in by.get((m, sample), {})})
    nulls = [s for s in sets if s.folder.startswith(de.NULL + '/')]
    contrasts = [s for s in sets if s not in nulls]
    features = {s.folder: reported(s.design) for s in nulls}
    if contrasts:
        features.update({s.folder: reported(contrasts[0].design) for s in contrasts})
    return features


def biological_de_texts(sets, methods, common, features, sizes=None):
    """(count(row, sample, feature), {folder: de.input_tables}) of a biological dataset's DE sets on the `common` rows,
    over `features` ({folder: features}). A sample's library size is every miRNA read the method reported in it (all its
    rows of `common`, those without an accession too), as the workflow scales counts per million; `sizes` ({(method,
    sample): size}) replaces them (the `fixed` arm of #488)."""
    by = biological.counts_by(common)
    count = lambda row, sample, feature: by.get((row, sample), {}).get(feature, 0.0)
    if sizes is None:
        sizes = collections.Counter()
        for method, sample, _, _, value in common:
            sizes[(method, sample)] += float(value)
    return count, {s.folder: de.input_tables(count, methods, features[s.folder], s.design, sizes, s.model) for s in sets}


def write_de_inputs(bench, texts):
    """Each DE set's inputs of de.smk ({folder: input texts}) and benchmark/de/levels.tsv, the sets besides the whole
    dataset's."""
    for folder, inputs in texts.items():
        de.write_inputs(bench / 'de' / folder / 'input', inputs)
    de.write_changed(bench / 'de' / de.LEVELS, de.levels_text([folder.split('/') for folder in texts if folder]))


def rounding_line(count, rows, features, samples, what):
    """A report line: the share of each row's `what` with reads that have a fractional count, rounded for the tests."""
    shares = ['{} {}'.format(row, de.text(de.fractional_share(count, row, features, samples))) for row in rows]
    return ('Counts are rounded to whole numbers for the tests (a method that divides a read between names reports fractions); the '
            'share of its {} with reads that have a fractional count: {}.'.format(what, ', '.join(shares)))


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


GENOME_KEYS = ['recovered', 'precision', 'exact', 'misassigned', 'unplaced', 'near_spill', 'deviation_correlation', 'de_tested',
               'de_sensitivity', 'de_fdr', 'de_fc_error', 'de_sensitivity_adjusted', 'de_fdr_adjusted', 'de_fc_error_adjusted',
               'structure_group_r2', 'structure_nuisance_r2', 'coexpression']
GENOME_CI_KEYS = ['recovered', 'precision', 'exact', 'misassigned', 'near_spill', 'deviation_correlation', 'de_sensitivity', 'de_fdr',
                  'de_fc_error', 'de_sensitivity_adjusted', 'de_fdr_adjusted', 'de_fc_error_adjusted']


def evaluate_synthetic_genome(args, project, bench, spec, ext_dir):
    """The genome-wide synthetic dataset (#412): genome_metrics on the collapsed common table against truth.tsv, with
    the design tables samples.tsv and mirnas.tsv written by prepare."""
    prefix = PREFIX[spec['organism']]
    mature_fa = synthetic_files(spec, ext_dir)[0]
    missing = [str(p) for p in [mature_fa] + [bench / n for n in ('truth.tsv', 'samples.tsv', 'mirnas.tsv')] if not p.is_file()]
    if missing:
        fail('missing: {}'.format(', '.join(missing)))
    mature = synthetic.read_mature(mature_fa, prefix)
    id_to_seq = {accession: seq for accession, seq in mature.values()}
    id_to_name = {accession: name for name, (accession, _) in mature.items()}
    methods = list(spec['methods'])
    gff = args.mirflowz_gff or (mirflowz_annotation(project) if 'mirflowz' in methods else None)
    pending = []  # tables written once the differential-expression checks (evaluate --de) have passed
    raw = common_table(args, project, bench, spec, ext_dir, methods, mature, gff, pending)
    common, _ = biological.collapse(raw, id_to_seq)
    truth, samples, mirnas = (read_tsv(bench / n) for n in ('truth.tsv', 'samples.tsv', 'mirnas.tsv'))
    ev = genome_metrics.Evaluation(common, truth, samples, mirnas, methods)
    label = lambda i: id_to_name.get(i) or ev.mirnas.get(i, {}).get('mature_name') or i
    rows = ['truth'] + methods
    point = {name: ev.mixed_metrics(name, ev.units, ev.pairs, ev.de_units, ev.dev_units) for name in rows}
    probes = {name: ev.probes(name) for name in rows}
    classes = sorted({c for v in probes.values() for c in v})
    for name in rows:
        point[name]['unplaced'] = ev.unplaced_share(name) if name != 'truth' else 0.0
        point[name]['structure_group_r2'], point[name]['structure_nuisance_r2'], point[name]['coexpression'] = ev.structure(name)
        for c in classes:
            point[name]['probe_' + c] = probes[name].get(c, (float('nan'),))[0]
    keys = GENOME_KEYS + ['probe_' + c for c in classes]
    decoy_scores = {name: ev.decoys(name) for name in rows}  # decoy reads counted as miRNAs (#412 PR 2)
    decoy_units = sorted({u for v in decoy_scores.values() for u in v})
    if decoy_units:
        for name in rows:
            scored = decoy_scores[name]
            reads = sum(scored[u][1] for u in decoy_units)
            point[name]['decoy_as_mirna'] = sum(scored[u][0] * scored[u][1] for u in decoy_units) / reads if reads else float('nan')
            point[name].update({'decoy_as_mirna_' + u: scored[u][0] for u in decoy_units})
        keys += ['decoy_as_mirna'] + ['decoy_as_mirna_' + u for u in decoy_units]
    overlap_scores = {name: ev.overlaps(name) for name in rows}  # fragments overlapping an embedded mature (#412 PR 2)
    overlap_classes = sorted({c for v in overlap_scores.values() for c in v})
    for name in rows:
        for c in overlap_classes:
            point[name]['overlap_on_mirna_' + c], point[name]['overlap_as_mirna_' + c] = overlap_scores[name][c][:2]
    keys += [k + c for c in overlap_classes for k in ('overlap_on_mirna_', 'overlap_as_mirna_')]
    iso_rows, no_isomirs = isomir_table(args, project, bench, spec, methods, mature, gff, raw, pending)
    with_isomirs = [m for m in methods if m not in no_isomirs]
    iso_tallies, iso_classes = isomirs.synthetic_tallies(truth, iso_rows, samples, with_isomirs, isomirs.representatives(id_to_seq))
    iso = isomirs.synthetic_values(iso_tallies, iso_classes)
    for name, values in iso.items():
        point[name].update(values)
    keys += isomirs.synthetic_keys(iso_classes)
    de_texts = de_results = None
    de_keys, de_lines, de_features = [], [], []
    de_way = 0
    if 'de' in spec:  # the tests of the ratiometric pools (de.R) on the mixed samples (#447)
        units = de_units(spec, samples)
        sets = de_sets(spec, units)
        column = {sample: k for k, sample in enumerate(ev.mixed)}
        count = lambda row, sample, unit: float(ev.counts(row)[ev.index[unit], column[sample]])
        sizes = {(row, sample): float(ev.t[:, k].sum()) if row == 'truth' else ev.total[(row, sample)]
                 for row in rows for sample, k in column.items()}
        set_sizes = filtered_sizes(project, sets, rows, sizes)
        de_texts = {s.folder: de.input_tables(count, rows, ev.units, s.design, set_sizes[s.folder], s.model) for s in sets}
        de_results, de_reports = differential_expression(args, bench, sets, de_texts, rows, ev.units)
    if de_results:
        de_features = ev.de_features()
        de_way = de.direction(de_features)  # the full design's, for every bootstrap draw (#453)
        for suffix, results in de_results.items():
            for name, values in de.scores(results, de_features, rows, suffix, 'truth').items():
                point[name].update(values)
        de_keys = de.metric_keys(list(de_results), 'truth')
        numerator, denominator = de_contrast(spec, units)
        changed = sum(1 for _, d in de_features if d)
        intro = de.SYNTHETIC_INTRO.format(numerator=numerator, denominator=denominator, features=len(de_features), changed=changed,
                                          min_mean=genome_metrics.DE_MIN, fdr=de.FDR)
        if not any(m.get('de') == '1' for m in ev.mirnas.values()):  # a null scenario (#453): the design, not only the scored
            intro += (' No miRNA is designed to change here, so every call is false: the false-call rate under the generator\'s '
                      'variance alone, not the workflow\'s on real samples (for that, the null comparisons of biological datasets, #449).')
        de_lines = de.report_lines(de_results[''], de_features, rows, *de_reports[''][:2], intro)
        de_lines += ['', rounding_line(count, methods, ev.units, ev.mixed, 'miRNAs')]
        de_lines += de.set_lines(
            'Adjusted for the nuisance factor',
            'The same samples with the nuisance factor\'s true sample scores as a covariate of the GLM test and voom (the exact test '
            'takes no covariates, so its values equal those above), scored as above. In metrics.tsv the metrics carry `.adjusted`, '
            'e.g. `voom_called.adjusted`; their intervals and the differences between methods are in metrics_ci.tsv and '
            'differences.tsv, not in the uncertainty section below. The metrics against the truth are empty for the truth.',
            'set', [(s.suffix[1:], s.folder, de_reports[s.suffix][2], de_reports[s.suffix][1], de_reports[s.suffix][0])
                    for s in sets if s.suffix and s.arm is None], point, rows, 'truth')
        de_lines += de.arm_lines([(s.suffix[1:], s.folder, de_reports[s.suffix][2], de_reports[s.suffix][1], de_reports[s.suffix][0])
                                  for s in sets if s.suffix == '' or s.arm], point, rows, 'truth')
    table_keys = keys + de.metric_keys(reference='truth') if de_results else keys  # the report's table: other sets have their section
    keys = keys + de_keys
    rnd = lambda v: round(v, 4) if isinstance(v, float) and not math.isnan(v) else ('' if isinstance(v, float) else v)
    metric_rows = [[k] + [rnd(point[name].get(k, float('nan'))) for name in rows] for k in keys]
    per = []
    for i in ev.units:
        k = ev.index[i]
        per.append([label(i), i, rnd(float(ev.t[k].sum()))] + [rnd(float(ev.o[name][k].sum())) for name in methods])
    for path, header, table in pending + [(bench / 'per_mirna.tsv', ['mature_name', 'mature_id', 'truth'] + methods, per)]:
        write_tsv(path, header, table)
    if de_texts:
        write_de_inputs(bench, de_texts)
    mixed = len(ev.mixed)
    lines = ['# {}'.format(spec['name']), '',
             '{} ({}), genome-wide synthetic reads (#412): {} mixed samples in two groups, {} probe samples; {} truth miRNAs, '
             'identical sequences counted as one. `truth` is the truth counts analysed as a method: the best the design allows.'.format(
                 spec['organism'], spec['libprepkit'], mixed, len(samples) - mixed,
                 sum(1 for m in mirnas if float(m['base_per_million']) > 0)), '',
             '## Metrics', '', 'Mixed samples (definitions in genome_metrics.py). DE, the built-in check (`de_*`): each row tests '
             'the miRNAs with at least {} reads on average in its own counts, as an analysis of them would (`de_tested`); {} miRNAs '
             'are truly DE, and a call on any other miRNA is false. {} `probe_*`: probe samples, recovery per read class.'.format(
                 genome_metrics.DE_MIN, len(ev.true_de),
                 'DE with the tests of the ratiometric pools (`edger_*`, `edger_glm_*`, `voom_*`): the section below.' if de_results else
                 'DE with the tests of the ratiometric pools: run de.smk, then `evaluate --de` (#447).' if 'de' in spec else ''), '']
    lines += md_table(['metric'] + rows, [row for row in metric_rows if row[0] in table_keys])
    strata = {name: ev.strata(name) for name in rows}
    names = sorted({k for v in strata.values() for k in v})
    lines += ['', '## Recovery by stratum', ''] + md_table(
        ['stratum', 'miRNAs'] + rows, [[n, strata['truth'].get(n, ('', ''))[1]] + [rnd(strata[r].get(n, (float('nan'),))[0]) for r in rows]
                                       for n in names])
    lines += ['', '## Probe samples', '', 'Per read class: the probe miRNAs\' observed reads (control subtracted) over their truth, '
              'and reads landing on other miRNAs (`misplaced`, relative to the same truth).', ''] + md_table(
        ['class', 'miRNAs'] + ['{}: recovery / misplaced'.format(r) for r in methods],
        [[c, probes['truth'][c][2]] + ['{} / {}'.format(rnd(probes[r][c][0]), rnd(probes[r][c][1])) for r in methods] for c in classes])
    if decoy_units:
        lines += ['', '## Decoys counted as miRNAs', '',
                  'Probe samples of one class and source of other RNA (#412 PR 2; `P_decoy_*`): the reads a method counts on '
                  'miRNAs other than the carrier, minus the control\'s, over the decoy reads (`decoy_as_mirna_<class>-<source>`; '
                  '`decoy_as_mirna` pools them). Decoys never share {} nt with a miRBase hairpin, so every such read is false.'.format(
                      decoys.OVERLAP), '']
        lines += md_table(['class-source', 'decoy reads'] + methods,
                          [[u, '{:.0f}'.format(decoy_scores['truth'][u][1])] + [rnd(decoy_scores[r][u][0]) for r in methods] for u in decoy_units])
        landed_lines = []
        for name in methods:
            total = collections.Counter()
            by = collections.defaultdict(collections.Counter)
            for u in decoy_units:
                for i, v in decoy_scores[name][u][2].items():
                    total[i] += v
                    by[i][u] += v
            top = [(i, v) for i, v in total.most_common(5) if v >= 1]
            if top:
                landed_lines.append('- {}: '.format(name) + '; '.join('{} {:.0f} ({})'.format(
                    label(i), v, ', '.join('{} {:.0f}'.format(u, c) for u, c in by[i].most_common(3))) for i, v in top))
        if landed_lines:
            lines += ['', 'The miRNAs with the most decoy reads (control subtracted), by class and source:', ''] + landed_lines
    embedded = [c for c in classes if c.startswith('embedded_')]
    if embedded:
        lines += ['', '## miRNAs inside other RNA classes', '',
                  'miRBase matures that lie exactly inside a molecule of another class (#412 PR 2; `P_embedded_*`, canonical reads): '
                  'recovery as in the probe samples (`probe_embedded_<class>`); 1 - recovery is the share lost to the other class.', '']
        inside = {i: m.get('inside', '') for i, m in ev.mirnas.items()}
        lines += md_table(['class', 'miRNAs'] + methods, [[c[len('embedded_'):], ', '.join(sorted(
            label(i) for i, m in ev.mirnas.items() if c in m['probe_classes'].split(',')))] + [rnd(probes[r][c][0]) for r in methods]
            for c in embedded])
        if overlap_classes:
            lines += ['', 'Fragments of their host molecules that cover at least half of an embedded mature but are no read a '
                      'miRNA could give (`P_overlap_*`; not false positives, as miRBase has the mature): the share of the fragment '
                      'reads a method counts on the mature they overlap (`overlap_on_mirna_<class>`) and on any miRNA '
                      '(`overlap_as_mirna_<class>`), control subtracted.', '']
            lines += md_table(['class', 'fragment reads'] + ['{}: on its miRNA / any'.format(r) for r in methods],
                              [[c, '{:.0f}'.format(overlap_scores['truth'][c][2])] + ['{} / {}'.format(rnd(overlap_scores[r][c][0]),
                               rnd(overlap_scores[r][c][1])) for r in methods] for c in overlap_classes])
        lines += ['', 'Their molecules: ' + '; '.join('{} ({})'.format(label(i), v) for i, v in sorted(inside.items(), key=lambda kv: label(kv[0]))
                                                       if v) + '.']
    lines += isomirs.synthetic_lines(iso, iso_classes, with_isomirs, no_isomirs)
    absent = []
    for name in methods:
        k = [(float(ev.o[name][ev.index[i]].sum()), i) for i in ev.units if float(ev.t[ev.index[i]].sum()) == 0]
        top = sorted(k, reverse=True)[:5]
        if top and top[0][0] > 0:
            absent.append('- {}: '.format(name) + ', '.join('{} ({:.0f})'.format(label(i), v) for v, i in top if v > 0))
    if absent:
        lines += ['', '## Reads on miRNAs absent from the truth (top 5, mixed samples)', ''] + absent
    false = [(name, ev.false_calls(name)) for name in rows]
    if any(calls for _, calls in false):
        lines += ['', '## False DE calls', '', 'Calls on miRNAs not designed DE: reads per mixed sample in the truth and in the '
                  'row\'s counts, and the estimated log2 fold change.', '']
        lines += ['- {}: '.format(name) + ', '.join('{} ({:.1f} → {:.1f} reads, {:+.2f})'.format(label(i), tm, om, est)
                                                    for i, tm, om, est in calls) for name, calls in false if calls]
    lines += de_lines

    def compute(units, pairs, de_units, dev_units, scored):
        result = {name: ev.mixed_metrics(name, units, pairs, de_units, dev_units) for name in rows}
        for name, values in isomirs.synthetic_values(iso_tallies, iso_classes, units).items():  # the isomiR table's tallies (#459)
            result[name].update(values)
        for suffix, results in (de_results or {}).items():  # the DE fits stay fixed: their scores resample the scored miRNAs
            for name, values in de.scores(results, scored, rows, suffix, 'truth', way=de_way).items():
                result[name].update(values)
        return result
    de_rng = random.Random(args.seed)  # its own generator: the other metrics' resamples are the same with and without --de
    extra, record, ci_lines = bootstrap(
        args, {name: {k: rnd(v) for k, v in point[name].items()} for name in rows}, compute,
        lambda rng: (uncertainty.resample(rng, ev.units), uncertainty.resample(rng, ev.pairs), uncertainty.resample(rng, ev.de_units),
                     uncertainty.resample(rng, ev.dev_units), uncertainty.resample(de_rng, de_features) if de_features else []),
        rows, GENOME_CI_KEYS + isomirs.synthetic_keys(iso_classes) + de_keys,
        'the miRNAs (each with all mixed samples), the near-identical pairs, the miRNAs any row tests '
        'for DE and the deviation miRNAs (DE calls made once; unplaced, structure, coexpression and probes without an interval){}'.format(
            '; the metrics of the pools\' tests resample the scored miRNAs, on the fixed fits' if de_results else ''),
        GENOME_CI_KEYS + isomirs.synthetic_keys(iso_classes) + [k for k in de_keys if k in table_keys])
    publish(project, bench, spec, ['metric'] + rows, metric_rows, lines + ci_lines, extra, record)
    print(bench / 'report.md')


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
    knockout = 'knockout' in spec and not lane_mode  # by library; lane samples are not compared (#452)
    if knockout and not filt:
        fail('`knockout` scores reads per million filtered reads, from the cutadapt route\'s reports; none found')

    gff = args.mirflowz_gff or (mirflowz_annotation(project) if 'mirflowz' in methods else None)
    pending = []  # tables written once the differential-expression checks (evaluate --de) have passed
    raw = common_table(args, project, bench, spec, ext_dir, methods, mature, gff, pending)
    iso_rows, no_isomirs = isomir_table(args, project, bench, spec, methods, mature, gff, raw, pending)
    common, groups = biological.collapse(raw, id_to_seq)
    unplaced = biological.unplaced(raw, methods, samples)
    de_texts, de_lines, null_values = None, [], {}
    if 'de' in spec and not lane_mode:  # the tests of the ratiometric pools (de.R) between two groups of samples (#447); not scored
        sets = de_sets(spec, de_units(spec, table))
        nulls = [s for s in sets if s.folder.startswith(de.NULL + '/')]  # null comparisons within a group (#449)
        contrasts = [s for s in sets if s not in nulls]
        features = biological_de_features(sets, methods, common)
        count, de_texts = biological_de_texts(sets, methods, common, features)
        de_results, de_reports = differential_expression(args, bench, sets, de_texts, methods, features)
        if de_results and contrasts:
            tested = [sample for sample, _, _ in contrasts[0].design]
            de_lines = de.call_lines([(s.suffix[1:], s.folder, s.samples, de_results[s.suffix], *de_reports[s.suffix][:2]) for s in contrasts],
                                     features[contrasts[0].folder], methods, *de_contrast(spec))
            de_lines += ['', rounding_line(count, methods, features[contrasts[0].folder], tested, 'miRNAs')]
        if de_results and nulls:
            null_values, spread, called = de.null_scores([(de_results[s.suffix], features[s.folder]) for s in nulls], methods)
            group_of = {row['sample']: row['group'] for row in table}
            splits = collections.Counter(group_of[s.design[0][0]] for s in nulls)
            in_splits = sorted({sample for s in nulls for sample, _, _ in s.design})
            de_lines += de.null_lines(null_values, spread, called, methods,
                                      [(g, sum(1 for r in table if r['group'] == g), splits[g]) for g in sorted(splits)],
                                      lambda feature: id_to_name.get(feature, feature),
                                      lambda row, feature: statistics.mean(count(row, sample, feature) for sample in in_splits))
    for path, header, rows in pending:
        write_tsv(path, header, rows)
    if de_texts:
        write_de_inputs(bench, de_texts)

    flag, molecules = None, {}
    if 'other_rna' in spec:  # the miRBase entries inside molecules of other RNA classes (#415)
        molecules = other_rna_molecules(spec, ext_dir)
        representative = {}
        for accession, seq in sorted(id_to_seq.items()):
            representative.setdefault(seq, accession)  # the first sorted accession, as biological.collapse
        found = decoys.embedded(molecules, {a: seq for seq, a in representative.items()})
        flag = collections.defaultdict(list)
        for cls, keys in sorted(found.items()):
            for accession, names in keys.items():
                flag[accession].append((cls, names))
        write_tsv(bench / 'other_rna.tsv', ['mature_id', 'mirna', 'class', 'molecules'],
                  [[a, id_to_name.get(a, a), cls, ','.join(names)] for a in sorted(flag) for cls, names in flag[a]])
    summary = biological.method_summary(common, methods, samples, frozenset(flag) if flag is not None else None)
    for m in methods:
        summary[m]['unplaced'] = unplaced[m]['share']
    keys = ['mirna_reads', 'unplaced', 'detected', 'detected_all_samples', 'top10_share']
    other_lines = []
    if flag is not None:
        other_lines = other_rna_section(flag, iso_rows, methods, no_isomirs, summary, common, molecules, id_to_name)
        keys += OTHER_RNA_KEYS
    if null_values:  # the null comparisons' false calls (#449), with their spread over the splits in the report
        for m in methods:
            summary[m].update(null_values[m])
        keys += de.null_keys()
    with_isomirs = [m for m in methods if m not in no_isomirs]
    iso_tallies, iso_detected, iso_kinds, iso_overlap = isomirs.biological_tallies(iso_rows, with_isomirs, id_to_seq,
                                                                                   isomirs.representatives(id_to_seq))
    iso = isomirs.biological_values(iso_tallies)
    for m, count in iso_detected.items():
        iso[m]['isomirs_detected'] = count
    for m in methods:
        summary[m].update(iso.get(m, dict.fromkeys(isomirs.BIOLOGICAL_KEYS, '')))
    keys += isomirs.BIOLOGICAL_KEYS
    metric_rows = [[k] + [summary[m][k] for m in methods] for k in keys]
    knockout_lines = []
    if knockout:
        expected = None
        if spec['knockout'].get('expected'):  # by the accession a method's reads count under, as biological.collapse
            representative, reasons = isomirs.representatives(id_to_seq), collections.defaultdict(set)
            listed = knockout_expected.read_list(bench / spec['knockout']['expected'])
            for accession, (_, reason, _) in listed.items():
                if accession in representative:  # a mature of another release is left out
                    reasons[representative[accession]].add(reason)
            for accession in flag or {}:
                reasons[accession].add(OTHER_CLASS)
            expected = {accession: '; '.join(sorted(r)) for accession, r in reasons.items()}
            in_release = (sum(accession in representative for accession in listed), len(listed))
        ko_rows, knockout_lines = knockout_section(bench, spec, table, common, methods, filt, id_to_name, flag, expected,
                                                   in_release if expected is not None else None)
        metric_rows += ko_rows
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
        ['metric'] + methods, metric_rows[:len(keys)])  # the knockout rows have their own section
    named = [(m, n, c) for m in methods for n, c in unplaced[m]['names']]
    if named:
        lines += ['', 'Most reads without an accession: ' + '; '.join(
            '{}: {}'.format(m, ', '.join('`{}` ({:.0f})'.format(n, c) for mm, n, c in named if mm == m)) for m in methods if unplaced[m]['names'])]
    lines += ['', '## Agreement between methods (medians over samples)', ''] + md_table(
        ['methods'] + pair_keys, [['{} / {}'.format(a, b)] + [v.get(k, '') for k in pair_keys] for (a, b), v in pairs.items()])
    lines += isomirs.biological_lines(iso, iso_kinds, iso_overlap, with_isomirs, no_isomirs)
    lines += other_lines
    lines += knockout_lines

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

    lines += de_lines
    stats = biological.sample_stats(common, raw, methods, samples)
    # The unit is the library: the lanes of one library are technical replicates, so a drawn library brings all its
    # lanes (a merged sample is its own library). detected_all_samples has no interval: a resample repeating a sample
    # has fewer distinct samples, so it is biased up.
    lanes_of = {}
    for row in table:
        lanes_of.setdefault(row['library'], []).append(row['sample'])
    libraries = sorted(lanes_of)
    def compute(drawn):
        result = biological.summary_of(stats, methods, drawn, all_samples=False, flagged=frozenset(flag) if flag is not None else None)
        for m, values in isomirs.biological_values(iso_tallies, drawn).items():  # the isomiR table's tallies (#459)
            result[m].update(values)
        return result
    extra, record, ci_lines = bootstrap(args, summary, compute,
                                        lambda rng: ([s for lib in uncertainty.resample(rng, libraries) for s in lanes_of[lib]],), methods,
                                        [k for k in keys if k not in ('detected_all_samples', 'isomirs_detected') + tuple(de.null_keys()) +
                                         tuple(k for k in OTHER_RNA_KEYS if k != 'detected_unflagged')],
                                        '{} (`detected_all_samples` and `isomirs_detected` without an interval: a resample repeating a '
                                        'sample has fewer distinct samples)'.format('the libraries, each with all its lanes' if lane_mode else 'the samples'))
    publish(project, bench, spec, ['metric'] + methods, metric_rows, lines + ci_lines, extra, record)
    print(bench / 'report.md')


KNOCKOUT_LISTED = 50  # surviving miRNAs shown in the report; knockout.tsv has every tested one
EXPECTED_NAMED = 10   # #478: expected miRNAs named per method in the report when tested but not surviving
OTHER_CLASS = 'inside another RNA class (#415)'  # a reason a miRNA is expected to survive a knockout
OTHER_RNA_KEYS = ['flagged_reads', 'flagged_detected', 'detected_unflagged', 'read_in_other_class']  # #415


def knockout_section(bench, spec, table, common, methods, filt, id_to_name, flag=None, expected=None, in_release=None):
    """The knockout comparison of a biological dataset (#452): its metrics.tsv rows and report lines; writes
    benchmark/knockout.tsv (every method, knockout group and tested miRNA). With `flag` (#415, {accession: [(class,
    molecules)]}), the surviving miRNAs inside molecules of other RNA classes too. With `expected` ({accession:
    reason}, #478: `knockout.expected` and the flagged miRNAs), the miRNAs expected to survive: knockout.tsv gives each one's
    reason, and each method's surviving miRNAs are split by it; `in_release`: (the list's entries in the project's miRBase
    release, its entries)."""
    declared = spec['knockout']
    control, survival = str(declared['control']), str(declared['survival'])
    units = knockout_units(spec, table)
    result = biological.knockout(common, methods, {s: filt[s]['kept'] for s, _, _ in units}, units, control)
    per, common_set = biological.survivors(result, methods, survival)
    groups = []
    for _, group, _ in units:
        if group != control and group not in groups:
            groups.append(group)
    size = collections.Counter(group for _, group, _ in units)
    rounded = lambda v, n=4: round(v, n) if v != '' else ''
    rows = [['knockout_mirna_share'] + [rounded(result[m]['share'][survival]) for m in methods],
            ['knockout_surviving'] + [len(per[m]['surviving']) for m in methods]]
    other = lambda i: ';'.join('{}:{}'.format(cls, ','.join(names[:3])) for cls, names in (flag or {}).get(i, []))
    tsv = [[m, g, i, id_to_name.get(i, i), round(c, 2), round(k, 2), round(r, 4), int(r > biological.SURVIVING)] +
           ([other(i)] if flag is not None else []) + ([expected.get(i, '') if g == survival else ''] if expected is not None else [])
           for m in methods for g in groups for i, (c, k, r) in result[m]['ratios'][g].items()]
    write_tsv(bench / 'knockout.tsv', ['method', 'group', 'mature_id', 'mirna', 'control_reads', 'knockout_reads', 'ratio', 'surviving'] +
              (['other_class'] if flag is not None else []) + (['expected'] if expected is not None else []), tsv)
    if flag is not None:
        rows.append(['knockout_surviving_flagged'] + [len(per[m]['surviving'] & set(flag)) for m in methods])
    split = biological.expected_split(result, per, methods, survival, set(expected)) if expected is not None else None
    if split:
        rows += [['knockout_surviving_expected'] + [len(split[m]['surviving_expected']) for m in methods],
                 ['knockout_surviving_unexpected'] + [len(per[m]['surviving'] - set(expected)) for m in methods]]  # whichever methods run
    within = declared.get('within')
    lines = ['', '## Knockouts (#452)', '',
             'Canonical miRNAs need DROSHA and DICER, so in those knockouts their reads should fall to near zero. The miRNA reads a',
             'method still counts there come from other RNAs it takes for miRNAs, or from miRNAs made by another pathway (e.g.',
             'mirtrons, without DROSHA): the surviving miRNAs common to every method are biology or annotation, those that differ',
             'are quantification. Reads per million filtered reads (those the cutadapt route kept), not per million miRNA reads,',
             'which the knockout changes; each knockout sample against the {} samples{}.'.format(
                 control, ' of its {}'.format(within) if within else ''), '',
             '**miRNA reads relative to {}** (first row: their miRNA reads per million filtered reads)'.format(control), ''] + md_table(
        ['group', 'samples'] + methods,
        [[control, size[control]] + [round(result[m]['control_rpm']) for m in methods]] +
        [[g, size[g]] + [rounded(result[m]['share'][g]) for m in methods] for g in groups])
    # Shares of the filtered reads: with the miRNAs gone, an RNA whose amount did not change takes a larger share
    inflation = {m: (1e6 - result[m]['control_rpm'] * result[m]['share'][survival]) / (1e6 - result[m]['control_rpm'])
                 for m in methods if result[m]['share'][survival] != '' and result[m]['control_rpm'] < 1e6}
    if inflation:
        lines += ['', 'The ratios are upper bounds: they are shares of the filtered reads, and with fewer miRNA reads the other RNAs take',
                  'a larger share, so an RNA whose amount did not change rises in {} by {:.1f}-{:.1f}x (the filtered reads not on miRNAs, knockout'.format(
                      survival, min(inflation.values()), max(inflation.values())),
                  'over controls). A miRNA kept at 25% has fallen to about {:.0%} of its amount; the routes share the denominator, so'.format(
                      0.25 / max(inflation.values())),
                  'they compare.']
    lines += ['', '**miRNAs surviving {}**: miRNAs with at least {} reads on average in the {} samples are tested; one survives'.format(
                  survival, biological.KNOCKOUT_MIN_READS, control),
              'when the knockout keeps more than {:.0%} of its reads per million filtered reads. {} survive in every method{}.'.format(
                  biological.SURVIVING, len(common_set),
                  ' ({})'.format(', '.join(sorted(id_to_name.get(i, i) for i in common_set))) if 0 < len(common_set) <= 20 else ''),
              '`only`: the surviving miRNAs no other method lists, and their reads (mean per {} sample).'.format(survival)]
    if split:
        reasons = collections.Counter(r for value in expected.values() for r in value.split('; '))
        lines += ['Expected to survive (`knockout.expected`, #478): {} miRNAs ({}); {} of the {} listed are in this miRBase release.'.format(
                      len(expected), ', '.join('{} {}'.format(n, r) for r, n in sorted(reasons.items(), key=lambda kv: (-kv[1], kv[0]))),
                      *in_release),
                  '`only` is split into the expected, those another method tests without listing them (it measures them lower) and',
                  'those no other method tests; `only_unexpected_reads`: the reads of the last two.']
    split_rows = [] if not split else [
        ['expected tested'] + [len(split[m]['tested_expected']) for m in methods],
        ['expected surviving'] + [len(split[m]['surviving_expected']) for m in methods]]
    only_rows = [] if not split else [
        ['only: expected'] + [len(split[m]['only_expected']) for m in methods],
        ['only: tested elsewhere'] + [len(split[m]['only_tested_elsewhere']) for m in methods],
        ['only: not tested elsewhere'] + [len(split[m]['only_untested_elsewhere']) for m in methods]]
    lines += [''] + md_table(
        ['metric'] + methods,
        [['tested'] + [per[m]['tested'] for m in methods], ['surviving'] + [len(per[m]['surviving']) for m in methods]] + split_rows +
        [['only'] + [len(per[m]['only']) for m in methods]] + only_rows + [['only_reads'] + [round(per[m]['only_reads'], 1) for m in methods]] +
        ([['only_unexpected_reads'] + [round(split[m]['only_unexpected_reads'], 1) for m in methods]] if split else []) +
        ([['flagged (#415)'] + [len(per[m]['surviving'] & set(flag)) for m in methods]] if flag is not None else []))
    if split:
        def missed(m):  # tested, expected, not surviving: by control reads
            ratios = result[m]['ratios'][survival]
            names = [id_to_name.get(i, i) for i in sorted(split[m]['tested_expected'] - per[m]['surviving'], key=lambda i: (-ratios[i][0], i))]
            return ', '.join(names[:EXPECTED_NAMED]) + (' and {} more'.format(len(names) - EXPECTED_NAMED) if len(names) > EXPECTED_NAMED else '')
        if any(split[m]['tested_expected'] - per[m]['surviving'] for m in methods):
            lines += ['', 'Expected miRNAs tested but not surviving, per method (most control reads first): ' + '; '.join(
                '{}: {}'.format(m, missed(m)) for m in methods if split[m]['tested_expected'] - per[m]['surviving']) + '.']
        candidates = sorted(id_to_name.get(i, i) for i in common_set - set(expected))
        if candidates and len(methods) > 1:  # methods agreeing is evidence about the miRNA, not about one method
            lines += ['', '{} of the {} miRNAs surviving in every method are not expected: candidates to check for the list (made without'
                      ' DROSHA, or from another RNA) rather than one method\'s false counts, though the methods may share a'
                      ' misassignment: {}.'.format(len(candidates), len(common_set), ', '.join(candidates))]
    if flag is not None and any(per[m]['surviving'] & set(flag) for m in methods):
        lines += ['', 'Surviving miRNAs that lie inside molecules of other RNA classes (#415), per method: ' + '; '.join(
            '{}: {}'.format(m, ', '.join(sorted(id_to_name.get(i, i) for i in per[m]['surviving'] & set(flag))))
            for m in methods if per[m]['surviving'] & set(flag)) + '.']
    listed = set().union(*(per[m]['surviving'] for m in methods))
    order = sorted(listed, key=lambda i: (-sum(i in per[m]['surviving'] for m in methods), id_to_name.get(i, i)))
    if order:
        def cell(m, i):
            ratio = result[m]['ratios'][survival].get(i, (None, None, None))[2]
            return '' if ratio is None else ('**{:.2f}**' if ratio > biological.SURVIVING else '{:.2f}').format(ratio)
        lines += ['', 'Surviving miRNAs, {} of {} by the number of methods (bold: surviving; empty: not tested); all tested miRNAs'
                  ' of every knockout are in `knockout.tsv`:'.format(min(len(order), KNOCKOUT_LISTED), len(order)), ''] + md_table(
            ['miRNA', 'methods'] + (['expected'] if split else []) + methods,
            [[id_to_name.get(i, i), sum(i in per[m]['surviving'] for m in methods)] + ([expected.get(i, '')] if split else []) +
             [cell(m, i) for m in methods] for i in order[:KNOCKOUT_LISTED]])
    return rows, lines


def mixture_de_sets(spec, design=None):
    """The DE sets of a mixture (#451): [(DeSet, numerator fraction, denominator fraction, between mixtures)] for
    mixture.contrasts, each paired by subject (`~ subject + group`) over the subjects whose series has both fractions,
    the groups named by the fraction in percent. The pure samples' contrast is the whole dataset's set (folder ''); a
    contrast between two mixtures' fractions is in `contrast/<numerator>_vs_<denominator>`. A row tests the miRNAs
    with at least biological.MIN_READS reads on average."""
    design = design or mixture.Design(spec['samples'])
    sets = []
    for numerator, denominator, between in mixture.contrasts(design):
        up, down = dict(design.at(numerator)), dict(design.at(denominator))
        subjects = sorted(set(up) & set(down))
        high, low = mixture.label(numerator), mixture.label(denominator)
        name = '{}_vs_{}'.format(high, low)
        if len(subjects) < 2:
            fail('the mixture contrast {} has {} subject(s) with both fractions; the paired design needs at least two'.format(name, len(subjects)))
        members = [(up[s], high, {'subject': s}) for s in subjects] + [(down[s], low, {'subject': s}) for s in subjects]
        sets.append((DeSet('.' + name if between else '', 'contrast/' + name if between else '', members,
                           de.model_text(high, low, (), biological.MIN_READS),
                           '{0} {1} and {0} {2} samples of {0} subjects'.format(len(subjects), high, low)), numerator, denominator, between))
    return sets


MIXTURE_DE_KEYS = ('called', 'wrong_direction', 'fc_bias', 'called_shared', 'unchanged_called')


def mixture_de_scores(results, expected, methods):
    """({method: {<test>_<metric>: value}}, {suffix: {method: (scored, share called)}}) of a mixture's contrasts between
    its mixtures' fractions (#451). `results`: {suffix: de.read_results}; `expected`: {method: {suffix: {feature:
    expected log2 fold change}}} (mixture.expected_changes), each method's own.

    A method is scored on the miRNAs it tests in a contrast whose expected change is at least mixture.MIN_EXPECTED
    (log2), all contrasts together, with de.scores: `called`, the share of those expected at least 2-fold apart that
    are significant in the expected direction; `wrong_direction`, the share of its calls against the expected
    direction; `fc_bias`, the median of sign(d) * (log2fc - d), negative when fold changes are compressed.
    `called_shared`: `called` over the miRNAs that every method expects at least 2-fold apart in the same direction in
    that contrast (an untested one is not called), so the methods are compared on the same miRNAs. `unchanged_called`:
    the share called of the miRNAs it tests with an expected change under mixture.MIN_EXPECTED, which the other metrics
    leave out; hardly any miRNA is truly unchanged between two mixtures, so it compares methods and is not a false-call
    rate. The second value gives, per contrast and method, the miRNAs scored for `called` and their share called by edgeR's GLM test."""
    strong = [{(suffix, feature): math.copysign(1, d) for suffix, changes in expected[m].items() for feature, d in changes.items()
               if abs(d) >= 1} for m in methods]
    shared = [key for key, sign in strong[0].items() if all(other.get(key) == sign for other in strong[1:])] if strong else []
    out, per_contrast = {}, collections.defaultdict(dict)
    for m in methods:
        fits = {(m, test, suffix + '|' + feature): fit for suffix in expected[m]
                for (row, test, feature), fit in results[suffix].items() if row == m}
        scored = [(suffix + '|' + feature, d) for suffix, changes in expected[m].items() for feature, d in changes.items()
                  if abs(d) >= mixture.MIN_EXPECTED and fits[(m, de.TESTS[0], suffix + '|' + feature)][1] is not None]
        values = de.scores(fits, scored, [m])[m]
        flat = [key for key in ((suffix + '|' + feature, d) for suffix, changes in expected[m].items() for feature, d in changes.items())
                if abs(key[1]) < mixture.MIN_EXPECTED and fits[(m, de.TESTS[0], key[0])][1] is not None]
        together = de.scores(fits, [(suffix + '|' + feature, expected[m][suffix][feature]) for suffix, feature in shared], [m])[m]
        out[m] = {}
        for test in de.TESTS:
            out[m].update({test + '_called': values[test + '_called'], test + '_wrong_direction': values[test + '_fdr'],
                           test + '_fc_bias': values[test + '_fc_bias'], test + '_called_shared': together[test + '_called'],
                           test + '_unchanged_called': de.share([de.significant(fits[(m, test, key)]) for key, _ in flat])})
        for suffix in expected[m]:
            strongly = [(key, d) for key, d in scored if key.startswith(suffix + '|') and abs(d) >= 1]
            called = [de.significant(fits[(m, 'edger_glm', key)]) and fits[(m, 'edger_glm', key)][0] * d > 0 for key, d in strongly]
            per_contrast[suffix][m] = (len(strongly), de.share(called))
    return out, dict(per_contrast)


def evaluate_mixture(args, project, bench, spec, ext_dir):
    """A mixture dataset (#451): mixture.py's linearity scores of each method against its own pure samples, and, with
    `de:`, de.R on the contrasts between fractions (mixture_de_sets), scored against each method's expected changes."""
    mature_fa, _ = biological_paths(spec, ext_dir)
    mature = synthetic.read_mature(mature_fa, PREFIX[spec['organism']])
    id_to_seq = {accession: seq for _, (accession, seq) in mature.items()}
    methods = list(spec['methods'])
    first, second = (str(c) for c in spec['mixture']['components'])
    design = mixture.Design(spec['samples'])
    samples = [row['sample'] for row in read_tsv(bench / 'samples.tsv')]
    if set(samples) != set(design.samples):
        fail('benchmark/samples.tsv does not hold the declaration\'s samples; prepare the project again (a mixture has no --lanes)')
    gff = args.mirflowz_gff or (mirflowz_annotation(project) if 'mirflowz' in methods else None)
    pending = []  # tables written once the differential-expression checks (evaluate --de) have passed
    raw = common_table(args, project, bench, spec, ext_dir, methods, mature, gff, pending)
    common, _ = biological.collapse(raw, id_to_seq)
    by = biological.counts_by(common)
    features = sorted({feature for counts in by.values() for feature in counts})
    counts = {m: mixture.Counts({s: by.get((m, s), {}) for s in design.samples}, design.samples, features) for m in methods}
    scores = {m: mixture.Scores(counts[m], design) for m in methods}
    described = biological.method_summary(common, methods, samples)
    point = {m: dict(mirna_reads=described[m]['mirna_reads'], detected=described[m]['detected'], **scores[m].values()) for m in methods}
    keys = ['mirna_reads', 'detected'] + list(mixture.KEYS)
    rnd = lambda v: '' if isinstance(v, float) and math.isnan(v) else v

    de_texts = de_results = None
    de_lines = []
    if 'de' in spec:
        sets = mixture_de_sets(spec, design)
        count = lambda row, sample, feature: by.get((row, sample), {}).get(feature, 0.0)
        # A method's reads on miRNAs with an accession: the total its proportions, and so its expected changes, use. With
        # every reported read, a share of unplaced reads that differs between two fractions would shift the test and not
        # the expectation.
        sizes = {(m, sample): float(counts[m].total[counts[m].column[sample]]) for m in methods for sample in design.samples}
        de_texts = {s.folder: de.input_tables(count, methods, features, s.design, sizes, s.model) for s, *_ in sets}
        de_results, de_reports = differential_expression(args, bench, [s for s, *_ in sets], de_texts, methods, features)
    if de_results:
        between = [(s, numerator, denominator) for s, numerator, denominator, inner in sets if inner]
        expected = {m: {s.suffix: mixture.expected_changes(counts[m], design, scores[m].ks, numerator, denominator)
                        for s, numerator, denominator in between} for m in methods}
        de_point, per_contrast = mixture_de_scores(de_results, expected, methods)
        for m in methods:
            point[m].update(de_point[m])
        keys += ['{}_{}'.format(test, key) for test in de.TESTS for key in MIXTURE_DE_KEYS]
        pure = sets[0][0]
        de_lines = de.call_lines([('', pure.folder, pure.samples, de_results[''], *de_reports[''][:2])], features, methods,
                                 'pure {} ({})'.format(second, mixture.label(1.0)), 'pure {} ({})'.format(first, mixture.label(0.0)),
                                 "a method's reads on miRNAs with an accession, the total its proportions use")
        de_lines += ['', '### Between the mixtures', '',
                     'The contrasts between two fractions of {second}, the higher fraction against the lower. The GLM test and voom '
                     'pair the samples by subject (`~ subject + group`); the exact test takes no covariates and compares the two '
                     'groups as they are. A method\'s expected change of a miRNA is the mean over the subjects of '
                     'log2(prediction at one fraction / prediction at the other), from its own pure samples, which are in none of '
                     'these contrasts, and the subject\'s content ratio, one number fitted on all its mixtures. Scored over the miRNAs a method tests with an expected change of at least {least} log2, all '
                     'contrasts together: `called`, those expected at least 2-fold apart that are called in the expected '
                     'direction; `wrong_direction`, the share of its calls against it; `fc_bias`, the median of sign(d) × (log2 '
                     'fold change − d), negative when fold changes are compressed; `called_shared`, `called` over the miRNAs every '
                     'method expects at least 2-fold apart in the same direction; `unchanged_called`, the share called of the miRNAs '
                     'with a smaller expected change than that least one, which compares the methods and is not a false-call rate: '
                     'hardly any miRNA is truly unchanged between two mixtures. The metrics have no bootstrap interval.'.format(
                         second=second, least=mixture.MIN_EXPECTED), '']
        de_lines += md_table(['metric'] + methods, [['{}_{}'.format(test, key)] + [rnd(point[m]['{}_{}'.format(test, key)]) for m in methods]
                                                    for test in de.TESTS for key in MIXTURE_DE_KEYS])
        de_lines += ['', 'Per contrast, edgeR\'s GLM test: the miRNAs a method tests that it expects at least 2-fold apart, and the share '
                     'called in the expected direction.', ''] + md_table(
            ['contrast', 'samples'] + methods,
            [[s.suffix[1:].replace('_vs_', ' against '), s.samples] + ['{} / {}'.format(n, rnd(share)) for n, share in
                                                                       (per_contrast[s.suffix][m] for m in methods)] for s, _, _ in between])
        de_lines += ['', rounding_line(count, methods, features, design.samples, 'miRNAs')]
        de_lines += [''] + ['{}: design `{}` (edgeR {}, limma {}); `de/{}/results.tsv` SHA-256 `{}`.'.format(
            s.suffix[1:], de_reports[s.suffix][0].get('formula', '?'), de_reports[s.suffix][0].get('edgeR', '?'),
            de_reports[s.suffix][0].get('limma', '?'), s.folder, de_reports[s.suffix][1]) for s, _, _ in between]
    for path, header, rows in pending:
        write_tsv(path, header, rows)
    if de_texts:
        write_de_inputs(bench, de_texts)

    metric_rows = [[k] + [rnd(point[m].get(k, float('nan'))) for m in methods] for k in keys]
    percent = lambda fraction: '{}%'.format(mixture.label(fraction))
    lines = ['# {}'.format(spec['name']), '',
             '{} ({}), {} samples: total RNA of {} and {} of {} subjects, mixed per subject at {} {} (#451). No known '
             'composition: each mixture is predicted from the method\'s own counts of the subject\'s two pure samples, so the '
             'scores show how linear and consistent a method counts across samples, not whether it assigns reads to the right '
             'miRNA.'.format(spec['organism'], spec['libprepkit'], len(samples), first, second, len(design.series),
                             ', '.join(percent(f) for f in design.fractions), second), '',
             '{} technical replicate pair(s){}.'.format(len(design.replicates), '; outside the series (e.g. a pool of several '
                                                        'subjects): ' + ', '.join(design.unused) if design.unused else ''), '',
             '## Linearity', '',
             'A mixture\'s share of miRNA reads from {second} is w = f k / (f k + 1 − f) for its fraction f, with a content ratio k '
             'per method and subject (the miRNA reads per unit of RNA of {second} over {first}, as the method counts them, fitted '
             'on the subject\'s mixtures; `content_ratio`: the median over the subjects; `content_ratio_spread`: the median '
             '|log2| deviation of the subjects from it). Over the miRNA-sample pairs '
             'with at least {reads} predicted reads: `linearity_error`, the median |log2(observed / predicted)|; `within_1_5fold`, '
             'the share within 1.5-fold. `monotonic`: of the tissue-specific miRNAs (at least {specific:g}-fold apart between the pure '
             'samples), the share of subject series that never move against the trend. `replicate_error`: the median |log2 ratio| '
             'between technical replicates, the error two libraries of one RNA have, so the floor for `linearity_error`.'.format(
                 first=first, second=second, reads=mixture.MIN_READS, specific=2 ** mixture.SPECIFIC), ''] + md_table(
        ['metric'] + methods, [row for row in metric_rows if row[0] in ('mirna_reads', 'detected') + mixture.KEYS])
    lines += ['', 'Per subject, `linearity_error` / content ratio (an outlying subject points to its samples, not to a method):', ''] + md_table(
        ['subject'] + methods, [[subject] + ['{} / {}'.format(rnd(scores[m].values(subjects=[subject])['linearity_error']),
                                                              round(scores[m].ks[subject], 4) if subject in scores[m].ks else '')
                                             for m in methods] for subject in sorted(design.series)])
    abundance = {m: scores[m].by_abundance() for m in methods}
    lines += ['', '`linearity_error` by predicted reads (pairs / error):', ''] + md_table(
        ['predicted reads'] + methods, [[abundance[methods[0]][i][0]] + ['{} / {}'.format(abundance[m][i][1], rnd(abundance[m][i][2]))
                                                                         for m in methods] for i in range(len(mixture.BINS))])
    lines += de_lines
    index = list(range(len(features)))
    subjects = sorted(design.series)
    extra, record, ci_lines = bootstrap(
        args, {m: {k: rnd(v) for k, v in point[m].items()} for m in methods},
        lambda drawn, people: {m: scores[m].values(drawn, people) for m in methods},
        lambda rng: (uncertainty.resample(rng, index), uncertainty.resample(rng, subjects)), methods, list(mixture.SCORE_KEYS),
        'the miRNAs and the subjects, on the fitted content ratios (`replicate_error` resamples the miRNAs only; the '
        'differential-expression metrics have no interval)')
    publish(project, bench, spec, ['metric'] + methods, metric_rows, lines + ci_lines, extra, record)
    print(bench / 'report.md')


def depth_de(args, out, spec, table, methods, tables, usable):
    """Differential expression by depth (#414) of a biological dataset that declares `de:`: every DE set of its
    evaluation (de_sets: the contrast, each level of `within`, the null comparisons of #449) on each depth's thinned
    `tables` (depth.thinned, the draw of the detection metrics) and on the full depth, over the full depth's features.
    The sets go to <out>/de/<depth>/<set>/ (depth.de_folder), listed in <out>/de/levels.tsv for de.smk. Without --de it
    writes their inputs and returns (depth.md's lines on running de.smk, None); with --de it checks de.smk's results
    like evaluate --de and returns (depth.md's section, depth_de.tsv's rows), writing nothing itself. Depths are scored
    one at a time against the full depth's results, so only those are kept."""
    sets = de_sets(spec, de_units(spec, table))
    features = biological_de_features(sets, methods, tables[depth_series.FULL][0])
    root = out / 'de'
    nulls = [s for s in sets if s.folder.startswith(de.NULL + '/')]
    contrasts = [s for s in sets if s not in nulls]
    scores, null_points, full_results, versions = {s.folder: [] for s in contrasts}, [], None, ''
    folders = []
    if not args.de:
        (out / 'depth_de.tsv').unlink(missing_ok=True)  # scores of an earlier draw, if any, before any input changes
    for depth in [depth_series.FULL] + usable:  # the full depth first: the other depths are scored against it
        at = {s.folder: s._replace(suffix=depth_series.de_folder(depth, s.folder), folder=depth_series.de_folder(depth, s.folder))
              for s in sets}
        _, texts = biological_de_texts(at.values(), methods, tables[depth][0], {at[f].folder: features[f] for f in at})
        folders += list(texts)
        if not args.de:
            for folder, inputs in texts.items():
                de.write_inputs(root / folder / 'input', inputs)
            continue
        results, reports = differential_expression(args, None, list(at.values()), texts, methods,
                                                   {at[f].folder: features[f] for f in at}, root=root, command='depth')
        if depth == depth_series.FULL:
            full_results = {s.folder: results[at[s.folder].suffix] for s in contrasts}
            record = reports[at[(contrasts or nulls)[0].folder].suffix][0]
            versions = 'edgeR {}, limma {}; {}. Results in `de/<depth>/<set>/results.tsv`.'.format(
                record.get('edgeR', '?'), record.get('limma', '?'), record.get('R', '?'))
        for s in contrasts:
            reference = None if depth == depth_series.FULL else full_results[s.folder]
            scores[s.folder].append((depth, depth_series.de_scores(results[at[s.folder].suffix], reference, features[s.folder], methods)))
        if nulls:
            null_points.append((depth, de.null_scores([(results[at[s.folder].suffix], features[s.folder]) for s in nulls], methods)[0]))
    if not args.de:
        de.write_changed(root / de.LEVELS, de.levels_text([folder.split('/') for folder in folders]))
        return ['', 'Differential expression by depth (#414): the inputs of {} DE sets at each of {} depths are in `de/`. Run '
                'de.smk on them (`snakemake -s <repo>/smallrna/benchmark/de.smk -d {} --config de=de --use-apptainer`), then '
                '`benchmark.py depth --de` with the same depths and seed.'.format(len(sets), len(usable) + 1, out)], None
    order = lambda points: sorted(points, key=lambda point: math.inf if point[0] == depth_series.FULL else point[0])
    titled = [(depth_series.de_set_name(s.folder), 'The contrast ({})'.format(s.samples) if not s.folder else '{} ({})'.format(s.folder, s.samples),
               order(scores[s.folder])) for s in contrasts]
    null_points = order(null_points)
    return depth_series.de_lines(titled, null_points, methods, versions, md_table), depth_series.de_tables(titled, null_points, methods)


def evaluated_filtering(project, bench, samples, why):
    """{sample: biological.filtering()} of an evaluated project's samples, from the cutadapt route's reports (kept with
    --notemp), checked to be those of the evaluated run: their reads sum to the evaluation's `reads_processed`. `why`
    says what needs them, in the message when they are missing."""
    stages = project / 'data/tmp/smallrna/filter/cutadapt_stages'
    reports = {s: {k: stages / '{}.{}.json'.format(s, k) for k in biological.CUTADAPT_PASSES} for s in samples}
    missing = [str(p) for paths in reports.values() for p in paths.values() if not p.is_file()]
    if missing:
        fail('cutadapt reports missing ({}): {}; run `snakemake quant_all --notemp`'.format(', '.join(missing[:3]), why))
    try:
        filtering = {s: biological.filtering(paths) for s, paths in reports.items()}
        recorded = evaluation.read_evaluation(bench)[1].get('reads_processed')
    except ValueError as error:
        fail(str(error))
    total = sum(f['reads_in'] for f in filtering.values())
    if recorded is not None and recorded != total:  # the count tables are the evaluated run's
        fail('the cutadapt reports give {:,} reads, the evaluation recorded {:,}: the filtering changed since; evaluate again'.format(
            total, recorded))
    return filtering


def evaluated_common(bench, spec, id_to_seq):
    """(the evaluation's common-table rows with identical sequences as one accession, biological.collapse; the declared
    methods that have rows)."""
    raw = [(r['method'], r['sample'], r['mature_name'], r['mature_id'], float(r['count'])) for r in read_tsv(bench / 'common.tsv')]
    common, _ = biological.collapse(raw, id_to_seq)
    return common, [m for m in spec['methods'] if any(row[0] == m for row in common)]


def depth(args):
    """`depth` (#409): the sequencing-depth series of an evaluated biological project (depth.py), from its common and
    isomiR tables and each sample's sequenced reads (the cutadapt route's reports), written to --output (default
    benchmark/depth/ in the project): depth.tsv, depth_agreement.tsv, depth.md and, with --plot, depth.png. A dataset
    that declares `de:` also gets differential expression by depth (#414, depth_de): the inputs of de.smk in de/, and
    with --de depth_de.tsv and its section of depth.md."""
    project = Path(args.project)
    bench = project / 'benchmark'
    spec = load_dataset(bench / 'dataset.yaml')
    if spec['kind'] != 'biological':
        fail('depth is for biological datasets, got a {} dataset'.format(spec['kind']))
    if args.de and 'de' not in spec:
        fail('--de: {} declares no `de:`, so there is no differential expression by depth'.format(spec['name']))
    table = read_tsv(bench / 'samples.tsv')
    if any(row['lane'] != 'merged' for row in table):
        fail('depth works on merged samples; this project has a sample per lane (prepare --lanes)')
    if not (bench / 'common.tsv').is_file():
        fail('{} missing: run `benchmark.py evaluate --project {}` first'.format(bench / 'common.tsv', project))
    if args.plot:
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            fail('--plot needs matplotlib')
    try:
        depths = depth_series.parse_depths(args.depths)
    except ValueError as error:
        fail(str(error))
    config = yaml.safe_load((project / 'config.yaml').read_text())
    mature_fa, _ = biological_paths(spec, Path(config['ext_dir']))
    id_to_seq = {accession: seq for _, (accession, seq) in synthetic.read_mature(mature_fa, PREFIX[spec['organism']]).items()}
    samples = [row['sample'] for row in table]
    groups = collections.defaultdict(list)
    for row in table:
        groups[row['group']].append(row['sample'])
    reads_in = {s: f['reads_in'] for s, f in evaluated_filtering(project, bench, samples, 'the depth axis is the reads each sample '
                                                                 'sequenced').items()}
    out = Path(args.output) if args.output else bench / 'depth'
    if 'de' in spec and (out / 'de').resolve() == (bench / 'de').resolve():
        fail('--output {}: its de/ would be the evaluation\'s benchmark/de; choose another folder (default benchmark/depth/)'.format(out))
    if args.de:  # before the thinning: every set at every depth needs de.smk's results
        sets = de_sets(spec, de_units(spec, table))
        missing = [str(Path(folder) / name) for depth in [depth_series.FULL] + depth_series.usable_depths(depths, samples, reads_in)
                   for folder in (depth_series.de_folder(depth, s.folder) for s in sets) for name in ('results.tsv', 'record.tsv')
                   if not (out / 'de' / folder / name).is_file()]
        if missing:
            fail('no {}{} in {}: run `benchmark.py depth` without --de (the same depths and seed), then `snakemake -s {} -d {} --config '
                 'de=de --use-apptainer`, then depth --de'.format(' or '.join(missing[:3]), ' ...' if len(missing) > 3 else '', out / 'de',
                                                                 Path(__file__).resolve().with_name('de.smk'), out))
    common, methods = evaluated_common(bench, spec, id_to_seq)
    isomir_rows = None
    if (bench / 'isomirs.tsv').is_file():
        isomir_rows = [tuple(r[k] for k in adapters.ISOMIR_COLUMNS[:-1]) + (float(r['count']),) for r in read_tsv(bench / 'isomirs.tsv')]
    tables, usable = depth_series.thinned(common, samples, reads_in, depths, args.seed, isomir_rows)
    skipped = [d for d in depths if d not in usable]
    points = depth_series.series_points(tables, usable, methods, samples, groups, isomir_rows)
    full = round(statistics.median(reads_in.values()))
    equal = len(set(reads_in.values())) == 1  # the full depth is a common measurement only then
    reached, rising = depth_series.saturation(points, methods, full, equal=equal), depth_series.rise(points, methods, full)
    de_lines, de_rows = depth_de(args, out, spec, table, methods, tables, usable) if 'de' in spec else ([], None)
    out.mkdir(parents=True, exist_ok=True)
    rows, pairs = depth_series.tables(points, methods)
    write_tsv(out / 'depth.tsv', ['method', 'depth', 'metric', 'value'], rows)
    write_tsv(out / 'depth_agreement.tsv', ['method_a', 'method_b', 'depth', 'metric', 'value'], pairs)
    if de_rows is not None:
        write_tsv(out / 'depth_de.tsv', ['method', 'depth', 'set', 'metric', 'value'], de_rows)
    lines = depth_series.report(spec['name'], points, methods, skipped, reached, rising, reads_in, args.seed, md_table) + de_lines
    (out / 'depth.md').write_text('\n'.join(lines) + '\n')
    if args.plot:
        depth_series.plot(points, methods, full, out / 'depth.png')
    print(out / 'depth.md')


def imbalance(args):
    """`imbalance` (#488): differential expression under controlled depth imbalance (imbalance.py) of an evaluated
    biological project that declares `de:`, written to --output (default benchmark/imbalance/ in the project). Without
    --de it writes every DE set under every design and arm to de/<scenario>_<arm>/<set>/ for de.smk (`--config de=de`)
    and imbalance.md saying so; with --de it checks de.smk's results like evaluate --de and writes imbalance.tsv and
    imbalance.md. Each scenario and arm is scored on its own, against the same arm's unthinned fits, so only the
    results of one are held at a time."""
    project = Path(args.project)
    bench = project / 'benchmark'
    spec = load_dataset(bench / 'dataset.yaml')
    if spec['kind'] != 'biological':
        fail('imbalance is for biological datasets, got a {} dataset'.format(spec['kind']))
    if 'de' not in spec:
        fail('{} declares no `de:`: imbalance tests the DE sets of an evaluation'.format(spec['name']))
    table = read_tsv(bench / 'samples.tsv')
    if any(row['lane'] != 'merged' for row in table):
        fail('imbalance works on merged samples; this project has a sample per lane (prepare --lanes)')
    if not (bench / 'common.tsv').is_file():
        fail('{} missing: run `benchmark.py evaluate --project {}` first'.format(bench / 'common.tsv', project))
    try:
        ratios = imbalance_series.parse_ratios(args.ratios)
    except ValueError as error:
        fail(str(error))
    out = Path(args.output) if args.output else bench / 'imbalance'
    if (out / 'de').resolve() == (bench / 'de').resolve():
        fail('--output {}: its de/ would be the evaluation\'s benchmark/de; choose another folder (default benchmark/imbalance/)'.format(out))
    samples = [row['sample'] for row in table]
    kept = {s: float(f['kept']) for s, f in evaluated_filtering(project, bench, samples, 'the `fixed` arm takes each sample\'s '
                                                                 'filtered reads').items()}
    sets = de_sets(spec, de_units(spec, table))
    scenarios = imbalance_series.scenarios(ratios)
    root = out / 'de'
    if args.de:  # before the thinning: every set in every scenario and arm needs de.smk's results
        missing = [str(Path(folder) / name) for scenario, _, _ in scenarios for arm in imbalance_series.ARMS for s in sets
                   for folder in [imbalance_series.folder(scenario, arm, s.folder)] for name in ('results.tsv', 'record.tsv')
                   if not (root / folder / name).is_file()]
        if missing:
            fail('no {}{} in {}: run `benchmark.py imbalance` without --de (the same ratios and seed), then `snakemake -s {} -d {} '
                 '--config de=de --use-apptainer`, then imbalance --de'.format(' or '.join(missing[:3]), ' ...' if len(missing) > 3 else '',
                                                                            root, Path(__file__).resolve().with_name('de.smk'), out))
    config = yaml.safe_load((project / 'config.yaml').read_text())
    mature_fa, _ = biological_paths(spec, Path(config['ext_dir']))
    id_to_seq = {accession: seq for _, (accession, seq) in synthetic.read_mature(mature_fa, PREFIX[spec['organism']]).items()}
    common, methods = evaluated_common(bench, spec, id_to_seq)
    features = biological_de_features(sets, methods, common)
    by_sample = collections.defaultdict(list)
    for row in common:
        by_sample[row[1]].append(row)

    def scenario_texts(scenario, design, ratio, arms):
        """({folder: input texts}, {(set folder, arm): the set's DeSet there}) of every set in a scenario and `arms`; one
        draw per set and scenario, the same in every arm."""
        texts, at = {}, {}
        for s in sets:
            members = [(sample, group) for sample, group, _ in s.design]
            rng = imbalance_series.generator(args.seed, scenario, s.folder)
            fraction = imbalance_series.fractions(design, ratio, members, imbalance_series.denominator(s.model), rng)
            rows = [row for sample, _ in members for row in by_sample[sample]]
            if design is not None:  # the samples kept whole stay as they are, unrounded, as in `none`
                rows = [row for row in rows if fraction[row[1]] == 1] + depth_series.thin([row for row in rows if fraction[row[1]] < 1], fraction, rng)
            for arm in arms:
                folder = imbalance_series.folder(scenario, arm, s.folder)
                at[(s.folder, arm)] = s._replace(suffix=folder, folder=folder, model=s.model + 'normalisation\ttmm\n' if arm == 'tmm' else s.model)
                sizes = {(m, sample): kept[sample] * fraction[sample] for m in methods for sample, _ in members} if arm == 'fixed' else None
                texts.update(biological_de_texts([at[(s.folder, arm)]], methods, rows, {folder: features[s.folder]}, sizes)[1])
        return texts, at

    realised = []  # (design, set folder, median filtered reads times fraction of the first group, of the second)
    for scenario, design, ratio in scenarios:
        for s in sets:
            members, second = [(sample, group) for sample, group, _ in s.design], imbalance_series.denominator(s.model)
            fraction = imbalance_series.fractions(design, ratio, members, second, imbalance_series.generator(args.seed, scenario, s.folder))
            reads = [statistics.median(kept[x] * fraction[x] for x, group in members if (group == second) == is_second) for is_second in (False, True)]
            realised.append((scenario, s.folder, *reads))
    header = imbalance_series.lines(spec['name'], [], [], methods, ratios, args.seed, '', md_table)[:3] + \
        imbalance_series.realised_lines(realised, md_table)
    if not args.de:
        for name in ('imbalance.tsv', 'imbalance.md'):  # scores of an earlier draw, if any, before any input changes
            (out / name).unlink(missing_ok=True)
        folders = []
        for scenario, design, ratio in scenarios:
            texts, _ = scenario_texts(scenario, design, ratio, imbalance_series.ARMS)
            for folder, inputs in texts.items():
                de.write_inputs(root / folder / 'input', inputs)
            folders += list(texts)
        de.write_changed(root / de.LEVELS, de.levels_text([folder.split('/') for folder in folders]))
        write_tsv(out / 'imbalance_depths.tsv', ['design', 'set', 'first', 'second', 'ratio'], imbalance_series.realised_rows(realised))
        header += ['', 'The inputs of {} DE sets in {} designs and {} arms are in `de/`. Run de.smk on them (`snakemake -s '
                   '<repo>/smallrna/benchmark/de.smk -d {} --config de=de --use-apptainer`), then `benchmark.py imbalance --de` with the '
                   'same ratios and seed.'.format(len(sets), len(scenarios), len(imbalance_series.ARMS), out)]
        (out / 'imbalance.md').write_text('\n'.join(header) + '\n')
        print(out / 'imbalance.md')
        return
    nulls = [s for s in sets if s.folder.startswith(de.NULL + '/')]
    contrasts = [s for s in sets if s not in nulls]
    points, null_points, versions = {s.folder: [] for s in contrasts}, [], ''
    for arm in imbalance_series.ARMS:
        reference = None
        for scenario, design, ratio in scenarios:  # `none` first: the others are scored against it
            texts, at = scenario_texts(scenario, design, ratio, [arm])
            found, reports = differential_expression(args, None, list(at.values()), texts, methods,
                                                     {at[(s.folder, arm)].folder: features[s.folder] for s in sets}, root=root, command='imbalance')
            result = {s.folder: found[at[(s.folder, arm)].suffix] for s in sets}
            if not versions:
                record = reports[at[(sets[0].folder, arm)].suffix][0]
                versions = 'edgeR {}, limma {}; {}. Results in `de/<design>_<arm>/<set>/results.tsv`.'.format(
                    record.get('edgeR', '?'), record.get('limma', '?'), record.get('R', '?'))
            if scenario == imbalance_series.NONE:
                reference = {s.folder: result[s.folder] for s in contrasts}
            for s in contrasts:
                values = depth_series.de_scores(result[s.folder], None if scenario == imbalance_series.NONE else reference[s.folder],
                                                features[s.folder], methods)
                if scenario != imbalance_series.NONE:
                    for m, shift in imbalance_series.shifts(result[s.folder], reference[s.folder], features[s.folder], methods).items():
                        values[m].update(shift)
                points[s.folder].append((scenario, arm, values))
            if nulls:
                splits = [(result[s.folder], features[s.folder]) for s in nulls]
                values = de.null_scores(splits, methods)[0]
                for m, offset in imbalance_series.null_offsets(splits, methods).items():
                    values[m].update(offset)
                null_points.append((scenario, arm, values))
    order = {(scenario, arm): i for i, (scenario, arm) in enumerate((sc[0], a) for sc in scenarios for a in imbalance_series.ARMS)}
    ordered = lambda items: sorted(items, key=lambda item: order[(item[0], item[1])])
    titled = [(depth_series.de_set_name(s.folder), 'The contrast ({})'.format(s.samples) if not s.folder else '{} ({})'.format(s.folder, s.samples),
               ordered(points[s.folder])) for s in contrasts]
    null_points = ordered(null_points)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'imbalance.tsv', ['method', 'design', 'arm', 'set', 'metric', 'value'], imbalance_series.rows(titled, null_points, methods))
    write_tsv(out / 'imbalance_depths.tsv', ['design', 'set', 'first', 'second', 'ratio'], imbalance_series.realised_rows(realised))
    body = imbalance_series.lines(spec['name'], titled, null_points, methods, ratios, args.seed, versions, md_table)
    (out / 'imbalance.md').write_text('\n'.join(header + body[3:]) + '\n')
    print(out / 'imbalance.md')


def fill(args):
    spec = load_dataset(args.dataset, args.mirbase_release)
    try:
        filler.fill(spec, Path(args.ext_dir), source=args.source, connections=args.connections, jobs=args.jobs,
                    fasterq_dump=args.fasterq_dump)
    except (OSError, ValueError) as error:
        fail(str(error))


def summary(args):
    try:
        output = summarizer.write(args.projects, args.output)
    except (OSError, ValueError) as error:
        fail(str(error))
    print(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--dataset', required=True)
    p.add_argument('--project', required=True)
    p.add_argument('--ext-dir', required=True, help='EXT_DIR holding the references and data named in the declaration')
    p.add_argument('--mirbase-release', help="a miRBase release the declaration lists under mirbase_releases (default: its mirbase)")
    p.add_argument('--methods', help="comma-separated subset of the declaration's methods, e.g. unitas,mirge3")
    p.add_argument('--reads', type=int, help='pool and biological datasets: copy the first N reads of each run instead of linking the runs; '
                   'genome-wide synthetic dataset: the mean reads of a mixed sample')
    p.add_argument('--lanes', action='store_true', help='biological datasets: each run (lane) its own sample, for the lane-replicate metric')
    p.add_argument('--per-group', type=int, help='biological datasets: the first N samples of each group, in declared order (a quick run)')
    p.add_argument('--synthetic-seed', type=int, help="synthetic datasets: the seed of the draw instead of the declaration's (#455); "
                   'each seed is a project of its own, and the summary groups the seeds of one declaration')
    p.add_argument('--set', action='append', default=[], metavar='quant.ROUTE.SETTING=VALUE',
                   help='a settings variant (#417): a route setting of smallrna.config changed for this project, e.g. '
                        'quant.tinyrna.mismatch=1 (repeatable); the summary compares it with a project of the same dataset run with the defaults')
    f = sub.add_parser('fill', help='download a public dataset into EXT_DIR (the only step that writes there)')
    f.add_argument('--dataset', required=True)
    f.add_argument('--ext-dir', required=True)
    f.add_argument('--mirbase-release', help="a miRBase release the declaration lists under mirbase_releases (default: its mirbase)")
    f.add_argument('--source', choices=['ena', 'sra'], default='ena',
                   help="where the runs come from: ena (the declared url) or sra (NCBI's copy on AWS, for runs declaring `sra`; needs fasterq-dump)")
    f.add_argument('--connections', type=int, default=8, help='--source sra: range requests at a time per run')
    f.add_argument('--jobs', type=int, default=4, help='--source sra: runs at a time')
    f.add_argument('--fasterq-dump', default='fasterq-dump', help="--source sra: the fasterq-dump command, e.g. 'apptainer exec <image> fasterq-dump'")
    e = sub.add_parser('evaluate')
    e.add_argument('--project', required=True)
    e.add_argument('--de', action='store_true', help="datasets that declare `de:`: score de.smk's differential "
                   'expression (#414, #447); its results must come from the inputs this evaluation writes')
    e.add_argument('--bootstrap', type=nonnegative, default=uncertainty.RESAMPLES, metavar='N',
                   help='bootstrap resamples for the metrics\' 95%% intervals (#413); 0 for none (default %(default)s)')
    e.add_argument('--seed', type=int, default=uncertainty.SEED, help='seed of the bootstrap (default %(default)s)')
    e.add_argument('--mirflowz-gff', type=Path, help="mirflowz's extended miRNA annotation; default the one the project's mirflowz run recorded")
    e.add_argument('--unitas-ref', type=Path, help="UNITAS reference dump (its miRBase hairpins place locus names); default the "
                   "dump the project's UNITAS run recorded")
    s = sub.add_parser('summary', help='one overview of several evaluated projects (Markdown or HTML by the extension)')
    s.add_argument('--projects', nargs='+', required=True)
    s.add_argument('--output', required=True, help='.md or .html; written outside the projects')
    d = sub.add_parser('depth', help='sequencing-depth curves of an evaluated biological project, by thinning its count tables (#409)')
    d.add_argument('--project', required=True)
    d.add_argument('--output', help='directory for depth.tsv, depth_agreement.tsv and depth.md (default benchmark/depth/ in the project)')
    d.add_argument('--plot', action='store_true', help='also depth.png, detection and replicate CV against depth (needs matplotlib)')
    d.add_argument('--depths', default=depth_series.DEPTHS, help='reads per sample, k and M suffixes (default %(default)s); '
                   'those above the shallowest sample are left out, and the full depth ends the series')
    d.add_argument('--seed', type=int, default=1, help='seed of the thinning (default %(default)s)')
    d.add_argument('--de', action='store_true', help="datasets that declare `de:`: score de.smk's differential expression by "
                   "depth (#414), after `depth` without it wrote the inputs to de/ and de.smk ran on them")
    m = sub.add_parser('imbalance', help='differential expression under controlled depth imbalance of an evaluated biological '
                       'project, by thinning its count tables unevenly (#488)')
    m.add_argument('--project', required=True)
    m.add_argument('--output', help='directory for de/, imbalance.tsv and imbalance.md (default benchmark/imbalance/ in the project)')
    m.add_argument('--ratios', default=imbalance_series.RATIOS, help='depth ratios between the deeper and the thinned samples, each '
                   'above 1 (default %(default)s)')
    m.add_argument('--seed', type=int, default=1, help='seed of the draws (default %(default)s)')
    m.add_argument('--de', action='store_true', help="score de.smk's results, after `imbalance` without it wrote the inputs to de/ "
                   "and de.smk ran on them")
    args = parser.parse_args()
    {'prepare': prepare, 'fill': fill, 'evaluate': evaluate, 'summary': summary, 'depth': depth, 'imbalance': imbalance}[args.command](args)


if __name__ == '__main__':
    main()
