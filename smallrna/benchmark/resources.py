"""Runtime and memory of each quantification method (#386), from the Snakemake benchmark files its rules write under
`data/tmp/smallrna/quant/benchmarks/` (QUANT_BENCHMARKS in rules/quant.smk): `<rule>/<sample>.tsv` for a rule run per
sample, `<rule>.tsv` for one run per project.

Per job Snakemake records the wall time (`s`), the CPU time (`cpu_time`) and the peak memory (`max_rss`, MB) of the
job's process tree, containers included. The wall time is exact, I/O waits included (e.g. a reference read over a
network mount). Memory and CPU are sampled, every 0.5 s for the first 15 s and every 30 s after: a short memory peak
can be missed, and the CPU time is a lower bound, since a child process's CPU after its last sample is lost. Routes
that start many short processes are undercounted most: mirflowz's nested Snakemake run showed 12 s of CPU where
/usr/bin/time measured 45 s (#386).

Each method's steps are listed in STEPS: `quant` steps make its count table from the filtered reads; `reference` steps
build, in the project, what it counts against (e.g. smrnaseq's hairpin index). A method is measured when every `quant`
step has its file; a run from before #386, or with a step left out, is reported as not measured.

A benchmark file is an output of its job, rewritten whenever the job runs, so each file describes the run that made
that job's current output. After a partial rerun a method's jobs can therefore come from different runs, possibly on
different machines; a benchmark project quantified from scratch has them from one run.

The shared reference builds under EXT_DIR (gcfdb/mirflowz.db's mirflowz_prepare and mirflowz_annotation, gcfdb/tesmall.db's tesmall_genome
and tesmall_db, gcfdb/biocore_htseq.db's biocore_htseq_bowtie2_index, gcfdb/biocore_miraligner.db's biocore_miraligner_db) are not reported here: their
benchmark files are in the reference they built (`<reference>/benchmarks/<rule>.tsv`), whichever project built it.
"""
import csv
import math
from pathlib import Path

# method: [(rule, per sample, kind)]
STEPS = {
    'unitas': [('unitas', True, 'quant'), ('unitas_mirtable', False, 'quant')],
    'mirge3': [('mirge3_project_library', False, 'reference'), ('mirge3_quant', False, 'quant'), ('mirge3_counts', False, 'quant')],
    'mirflowz': [('mirflowz_quant', False, 'quant')],
    'tesmall': [('tesmall_quant', True, 'quant'), ('tesmall_counts', False, 'quant')],
    'tinyrna': [('tinyrna_features', False, 'reference'), ('tinyrna_collapse', True, 'quant'),
                ('tinyrna_align', True, 'quant'), ('tinyrna_count', True, 'quant'), ('tinyrna_counts', False, 'quant')],
    'smrnaseq': [('smrnaseq_hairpin', False, 'reference'), ('smrnaseq_hairpin_index', False, 'reference'),
                 ('smrnaseq_collapse', True, 'quant'), ('smrnaseq_align', True, 'quant'), ('smrnaseq_mirtop', True, 'quant'),
                 ('smrnaseq_counts', False, 'quant')],
    'biocore_htseq': [('biocore_htseq_gff', False, 'reference'), ('biocore_htseq_collapse', True, 'quant'),
                      ('biocore_htseq_align', True, 'quant'), ('biocore_htseq_annotate', True, 'quant'),
                      ('biocore_htseq_count', True, 'quant'), ('biocore_htseq_counts', False, 'quant')],
    'biocore_miraligner': [('biocore_miraligner_collapse', True, 'quant'), ('biocore_miraligner_align', True, 'quant'),
                           ('biocore_miraligner_count', True, 'quant'), ('biocore_miraligner_counts', False, 'quant')],
}
HEADER = ['method', 'kind', 'step', 'sample', 'wall_s', 'cpu_s', 'max_rss_mb']
FIELDS = {'s': 'wall_s', 'cpu_time': 'cpu_s', 'max_rss': 'max_rss_mb'}


def read_benchmark(path):
    """{wall_s, cpu_s, max_rss_mb} of one Snakemake benchmark file: one row (no repeats); `NA` (Snakemake's value when
    it collected no process data, e.g. a job of a few milliseconds) becomes None. Anything else raises ValueError."""
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    if len(rows) != 2 or len(rows[0]) != len(rows[1]) or not set(FIELDS) <= set(rows[0]):
        raise ValueError('{}: not a Snakemake benchmark file with one row and the columns {}'.format(path, ', '.join(FIELDS)))
    record = dict(zip(rows[0], rows[1]))
    out = {}
    for field, name in FIELDS.items():
        if record[field] == 'NA':
            out[name] = None
            continue
        try:
            value = float(record[field])
        except ValueError:
            raise ValueError('{}: {} is {!r}, not a number'.format(path, field, record[field])) from None
        if not math.isfinite(value) or value < 0:
            raise ValueError('{}: {} is {!r}'.format(path, field, record[field]))
        out[name] = value
    if out['wall_s'] is None:
        raise ValueError('{}: no wall time'.format(path))
    return out


def jobs(quant, method, samples):
    """([(kind, step, sample or '', measurement)], [missing quant files]) of a method's benchmark files."""
    base = Path(quant) / 'benchmarks'
    found, missing = [], []
    for step, per_sample, kind in STEPS[method]:
        for sample in (samples if per_sample else ['']):
            path = base / step / '{}.tsv'.format(sample) if per_sample else base / '{}.tsv'.format(step)
            if path.is_file():
                found.append((kind, step, sample, read_benchmark(path)))
            elif kind == 'quant':
                missing.append(str(path.relative_to(base)))
    return found, missing


def total(values):
    values = [v for v in values if v is not None]
    return sum(values) if values else None


def summarize(quant, methods, samples, reads):
    """(resources.tsv text, report lines): every measured job, and per method its quantification's wall and CPU time
    summed over the jobs, the largest job's peak memory and, with `reads` (the run's input reads), the reads per second
    of summed wall time; then the reference steps measured in this run."""
    rows, summary, references, unmeasured = [], [], [], []
    for method in methods:
        found, missing = jobs(quant, method, samples)
        rows += [[method, kind, step, sample] + ['' if m[k] is None else '{:g}'.format(m[k]) for k in HEADER[4:]] for kind, step, sample, m in found]
        if missing:
            unmeasured.append('{} ({} missing, e.g. `{}`)'.format(method, len(missing), missing[0]))
            summary.append([method] + ['not measured'] + [''] * 4)
        else:
            quant_jobs = [m for kind, _, _, m in found if kind == 'quant']
            wall, cpu = total(m['wall_s'] for m in quant_jobs), total(m['cpu_s'] for m in quant_jobs)
            peak = max(((m['max_rss_mb'], step) for kind, step, _, m in found if kind == 'quant' and m['max_rss_mb'] is not None), default=None)
            summary.append([method, len(quant_jobs), fmt(wall), fmt(cpu), '{} ({})'.format(fmt(peak[0]), peak[1]) if peak else '', fmt(reads / wall) if reads and wall else ''])
        by_step = {}
        for kind, step, _, m in found:
            if kind == 'reference':
                by_step.setdefault(step, []).append(m)
        for step, measured in by_step.items():
            references.append([method, step, fmt(total(m['wall_s'] for m in measured)), fmt(total(m['cpu_s'] for m in measured)),
                               fmt(max((m['max_rss_mb'] for m in measured if m['max_rss_mb'] is not None), default=None))])
    text = '\t'.join(HEADER) + '\n' + ''.join('\t'.join(map(str, row)) + '\n' for row in rows)
    lines = ['', '## Runtime and memory', '',
             'Quantification steps from the filtered reads to the count table (resources.py), from the Snakemake benchmark files '
             'of their jobs: wall time summed over the jobs (jobs may have run in parallel; I/O waits included), CPU time summed '
             '(sampled, a lower bound: the CPU of short-lived child processes is partly lost), peak memory the largest job\'s '
             '(its process tree, containers included; sampled every 30 s after the first 15 s, so short peaks can be missed){}. '
             'The numbers depend on the machine, its storage and the cores and memory the run was given; every job is in '
             '`resources.tsv`.'.format(
                 ', and the run\'s input reads ({}) per second of summed wall time'.format(reads) if reads else ''), '']
    lines += md_table(['method', 'jobs', 'wall s', 'CPU s (at least)', 'peak MB (step)', 'reads/s'], summary)
    if unmeasured:
        lines += ['', 'Not measured (run before benchmark files were written, or steps left out): ' + '; '.join(unmeasured) + '.']
    lines += ['', 'Reference steps built in the project (the miRge3 project library, smrnaseq\'s hairpin index, tiny-count\'s '
              'feature sheet). Shared references under EXT_DIR record their build in their own `benchmarks/` folder.', '']
    lines += md_table(['method', 'step', 'wall s', 'CPU s (at least)', 'peak MB'], references) if references else ['None.']
    return text, lines


def fmt(value):
    if value is None:
        return ''
    return '{:.0f}'.format(value) if value >= 100 else '{:.3g}'.format(value)


def md_table(header, rows):
    return ['| ' + ' | '.join(header) + ' |', '|---' * len(header) + '|'] + ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows]
