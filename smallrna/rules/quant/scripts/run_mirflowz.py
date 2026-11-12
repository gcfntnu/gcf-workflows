#!/usr/bin/env python
"""Run mirflowz 0.10.0's mapping and quantification on filtered small-RNA reads, against a prepared reference.

The script writes a small Snakefile that includes `$MIRFLOWZ_DIR/workflow/rules/map.smk` and `quantify.smk` (not
mirflowz's top-level Snakefile, which deletes its intermediates directory, and not its preparation rules), a
configuration and a sample table, and runs it in a work directory next to `--output`. The work directory's
intermediates directory links the files of the prepared reference, so they are read, never rebuilt or removed.

mirflowz's read preprocessing is set so it does not filter the already filtered reads again: reads are passed as
FASTA, which skips mirflowz's quality filter (it only filters FASTQ, then converts to FASTA itself), `overlap` equal to
the adapter length (only a complete adapter is trimmed, so 3' ends matching
the start of the adapter are kept), `minimum_length` the kit's final minimum, and no limit on N bases. mirflowz's
cutadapt call always trims N bases at read ends (`--trim-n`).

mirflowz reads its sample table with pandas type inference, so names such as `1` or `01` would become numbers (and
collide). Samples are therefore passed to mirflowz as `gcf_sample_<n>` and given their own names back in the count
table headers and alignment directories; `mirflowz_samples.tsv` keeps both.

mirflowz declares its threads and memory (GB) as custom resources (`resources: threads=8, mem=50`), which `--cores`
does not limit. Each rule's `threads` and `mem` are capped at `--threads` and `--mem-mb` (the outer job's allocation, in whole GB), Snakemake's core count is
set to the same threads, and `--resources` limits the totals, so concurrent jobs stay within the outer rule's budget.

Outputs in `--output`: `all_mirna_counts.tab` (miRNAs and isomiRs), `all_pri-mir_counts.tab`, `alignments/<sample>/`
(mirflowz's per-sample SAM and BAM files) and the generated `mirflowz_config.yaml` and `mirflowz_samples.tsv`.

Usage: run_mirflowz.py --prepared DIR --adapter SEQ --minimum-length N --threads N --mem-mb N --output DIR SAMPLE=FASTQ ...
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

EXTENSION = 6
PREPARED = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
            'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv',
            'extended_mirna_annotation_{}_nt.gff3'.format(EXTENSION), 'extended_primir_annotation_{}_nt.gff3'.format(EXTENSION))
TABLES = ('all_mirna_counts.tab', 'all_pri-mir_counts.tab')
# Resources (threads, mem in GB) that mirflowz 0.10.0's mapping and quantification rules declare.
RESOURCES = {'remove_adapters': (8, None), 'map_genome_segemehl': (8, 50), 'map_transcriptome_segemehl': (8, 10),
             'map_genome_oligomap': (8, 50), 'map_transcriptome_oligomap': (8, 10), 'sort_genome_oligomap': (8, None),
             'sort_transcriptome_oligomap': (8, None), 'remove_inferiors': (4, 15), 'filter_by_indels': (4, 15)}
TARGET = '''
rule gcf_mirflowz_all:
    input:
        expand(OUT_DIR / "TABLES" / "all_{mir}_counts.tab", mir=[m for m in config["mir_list"] if m != "isomir"]),
        expand(OUT_DIR / "{sample}" / "alignments_intersecting_mirna_uncollapsed_sorted.bam.bai",
               sample=pd.unique(samples_table.index.values)),
'''


def fail(message):
    sys.exit('run_mirflowz.py: ' + message)


def sample_arg(text):
    name, sep, path = text.partition('=')
    if not sep or not name or not path:
        raise argparse.ArgumentTypeError('expected SAMPLE=FASTQ, got {!r}'.format(text))
    return name, Path(path)


def fastq_to_fasta(source, target):
    """Write the FASTQ reads (plain or gzipped) as gzipped FASTA, names and sequences unchanged; return the count."""
    with open(source, 'rb') as fh:
        gzipped = fh.read(2) == b'\x1f\x8b'
    count = 0
    with (gzip.open(source, 'rt') if gzipped else open(source)) as src, gzip.open(target, 'wt', compresslevel=1) as out:
        while True:
            header, sequence, plus, quality = (src.readline() for _ in range(4))
            if not header:
                return count
            if not header.startswith('@') or not plus.startswith('+'):
                fail('{} is not a four-line FASTQ file (read {})'.format(source, count + 1))
            out.write('>{}\n{}'.format(header[1:].split()[0], sequence))
            count += 1


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


def mirflowz_config(work, prepared, mirflowz_dir, adapter, minimum_length):
    """mirflowz's configuration for mapping and quantification; every key of its template must be present."""
    return {'samples': str(work / 'samples.tsv'),
            # Read only by the preparation rules, which are not included; the prepared reference is linked instead.
            'genome_file': str(prepared / 'genome_processed.fa'), 'gtf_file': 'NA', 'mirna_file': 'NA', 'map_chr_file': 'NA',
            'bed_file': str(work / 'empty.bed'), 'output_dir': str(work / 'results'),
            'pileups_dir': str(work / 'results/pileups'), 'intermediates_dir': str(work / 'intermediates'),
            'local_log': str(work / 'logs/local'), 'cluster_log': str(work / 'logs/cluster'),
            'scripts_dir': str(mirflowz_dir / 'scripts') + '/', 'extension': EXTENSION,
            'q_value': 10, 'p_value': 50, 'error_rate': 0.1, 'minimum_length': minimum_length, 'overlap': len(adapter),
            'max_n': 1000, 'max_length_reads': 30, 'nh': 100, 'mir_list': ['isomir', 'mirna', 'pri-mir'], 'lib_dict': {}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--prepared', required=True, type=Path, help='prepared mirflowz reference (gcfdb/mirflowz.db)')
    parser.add_argument('--adapter', required=True, help="the kit's 3' adapter")
    parser.add_argument('--minimum-length', required=True, type=int)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--mem-mb', type=int, required=True, help="the outer job's memory allocation (MB)")
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('samples', nargs='+', type=sample_arg, metavar='SAMPLE=FASTQ')
    args = parser.parse_args()

    mirflowz_dir = Path(os.environ.get('MIRFLOWZ_DIR', ''))
    if not (mirflowz_dir / 'workflow/rules/map.smk').is_file():
        fail('MIRFLOWZ_DIR does not hold mirflowz ({}); run in the mirflowz image'.format(mirflowz_dir))
    missing = [name for name in PREPARED if not (args.prepared / name).is_file()]
    if missing:
        fail('prepared mirflowz reference {} lacks {}'.format(args.prepared, ', '.join(missing)))
    names = [name for name, _ in args.samples]
    if len(set(names)) != len(names):
        fail('sample names are not unique')
    if not args.adapter.isalpha():
        fail('adapter must be a nucleotide sequence, got {!r}'.format(args.adapter))

    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.mirflowz-run-', dir=args.output) as work:
        work = Path(work).resolve()
        (work / 'intermediates').mkdir()
        for name in PREPARED:
            os.symlink((args.prepared / name).resolve(), work / 'intermediates' / name)
        (work / 'reads').mkdir()
        rows, internal = [], {}
        for i, (name, path) in enumerate(args.samples, 1):
            internal[name] = 'gcf_sample_{}'.format(i)
            target = work / 'reads' / (internal[name] + '.fa.gz')
            reads = fastq_to_fasta(path, target)
            if not reads:
                fail('{} ({}) holds no reads'.format(path, name))
            rows.append('{}\t{}\t{}\tfa\n'.format(internal[name], target, args.adapter))
        (work / 'samples.tsv').write_text('sample\tsample_file\tadapter\tformat\n' + ''.join(rows))
        (work / 'empty.bed').touch()
        config = mirflowz_config(work, args.prepared.resolve(), mirflowz_dir, args.adapter, args.minimum_length)
        (work / 'config.yaml').write_text(''.join('{}: {}\n'.format(k, json.dumps(v)) for k, v in config.items()))
        rules = mirflowz_dir / 'workflow/rules'
        (work / 'Snakefile').write_text('include: {!r}\ninclude: {!r}\n'.format(str(rules / 'map.smk'), str(rules / 'quantify.smk')) + TARGET)
        env = dict(os.environ, XDG_CACHE_HOME=str(work / '.cache'))
        result = subprocess.run(['snakemake', '-s', str(work / 'Snakefile'), '--configfile', str(work / 'config.yaml'),
                                 '--directory', str(work), '--cores', str(args.threads), '--nolock', '--workflow-profile', 'none', 'gcf_mirflowz_all',
                                 # The work directory is removed afterwards, so no metadata is kept; Snakemake then also
                                 # skips reading the conda files that mirflowz's rules name (conda is not used).
                                 '--drop-metadata',
                                 # Last: these options take several values each, so a target after them would be read as one.
                                 *resource_args(RESOURCES, args.threads, mem_budget_gb(args.mem_mb))],
                                capture_output=True, text=True, env=env)
        if result.returncode != 0:
            logs = sorted((work / 'logs/local').glob('*.log'), key=os.path.getmtime)
            tail = ''.join('\n--- {}\n{}'.format(p.name, p.read_text()[-1500:]) for p in logs[-3:])
            fail('mirflowz failed:\n{}{}'.format(result.stderr[-3000:], tail))

        results = work / 'results'
        missing = [name for name in TABLES if not (results / 'TABLES' / name).is_file()]
        if missing:
            fail('mirflowz wrote no {}'.format(', '.join(missing)))
        alignments = args.output / 'alignments'
        if alignments.exists():
            shutil.rmtree(alignments)
        alignments.mkdir()
        for name in names:
            os.replace(results / internal[name], alignments / name)
        original = {v: k for k, v in internal.items()}
        for name in TABLES:
            lines = (results / 'TABLES' / name).read_text().split('\n', 1)
            header = lines[0].split('\t')
            if sorted(header[1:]) != sorted(original):
                fail('{} has columns {}, expected the samples {}'.format(name, ', '.join(header[1:]), ', '.join(sorted(original))))
            lines[0] = '\t'.join(header[:1] + [original[c] for c in header[1:]])
            (args.output / name).write_text('\n'.join(lines))
        shutil.copy(work / 'config.yaml', args.output / 'mirflowz_config.yaml')
        (args.output / 'mirflowz_samples.tsv').write_text('sample\tmirflowz_sample\tsample_file\tadapter\tformat\n' + ''.join(
            '{}\t{}'.format(original[row.split('\t', 1)[0]], row) for row in rows))


if __name__ == '__main__':
    main()
