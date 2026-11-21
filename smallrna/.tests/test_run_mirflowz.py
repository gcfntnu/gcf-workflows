"""run_mirflowz.py: mirflowz mapping and quantification on filtered reads with a linked reference; snakemake is mocked."""
import gzip
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/run_mirflowz.py'
ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'
PREPARED = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
            'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv',
            'extended_mirna_annotation_6_nt.gff3', 'extended_primir_annotation_6_nt.gff3')

# Records the call, the sample table and the reads as mirflowz sees them, and writes mirflowz's outputs. Like mirflowz,
# it reads the sample table with pandas type inference and looks samples up by their string wildcard.
FAKE_SNAKEMAKE = '''#!{python}
import gzip, json, os, sys
from pathlib import Path
import pandas as pd
args = sys.argv[1:]
config = {k: json.loads(v) for k, v in (l.split(': ', 1) for l in open(args[args.index('--configfile') + 1]).read().splitlines())}
work = Path(args[args.index('--directory') + 1])
inter = Path(config['intermediates_dir'])
rows = [l.split('\\t') for l in open(config['samples']).read().splitlines()]
table = pd.read_csv(config['samples'], header=0, index_col=0, comment='#', engine='python', sep='\\t')
assert table.index.is_unique, 'duplicate sample index'
for r in rows[1:]:
    table.loc[r[0]]  # KeyError when pandas turned the ID into a number
call = {'args': args, 'config': config, 'snakefile': Path(args[args.index('-s') + 1]).read_text(), 'rows': rows,
        'links': {p.name: os.readlink(p) for p in inter.iterdir() if p.is_symlink()},
        'reads': {r[0]: gzip.open(r[1], 'rt').read() for r in rows[1:]},
        'env': {k: os.environ.get(k) for k in ('OLIGOMAP_JOBS', 'OLIGOMAP_SPLIT_DIR')}}
Path(os.environ['CALL_COPY']).write_text(json.dumps(call))
if os.environ.get('FAIL'):
    Path(config['local_log']).mkdir(parents=True); (Path(config['local_log']) / 'remove_adapters_S1.log').write_text('cutadapt broke')
    sys.exit(1)
out = Path(config['output_dir'])
(out / 'TABLES').mkdir(parents=True)
for name in ('mirna', 'pri-mir'):
    (out / 'TABLES' / 'all_{}_counts.tab'.format(name)).write_text('ID\\t' + '\\t'.join(r[0] for r in rows[1:]) + '\\n')
for r in rows[1:]:
    (out / r[0]).mkdir()
    (out / r[0] / 'alignments_intersecting_mirna_uncollapsed_sorted.bam').write_text('bam ' + r[0])
'''.replace('{python}', sys.executable)


@pytest.fixture
def setup(tmp_path):
    mirflowz = tmp_path / 'opt/mirflowz'
    (mirflowz / 'workflow/rules').mkdir(parents=True)
    for name in ('map.smk', 'quantify.smk', 'prepare.smk'):
        (mirflowz / 'workflow/rules' / name).write_text('# ' + name)
    (tmp_path / 'bin').mkdir()
    (tmp_path / 'bin/snakemake').write_text(FAKE_SNAKEMAKE)
    (tmp_path / 'bin/snakemake').chmod(0o755)
    prepared = tmp_path / 'ext/mirflowz/ref/mus_musculus/prepared'
    prepared.mkdir(parents=True)
    for name in PREPARED:
        (prepared / name).write_text(name)
    (tmp_path / 'S1_R1.fastq').write_text('@r1\nTGAGGTAGTAGGTTGTATAGTT\n+\nIIIIIIIIIIIIIIIIIIIIII\n')
    (tmp_path / 'S2_R1.fastq.gz').write_bytes(gzip.compress(b'@r2\nTAGCTTATCAGACTGATGTTGA\n+\nIIIIIIIIIIIIIIIIIIIIII\n'))
    env = dict(os.environ, PATH='{}:{}'.format(tmp_path / 'bin', os.environ['PATH']), MIRFLOWZ_DIR=str(mirflowz),
               CALL_COPY=str(tmp_path / 'call.json'))
    return tmp_path, mirflowz, prepared, env


def run(setup, fail=False, prepared=None, samples=None, threads='3', mem_mb='12800'):
    tmp_path, _, default_prepared, env = setup
    env = dict(env, FAIL='1' if fail else '')
    out = tmp_path / 'quant/mirflowz'
    samples = samples or ['S1={}'.format(tmp_path / 'S1_R1.fastq'), 'S2={}'.format(tmp_path / 'S2_R1.fastq.gz')]
    result = subprocess.run([sys.executable, str(SCRIPT), '--prepared', str(prepared or default_prepared), '--adapter', ADAPTER,
                             '--minimum-length', '16', '--threads', threads, '--mem-mb', mem_mb, '--output', str(out), *samples],
                            capture_output=True, text=True, env=env, timeout=60)
    call = json.loads((tmp_path / 'call.json').read_text()) if (tmp_path / 'call.json').exists() else None
    return result, out, call


def test_mapping_and_quantification_use_the_linked_reference_and_keep_the_outputs(setup):
    _, mirflowz, prepared, _ = setup
    result, out, call = run(setup)
    assert result.returncode == 0, result.stderr
    rules = mirflowz / 'workflow/rules'
    assert call['snakefile'].startswith("include: '{}'\ninclude: '{}'\n".format(rules / 'map.smk', rules / 'quantify.smk'))
    assert 'prepare.smk' not in call['snakefile'] and 'rule gcf_mirflowz_all:' in call['snakefile']
    assert call['args'][call['args'].index('--nolock') + 3] == 'gcf_mirflowz_all' and call['args'][call['args'].index('--cores') + 1] == '3'
    assert call['links'] == {name: str((prepared / name).resolve()) for name in PREPARED}
    assert [row[0] for row in call['rows']] == ['sample', 'gcf_sample_1', 'gcf_sample_2']
    assert all(row[2:] == [ADAPTER, 'fa'] for row in call['rows'][1:])
    # FASTA input skips mirflowz's quality filter; names and sequences are kept.
    assert call['reads'] == {'gcf_sample_1': '>r1\nTGAGGTAGTAGGTTGTATAGTT\n', 'gcf_sample_2': '>r2\nTAGCTTATCAGACTGATGTTGA\n'}
    assert sorted(p.name for p in out.iterdir()) == ['alignments', 'all_mirna_counts.tab', 'all_pri-mir_counts.tab',
                                                    'mirflowz_config.yaml', 'mirflowz_samples.tsv']
    assert (out / 'all_mirna_counts.tab').read_text() == 'ID\tS1\tS2\n'
    assert (out / 'alignments/S2/alignments_intersecting_mirna_uncollapsed_sorted.bam').read_text() == 'bam gcf_sample_2'
    assert all((prepared / name).read_text() == name for name in PREPARED)


def test_mirflowz_does_not_filter_the_filtered_reads_again(setup):
    """FASTA input (no quality filter), adapter trimmed only where it occurs in full, the kit's minimum length, no N limit."""
    _, _, call = run(setup)
    config = call['config']
    assert all(row[3] == 'fa' for row in call['rows'][1:])
    assert config['overlap'] == len(ADAPTER)
    assert config['minimum_length'] == 16
    assert config['max_n'] >= 1000
    assert config['extension'] == 6 and config['lib_dict'] == {}


@pytest.mark.parametrize('names', [('1', '2'), ('01', '1'), ('S1', 'NA')])
def test_numeric_and_leading_zero_sample_names_keep_their_identity(setup, names):
    """mirflowz's pandas sample table would turn 1 into a number and collapse 01 and 1; internal IDs avoid that."""
    tmp_path = setup[0]
    result, out, call = run(setup, samples=['{}={}'.format(names[0], tmp_path / 'S1_R1.fastq'),
                                            '{}={}'.format(names[1], tmp_path / 'S2_R1.fastq.gz')])
    assert result.returncode == 0, result.stderr
    assert (out / 'all_mirna_counts.tab').read_text() == 'ID\t{}\t{}\n'.format(*names)
    assert (out / 'all_pri-mir_counts.tab').read_text().splitlines()[0] == 'ID\t{}\t{}'.format(*names)
    assert sorted(p.name for p in (out / 'alignments').iterdir()) == sorted(names)
    assert (out / 'alignments' / names[1] / 'alignments_intersecting_mirna_uncollapsed_sorted.bam').read_text() == 'bam gcf_sample_2'
    samples = (out / 'mirflowz_samples.tsv').read_text().splitlines()
    assert [line.split('\t')[:2] for line in samples[1:]] == [[names[0], 'gcf_sample_1'], [names[1], 'gcf_sample_2']]


def test_nested_resources_stay_within_the_outer_budget(setup):
    """mirflowz's custom threads/mem resources are capped per rule and limited in total."""
    _, _, call = run(setup)
    args = call['args']
    set_threads = args[args.index('--set-threads') + 1:args.index('--set-resources')]
    set_resources = args[args.index('--set-resources') + 1:args.index('--resources')]
    assert 'remove_adapters=3' in set_threads and 'map_genome_segemehl=3' in set_threads
    assert all(int(item.split('=')[1]) <= 3 for item in set_threads)
    assert 'map_genome_segemehl:threads=3' in set_resources and 'map_genome_segemehl:mem=12' in set_resources
    assert 'map_transcriptome_segemehl:mem=10' in set_resources
    assert args[args.index('--resources') + 1:args.index('--resources') + 3] == ['threads=3', 'mem=12']
    assert args[args.index('--cores') + 1] == '3'


@pytest.mark.parametrize('threads, parts', [('3', '3'), ('16', '16'), ('32', '16')])
def test_oligomap_maps_in_parallel_parts(setup, threads, parts):
    """The image's oligomap runs OLIGOMAP_JOBS parts in parallel: the nested run sets it to the job's threads, at most
    16, gives both oligomap rules as many (32 threads: two oligomap steps at once) and memory for their queries, keeps the parts in the
    work directory, and gives the single-threaded remove_inferiors and filter_by_indels one thread."""
    _, _, call = run(setup, threads=threads, mem_mb='61440')
    args = call['args']
    set_threads = args[args.index('--set-threads') + 1:args.index('--set-resources')]
    set_resources = args[args.index('--set-resources') + 1:args.index('--resources')]
    assert call['env']['OLIGOMAP_JOBS'] == parts
    split_dir = Path(call['env']['OLIGOMAP_SPLIT_DIR'])
    assert split_dir.name == 'oligomap_parts' and split_dir.parent.name.startswith('.mirflowz-run-')
    assert split_dir.parent == Path(call['config']['intermediates_dir']).parent
    for rule in ('map_genome_oligomap', 'map_transcriptome_oligomap'):
        assert '{}={}'.format(rule, parts) in set_threads and '{}:threads={}'.format(rule, parts) in set_resources
    # the fixture's samples have a handful of queries: one process is about OLIGOMAP_BASE_GB
    assert 'map_genome_oligomap:mem={}'.format(math.ceil(int(parts) * 0.1)) in set_resources
    assert 'map_transcriptome_oligomap:mem={}'.format(math.ceil(int(parts) * 0.1)) in set_resources
    assert args[args.index('--resources') + 1] == 'threads={}'.format(threads)
    assert 'remove_inferiors=1' in set_threads and 'filter_by_indels=1' in set_threads
    assert 'map_genome_segemehl={}'.format(min(8, int(threads))) in set_threads


def load_script():
    spec = importlib.util.spec_from_file_location('run_mirflowz', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_oligomap_memory_follows_the_queries():
    """Every part holds all of a sample's queries, so a step's memory is the parts times one process, and a process
    grows with the queries. Measured at 16 parts on GRCh38 (#446): genome 0.45 GB per process for 45.7 thousand
    queries, 1.07 GB for 124.5 thousand, over 5.3 GB for 716 thousand; transcriptome 0.23 and 0.53 GB. The model is
    above each, and the parts are capped so that a genome step fits the budget."""
    script = load_script()
    process = lambda rule, queries: script.OLIGOMAP_BASE_GB + script.QUERY_GB[rule] * queries
    for queries, genome, transcriptome in ((45_717, 0.446, 0.233), (124_493, 1.066, 0.532), (716_066, 5.3, None)):
        assert process('map_genome_oligomap', queries) > genome
        assert transcriptome is None or process('map_transcriptome_oligomap', queries) > transcriptome
    assert script.oligomap_resources(45_717, 8, 60) == (8, {'map_genome_oligomap': (8, 5), 'map_transcriptome_oligomap': (8, 3)})
    assert script.oligomap_resources(124_493, 32, 60) == (16, {'map_genome_oligomap': (16, 22), 'map_transcriptome_oligomap': (16, 12)})
    assert script.oligomap_resources(716_066, 32, 60) == (8, {'map_genome_oligomap': (8, 59), 'map_transcriptome_oligomap': (8, 30)})
    with pytest.raises(SystemExit, match='raise quant.mirflowz.mem_gb'):
        script.oligomap_resources(716_066, 8, 6)


def test_queries_are_distinct_reads_of_at_most_30_nt(tmp_path):
    script = load_script()
    reads = ['ACGTACGTACGTACGTACGTAC', 'ACGTACGTACGTACGTACGTAC', 'TTTTACGTACGTACGTACGTAC', 'A' * 30, 'A' * 31]
    source = tmp_path / 'S1.fastq'
    source.write_text(''.join('@r{}\n{}\n+\n{}\n'.format(i, r, 'I' * len(r)) for i, r in enumerate(reads)))
    assert script.fastq_to_fasta(source, tmp_path / 'S1.fa.gz') == (5, 3)   # the 31-nt read is not a query


def test_memory_budget_below_one_gb_is_rejected(setup):
    tmp_path, _, prepared, env = setup
    result = subprocess.run([sys.executable, str(SCRIPT), '--prepared', str(prepared), '--adapter', ADAPTER, '--minimum-length', '16',
                             '--threads', '2', '--mem-mb', '512', '--output', str(tmp_path / 'out'), 'S1={}'.format(tmp_path / 'S1_R1.fastq')],
                            capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode != 0 and 'below 1 GB' in result.stderr


def test_incomplete_prepared_reference_fails_before_running(setup):
    tmp_path, _, prepared, _ = setup
    (prepared / 'segemehl_genome_index.idx').unlink()
    result, out, call = run(setup)
    assert result.returncode != 0
    assert 'lacks segemehl_genome_index.idx' in result.stderr
    assert call is None


def test_failed_run_shows_the_tool_log_and_leaves_no_outputs(setup):
    result, out, _ = run(setup, fail=True)
    assert result.returncode != 0
    assert 'mirflowz failed' in result.stderr and 'cutadapt broke' in result.stderr
    assert sorted(p.name for p in out.iterdir()) == []


def test_duplicate_sample_names_fail(setup):
    tmp_path = setup[0]
    result, _, call = run(setup, samples=['S1={}'.format(tmp_path / 'S1_R1.fastq'), 'S1={}'.format(tmp_path / 'S2_R1.fastq.gz')])
    assert result.returncode != 0 and 'not unique' in result.stderr and call is None
