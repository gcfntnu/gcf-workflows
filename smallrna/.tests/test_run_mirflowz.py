"""run_mirflowz.py: mirflowz mapping and quantification on filtered reads with a linked reference; snakemake is mocked."""
import gzip
import json
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

# Records the call, the sample table and the reads as mirflowz sees them, and writes mirflowz's outputs.
FAKE_SNAKEMAKE = '''#!/usr/bin/env python3
import gzip, json, os, sys
from pathlib import Path
args = sys.argv[1:]
config = {k: json.loads(v) for k, v in (l.split(': ', 1) for l in open(args[args.index('--configfile') + 1]).read().splitlines())}
work = Path(args[args.index('--directory') + 1])
inter = Path(config['intermediates_dir'])
rows = [l.split('\\t') for l in open(config['samples']).read().splitlines()]
call = {'args': args, 'config': config, 'snakefile': Path(args[args.index('-s') + 1]).read_text(), 'rows': rows,
        'links': {p.name: os.readlink(p) for p in inter.iterdir() if p.is_symlink()},
        'reads': {r[0]: gzip.open(r[1], 'rt').read() for r in rows[1:]}}
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
'''


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


def run(setup, fail=False, prepared=None, samples=None):
    tmp_path, _, default_prepared, env = setup
    env = dict(env, FAIL='1' if fail else '')
    out = tmp_path / 'quant/mirflowz'
    samples = samples or ['S1={}'.format(tmp_path / 'S1_R1.fastq'), 'S2={}'.format(tmp_path / 'S2_R1.fastq.gz')]
    result = subprocess.run([sys.executable, str(SCRIPT), '--prepared', str(prepared or default_prepared), '--adapter', ADAPTER,
                             '--minimum-length', '16', '--threads', '3', '--output', str(out), *samples],
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
    assert call['args'][-1] == 'gcf_mirflowz_all' and call['args'][call['args'].index('--cores') + 1] == '3'
    assert call['links'] == {name: str((prepared / name).resolve()) for name in PREPARED}
    assert [row[0] for row in call['rows']] == ['sample', 'S1', 'S2'] and all(row[2:] == [ADAPTER, 'fa'] for row in call['rows'][1:])
    # FASTA input skips mirflowz's quality filter; names and sequences are kept.
    assert call['reads'] == {'S1': '>r1\nTGAGGTAGTAGGTTGTATAGTT\n', 'S2': '>r2\nTAGCTTATCAGACTGATGTTGA\n'}
    assert sorted(p.name for p in out.iterdir()) == ['alignments', 'all_mirna_counts.tab', 'all_pri-mir_counts.tab',
                                                    'mirflowz_config.yaml', 'mirflowz_samples.tsv']
    assert (out / 'all_mirna_counts.tab').read_text() == 'ID\tS1\tS2\n'
    assert (out / 'alignments/S2/alignments_intersecting_mirna_uncollapsed_sorted.bam').read_text() == 'bam S2'
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
