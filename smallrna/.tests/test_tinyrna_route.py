"""tiny-count route (#373): selection, the reference on the mirflowz reference (TEsmall's shared genome index and the
renamed miRBase GFF), alignment and counting settings, guards, the Features Sheet, the tiny-count wrapper and the
joined table."""
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/tinyrna.py'
GENOME = 'data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus/genome'
OUT = 'data/tmp/smallrna/quant/tinyrna'

spec = importlib.util.spec_from_file_location('tinyrna', SCRIPT)
tinyrna = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tinyrna)


def test_route_plans_reference_and_counting_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,tinyrna', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule tinyrna_count:') == 4
    gff = rule_block(output, 'tinyrna_gff')
    assert ('--mirna-gff data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3 '
            '--chr-map data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus/chr_map.tsv') in gff
    assert '--output data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/mirna.gff3' in gff
    assert 'input: data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa' in rule_block(output, 'tesmall_genome')
    assert 'rule tesmall_db:' not in output  # only the shared genome index, not the TEsmall database
    assert 'features --offsets="-1,0,1" --extension 4 --length "16-28" --strand sense' in rule_block(output, 'tinyrna_features')
    assert 'tiny-collapse -i data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq -o {}/SR01/SR01'.format(OUT) in output
    align = rule_block(output, 'tinyrna_align')
    assert '-v 0 --all --best --strata --no-unal --sam --seed 0' in align and '| samtools view -b --no-PG' in align
    assert 'bowtie -x {}/bowtie_index/genome -f {}/SR0'.format(GENOME, OUT) in align
    count = rule_block(output, 'tinyrna_count')
    assert '--gff data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/mirna.gff3 --features {}/features.csv'.format(OUT) in count
    assert '--normalize-by-genomic-hits T --normalize-by-feature-hits T' in count
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'tinyrna_counts')
    for name in ('assigned_nt_len_dist.csv', 'mapped_nt_len_dist.csv', 'config'):  # declared, so a rerun restores them
        assert '{}/SR01/{}'.format(OUT, name) in count
    returncode, output = dry_run(tmp_path, 'unitas,tinyrna', target='bfq_all')
    assert returncode == 0, output
    assert 'tinyrna' not in output and 'tiny-count' not in output


def test_settings_reach_the_rules(tmp_path):
    def settings(s):
        s['quant'].setdefault('tinyrna', {}).update(mismatch=1, five_prime_offsets='0', three_prime_extension=2, length='18-25',
                                                     strand='both', normalize_by_genomic_hits=False)
    returncode, output = dry_run(tmp_path, 'tinyrna', settings, printshell=True)
    assert returncode == 0, output
    assert '-v 1 --all' in rule_block(output, 'tinyrna_align')
    assert '--offsets="0" --extension 2 --length "18-25" --strand both' in rule_block(output, 'tinyrna_features')
    assert '--normalize-by-genomic-hits F --normalize-by-feature-hits T' in rule_block(output, 'tinyrna_count')


def test_rat_uses_its_mirflowz_reference(tmp_path):
    returncode, output = dry_run(tmp_path, 'tinyrna', lambda s: s.update(organism='rattus_norvegicus'), printshell=True)
    assert returncode == 0, output
    assert '--output data/ext/tinyrna/mirbase-23.1_ensembl-104/rattus_norvegicus/mirna.gff3' in rule_block(output, 'tinyrna_gff')


def test_guards_stop_at_load(tmp_path):
    returncode, output = dry_run(tmp_path, 'tinyrna', lambda s: s.update(organism='danio_rerio'))
    assert returncode != 0 and 'organism' in output  # the mirflowz reference checks apply to this route
    returncode, output = dry_run(tmp_path / 'b', 'tinyrna', lambda s: s['quant'].setdefault('tinyrna', {}).update(mismatch=4))
    assert returncode != 0 and 'quant.tinyrna.mismatch must be 0-3' in output
    returncode, output = dry_run(tmp_path / 'c', 'tinyrna',
                                 lambda s: s['quant'].setdefault('tinyrna', {}).update(normalize_by_feature_hits='no'))
    assert returncode != 0 and 'quant.tinyrna.normalize_by_feature_hits must be true or false' in output


def test_features_sheet_has_one_rule_per_five_prime_offset():
    rows = tinyrna.features_sheet('1,-1,0', 4, '16-28', 'sense')
    assert rows[0] == tinyrna.FEATURES_HEADER
    assert [r[6] for r in rows[1:]] == ["5' anchored, -1, 4", "5' anchored, 0, 4", "5' anchored, 1, 4"]
    assert {(r[0], r[1], r[2], r[4], r[5], r[7], r[8], r[9], r[10]) for r in rows[1:]} == {
        ('all', 'all', 'miRNA', 'miRNA', '1', '', 'sense', 'all', '16-28')}
    assert tinyrna.features_sheet('0', 0, '18, 20-22', 'both')[1][10] == '18,20-22'


@pytest.mark.parametrize('args, message', [
    (('', 4, '16-28', 'sense'), 'one or more distinct'),
    (('0,0', 4, '16-28', 'sense'), 'one or more distinct'),
    (('a', 4, '16-28', 'sense'), 'integers'),
    (('0', -1, '16-28', 'sense'), "3' extension"),
    (('0', 4, '16-28', 'plus'), 'strand'),
    (('0', 4, '16 to 28', 'sense'), 'length'),
])
def test_features_sheet_rejects_bad_settings(args, message):
    with pytest.raises(SystemExit, match=message):
        tinyrna.features_sheet(*args)


def test_features_command_writes_csv(tmp_path):
    result = subprocess.run([sys.executable, str(SCRIPT), 'features', '--offsets=-1,0,1', '--extension', '4', '--length', '16-28',
                             '--strand', 'sense'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    rows = list(csv.reader(result.stdout.splitlines()))
    assert rows[0][6] == 'Overlap' and rows[1][6] == "5' anchored, -1, 4" and len(rows) == 4


def fake_tiny_count(bin_dir, body):
    """A tiny-count stand-in on PATH; body is shell run in tiny-count's working directory."""
    bin_dir.mkdir(exist_ok=True)
    tool = bin_dir / 'tiny-count'
    tool.write_text('#!/bin/sh\n' + body)
    tool.chmod(0o755)
    return dict(os.environ, PATH='{}:{}'.format(bin_dir, os.environ['PATH']))


def run_count(tmp_path, env):
    (tmp_path / 'S1.bam').touch()
    (tmp_path / 'mirna.gff3').touch()
    (tmp_path / 'features.csv').write_text('x\n')
    return subprocess.run([sys.executable, str(SCRIPT), 'count', '--sample', 'S1', '--alignments', str(tmp_path / 'S1.bam'),
                           '--gff', str(tmp_path / 'mirna.gff3'), '--features', str(tmp_path / 'features.csv'),
                           '--output', str(tmp_path / 'out')], capture_output=True, text=True, env=env)


def test_count_keeps_tiny_count_tables_under_fixed_names(tmp_path):
    tables = ' '.join('tinyrna_2026-01-01_00-00-00_{}.csv'.format(n) for n in tinyrna.RUN_TABLES)
    tables += ' ' + ' '.join('tinyrna_2026-01-01_00-00-00_S1_rep_1_{}.csv'.format(n) for n in tinyrna.LIBRARY_TABLES)
    env = fake_tiny_count(tmp_path / 'bin', 'cat paths.yml samples.csv > args.txt; echo "$@" >> args.txt\n'
                          'd=tiny-count_2026-01-01_00-00-00; mkdir -p $d/config && cp args.txt features.csv $d/config/\n'
                          'cd $d && touch {}\n'.format(tables))
    result = run_count(tmp_path, env)
    assert result.returncode == 0, result.stderr
    out = tmp_path / 'out'
    assert sorted(p.name for p in out.glob('*.csv')) == sorted('{}.csv'.format(n) for n in tinyrna.RUN_TABLES + tinyrna.LIBRARY_TABLES)
    record = (out / 'config/args.txt').read_text()
    for key in ('samples_csv: samples.csv', 'features_csv: features.csv', 'alias: [Name]', 'ebwt:', 'plot_style_sheet:'):
        assert key in record  # every Paths File key, as tiny-count requires
    assert '{},S1,1,,'.format(tmp_path / 'S1.bam') in record and '-pf paths.yml -o tinyrna -ng T -nf T' in record
    assert (out / 'config/features.csv').exists() and not (out / 'work').exists()


def test_count_fails_without_tiny_count_config_folder(tmp_path):
    tables = ' '.join('tinyrna_2026-01-01_00-00-00_{}.csv'.format(n) for n in tinyrna.RUN_TABLES)
    tables += ' ' + ' '.join('tinyrna_2026-01-01_00-00-00_S1_rep_1_{}.csv'.format(n) for n in tinyrna.LIBRARY_TABLES)
    env = fake_tiny_count(tmp_path / 'bin', 'd=tiny-count_2026-01-01_00-00-00; mkdir -p $d && cd $d && touch {}\n'.format(tables))
    result = run_count(tmp_path, env)
    assert result.returncode != 0 and 'tiny-count wrote no config folder' in result.stderr


def test_count_fails_when_tiny_count_reports_an_error_but_exits_0(tmp_path):
    env = fake_tiny_count(tmp_path / 'bin', 'echo "Traceback (most recent call last):"; echo "AssertionError: bad"; exit 0\n')
    result = run_count(tmp_path, env)
    assert result.returncode != 0 and 'tiny-count failed for S1' in result.stderr


def test_count_fails_when_a_table_is_missing(tmp_path):
    env = fake_tiny_count(tmp_path / 'bin', 'mkdir tiny-count_x && touch tiny-count_x/tinyrna_x_feature_counts.csv\n')
    result = run_count(tmp_path, env)
    assert result.returncode != 0 and 'wrote no counts_by_rule' in result.stderr


def test_join_combines_samples_and_checks_the_library_column(tmp_path):
    (tmp_path / 'a.csv').write_text('Feature ID,Classifier,Feature Name,S1_rep_1\nMIMAT1,miRNA,mmu-miR-1,5.0\nMIMAT1_1,miRNA,mmu-miR-1,2.5\n')
    (tmp_path / 'b.csv').write_text('Feature ID,Classifier,Feature Name,S2_rep_1\nMIMAT1,miRNA,mmu-miR-1,7.0\nMIMAT2,miRNA,mmu-miR-2,1.0\n')
    command = [sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(tmp_path / 'a.csv'), str(tmp_path / 'b.csv')]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ['Feature ID\tClassifier\tFeature Name\tS1\tS2', 'MIMAT1\tmiRNA\tmmu-miR-1\t5.0\t7.0',
                                          'MIMAT1_1\tmiRNA\tmmu-miR-1\t2.5\t0', 'MIMAT2\tmiRNA\tmmu-miR-2\t0\t1.0']
    command[4:6] = ['S2', 'S1']
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0 and 'library column S1_rep_1 is not sample S2' in result.stderr
