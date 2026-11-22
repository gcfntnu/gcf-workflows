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
import yaml

from test_quant_methods import MIRBASE_GFF, MIRBASE_GFF_ID, dry_run, rule_block

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
    assert (f'--mirna-gff {MIRBASE_GFF} '
            f'--chr-map data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus/gff-{MIRBASE_GFF_ID}/chr_map.tsv') in gff
    assert f'--output data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/gff-{MIRBASE_GFF_ID}/mirna.v2.gff3' in gff  # named by the GFF3 (#440)
    assert 'input: data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa' in rule_block(output, 'tesmall_genome')
    assert 'rule tesmall_db:' not in output  # only the shared genome index, not the TEsmall database
    assert 'features --offsets=-1,0,1 --extension 4 --length 16-28 --strand sense' in rule_block(output, 'tinyrna_features')
    assert 'tiny-collapse -i data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq -o {}/SR01/SR01'.format(OUT) in output
    align = rule_block(output, 'tinyrna_align')
    assert '-v 0 --all --best --strata --no-unal --sam --seed 0' in align and '| samtools view -b --no-PG' in align
    assert '--reorder' not in align and ' -p ' not in align   # one thread: --reorder with threads deadlocks Bowtie 1.3.1
    assert 'mem_mb=2000' in align.replace(' ', '')             # the index in memory, so --resources caps parallel jobs
    assert 'bowtie -x {}/bowtie_index/genome -f {}/SR0'.format(GENOME, OUT) in align
    count = rule_block(output, 'tinyrna_count')
    assert f'--gff data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/gff-{MIRBASE_GFF_ID}/mirna.v2.gff3 --features {OUT}/features.csv' in count
    assert '--normalize-by-genomic-hits T --normalize-by-feature-hits T' in count
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'tinyrna_counts')
    for name in ('assigned_nt_len_dist.csv', 'mapped_nt_len_dist.csv', 'config'):  # declared, so a rerun restores them
        assert all(f'{OUT}/{sample}/{name}' in output for sample in ('SR01', 'SR02', 'SR03', 'SR04'))
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
    assert '--offsets=0 --extension 2 --length 18-25 --strand both' in rule_block(output, 'tinyrna_features')
    assert '--normalize-by-genomic-hits F --normalize-by-feature-hits T' in rule_block(output, 'tinyrna_count')


def test_settings_reach_the_script_as_written(tmp_path):
    """Each Features Sheet setting is one shell word, so quotes or spaces in a project value reach the validator."""
    def settings(s):
        s['quant'].setdefault('tinyrna', {}).update(strand='"sense"', length='16 to 28')
    returncode, output = dry_run(tmp_path, 'tinyrna', settings, printshell=True)
    assert returncode == 0, output
    command = rule_block(output, 'tinyrna_features')
    assert """--length '16 to 28' --strand '"sense"'""" in command
    result = subprocess.run([sys.executable, str(SCRIPT), 'features', '--offsets=0', '--extension', '4', '--length', '16-28',
                             '--strand', '"sense"'], capture_output=True, text=True)
    assert result.returncode != 0 and """strand must be one of sense, antisense, both, got '"sense"'""" in result.stderr


def test_rat_uses_its_mirflowz_reference(tmp_path):
    returncode, output = dry_run(tmp_path, 'tinyrna', lambda s: s.update(organism='rattus_norvegicus'), printshell=True)
    assert returncode == 0, output
    assert f'--output data/ext/tinyrna/mirbase-23.1_ensembl-104/rattus_norvegicus/gff-{MIRBASE_GFF_ID}/mirna.v2.gff3' in rule_block(output, 'tinyrna_gff')


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


@pytest.mark.parametrize('folder', ['project # one', 'project \U0001F9EC'])
def test_count_keeps_tiny_count_tables_under_fixed_names(tmp_path, folder):
    tables = ' '.join('tinyrna_2026-01-01_00-00-00_{}.csv'.format(n) for n in tinyrna.RUN_TABLES)
    tables += ' ' + ' '.join('tinyrna_2026-01-01_00-00-00_S1_rep_1_{}.csv'.format(n) for n in tinyrna.LIBRARY_TABLES)
    tmp_path = tmp_path / folder  # a path outside the Basic Multilingual Plane must survive paths.yml
    tmp_path.mkdir()
    env = fake_tiny_count(tmp_path / 'bin', 'cat samples.csv > args.txt; echo "$@" >> args.txt\n'
                          'd=tiny-count_2026-01-01_00-00-00; mkdir -p $d/config && cp args.txt paths.yml features.csv $d/config/\n'
                          'cd $d && touch {}\n'.format(tables))
    result = run_count(tmp_path, env)
    assert result.returncode == 0, result.stderr
    out = tmp_path / 'out'
    assert sorted(p.name for p in out.glob('*.csv')) == sorted('{}.csv'.format(n) for n in tinyrna.RUN_TABLES + tinyrna.LIBRARY_TABLES)
    record = (out / 'config/args.txt').read_text()
    paths = yaml.safe_load((out / 'config/paths.yml').read_text(encoding='utf-8'))
    assert paths == {'samples_csv': 'samples.csv', 'features_csv': 'features.csv',
                     'gff_files': [{'path': str(tmp_path / 'mirna.gff3'), 'alias': ['Name']}], 'run_directory': 'run_directory',
                     'tmp_directory': None, 'ebwt': None, 'reference_genome_files': None, 'adapter_fasta': None, 'plot_style_sheet': None}
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


def test_direct_count_target_respects_normalization_and_validates_settings(tmp_path):
    def settings(config):
        config['quant'].setdefault('tinyrna', {}).update(normalize_by_genomic_hits=False, normalize_by_feature_hits=False)
    target = f'{OUT}/SR01/feature_counts.csv'
    code, output = dry_run(tmp_path, 'unitas', settings, printshell=True, target=target)
    assert code == 0, output
    assert '--normalize-by-genomic-hits F --normalize-by-feature-hits F' in rule_block(output, 'tinyrna_count')
    code, output = dry_run(tmp_path, 'unitas', lambda c: c['quant']['tinyrna'].update(normalize_by_genomic_hits='no'), target=target)
    assert code != 0 and 'normalize_by_genomic_hits must be true or false' in output


def test_join_rejects_duplicate_keys_in_a_native_table(tmp_path):
    native = tmp_path / 'counts.csv'
    native.write_text('Feature ID,Classifier,Feature Name,S1_rep_1\nMIMAT1,miRNA,mir-1,2\nMIMAT1,miRNA,mir-1,3\n')
    with pytest.raises(SystemExit, match='duplicate feature ID and classifier'):
        tinyrna.read_counts(native, 'S1')


def test_seed_alignment_mode(tmp_path):
    """quant.tinyrna.alignment: seed (#417): Bowtie's seed mode, -n mismatch -l seed_length -e 40 per mismatch allowed in
    the read (FASTA reads score quality 40 per base); the default keeps the end-to-end command; bad values stop at load."""
    seed = lambda **kw: (lambda s: s['quant'].setdefault('tinyrna', {}).update(alignment='seed', **kw))
    returncode, output = dry_run(tmp_path, 'tinyrna', seed(seed_length=16, seed_total=2), printshell=True)
    assert returncode == 0, output
    align = rule_block(output, 'tinyrna_align')
    assert '-n 0 -l 16 -e 80 --nomaqround --all --best --strata --no-unal --sam --seed 0' in align and ' -v ' not in align
    returncode, output = dry_run(tmp_path / 'default', 'tinyrna', printshell=True)
    assert returncode == 0 and '-f {}/SR0'.format(OUT) in output and '-v 0 --all --best' in rule_block(output, 'tinyrna_align')
    # requested directly, without quant.method selecting tiny-count: the same options
    returncode, output = dry_run(tmp_path / 'direct', 'unitas', seed(), printshell=True, target=OUT + '/tinyrna_counts.tsv')
    assert returncode == 0 and '-n 0 -l 18 -e 40 --nomaqround' in rule_block(output, 'tinyrna_align'), output
    for value, message in ((dict(seed_length=4), 'quant.tinyrna.seed_length must be a whole number from 5 to 32'),
                           (dict(seed_total=0, mismatch=1), 'quant.tinyrna.seed_total must be a whole number from quant.tinyrna.mismatch to 3')):
        returncode, output = dry_run(tmp_path / str(len(message)), 'tinyrna', seed(**value))
        assert returncode != 0 and message in output, output
    returncode, output = dry_run(tmp_path / 'mode', 'tinyrna', lambda s: s['quant'].setdefault('tinyrna', {}).update(alignment='local'))
    assert returncode != 0 and "quant.tinyrna.alignment must be end_to_end or seed, got 'local'" in output
