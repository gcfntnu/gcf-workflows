"""TEsmall route (#374): selection, database planning on the mirflowz reference, run settings, guards and the joined table."""
from pathlib import Path
import re
import subprocess
import sys

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
REF = 'data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus'
DB = REF + '/db/genomes/gcf'


def test_route_plans_database_and_runs_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule tesmall_quant:') == 4
    genome = rule_block(output, 'tesmall_genome')
    assert 'input: data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa' in genome  # the mirflowz reference genome
    assert re.search(rf'bowtie-build -q --threads \d+ {DB}/sequence/genome.fa {DB}/sequence/bowtie_index/genome', genome)
    db = rule_block(output, 'tesmall_db')
    assert '--chr-map data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus/chr_map.tsv' in db
    assert '--pirna-clusters data/ext/unitas/05-11-2026/UNITAS_refdump_mus_musculus/piRclDB' in db
    assert f'--trna-archive {REF}/sources/mm39-tRNAs.tar.gz --rmsk {REF}/sources/rmsk.txt.gz' in db
    assert 'https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz' in rule_block(output, 'tesmall_gtrnadb')
    assert 'goldenPath/mm39/database/rmsk.txt.gz' in rule_block(output, 'tesmall_rmsk')
    assert 'id=BK000964.3' in rule_block(output, 'tesmall_rdna')
    run = rule_block(output, 'tesmall_quant')
    assert 'TEsmall --no-trim -f "$fastq" -l SR0' in run and f'-g gcf --dbfolder "$db"' in run
    assert '-m 15 -M 36 --mismatch 0 --maxaln 100' in run  # the kit's filtered minimum (NEXTflex v4 fixture)
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in run
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', target='bfq_all')
    assert returncode == 0, output
    assert 'tesmall' not in output


def test_rat_database_uses_rnor6_sources(tmp_path):
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s.update(organism='rattus_norvegicus'), printshell=True)
    assert returncode == 0, output
    assert 'Rnorv6/rn6-tRNAs.tar.gz' in rule_block(output, 'tesmall_gtrnadb') and 'id=NR_046239.1' in rule_block(output, 'tesmall_rdna')
    assert 'data/ext/tesmall/mirbase-23.1_ensembl-104/rattus_norvegicus/' in rule_block(output, 'tesmall_db')


def test_guards_stop_at_load(tmp_path):
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s.update(organism='danio_rerio'))
    assert returncode != 0 and 'organism' in output
    returncode, output = dry_run(tmp_path / 'b', 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(max_length=10))
    assert returncode != 0 and 'quant.tesmall.max_length 10 is below the minimum read length' in output


def test_dotted_sample_names_stop_at_load(tmp_path):
    def rename(settings):
        settings['samples'] = {'S.1' if k == 'SR01' else k: v for k, v in settings['samples'].items()}
    returncode, output = dry_run(tmp_path, 'tesmall', rename)
    assert returncode != 0 and 'TEsmall cuts sample labels at dots' in output


def test_counts_join_samples_and_reject_unexpected_tables(tmp_path):
    (tmp_path / 'a.txt').write_text('fid\tftype\tS1\nmmu-miR-1a-3p\tmiRNA\t5\nL1:L1:LINE\tsense_TE\t2\n')
    (tmp_path / 'b.txt').write_text('fid\tftype\tS2\nmmu-miR-1a-3p\tmiRNA\t7\nmmu-mir-1a-1\thairpin\t1\n')
    script = [sys.executable, str(ROOT / 'smallrna/rules/quant/scripts/tesmall_counts.py')]
    result = subprocess.run(script + ['--samples', 'S1', 'S2', '--tables', str(tmp_path / 'a.txt'), str(tmp_path / 'b.txt')],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ['ftype\tfid\tS1\tS2', 'hairpin\tmmu-mir-1a-1\t0\t1', 'miRNA\tmmu-miR-1a-3p\t5\t7',
                                          'sense_TE\tL1:L1:LINE\t2\t0']
    result = subprocess.run(script + ['--samples', 'S2', '--tables', str(tmp_path / 'a.txt')], capture_output=True, text=True)
    assert result.returncode != 0 and 'expected columns fid, ftype, S2' in result.stderr
