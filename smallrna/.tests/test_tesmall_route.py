"""TEsmall route (#374): selection, database planning on the mirflowz reference (one database per UNITAS dump, the
genome shared), run settings, guards and the joined table."""
from pathlib import Path
import re
import subprocess
import sys

from test_quant_methods import configure, dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
REF = 'data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus'
GENOME = REF + '/genome'
DUMP = REF + '/unitas-05-11-2026'
DB = DUMP + '/db/genomes/gcf'
CLASSES = ('structural_RNA', 'miRNA', 'hairpin', 'exon', 'TE', 'intron', 'piRNA_cluster')
INDEX = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')


def test_route_plans_database_and_runs_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule tesmall_quant:') == 4
    genome = rule_block(output, 'tesmall_genome')
    assert 'input: data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa' in genome  # the mirflowz reference genome
    assert re.search(rf'bowtie-build -q --threads \d+ {GENOME}/genome.fa {GENOME}/bowtie_index/genome', genome)
    assert f'ln -sfr' in rule_block(output, 'tesmall_genome_links') and f'link={DB}/sequence/${{target#{GENOME}/}}' in output
    db = rule_block(output, 'tesmall_db')
    assert '--chr-map data/ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus/chr_map.tsv' in db
    assert '--pirna-clusters data/ext/unitas/05-11-2026/UNITAS_refdump_mus_musculus/piRclDB' in db
    assert f'--trna-archive {REF}/sources/mm39-tRNAs.tar.gz --rmsk {REF}/sources/rmsk.txt.gz' in db
    assert f'--genome-index {GENOME}/bowtie_index/genome' in db and f'--record {DUMP}/tesmall_reference.tsv' in db
    assert 'https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz' in rule_block(output, 'tesmall_gtrnadb')
    assert 'goldenPath/mm39/database/rmsk.txt.gz' in rule_block(output, 'tesmall_rmsk')
    assert 'id=BK000964.3' in rule_block(output, 'tesmall_rdna')
    run = rule_block(output, 'tesmall_quant')
    assert 'TEsmall --no-trim -f "$fastq" -l SR0' in run and f'-g gcf --dbfolder "$db"' in run
    assert '-m 15 -M 36 --mismatch 0 --maxaln 100' in run  # the kit's filtered minimum (NEXTflex v4 fixture)
    assert '-o miRNA hairpin structural_RNA exon TE intron piRNA_cluster' in run  # miRNA first by default
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in run
    assert f'db=$(realpath {DUMP}/db)' in run
    for path in (f'{DB}/annotation/piRNA_cluster.bed', f'{DB}/sequence/bowtie_index/tDNA.rev.2.ebwt', f'{DB}/sequence/genome.fa.fai'):
        assert path in run.split('output:', 1)[0]  # every file TEsmall checks at startup is an input
    returncode, output = dry_run(tmp_path, 'unitas,tesmall', target='bfq_all')
    assert returncode == 0, output
    assert 'tesmall' not in output


def test_each_unitas_dump_has_its_own_database_on_the_shared_genome(tmp_path):
    release = lambda s: s.setdefault('db', {}).setdefault('unitas', {}).update(release='01-01-2027')
    returncode, output = dry_run(tmp_path, 'tesmall', release, printshell=True)
    assert returncode == 0, output
    db = rule_block(output, 'tesmall_db')
    assert '--pirna-clusters data/ext/unitas/01-01-2027/UNITAS_refdump_mus_musculus/piRclDB' in db
    assert f'--record {REF}/unitas-01-01-2027/tesmall_reference.tsv' in db and f'{DUMP}/' not in output
    assert f'--genome-index {GENOME}/bowtie_index/genome' in db


def test_a_partial_database_is_completed_before_tesmall_runs(tmp_path):
    """With the genome, its links and the database in place except one BED, quantification first rebuilds the database."""
    configure(tmp_path, 'tesmall')
    present = ([f'{GENOME}/genome.fa', f'{GENOME}/genome.fa.fai'] + [f'{GENOME}/bowtie_index/genome.{s}' for s in INDEX] +
               [f'{DB}/sequence/genome.fa', f'{DB}/sequence/genome.fa.fai'] + [f'{DB}/sequence/bowtie_index/genome.{s}' for s in INDEX] +
               [f'{DB}/annotation/{c}.bed' for c in CLASSES if c != 'TE'] +
               [f'{DB}/sequence/{n}.fa{e}' for n in ('rDNA', 'tDNA') for e in ('', '.fai')] +
               [f'{DB}/sequence/bowtie_index/{n}.{s}' for n in ('rDNA', 'tDNA') for s in INDEX] + [f'{DUMP}/tesmall_reference.tsv'])
    for name in present:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).touch()
    returncode, output = dry_run(tmp_path, 'tesmall')
    assert returncode == 0, output
    assert 'rule tesmall_db:' in output and f'{DB}/annotation/TE.bed' in output


def test_no_shared_minimum_length_gives_tesmall_none(tmp_path):
    fastp = lambda s: s.update(filter={'trim': {'quantifier': 'fastp', 'fastp': {'args': '-L'}}})
    returncode, output = dry_run(tmp_path, 'tesmall', fastp, printshell=True)
    assert returncode == 0, output
    assert '-m 1 -M 36 ' in rule_block(output, 'tesmall_quant')


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
    returncode, output = dry_run(tmp_path / 'c', 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(order='miRNA,hairpin'))
    assert returncode != 0 and 'quant.tesmall.order must list each TEsmall class once' in output
    # TEsmall gets the resolved database path: a symbolic link to a folder named genomes... is caught at load.
    (tmp_path / 'genomes-store').mkdir()
    (tmp_path / 'ext-link').symlink_to(tmp_path / 'genomes-store')
    returncode, output = dry_run(tmp_path / 'd', 'tesmall', lambda s: s.update(ext_dir=str(tmp_path / 'ext-link')))
    assert returncode != 0 and 'contains "genomes", which TEsmall strips from it' in output


def test_tesmall_own_order_can_be_set(tmp_path):
    order = 'structural_RNA,miRNA,hairpin,exon,TE,intron,piRNA_cluster'
    returncode, output = dry_run(tmp_path, 'tesmall', lambda s: s['quant'].setdefault('tesmall', {}).update(order=order), printshell=True)
    assert returncode == 0, output
    assert '-o structural_RNA miRNA hairpin exon TE intron piRNA_cluster' in rule_block(output, 'tesmall_quant')


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
