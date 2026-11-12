"""Dry-run contract for the mirflowz route: selection, reference planning, assembly checks and read settings."""
import pytest

from test_quant_methods import dry_run, rule_block

REF_116 = 'data/ext/mirflowz/mirbase-23.1_ensembl-116'


def organism(name):
    return lambda s: s.update(organism=name)


def test_combined_route_plans_both_quantifiers_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,mirflowz')
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and 'rule mirflowz_quant:' in output
    returncode, output = dry_run(tmp_path, 'unitas,mirflowz', target='bfq_all')
    assert returncode == 0, output
    assert 'mirflowz' not in output


def test_mouse_reference_reuses_the_project_ensembl_files(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirflowz', printshell=True)
    assert returncode == 0, output
    prepare = rule_block(output, 'mirflowz_prepare')
    assert ('--genome data/ext/ensembl/release-116/mus_musculus/GRCm39/fasta/genome.fa '
            '--gtf data/ext/ensembl/release-116/mus_musculus/GRCm39/anno/mus_musculus.GRCm39.116.gtf '
            '--mirna-gff data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3') in prepare
    assert f'--output {REF_116}/mus_musculus/prepared --record {REF_116}/mus_musculus/mirflowz_reference.tsv' in prepare
    assert '--assembly GRCm39 --output' in rule_block(output, 'mirflowz_chr_map')
    assert 'rule mirflowz_ensembl_genome:' not in output


def test_rat_reference_downloads_ensembl_104_and_removes_the_misplaced_contig(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirflowz', organism('rattus_norvegicus'), printshell=True)
    assert returncode == 0, output
    assert ('https://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz'
            in rule_block(output, 'mirflowz_ensembl_genome'))
    assert 'Rattus_norvegicus.Rnor_6.0.104.gtf.gz' in rule_block(output, 'mirflowz_ensembl_gtf')
    chr_map = rule_block(output, 'mirflowz_chr_map')
    assert '--assembly Rnor_6.0 --remove 901866519 --output data/ext/mirflowz/mirbase-23.1_ensembl-104/rattus_norvegicus/' in chr_map
    assert 'rule ensembl_genome:' not in output


def test_custom_reference_replaces_the_downloads(tmp_path):
    refs = tmp_path / 'refs'
    refs.mkdir(parents=True)
    custom = {'name': 'mm-chr19'}
    for key, name in (('genome', 'chr19.fa.gz'), ('gtf', 'chr19.gtf'), ('mirna_gff', 'chr19.gff3')):
        (refs / name).touch()
        custom[key] = str(refs / name)
    returncode, output = dry_run(tmp_path, 'mirflowz', lambda s: s.setdefault('db', {}).setdefault('mirflowz', {}).update(custom=custom),
                                 printshell=True)
    assert returncode == 0, output
    prepare = rule_block(output, 'mirflowz_prepare')
    assert f'--genome {refs}/chr19.fa.gz --gtf {refs}/chr19.gtf --mirna-gff {refs}/chr19.gff3' in prepare
    assert '--output data/ext/mirflowz/custom-mm-chr19/mus_musculus/prepared' in prepare
    assert '--assembly' not in rule_block(output, 'mirflowz_chr_map')
    assert 'rule ensembl_genome:' not in output and 'rule mirbase_gff:' not in output


def test_read_settings_follow_the_kit(tmp_path):
    kit = 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4'
    returncode, output = dry_run(tmp_path, 'mirflowz', lambda s: s.update(libprepkit=kit), printshell=True)
    assert returncode == 0, output
    quant = rule_block(output, 'mirflowz_quant')
    assert '--adapter TGGAATTCTCGGGTGCCAAGG --minimum-length 15 --threads' in quant
    assert 'SR01=data/tmp/smallrna/filter/fastp/SR01_R1.fastq' in quant


@pytest.mark.parametrize('update, message', [
    (lambda s: s.update(organism='rattus_norvegicus') or
     s.setdefault('db', {}).setdefault('mirflowz', {}).setdefault('ensembl_release', {}).update(rattus_norvegicus='116'),
     'db.mirflowz.ensembl_release 116 for rattus_norvegicus is on GRCr8, but miRBase 23.1 GFF coordinates are on Rnor_6.0'),
    (lambda s: s.setdefault('db', {}).setdefault('mirflowz', {}).update(custom={'name': 'x', 'genome': '/g.fa'}),
     'db.mirflowz.custom: set name, genome, gtf and mirna_gff together; missing gtf, mirna_gff'),
    (lambda s: s.setdefault('db', {}).setdefault('mirbase', {}).update(release='22'),
     'mirflowz: no GFF assembly known for db.mirbase.release 22')])
def test_invalid_mirflowz_settings_fail_at_load_only_when_selected(tmp_path, update, message):
    returncode, output = dry_run(tmp_path, 'mirflowz', update)
    assert returncode != 0
    assert message in output
    returncode, output = dry_run(tmp_path / 'unitas', 'unitas', update)
    assert message not in output
