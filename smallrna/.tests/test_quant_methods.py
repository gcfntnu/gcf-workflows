"""Dry-run contract for the small-RNA `quant.method` selection of `quant_all`."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def dry_run(tmp_path, method, update=None, printshell=False, target='quant_all'):
    shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    settings.setdefault('quant', {})['method'] = method
    if update:
        update(settings)
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                      f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', *(['-p'] if printshell else []), '--cores', '1', target],
                            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    return result.returncode, result.stdout + result.stderr


@pytest.mark.parametrize('method, rule', [('unitas', 'unitas_mirtable'), ('mirge3', 'mirge3_quant'),
                                          ('unitas,mirge3', 'mirge3_quant')])
def test_supported_methods_plan_their_quantifier(tmp_path, method, rule):
    returncode, output = dry_run(tmp_path, method)
    assert returncode == 0, output
    assert f'rule {rule}:' in output


@pytest.mark.parametrize('method', ['mirge', 'quickmirseq', 'seqbuster', 'unitass'])
def test_removed_or_unknown_methods_fail(tmp_path, method):
    returncode, output = dry_run(tmp_path, method)
    assert returncode != 0
    assert f'quant.method: unsupported value(s) {method}' in output


def mirge3_command(output):
    return next(line for line in output.splitlines() if 'miRge3.0 ' in line)


@pytest.mark.parametrize('kit, umi', [('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', '--uniq-mol-ids 4,4'),
                                      ('Bioo Scientific NEXTflex Small RNA-Seq Kit v4', None),
                                      ('TruSeq Small RNA Library Prep Kit', None)])
def test_mirge3_unique_molecule_ids_follow_the_kit(tmp_path, kit, umi):
    returncode, output = dry_run(tmp_path, 'mirge3', lambda s: s.update(libprepkit=kit), printshell=True)
    assert returncode == 0, output
    command = ' '.join(mirge3_command(output).split())
    if umi:
        assert umi in command
    else:
        assert '--uniq-mol-ids' not in command


def test_mirge3_reads_merged_fastq_from_every_lane(tmp_path):
    def two_lanes(settings):
        sample = next(iter(settings['samples']))
        first = settings['samples'][sample]['R1']
        settings['samples'][sample]['R1'] = first + ',' + first.replace('_R1', '_R1_lane2')
        lane2 = tmp_path / settings['fastq_dir'] / first.replace('_R1', '_R1_lane2')
        shutil.copy(tmp_path / settings['fastq_dir'] / first, lane2)
    returncode, output = dry_run(tmp_path, 'mirge3', two_lanes, printshell=True)
    assert returncode == 0, output
    links = [line for line in output.splitlines() if 'ln -sr ' in line and '/mirge3/fastq/' in line]
    assert links and all('/filter/fastq/' in line for line in links)
    merge = [line for line in output.splitlines() if '_R1_lane2' in line]
    assert merge, 'the second lane was not merged'


def test_mirge3_image_and_spikein_come_from_configuration(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirge3', printshell=True)
    assert returncode == 0, output
    assert '--spikeIn' not in mirge3_command(output)
    assert 'quay.io/biocontainers/mirge3' not in (ROOT / 'smallrna/rules/quant/mirge3.smk').read_text()


def rule_block(output, rule):
    return output.split(f'rule {rule}:', 1)[1].split('\nrule ', 1)[0]


@pytest.mark.parametrize('ref, fasta', [(None, []), ('ERCC92', ['ERCC92']), ('smallrna_calibrators', ['smallrna_calibrators']),
                                        ('smallrna_calibrators, ERCC92', ['smallrna_calibrators', 'ERCC92'])])
def test_mirge3_spikein_index_follows_the_configured_sets(tmp_path, ref, fasta):
    """The spike-in index is built from the gcfdb sets of filter.spikein.ref before miRge3 runs, or removed."""
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref=ref)
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    spikein = rule_block(output, 'mirge3_spikein')
    command = next(line for line in spikein.splitlines() if 'mirge3_spikein.py' in line)
    assert '--index-dir data/ext/mirge3/09-11-2026/mouse/index.Libs --organism-name mouse' in command
    assert command.split('--output data/tmp/smallrna/quant/mirge3/spikein.tsv', 1)[1].split() == \
        [f'data/ext/spikein/fasta/{name}.fa' for name in fasta]
    assert 'data/tmp/smallrna/quant/mirge3/spikein.tsv' in rule_block(output, 'mirge3_quant').split('output:', 1)[0]
    assert ('--spikeIn' in mirge3_command(output)) == bool(fasta)


def test_unknown_spikein_set_fails_at_load(tmp_path):
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref='ERCC')
    returncode, output = dry_run(tmp_path, 'mirge3', update)
    assert returncode != 0
    assert "filter.spikein.ref: unsupported value 'ERCC'" in output


@pytest.mark.parametrize('db, name', [('mirbase', 'miRBase'), ('mirgenedb', 'MirGeneDB'), ('MirGeneDB', 'MirGeneDB')])
def test_mirge3_library_release_and_database_selection(tmp_path, db, name):
    update = lambda s: s['quant'].setdefault('mirge3', {}).update(db=db)
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    library = rule_block(output, 'mirge3_library')
    assert '--release-url mirge3/09-11-2026 --release 09-11-2026 --species mus_musculus --organism-name mouse' in library
    assert '--output data/ext/mirge3/09-11-2026/verify/miRge3_lib_mus_musculus.tsv' in library
    reference = rule_block(output, 'mirge3_reference')
    assert (f'--library data/ext/mirge3/09-11-2026/mouse --mirna-db {db.lower()} --release 09-11-2026 '
            '--verification data/ext/mirge3/09-11-2026/verify/miRge3_lib_mus_musculus.tsv') in reference
    command = ' '.join(mirge3_command(output).split())
    assert f'--mir-DB {name} --libraries-path data/ext/mirge3/09-11-2026 --organism-name mouse' in command
    assert 'data/tmp/smallrna/quant/mirge3/reference.tsv' in rule_block(output, 'quant_all')


def test_explicit_mirge3_release_selects_its_own_directory(tmp_path):
    update = lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release='01-12-2026')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert '--release-url mirge3/01-12-2026' in rule_block(output, 'mirge3_library')
    assert '--libraries-path data/ext/mirge3/01-12-2026 ' in mirge3_command(output)


def test_sourceforge_library_stays_at_its_old_place(tmp_path):
    update = lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release='sourceforge')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert 'rule mirge3_library:' not in output
    assert 'data/ext/mirge3/mus_musculus/mouse/index.Libs/mouse_snorna.3.ebwt' in rule_block(output, 'mirge3_prebuild')
    reference = rule_block(output, 'mirge3_reference')
    assert '--library data/ext/mirge3/mus_musculus/mouse --mirna-db mirbase --release sourceforge  --output' in reference
    assert '--libraries-path data/ext/mirge3/mus_musculus --organism-name mouse' in mirge3_command(output)
    assert '--index-dir data/ext/mirge3/mus_musculus/mouse/index.Libs' in rule_block(output, 'mirge3_spikein')


@pytest.mark.parametrize('update, message', [
    (lambda s: s['quant'].setdefault('mirge3', {}).update(db='mirbase22'), "quant.mirge3.db: mirbase or mirgenedb, got 'mirbase22'"),
    (lambda s: s.setdefault('db', {}).setdefault('mirge3', {}).update(release=None), 'db.mirge3.release is not set')])
def test_invalid_mirge3_selection_fails_at_load(tmp_path, update, message):
    returncode, output = dry_run(tmp_path, 'mirge3', update)
    assert returncode != 0
    assert message in output


@pytest.mark.parametrize('mirna_db, feature_info', [
    ('mirbase', 'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.tsv'),
    ('mirgenedb', 'data/ext/unitas/06-11-2026-mirgenedb/anno/mirgenedb_genes_mus_musculus.tsv')])
def test_unitas_mirna_feature_info_follows_the_dump(tmp_path, mirna_db, feature_info):
    """A MirGeneDB dump names miRNAs as MirGeneDB does, so the AnnData annotation comes from the dump's GFF."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata')
    assert returncode == 0, output
    rule = output.split('rule unitas_mirna_anndata:', 1)[1]
    assert feature_info in rule.split('output:', 1)[0], output


@pytest.mark.parametrize('mirna_db, annotation', [
    ('mirbase', 'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3'),
    ('mirgenedb', 'data/ext/unitas/06-11-2026-mirgenedb/UNITAS_refdump_mus_musculus/gcf_mirgenedb.gff3')])
def test_unitas_reference_is_checked_before_counts_meet_annotations(tmp_path, mirna_db, annotation):
    """The dump is verified once in EXT_DIR, then checked against the annotation the AnnData joins it to."""
    def select(settings):
        settings.setdefault('db', {}).setdefault('unitas', {})['mirna_db'] = mirna_db
    returncode, output = dry_run(tmp_path, 'unitas', update=select, target='unitas_mirna_anndata', printshell=True)
    assert returncode == 0, output
    release = '06-11-2026-mirgenedb' if mirna_db == 'mirgenedb' else '05-11-2026'
    assert f'--output data/ext/unitas/{release}/verify/UNITAS_refdump_mus_musculus.tsv' in output
    check = output.split('rule unitas_reference:', 1)[1].split('rule ', 1)[0]
    assert f'--mirna-db {mirna_db} --annotation {annotation} --release {release}' in check
    assert 'data/tmp/smallrna/quant/unitas/reference.tsv' in output.split('rule unitas_mirna_anndata:', 1)[1].split('output:', 1)[0]
