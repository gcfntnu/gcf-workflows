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


def test_mirge3_spikein_flag_when_sets_are_listed(tmp_path):
    update = lambda s: s.setdefault('filter', {}).setdefault('spikein', {}).update(ref='ERCC92')
    returncode, output = dry_run(tmp_path, 'mirge3', update, printshell=True)
    assert returncode == 0, output
    assert '--spikeIn' in mirge3_command(output)


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
