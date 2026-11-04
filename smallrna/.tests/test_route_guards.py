"""Unsupported small-RNA routes stop at workflow load instead of failing mid-run."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def dry_run(tmp_path, update):
    shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    update(settings)
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                      f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '--cores', '1', 'multiqc_report'],
                            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    return result.returncode, result.stdout + result.stderr


def set_aligner(name):
    return lambda settings: settings.setdefault('align', {}).setdefault('genome', {}).update(aligner=name)


def set_dedup(value):
    return lambda settings: settings.setdefault('filter', {}).update(dedup=value)


@pytest.mark.parametrize('aligner', ['bowtie2', 'bwa'])
def test_genome_aligner_other_than_star_fails(tmp_path, aligner):
    returncode, output = dry_run(tmp_path, set_aligner(aligner))
    assert returncode != 0
    assert f"align.genome.aligner: smallrna supports only star, got '{aligner}'" in output


def test_dedup_true_fails(tmp_path):
    returncode, output = dry_run(tmp_path, set_dedup(True))
    assert returncode != 0
    assert 'filter.dedup: read deduplication is not supported in smallrna' in output


@pytest.mark.parametrize('update', [set_aligner('star'), set_dedup(False)], ids=['star', 'dedup-false'])
def test_supported_settings_plan_delivery(tmp_path, update):
    returncode, output = dry_run(tmp_path, update)
    assert returncode == 0, output
    assert 'rule star_align:' in output


def set_spikein(ref):
    return lambda settings: settings.setdefault('filter', {}).setdefault('spikein', {}).update(quantifier='unitas', ref=ref)


@pytest.mark.parametrize('ref, fastas', [('ERCC92,smallrna_calibrators', ['ERCC92.fa', 'smallrna_calibrators.fa']),
                                         ('ERCC92', ['ERCC92.fa']), ('smallrna_calibrators', ['smallrna_calibrators.fa'])])
def test_spikein_sets_select_their_fasta_files(tmp_path, ref, fastas):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode == 0, output
    inputs = output.split('rule spikein_unitas:')[1].split('output:')[0]
    assert sorted(path.rsplit('/', 1)[1] for path in inputs.split('input:')[1].strip().split(', ')) == sorted(fastas)


@pytest.mark.parametrize('ref', ['ercc', 'ERCC92,mirbase'])
def test_unknown_spikein_set_fails(tmp_path, ref):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode != 0
    assert 'filter.spikein.ref: unsupported value' in output


@pytest.mark.parametrize('ref', ['', None])
def test_spikein_quantification_without_sets_fails(tmp_path, ref):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode != 0
    assert 'filter.spikein.ref lists no spike-in sets' in output


def test_default_has_no_spikeins(tmp_path):
    returncode, output = dry_run(tmp_path, lambda settings: None)
    assert returncode == 0, output
    assert 'rule spikein_unitas:' not in output and 'spikein/fasta' not in output
