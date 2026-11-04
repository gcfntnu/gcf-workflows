"""Dry-run contract for the small-RNA `quant.method` selection of `quant_all`."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def dry_run(tmp_path, method):
    shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    settings.setdefault('quant', {})['method'] = method
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                      f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '--cores', '1', 'quant_all'],
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
