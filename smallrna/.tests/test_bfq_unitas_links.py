"""bfq/logs/unitas/<sample> links are relative, so delivery survives a moved work directory (#263)."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_unitas_delivery_links_are_relative(tmp_path):
    work = tmp_path / 'project'
    shutil.copytree(ROOT / 'smallrna/.tests', work, ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    (work / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                  f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    samples = list(yaml.safe_load((work / 'config.yaml').read_text())['samples'])
    unitas = work / 'data/tmp/smallrna/quant/unitas/unitas'
    for sample in samples:
        (unitas / sample).mkdir(parents=True)
        (unitas / sample / 'results.html').write_text(sample + '\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    # Only bfq_level2_unitas is planned: its results.html inputs exist and need nothing upstream.
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '--cores', '1', 'bfq_level2_unitas'],
                            cwd=work, env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr

    moved = tmp_path / 'moved'
    work.rename(moved)
    for sample in samples:
        link = moved / 'data/tmp/smallrna/bfq/logs/unitas' / sample
        assert link.is_symlink() and not os.path.isabs(os.readlink(link))
        assert (link / 'results.html').read_text() == sample + '\n'
