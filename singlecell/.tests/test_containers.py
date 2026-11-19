"""Resolved rule containers for 10x and Parse configurations (#285); the workflow is parsed, no jobs run."""
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from test_velocity import PARSE_KIT, TENX_KIT, protected_env, workflow_fixture

ROOT = Path(__file__).resolve().parents[2]
DOCKER = yaml.safe_load((ROOT / 'docker.config').read_text())['docker']
DOUBLET_RULES = {'dbl_skip_doubletdetection', 'dbl_doubletdetection', 'dbl_scdblfinder', 'dbl_scds', 'dbl_scrublet',
                 'dbl_doublet_rank_aggr'}
SCANPY_RULES = {'tmp_lightweight_raw', 'tmp_lightweight_filtered', 'scanpy_aggr_filtered'}
INSPECT = '''
import json
from pathlib import Path
Path('rules.json').write_text(json.dumps({r.name: {'image': r.container_img, 'shell': r.is_shell} for r in workflow.rules}))
'''


def resolved_containers(tmp_path, kit, methods):
    workflow_fixture(tmp_path, methods, kit=kit)
    (tmp_path / 'Snakefile').write_text((tmp_path / 'Snakefile').read_text() + INSPECT)
    command = [str(Path(sys.executable).with_name('snakemake')), '--list-rules']
    result = subprocess.run(command, cwd=tmp_path, env=protected_env(tmp_path), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads((tmp_path / 'rules.json').read_text())


@pytest.mark.parametrize('kit, methods', [(TENX_KIT, 'cellranger'), (TENX_KIT, '10x_starsolo'),
                                          (PARSE_KIT, 'splitpipe'), (PARSE_KIT, 'parsebio_starsolo')])
def test_doublet_image_is_set_per_rule_only(tmp_path, kit, methods):
    rules = resolved_containers(tmp_path, kit, methods)
    doublet = 'docker://' + DOCKER['doublet-detection']
    assert {name for name, rule in rules.items() if rule['image'] == doublet} == DOUBLET_RULES
    assert all(rules[name]['shell'] for name in DOUBLET_RULES)
    for name in SCANPY_RULES:
        assert rules[name]['image'] == 'docker://' + DOCKER['scanpy'], name
    if 'bfq_parse_metadata' in rules:
        assert rules['bfq_parse_metadata']['image'] is None


@pytest.mark.parametrize('method', ['souporcell_ref', 'souporcell_noref'])
def test_multiplex_images_and_souporcell_dag(tmp_path, method):
    config = workflow_fixture(tmp_path, 'cellranger', kit=TENX_KIT)
    config['organism'] = 'homo_sapiens'
    config['quant']['demultiplex'] = {'method': method, 'n_individuals': 2, 'donor_dir': 'donors'}
    (tmp_path / 'donors/libA').mkdir(parents=True)
    for name in ['cellSNP.cells.vcf', 'cellSNP.cells.chr.vcf']:
        (tmp_path / 'donors/libA' / name).touch()
    if method == 'souporcell_noref':
        config['docker'] = {'souporcell': 'example/souporcell:override', 'bcftools': 'example/bcftools:override',
                            'popscle': 'example/popscle:override'}
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(config))
    with (tmp_path / 'Snakefile').open('a') as snakefile:
        snakefile.write(INSPECT)
    target = f'data/tmp/singlecell/quant/cellranger/libA/demultiplexing/{method}/clusters.tsv'
    command = [str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '1', target]
    result = subprocess.run(command, cwd=tmp_path, env=protected_env(tmp_path), capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f'rule {method}:' in result.stdout
    resolved = json.loads((tmp_path / 'rules.json').read_text())
    images = dict(DOCKER, **config.get('docker', {}))
    assert DOCKER['souporcell'] == 'gcfntnu/souporcell:2.5'
    assert '@sha256:' in DOCKER['popscle']
    for rule, key in [('souporcell_ref', 'souporcell'), ('souporcell_noref', 'souporcell'),
                      ('vireo_donor_subset', 'bcftools'), ('freemuxlet_noref', 'popscle'), ('freemuxlet_pileup', 'popscle')]:
        assert resolved[rule]['image'] == 'docker://' + images[key]


def test_every_multiplex_container_uses_the_configured_docker_image():
    import re
    source = (ROOT / 'singlecell/rules/quant/multiplex.smk').read_text()
    directives = re.findall(r'^    container:\s*([^\n]+)', source, re.MULTILINE)
    assert directives
    for expression in directives:
        match = re.fullmatch(r"""['"]docker://['"]\s*\+\s*config\[['"]docker['"]\]\[['"]([\w-]+)['"]\]""", expression)
        assert match and match.group(1) in DOCKER, expression
