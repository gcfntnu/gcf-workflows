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
