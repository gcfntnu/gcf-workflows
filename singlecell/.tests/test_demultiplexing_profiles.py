"""Workflow-owned conversion settings remain independent of scientific presets."""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_conversion_matrix_is_bounded_and_preserves_legacy_images():
    profiles = yaml.safe_load((ROOT / 'demultiplexing.config').read_text())
    images = yaml.safe_load((ROOT / 'docker.config').read_text())['docker']
    assert profiles['schema_version'] == 1
    assert profiles['bcl_convert_versions'] == ['4.2.4']
    assert images['bcl-convert'] == 'gcfntnu/bcl-convert:4.2.4'
    names = []
    for profile in profiles['profiles'].values():
        names.extend(profile['kits'])
        image = images[profile['legacy_tool']]
        assert ':' in image and not image.endswith(':latest')
        if profile['default_backend'] == 'bcl-convert':
            assert profile['assay'] in {'gex', 'multiome-gex'}
            assert profile['minimum_biological_cycles'][0] >= 28
        else:
            assert profile['default_backend'] == 'mkfastq' and profile['assay'] == 'legacy'
    assert len(names) == len(set(names))
    for label in ('v1', 'v2', 'ATAC', 'Visium'):
        assert any(label in kit for kit in names)
    assert images['cellranger-mkfastq'] == 'gcfntnu/cellranger:7.2.0'
    assert images['cellranger'] == 'gcfntnu/cellranger:10.0.0'
