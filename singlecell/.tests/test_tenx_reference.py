"""Reference archive packaging checks use tiny text files, with no workflow or network."""
import gzip
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'gcfdb/scripts/copy_10x_gtf.py'
spec = importlib.util.spec_from_file_location('copy_10x_gtf', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize('compressed', [False, True])
def test_reference_annotation_cli_preserves_plain_and_compressed_bytes(tmp_path, compressed):
    genes = tmp_path / 'genes'
    genes.mkdir()
    expected = b'# reference annotation\nchr1\tsynthetic\tgene\t1\t9\t.\t+\t.\tgene_id "example";\n'
    source = genes / ('genes.gtf.gz' if compressed else 'genes.gtf')
    source.write_bytes(gzip.compress(expected) if compressed else expected)
    destination = tmp_path / 'anno/genes.gtf'
    subprocess.run([sys.executable, str(SCRIPT), str(genes), str(destination)], check=True)
    assert destination.read_bytes() == expected


@pytest.mark.parametrize('ambiguous', [False, True])
def test_missing_or_ambiguous_annotation_preserves_existing_destination(tmp_path, ambiguous):
    genes = tmp_path / 'genes'
    genes.mkdir()
    if ambiguous:
        (genes / 'genes.gtf').write_text('plain')
        (genes / 'genes.gtf.gz').write_bytes(gzip.compress(b'compressed'))
    destination = tmp_path / 'retained.gtf'
    destination.write_bytes(b'previous annotation')
    with pytest.raises(ValueError, match='Expected exactly one'):
        module.copy_gtf(genes, destination)
    assert destination.read_bytes() == b'previous annotation'


def test_corrupt_gzip_preserves_existing_annotation(tmp_path):
    genes = tmp_path / 'genes'
    genes.mkdir()
    (genes / 'genes.gtf.gz').write_bytes(gzip.compress(b'new annotation')[:-4])
    destination = tmp_path / 'retained.gtf'
    destination.write_bytes(b'previous annotation')
    with pytest.raises(EOFError):
        module.copy_gtf(genes, destination)
    assert destination.read_bytes() == b'previous annotation'
    assert {path.name for path in tmp_path.iterdir()} == {'genes', 'retained.gtf'}
