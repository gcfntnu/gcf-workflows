"""CellBender expression presence decodes the posterior on CellBender's own input matrix (#303).

CellBender flattens (cell, feature) as m = cell * n_features + feature over its full input matrix, including
non-expression features. A matrix with a different feature axis (for example a gene-expression-only copy)
must be rejected, not silently decoded.
"""
import gzip
from types import SimpleNamespace
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'singlecell/rules/quant/scripts'))
import cellbender_feature_presence as presence  # noqa: E402
from test_velocity import TENX_KIT, protected_env, workflow_fixture  # noqa: E402

FEATURES = [('ENSG1', 'GeneA', 'Gene Expression'), ('AB1', 'AbX', 'Antibody Capture'), ('ENSG2', 'GeneB', 'Gene Expression')]
BARCODES = ['AAAC-1', 'AAAG-1']


def write_mtx(path, features, counts, gz):
    path.mkdir(parents=True)
    suffix = '.gz' if gz else ''
    opener = gzip.open if gz else open
    with opener(path / f'matrix.mtx{suffix}', 'wb') as fh:
        scipy.io.mmwrite(fh, sp.coo_matrix(counts.T))
    with opener(path / f'barcodes.tsv{suffix}', 'wt') as fh:
        fh.write(''.join(b + '\n' for b in BARCODES))
    with opener(path / (f'features.tsv{suffix}' if gz else 'genes.tsv'), 'wt') as fh:
        fh.write(''.join('\t'.join(f[:3] if gz else f[:2]) + '\n' for f in features))


def posterior_axes():
    """Entries for (cell 1, GeneB) and (cell 0, AbX), flattened over the full 3-feature input."""
    cells, feats = np.array([1, 0]), np.array([2, 1])
    return cells * len(FEATURES) + feats, cells, feats


def test_decodes_on_cellbender_input_with_all_features(tmp_path):
    counts = np.array([[1, 5, 0], [0, 2, 3]])
    write_mtx(tmp_path / 'raw', FEATURES, counts, gz=True)
    X, barcodes, genes = presence.read_matrix(tmp_path / 'raw')
    m, cells, feats = posterior_axes()
    n, g = presence.decode_and_validate(m, *X.shape, cells, feats)
    assert n.tolist() == [1, 0] and g.tolist() == [2, 1]
    assert [genes.iloc[i, 1] for i in g] == ['GeneB', 'AbX'] and X[1, 2] == 3


def test_rejects_gene_expression_only_copy(tmp_path):
    """A gex-only copy can decode some indexes without error (m=1 is (0, 1) for 2 or 3 features); identity must fail."""
    gex = [f for f in FEATURES if f[2] == 'Gene Expression']
    write_mtx(tmp_path / 'gex', gex, np.array([[1, 0], [0, 3]]), gz=False)
    X, barcodes, genes = presence.read_matrix(tmp_path / 'gex')
    n, g = presence.decode_and_validate(np.array([1]), *X.shape, np.array([0]), np.array([1]))
    assert genes.iloc[g[0], 1] == 'GeneB'  # CellBender meant AbX: decoding alone cannot tell
    with pytest.raises(ValueError, match='exact matrix'):
        presence.validate_identity(genes.iloc[:, 0].tolist(), barcodes, [f[0] for f in FEATURES], BARCODES)
    m, cells, feats = posterior_axes()
    with pytest.raises(ValueError):
        presence.decode_and_validate(m, *X.shape, cells, feats)


def test_identity_accepts_cellbender_input(tmp_path):
    write_mtx(tmp_path / 'raw', FEATURES, np.array([[1, 5, 0], [0, 2, 3]]), gz=True)
    _, barcodes, genes = presence.read_matrix(tmp_path / 'raw')
    presence.validate_identity(genes.iloc[:, 0].tolist(), barcodes, [f[0] for f in FEATURES], BARCODES[1:])
    with pytest.raises(ValueError, match='not in the counts'):
        presence.validate_identity(genes.iloc[:, 0].tolist(), barcodes, [f[0] for f in FEATURES], ['TTTT-1'])


def test_rejects_inconsistent_mtx_directory(tmp_path):
    write_mtx(tmp_path / 'raw', FEATURES[:2], np.array([[1, 5, 0], [0, 2, 3]]), gz=True)
    with pytest.raises(ValueError, match='does not match'):
        presence.read_matrix(tmp_path / 'raw')


def test_reads_cellbender_h5_layouts(tmp_path):
    """Synthetic files in the layout CellBender writes (write_matrix_to_cellranger_h5, write_posterior_coo_to_h5)."""
    with h5py.File(tmp_path / 'out.h5', 'w') as f:
        f['matrix/features/id'] = np.array([x[0] for x in FEATURES], dtype='S')
        f['matrix/barcodes'] = np.array(BARCODES[1:], dtype='S')
    with h5py.File(tmp_path / 'out_v2.h5', 'w') as f:
        f['matrix_v2/genes'] = np.array([x[0] for x in FEATURES], dtype='S')
        f['matrix_v2/barcodes'] = np.array(BARCODES, dtype='S')
    with h5py.File(tmp_path / 'posterior.h5', 'w') as f:
        f['metadata/barcode_inds'] = np.array([1, 0], dtype=np.uint32)
        f['metadata/feature_inds'] = np.array([2, 1], dtype=np.uint32)
    assert presence.read_cellbender_axes(tmp_path / 'out.h5') == ([x[0] for x in FEATURES], BARCODES[1:])
    assert presence.read_cellbender_axes(tmp_path / 'out_v2.h5')[1] == BARCODES
    cells, feats = presence.read_posterior_axes(tmp_path / 'posterior.h5')
    assert cells.tolist() == [1, 0] and feats.tolist() == [2, 1]


def test_main_end_to_end(tmp_path, monkeypatch):
    """Whole script on synthetic CellBender files; only CellBender's posterior loader is replaced."""
    counts = np.array([[1, 5, 0], [0, 2, 3]])
    write_mtx(tmp_path / 'raw', FEATURES, counts, gz=True)
    m, cells, feats = posterior_axes()
    with h5py.File(tmp_path / 'out.h5', 'w') as f:
        f['matrix/features/id'] = np.array([x[0] for x in FEATURES], dtype='S')
        f['matrix/barcodes'] = np.array(BARCODES, dtype='S')
    with h5py.File(tmp_path / 'posterior.h5', 'w') as f:
        f['metadata/barcode_inds'], f['metadata/feature_inds'] = cells, feats
    # noise count equal to the measured count, with probability 0.2 -> P(real) reported as exp(log_prob) of that entry
    coo = sp.coo_matrix((np.log([0.2, 0.4]), (m, [3, 5])), shape=(len(FEATURES) * len(BARCODES), 6))
    monkeypatch.setattr(presence, 'load_posterior', lambda path: SimpleNamespace(_noise_count_posterior_coo=coo))
    out = tmp_path / 'presence.tsv'
    monkeypatch.setattr(sys, 'argv', ['x', '--counts', str(tmp_path / 'raw'), '--cb-out-h5', str(tmp_path / 'out.h5'),
                                      '--posterior', str(tmp_path / 'posterior.h5'), '--subset', 'GeneB', '--output', str(out)])
    presence.main()
    table = out.read_text().splitlines()
    assert table[0].split('\t')[:2] == ['Barcode', 'P_gene_GeneB'] and table[1].startswith('AAAG-1\t0.2')


@pytest.mark.parametrize('methods', ['cellranger', '10x_starsolo'])
def test_presence_reads_the_cellbender_input_directory(tmp_path, methods):
    source = workflow_fixture(tmp_path, methods, kit=TENX_KIT)
    source['quant']['cellbender'] = {'enabled': True}
    source['quant']['cellbender_call'] = {'subset': 'GeneA,GeneB'}
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(source))
    target = f'data/tmp/singlecell/quant/aggregate/{methods}/cellbender/all_samples_expression_presence.tsv'
    command = [str(Path(sys.executable).with_name('snakemake')), '--dry-run', '--printshellcmds', '--cores', '1', target]
    result = subprocess.run(command, cwd=tmp_path, env=protected_env(tmp_path), capture_output=True, text=True)
    log = result.stdout + result.stderr
    assert result.returncode == 0, log
    inputs = {s: d for d, s in re.findall(r'remove-background --input (\S+) --output \S+/(\w+)\.h5 ', log)}
    counts = {s: d for d, s in re.findall(r'cellbender_feature_presence\.py --counts (\S+) --cb-out-h5 \S+/(\w+)\.h5 ', log)}
    assert set(inputs) == {'libA', 'libB'}, log
    assert counts == inputs
