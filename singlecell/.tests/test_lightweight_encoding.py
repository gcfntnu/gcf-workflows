"""Lightweight AnnData stays readable by older anndata readers (#285).

solo-sc (anndata 0.10), socube and cellbender (anndata 0.8) read `anndata.light.h5ad`; they cannot read
`nullable-string-array`, which anndata >= 0.11 writes for pandas string indexes (the pandas 3 default).
"""
import sys
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'rules' / 'quant' / 'scripts'))
import convert_scanpy  # noqa: E402

KEYS = ['obs/_index', 'var/gene_id', 'var/gene_symbol']


@pytest.mark.parametrize('dtype', [None, 'string', 'infer_string'])
def test_lightweight_writes_string_array(tmp_path, dtype):
    with pd.option_context('future.infer_string', dtype == 'infer_string'):
        dt = 'string' if dtype == 'string' else None
        obs = pd.DataFrame(index=pd.Index([f'BC{i}-1' for i in range(5)], dtype=dt))
        var = pd.DataFrame({'gene_symbol': pd.array([f'G{i}' for i in range(4)], dtype=dt)},
                           index=pd.Index([f'ENSG{i}' for i in range(4)], dtype=dt))
        lw = convert_scanpy.to_anndata_lightweight(anndata.AnnData(X=sp.csr_matrix(np.ones((5, 4), np.float32)), obs=obs, var=var))
        lw.write(tmp_path / 'anndata.light.h5ad')
    with h5py.File(tmp_path / 'anndata.light.h5ad') as f:
        assert {k: f[k].attrs['encoding-type'] for k in KEYS} == dict.fromkeys(KEYS, 'string-array')
    back = anndata.read_h5ad(tmp_path / 'anndata.light.h5ad')
    assert list(back.obs_names) == [f'BC{i}-1' for i in range(5)] and list(back.var['gene_symbol']) == [f'G{i}' for i in range(4)]
