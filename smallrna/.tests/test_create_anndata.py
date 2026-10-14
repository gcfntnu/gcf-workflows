"""Run with pytest in an environment containing anndata, numpy and pandas.

CLI tests exercise real H5AD serialization. The VSN boundary test substitutes R;
actual DESeq2 normalization is covered by the deployed small-RNA integration run.
"""

import contextlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import types

import anndata as ad
import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "rules/quant/scripts/create_anndata.py"


@pytest.fixture
def inputs(tmp_path):
    counts = pd.DataFrame({"sample_b": [12, 24, 36], "sample_a": [11, 23, 35]}, index=["mir_c", "mir_a", "mir_b"])
    samples = pd.DataFrame({"group": ["A", "B", "unused"]}, index=["sample_a", "sample_b", "sample_extra"])
    features = pd.DataFrame({"annotation": ["alpha", "beta", "gamma", "unused"]}, index=["mir_a", "mir_b", "mir_c", "mir_extra"])
    for name, frame in (("counts", counts), ("samples", samples), ("features", features)):
        frame.to_csv(tmp_path / f"{name}.tsv", sep="\t")
    return tmp_path, counts, samples, features


def run_cli(root, *extra):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(root / "counts.tsv"), "--output", str(root / "result.h5ad"), *extra],
        capture_output=True, text=True, check=False,
    )


def test_cli_preserves_float32_counts_and_reorders_metadata(inputs):
    root, counts, _, _ = inputs
    result = run_cli(root, "--sample-info", str(root / "samples.tsv"), "--feature-info", str(root / "features.tsv"))
    assert result.returncode == 0, result.stderr
    actual = ad.read_h5ad(root / "result.h5ad")
    assert actual.X.dtype == np.dtype("float32")
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    assert actual.obs_names.tolist() == ["sample_b", "sample_a"]
    assert actual.var_names.tolist() == ["mir_c", "mir_a", "mir_b"]
    assert actual.obs["group"].tolist() == ["B", "A"]
    assert actual.var["annotation"].tolist() == ["gamma", "alpha", "beta"]


def test_cli_without_optional_metadata(inputs):
    root, counts, _, _ = inputs
    result = run_cli(root)
    assert result.returncode == 0, result.stderr
    actual = ad.read_h5ad(root / "result.h5ad")
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    assert actual.X.dtype == np.dtype("float32")
    assert actual.obs_names.tolist() == counts.columns.tolist()
    assert actual.var_names.tolist() == counts.index.tolist()


def test_cli_missing_feature_annotation_warns_without_dropping_counts(inputs):
    root, counts, _, features = inputs
    features.drop(index="mir_c").to_csv(root / "features.tsv", sep="\t")
    result = run_cli(root, "--feature-info", str(root / "features.tsv"))
    assert result.returncode == 0, result.stderr
    assert "missing annotations in feature info" in result.stderr
    actual = ad.read_h5ad(root / "result.h5ad")
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    assert actual.var_names.tolist() == counts.index.tolist()
    assert pd.isna(actual.var.loc["mir_c", "annotation"])
    assert actual.var.loc["mir_a", "annotation"] == "alpha"


def test_cli_missing_sample_metadata_fails(inputs):
    root, _, samples, _ = inputs
    samples.drop(index="sample_b").to_csv(root / "samples.tsv", sep="\t")
    result = run_cli(root, "--sample-info", str(root / "samples.tsv"))
    assert result.returncode != 0
    assert "missing samples in sample info" in result.stderr
    assert not (root / "result.h5ad").exists()


@pytest.mark.parametrize("sample_count,method,r_function", [(2, "deseq2_rlog", "rlogTransformation"), (51, "deseq2_vst", "varianceStabilizingTransformation")])
def test_vsn_scoped_conversion_preserves_axes_and_h5ad_layer(monkeypatch, tmp_path, sample_count, method, r_function):
    """Check the Python/R boundary, not the numerical behavior of DESeq2."""
    spec = importlib.util.spec_from_file_location("smallrna_create_anndata", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    active = False
    environment = {}
    expressions = []
    numpy_converter, combined = object(), object()

    class Converter:
        def __add__(self, other):
            assert other is numpy_converter
            return combined

    @contextlib.contextmanager
    def localconverter(converter):
        nonlocal active
        assert converter is combined
        active = True
        try:
            yield
        finally:
            active = False

    class GlobalEnv(dict):
        def __setitem__(self, name, value):
            assert active, "NumPy input conversion must be scoped"
            super().__setitem__(name, value)

    class R:
        def __call__(self, expression):
            assert active
            expressions.append(expression)
            # Distinct values expose accidental transposition or count replacement.
            environment["vsd"] = np.asarray(globalenv["adata"], dtype=float) + 1000

        def __getitem__(self, name):
            assert active, "R result conversion must be scoped"
            return environment[name]

    globalenv = GlobalEnv()
    imported_packages = []
    robjects = types.ModuleType("rpy2.robjects")
    robjects.r = R()
    robjects.globalenv = globalenv
    robjects.default_converter = Converter()
    robjects.numpy2ri = types.SimpleNamespace(converter=numpy_converter)
    packages = types.ModuleType("rpy2.robjects.packages")
    packages.importr = imported_packages.append
    conversion = types.ModuleType("rpy2.robjects.conversion")
    conversion.localconverter = localconverter
    for name, value in (("rpy2", types.ModuleType("rpy2")), ("rpy2.robjects", robjects), ("rpy2.robjects.packages", packages), ("rpy2.robjects.conversion", conversion)):
        monkeypatch.setitem(sys.modules, name, value)

    counts = np.arange(sample_count * 3, dtype=np.float32).reshape(sample_count, 3)
    actual = ad.AnnData(counts.copy(), obs=pd.DataFrame(index=[f"sample_{i}" for i in range(sample_count)]), var=pd.DataFrame(index=["mir_c", "mir_a", "mir_b"]))
    transformed = module.vsn(actual)
    assert not active
    assert imported_packages == ["DESeq2"]
    np.testing.assert_array_equal(globalenv["adata"], counts.T)
    assert expressions == [f'vsd <- {r_function}(round(adata+1), fitType="mean")']
    transformed.write_h5ad(tmp_path / "vsn.h5ad")
    restored = ad.read_h5ad(tmp_path / "vsn.h5ad")
    np.testing.assert_array_equal(restored.X, counts)
    np.testing.assert_array_equal(restored.layers["vsd"], counts + 1000)
    assert restored.obs_names.tolist() == actual.obs_names.tolist()
    assert restored.var_names.tolist() == actual.var_names.tolist()
    assert restored.uns["varians_transform"] == method


def test_cli_precomputed_vsn_roundtrip(inputs):
    root, counts, _, _ = inputs
    transformed = np.log2(counts + 1)
    transformed.to_csv(root / "vsn.tsv", sep="\t")
    result = run_cli(root, "--add-vsn", "--vsn-input", str(root / "vsn.tsv"))
    assert result.returncode == 0, result.stderr
    actual = ad.read_h5ad(root / "result.h5ad")
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    np.testing.assert_allclose(actual.layers["vsd"], transformed.to_numpy().T)
    assert actual.X.dtype == np.dtype("float32")
    assert actual.uns["varians_transform"] == "deseq2_rlog"
    assert actual.obs_names.tolist() == counts.columns.tolist()
    assert actual.var_names.tolist() == counts.index.tolist()


@pytest.mark.parametrize("invalid", ["reordered_samples", "reordered_features", "missing_sample", "missing_feature", "nonfinite"])
def test_cli_precomputed_vsn_rejects_invalid_input(inputs, invalid):
    root, counts, _, _ = inputs
    transformed = counts.astype(float)
    if invalid == "reordered_samples":
        transformed = transformed.loc[:, transformed.columns[::-1]]
    elif invalid == "reordered_features":
        transformed = transformed.iloc[::-1]
    elif invalid == "missing_sample":
        transformed = transformed.iloc[:, :1]
    elif invalid == "missing_feature":
        transformed = transformed.iloc[:2]
    else:
        transformed.iloc[0, 0] = np.inf
    transformed.to_csv(root / "vsn.tsv", sep="\t")
    result = run_cli(root, "--add-vsn", "--vsn-input", str(root / "vsn.tsv"))
    assert result.returncode != 0
    assert not (root / "result.h5ad").exists()


def test_cli_precomputed_vsn_requires_add_vsn(inputs):
    root, counts, _, _ = inputs
    counts.to_csv(root / "vsn.tsv", sep="\t")
    result = run_cli(root, "--vsn-input", str(root / "vsn.tsv"))
    assert result.returncode != 0
    assert "--add-vsn" in result.stderr
    assert not (root / "result.h5ad").exists()


def test_cli_precomputed_vsn_large_cohort_records_vst(tmp_path):
    counts = pd.DataFrame(
        np.arange(3 * 51).reshape(3, 51),
        index=["mir_c", "mir_a", "mir_b"],
        columns=[f"sample_{i}" for i in range(51)],
    )
    transformed = np.log2(counts + 1)
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    transformed.to_csv(tmp_path / "vsn.tsv", sep="\t")
    result = run_cli(tmp_path, "--add-vsn", "--vsn-input", str(tmp_path / "vsn.tsv"))
    assert result.returncode == 0, result.stderr
    actual = ad.read_h5ad(tmp_path / "result.h5ad")
    assert actual.uns["varians_transform"] == "deseq2_vst"
    assert actual.X.dtype == np.dtype("float32")
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    np.testing.assert_allclose(actual.layers["vsd"], transformed.to_numpy().T)
    assert actual.obs_names.tolist() == counts.columns.tolist()
    assert actual.var_names.tolist() == counts.index.tolist()


def test_cli_metadata_column_names_escape_hdf5_separators_without_collisions(inputs):
    root, counts, samples, features = inputs
    samples["260/280"] = [1.8, 1.9, 2.0]
    samples["260/230"] = [2.1, 2.2, 2.3]
    samples["260%2F280"] = [81, 82, 83]
    features["reference/accession"] = ["ref_a", "ref_b", "ref_c", "ref_unused"]
    features["reference%2Faccession"] = [11, 12, 13, 14]
    samples.to_csv(root / "samples.tsv", sep="\t")
    features.to_csv(root / "features.tsv", sep="\t")
    result = run_cli(root, "--sample-info", str(root / "samples.tsv"), "--feature-info", str(root / "features.tsv"))
    assert result.returncode == 0, result.stderr
    actual = ad.read_h5ad(root / "result.h5ad")
    assert actual.obs.columns.tolist() == ["group", "260%2F280", "260%2F230", "260%252F280"]
    assert actual.var.columns.tolist() == ["annotation", "reference%2Faccession", "reference%252Faccession"]
    np.testing.assert_allclose(actual.obs["260%2F280"], [1.9, 1.8])
    np.testing.assert_allclose(actual.obs["260%2F230"], [2.2, 2.1])
    assert actual.obs["260%252F280"].tolist() == [82, 81]
    assert actual.obs["group"].tolist() == ["B", "A"]
    assert actual.var["reference%2Faccession"].tolist() == ["ref_c", "ref_a", "ref_b"]
    assert actual.var["reference%252Faccession"].tolist() == [13, 11, 12]
    assert actual.var["annotation"].tolist() == ["gamma", "alpha", "beta"]
    assert actual.uns["metadata_column_names"] == {
        "obs": {"260%2F280": "260/280", "260%2F230": "260/230", "260%252F280": "260%2F280"},
        "var": {"reference%2Faccession": "reference/accession", "reference%252Faccession": "reference%2Faccession"},
    }
    np.testing.assert_array_equal(actual.X, counts.to_numpy().T)
    assert actual.X.dtype == np.dtype("float32")
    assert actual.obs_names.tolist() == counts.columns.tolist()
    assert actual.var_names.tolist() == counts.index.tolist()
