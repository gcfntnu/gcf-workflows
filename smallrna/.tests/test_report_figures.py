"""Small-RNA report figures (#354): the PCA on log2(CPM + 1) and the top-miRNA plot in counts per million.

The figure tests need the workflow's scanpy and anndata (the `scanpy` image) and skip without them; the dry runs
need Snakemake and skip without it.
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from test_quant_methods import ROOT, configure, dry_run, record_as_done

SCRIPTS = ROOT / 'smallrna/rules/bfq/scripts'


def mirna_adata(counts, samples):
    ad = pytest.importorskip('anndata')
    pytest.importorskip('scanpy')
    return ad.AnnData(np.asarray(counts, dtype=float), obs={'sample': samples},
                      var={'mirna': ['miR-{}'.format(i) for i in range(len(counts[0]))]}).copy()


def write_adata(tmp_path, counts):
    adata = mirna_adata(counts, ['S{}'.format(i + 1) for i in range(len(counts))])
    adata.obs_names = list(adata.obs['sample'])
    adata.var_names = list(adata.var['mirna'])
    path = tmp_path / 'adata.h5ad'
    adata.write_h5ad(path)
    return path


def run_pca(tmp_path, counts):
    path = write_adata(tmp_path, counts)
    output = tmp_path / 'pca_mqc.yaml'
    result = subprocess.run([sys.executable, str(SCRIPTS / 'plotpca.py'), '--recipe', 'smallrna', '--input', str(path),
                             '--output', str(output)], capture_output=True, text=True, timeout=300, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    return yaml.safe_load(output.read_text())


def test_pca_uses_log_cpm_so_depth_alone_does_not_move_a_sample(tmp_path):
    rng = np.random.default_rng(1)
    base = rng.negative_binomial(5, 0.01, size=(4, 60))
    counts = np.vstack([base, base[0] * 10])  # S5: S1 sequenced ten times deeper
    counts[1, 0] = 0  # miR-0 detected in four samples: kept; miR-59 only in S1 and S5: below three samples, left out
    counts[1:4, 59] = 0
    section = run_pca(tmp_path, counts)
    points = section['data'][0]
    assert points['S1']['x'] == pytest.approx(points['S5']['x'], abs=1e-6)
    assert points['S1']['y'] == pytest.approx(points['S5']['y'], abs=1e-6)
    assert section['pconfig']['xlab'].startswith('PC1 (') and section['pconfig']['xlab'].endswith('%)')
    assert section['pconfig']['ylab'].startswith('PC2 (')
    assert 'log2(CPM + 1)' in section['description'] and 'at least 3 samples' in section['description']


def test_pca_of_three_samples(tmp_path):
    counts = np.random.default_rng(2).negative_binomial(5, 0.01, size=(3, 40))
    section = run_pca(tmp_path, counts)
    assert sorted(section['data'][0]) == ['S1', 'S2', 'S3']


def test_pca_components_come_from_the_filtered_matrix(tmp_path):
    """Three of six miRNAs are detected in at least three samples: two components from the filtered matrix."""
    counts = [[100, 20, 10, 5, 0, 0], [20, 100, 10, 0, 5, 0], [40, 20, 100, 0, 0, 5], [60, 40, 20, 0, 0, 0]]
    section = run_pca(tmp_path, counts)
    assert section['plot_type'] == 'scatter' and sorted(section['data'][0]) == ['S1', 'S2', 'S3', 'S4']
    assert section['pconfig']['xlab'].startswith('PC1 (')


@pytest.mark.parametrize('counts, reason', [
    ([[5, 0, 3, 8], [0, 5, 3, 8], [5, 5, 0, 0]], 'fewer than three miRNAs are detected in at least three samples'),
    ([[5, 0, 3], [0, 5, 3]], 'the project has fewer than three samples')])
def test_pca_is_a_note_when_two_components_are_impossible(tmp_path, counts, reason):
    section = run_pca(tmp_path, counts)
    assert section['plot_type'] == 'html' and section['data'] == '<p>PCA not shown: {}.</p>'.format(reason)
    assert section['id'] == 'pca' and section['parent_id'] == 'exploratory'


def load_highly_expressed():
    pytest.importorskip('scanpy')
    spec = importlib.util.spec_from_file_location('plot_highly_expressed', SCRIPTS / 'plot_highly_expressed.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_top_mirnas_are_counts_per_million_mirna_reads():
    module = load_highly_expressed()
    counts = np.array([[700, 200, 100] + [0] * 7, [10, 80, 10] + [0] * 7], dtype=float)
    adata = mirna_adata(counts, ['S1', 'S2'])
    adata.obs_names, adata.var_names = ['S1', 'S2'], list(adata.var['mirna'])
    section = module.highly_expressed_yaml(adata)
    cpm, reads = section['data']
    assert cpm['S1']['miR-0'] == pytest.approx(7e5) and cpm['S2']['miR-1'] == pytest.approx(8e5)
    assert all(sum(values.values()) == pytest.approx(1e6) for values in cpm.values())  # all ten miRNAs shown
    assert reads['S1']['miR-0'] == 700
    assert list(section['categories'])[:3] == ['miR-1', 'miR-0', 'miR-2']  # ranked by mean CPM: 5e5, 4e5, 1e5
    assert section['section_name'] == section['pconfig']['title'] == 'Highly expressed miRNAs'
    assert [label['name'] for label in section['pconfig']['data_labels']] == ['CPM', 'counts']


def keep_samples(names):
    def update(settings):
        settings['samples'] = {name: settings['samples'][name] for name in names}
    return update


def cut_samples(tmp_path, names):
    configure(tmp_path, 'unitas', keep_samples(names))
    table = (tmp_path / 'pep/sample_table.csv').read_text().splitlines()
    (tmp_path / 'pep/sample_table.csv').write_text('\n'.join([table[0]] + [r for r in table[1:] if r.split(',')[0] in names]) + '\n')


@pytest.mark.parametrize('names', [['SR01', 'SR02'], ['SR01', 'SR02', 'SR03']])
def test_pca_section_is_always_written(tmp_path, names):
    pytest.importorskip('snakemake')  # dry runs; the figure tests also run in a scanpy-only environment
    configure(tmp_path, 'unitas')  # copy the fixture
    cut_samples(tmp_path, names)
    returncode, output = dry_run(tmp_path, 'unitas', keep_samples(names), target='multiqc_report')
    assert returncode == 0, output
    assert 'rule bfq_level2_pca:' in output and 'rule bfq_level2_mirna_high:' in output


def test_earlier_pca_is_rewritten_when_a_project_drops_below_three_samples(tmp_path):
    """A PCA made with four samples is not left in figs/ for the report when two samples remain."""
    pytest.importorskip('snakemake')
    configure(tmp_path, 'unitas')
    record_as_done(tmp_path, 'data/tmp/smallrna/quant/unitas/adata.h5ad', 'data/tmp/smallrna/bfq/figs/pca_mqc.yaml')
    cut_samples(tmp_path, ['SR01', 'SR02'])
    returncode, output = dry_run(tmp_path, 'unitas', keep_samples(['SR01', 'SR02']), target='multiqc_report')
    assert returncode == 0, output
    pca = output.split('rule bfq_level2_pca:', 1)[1].split('\nrule ', 1)[0]
    assert 'output: data/tmp/smallrna/bfq/figs/pca_mqc.yaml' in pca
