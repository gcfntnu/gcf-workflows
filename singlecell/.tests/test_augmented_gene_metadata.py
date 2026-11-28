"""Join an augmented reference's gene table to a synthetic Cell Ranger count matrix."""
import gzip
from pathlib import Path
import subprocess
import sys

import anndata
import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[2]


def write_matrix(tmp_path, host_id, reporter_id):
    matrix = tmp_path / 'libA/outs/filtered_feature_bc_matrix'
    matrix.mkdir(parents=True)
    reporter_symbol = reporter_id or 'ReporterSymbol'
    for name, contents in {
        'matrix.mtx.gz': '%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 4\n2 1 3\n2 2 9\n',
        'features.tsv.gz': f'{host_id}\tHostGene\tGene Expression\n{reporter_id}\t{reporter_symbol}\tGene Expression\n',
        'barcodes.tsv.gz': 'AAAA-1\nCCCC-1\n',
    }.items():
        with gzip.open(matrix / name, 'wt') as output:
            output.write(contents)
    return matrix / 'matrix.mtx.gz'


def convert_matrix(tmp_path, matrix, genes, converted):
    command = [sys.executable, ROOT / 'singlecell/rules/quant/scripts/convert_scanpy.py', matrix, '-f', 'cellranger',
               '--barcode-rename', 'skip', '--canonical-filtered', '--feature-info', genes, '-F', 'anndata', '-o', converted]
    return subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)


@pytest.mark.parametrize('host_id,reporter_id', [('HOST', 'REPORTER'), ('HOST', 'NA'), ('1', '001')])
def test_exon_only_reporter_counts_and_annotation(tmp_path, host_id, reporter_id):
    gtf = tmp_path / 'genes.gtf'
    gtf.write_text(
        f'chrHost\tfixture\tgene\t1\t9\t.\t+\t.\tgene_id "{host_id}"; gene_name "HostGene"; gene_biotype "protein_coding";\n'
        f'chrHost\tfixture\texon\t3\t6\t.\t+\t.\tgene_id "{host_id}"; transcript_id "HOST-T";\n'
        f'chrReporter\tfixture\texon\t2\t5\t.\t-\t.\tgene_id "{reporter_id}"; transcript_id "REPORTER-T";\n'
        f'chrReporter\tfixture\texon\t10\t11\t.\t-\t.\tgene_id "{reporter_id}"; transcript_id "REPORTER-T";\n'
    )
    exons = tmp_path / 'exons.fa'
    exons.write_text('>chrHost:2-6\nAGCT\n>chrReporter:1-5\nGGCC\n>chrReporter:9-11\nAT\n')
    genes = tmp_path / 'genes.tsv'
    subprocess.run([sys.executable, ROOT / 'gcfdb/scripts/gtf2tsv.py', '--gtf', gtf, '--fasta', exons,
                    '--output', genes, '--feature', 'gene', '--add-gc-content'], cwd=tmp_path, check=True)

    matrix = write_matrix(tmp_path, host_id, reporter_id)
    converted = tmp_path / 'filtered.h5ad'
    result = convert_matrix(tmp_path, matrix, genes, converted)
    assert result.returncode == 0, result.stderr
    data = anndata.read_h5ad(converted)
    assert list(data.obs_names) == ['AAAA-1', 'CCCC-1']
    assert list(data.var_names) == [host_id, reporter_id]
    np.testing.assert_array_equal(data.X.toarray(), [[4, 3], [0, 9]])
    assert data.var.loc[reporter_id, 'gene_symbols'] == reporter_id
    assert data.var.loc[reporter_id, 'chrom'] == 'chrReporter'
    assert data.var.loc[reporter_id, 'gc_content'] == pytest.approx(100 * 4 / 6)
    assert data.var.loc[host_id, 'gene_biotype'] == 'protein_coding'
    assert pd.isna(data.var.loc[reporter_id, 'gene_biotype'])

    table = pd.read_csv(genes, sep='\t', dtype={'gene_id': str}, keep_default_na=False).set_index('gene_id', verify_integrity=True)
    assert set(table.index) == {host_id, reporter_id}
    fields = ['gene_name', 'chrom', 'strand', 'start', 'end']
    assert table.loc[host_id, fields].tolist() == ['HostGene', 'chrHost', '+', 1, 9]
    assert table.loc[reporter_id, fields].tolist() == [reporter_id, 'chrReporter', '-', 2, 11]


@pytest.mark.parametrize('missing_from', ['feature_info', 'matrix'])
def test_missing_gene_id_is_rejected(tmp_path, missing_from):
    genes = tmp_path / 'genes.tsv'
    reporter = '' if missing_from == 'feature_info' else 'REPORTER'
    genes.write_text(f'gene_id\tchrom\nHOST\tchrHost\n{reporter}\tchrReporter\n')
    matrix = write_matrix(tmp_path, 'HOST', '' if missing_from == 'matrix' else 'REPORTER')
    converted = tmp_path / 'filtered.h5ad'
    result = convert_matrix(tmp_path, matrix, genes, converted)
    assert result.returncode != 0 and not converted.exists()
