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


def test_exon_only_reporter_counts_and_annotation(tmp_path):
    gtf = tmp_path / 'genes.gtf'
    gtf.write_text(
        'chrHost\tfixture\tgene\t1\t9\t.\t+\t.\tgene_id "HOST"; gene_name "HostGene"; gene_biotype "protein_coding";\n'
        'chrHost\tfixture\texon\t3\t6\t.\t+\t.\tgene_id "HOST"; transcript_id "HOST-T";\n'
        'chrReporter\tfixture\texon\t2\t5\t.\t-\t.\tgene_id "REPORTER"; transcript_id "REPORTER-T";\n'
        'chrReporter\tfixture\texon\t10\t11\t.\t-\t.\tgene_id "REPORTER"; transcript_id "REPORTER-T";\n'
    )
    exons = tmp_path / 'exons.fa'
    exons.write_text('>chrHost:2-6\nAGCT\n>chrReporter:1-5\nGGCC\n>chrReporter:9-11\nAT\n')
    genes = tmp_path / 'genes.tsv'
    subprocess.run([sys.executable, ROOT / 'gcfdb/scripts/gtf2tsv.py', '--gtf', gtf, '--fasta', exons,
                    '--output', genes, '--feature', 'gene', '--add-gc-content'], cwd=tmp_path, check=True)

    matrix = tmp_path / 'libA/outs/filtered_feature_bc_matrix'
    matrix.mkdir(parents=True)
    for name, contents in {
        'matrix.mtx.gz': '%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 4\n2 1 3\n2 2 9\n',
        'features.tsv.gz': 'HOST\tHostGene\tGene Expression\nREPORTER\tREPORTER\tGene Expression\n',
        'barcodes.tsv.gz': 'AAAA-1\nCCCC-1\n',
    }.items():
        with gzip.open(matrix / name, 'wt') as output:
            output.write(contents)
    converted = tmp_path / 'filtered.h5ad'
    subprocess.run([sys.executable, ROOT / 'singlecell/rules/quant/scripts/convert_scanpy.py', matrix / 'matrix.mtx.gz',
                    '-f', 'cellranger', '--barcode-rename', 'skip', '--canonical-filtered', '--feature-info', genes,
                    '-F', 'anndata', '-o', converted], cwd=tmp_path, check=True)
    data = anndata.read_h5ad(converted)
    assert list(data.obs_names) == ['AAAA-1', 'CCCC-1']
    assert list(data.var_names) == ['HOST', 'REPORTER']
    np.testing.assert_array_equal(data.X.toarray(), [[4, 3], [0, 9]])
    assert data.var.loc['REPORTER', 'chrom'] == 'chrReporter'
    assert data.var.loc['REPORTER', 'gc_content'] == pytest.approx(100 * 4 / 6)
    assert data.var.loc['HOST', 'gene_biotype'] == 'protein_coding'

    table = pd.read_csv(genes, sep='\t').set_index('gene_id', verify_integrity=True)
    assert set(table.index) == {'HOST', 'REPORTER'}
    fields = ['gene_name', 'chrom', 'strand', 'start', 'end']
    assert table.loc['HOST', fields].tolist() == ['HostGene', 'chrHost', '+', 1, 9]
    assert table.loc['REPORTER', fields].tolist() == ['REPORTER', 'chrReporter', '-', 2, 11]
