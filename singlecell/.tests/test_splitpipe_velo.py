"""splitpipe_velo.py on a synthetic transcript-assignment table (#288).

Runs the script with DeprecationWarning as an error: polars deprecates an API before a later release removes it
(`collect(streaming=True)`, deprecated in 1.25 and removed in 2.0, broke this script), so this catches such calls
before the default container moves to that release.
"""
from pathlib import Path
import gzip
import subprocess
import sys

import numpy as np
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'singlecell/rules/quant/scripts/splitpipe_velo.py'
HEADER = 'bc_wells,genome,gene,gene_name,count,exonic,rt_type,bc_bcis,cell_barcode,polyN'
ROWS = [
    # bc_wells, gene, count, exonic, rt_type
    ('01_01_01__s1', 'GENE1', 1, 'True', 'R'),
    ('01_01_01__s1', 'GENE1', 2, 'True', 'T'),   # summed with the row above
    ('02_02_02__s1', 'GENE2', 4, 'True', 'R'),
    ('01_01_01__s1', 'GENE2', 5, 'False', 'T'),
    ('03_03_03__s1', 'GENE1', 6, 'False', 'R'),  # barcode with intronic reads only
    ('02_02_02__s1', 'GENE3', 7, 'True', 'X'),   # rt_type outside R/T is dropped; GENE3 is kept as a zero column
    ('02_02_02__s1', 'GENE1', 8, 'False', 'R'),
]
BARCODES = ['01_01_01__s1', '02_02_02__s1', '03_03_03__s1']
GENES = ['GENE1', 'GENE2', 'GENE3']
# barcodes x genes
SPLICED = np.array([[3, 0, 0], [0, 4, 0], [0, 0, 0]])
UNSPLICED = np.array([[0, 5, 0], [8, 0, 0], [6, 0, 0]])


def test_splitpipe_velo_writes_spliced_and_unspliced_matrices(tmp_path):
    table = tmp_path / 'tscp_assignment.csv.gz'
    with gzip.open(table, 'wt') as fh:
        fh.write(HEADER + '\n')
        for bc, gene, count, exonic, rt_type in ROWS:
            fh.write(f'{bc},synthetic_genome,{gene},{gene}_name,{count},{exonic},{rt_type},bcis,cell,A\n')
    outdir = tmp_path / 'velo'

    result = subprocess.run([sys.executable, '-W', 'error::DeprecationWarning', str(SCRIPT),
        '--input', str(table), '--output', str(outdir)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    assert (outdir / 'barcodes.tsv').read_text().splitlines() == BARCODES
    assert (outdir / 'genes.tsv').read_text().splitlines() == GENES
    for name, expected in [('spliced.mtx', SPLICED), ('unspliced.mtx', UNSPLICED)]:
        matrix = mmread(outdir / name)
        assert np.issubdtype(matrix.dtype, np.integer), name
        np.testing.assert_array_equal(matrix.toarray(), expected.T, err_msg=name)  # written genes x barcodes
