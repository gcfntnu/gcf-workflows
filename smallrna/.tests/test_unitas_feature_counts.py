"""UNITAS feature counts (`unitas_allfeatures`) count each Ensembl ncRNA gene once (gcfntnu/gcf-docker#56).

Dumps from 20-11-2026 on name a gene's further ncRNA transcripts `<symbol>-<gene>_(n)`; older dumps wrote
`<symbol>-(n)` (and kept the FASTA `>` in names without a gene symbol).
"""
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/unitas_feature_counts.py'

# sequence, reads, annotations: rows as in unitas.full_annotation_matrix.txt
S1 = [
    ('AAAAAAAAAAAAAAAAAAAA', 10, ['lncRNA|LINCX-ENSG00000000001.1 (1-20)']),
    ('CCCCCCCCCCCCCCCCCCCC', 4, ['lncRNA|LINCX-ENSG00000000001.1_(2) (5-24)']),
    ('CACACACACACACACACACA', 4, ['lncRNA|LINCX-ENSG00000000001.1_(12) (1-20)', 'processed_pseudogene|PSX (2-21)']),
    ('GGGGGGGGGGGGGGGGGGGG', 6, ['lncRNA|LINCX-ENSG00000000002.1_(2) (1-20)']),
    ('TTTTTTTTTTTTTTTTTTTT', 3, ['lncRNA|ENST00000000005.1-ENSG00000000005.1 (3-22)']),
    ('ACACACACACACACACACAC', 2, ['lncRNA|OLDX-(2) (1-20)']),
    ('TAGCTTATCAGACTGATGTTGA', 8, ['miR-21-5p']),
    ('GCATTGGTGGTTCAGTGGTAGA', 7, ['tRNA|tRNA-Gly-GCC-1-1 (1-22)']),
    ('TGTGTGTGTGTGTGTGTGTG', 5, []),
]
S2 = [('AAAAAAAAAAAAAAAAAAAA', 1, ['lncRNA|LINCX-ENSG00000000001.1 (1-20)'])]
EXPECTED = {
    'lncRNA|LINCX-ENSG00000000001.1': (16.0, 1.0),  # 10 + 4 + 4 / 2: all its transcripts
    'lncRNA|LINCX-ENSG00000000002.1': (6.0, 0.0),   # the other gene with the same symbol stays apart
    'lncRNA|ENST00000000005.1-ENSG00000000005.1': (3.0, 0.0),
    'lncRNA|OLDX': (2.0, 0.0),                      # an older dump's '-(n)'
    'miR-21-5p': (8.0, 0.0),
    'processed_pseudogene|PSX': (2.0, 0.0),
    'tRNA-Gly-GCC': (7.0, 0.0),
    'no-annotation': (5.0, 0.0),
}


def write_matrix(path, rows):
    path.parent.mkdir(parents=True)
    lines = ['sequence\treads\tannotation(s)\n']
    lines += [f'{seq}\t{reads}\t' + ''.join(a + '\t' for a in annotations) + '\n' for seq, reads, annotations in rows]
    path.write_text(''.join(lines))


def test_transcripts_of_a_gene_count_as_one_feature(tmp_path):
    inputs = []
    for sample, rows in (('S1', S1), ('S2', S2)):
        inputs.append(tmp_path / sample / 'unitas.full_annotation_matrix.txt')
        write_matrix(inputs[-1], rows)
    output = tmp_path / 'out' / 'allfeatures.tsv'
    subprocess.run([sys.executable, str(SCRIPT), '--frac-counts', '--output', str(output)] + [str(p) for p in inputs],
                   check=True, capture_output=True, text=True)

    table = pd.read_csv(output, sep='\t', index_col='feature_id')
    assert list(table.columns) == ['S1', 'S2']
    assert {name: (row.S1, row.S2) for name, row in table.iterrows()} == EXPECTED

    annotations = pd.read_csv(str(output) + '.annotations', sep='\t').set_index('feature_name')['feature_class']
    assert annotations['LINCX-ENSG00000000001.1'] == 'lncRNA'
    assert annotations['ENST00000000005.1-ENSG00000000005.1'] == 'lncRNA'
    assert annotations['miR-21-5p'] == 'miRNA'
    assert annotations['tRNA-Gly-GCC'] == 'tRNA'
