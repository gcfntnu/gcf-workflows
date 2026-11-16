"""miRge3 miRNA counts from mapped.csv, with a project-level exact-read rule instead of miRge3's per-sample one (#349):
synthetic mapped.csv, merge file and miRge3 table."""
import importlib.util
import io
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'rules/quant/scripts/mirge3_counts.py'
spec = importlib.util.spec_from_file_location('mirge3_counts', SCRIPT)
counts_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(counts_module)

HEADER = 'Sequence,annotFlag,exact miRNA,hairpin miRNA,mature tRNA,primary tRNA,snoRNA,rRNA,ncrna others,mRNA,isomiR miRNA,S1,S2\n'


def mapped_row(seq, exact='', isomir='', s1=0, s2=0, other=''):
    return '{},1,{},,,,{},,,,{},{},{}\n'.format(seq, exact, other, isomir, s1, s2)


MAPPED = HEADER + ''.join([
    mapped_row('AAAA', exact='mmu-miR-1a-3p', s1=10, s2=1),        # exact; S2 below miRge3's 2-exact-read rule
    mapped_row('AAAAC', isomir='mmu-miR-1a-3p', s1=5, s2=4),      # isomiR of the same miRNA
    mapped_row('CCCC', isomir='mmu-miR-2-5p', s1=3, s2=0),        # isomiR only: dropped
    mapped_row('GGGG', exact='mmu-miR-3a-5p', s1=2, s2=0),        # merged with mmu-miR-3b-5p (identical sequence)
    mapped_row('GGGGT', isomir='mmu-miR-3b-5p', s1=1, s2=1),
    mapped_row('TTTT', exact='mmu-miR-4-3p.SNPC', s1=0, s2=6),    # SNP variant, folded into mmu-miR-4-3p
    mapped_row('TTTA', exact='mmu-miR-4-3p', s1=1, s2=0),
    mapped_row('CACA', exact='mmu-miR-6-5p', s1=1, s2=1),         # 1 exact read per sample, 2 in the project: kept
    mapped_row('GTGT', exact='mmu-miR-7-3p', s1=1, s2=0),         # 1 exact read in the project: dropped
    mapped_row('GTGTA', isomir='mmu-miR-7-3p', s1=0, s2=3),
    mapped_row('ACGT', other='ENSMUST00000000001', s1=50, s2=50),  # not a miRNA
])
MERGES = 'mmu-miR-3a-5p/3b-5p,mmu-miR-3a-5p,mmu-miR-3b-5p\nmmu-miR-4-3p,mmu-miR-4-3p.SNPC,mmu-miR-4-3p.SNPA\n'
NAMES = ['mmu-miR-1a-3p', 'mmu-miR-2-5p', 'mmu-miR-3a-5p/3b-5p', 'mmu-miR-4-3p', 'mmu-miR-5-5p', 'mmu-miR-6-5p',
         'mmu-miR-7-3p']


def merges(tmp_path):
    path = tmp_path / 'merges.csv'
    path.write_text(MERGES)
    return counts_module.read_merges(path)


def test_mirnas_with_two_exact_reads_in_the_project_keep_all_their_reads(tmp_path):
    counts = counts_module.count(io.StringIO(MAPPED), merges(tmp_path), NAMES, ['S1', 'S2'])
    assert counts['mmu-miR-1a-3p'] == {'S1': 15, 'S2': 5}  # 1 exact read in S2: miRge3 gives 0
    assert counts['mmu-miR-3a-5p/3b-5p'] == {'S1': 3, 'S2': 1}  # merged names
    assert counts['mmu-miR-4-3p'] == {'S1': 1, 'S2': 6}  # SNP variant folded in
    assert counts['mmu-miR-6-5p'] == {'S1': 1, 'S2': 1}  # 2 exact reads in the project: miRge3 gives 0 in both
    assert counts['mmu-miR-2-5p'] == {'S1': 0, 'S2': 0}  # isomiRs only
    assert counts['mmu-miR-7-3p'] == {'S1': 0, 'S2': 0}  # 1 exact read in the project
    assert counts['mmu-miR-5-5p'] == {'S1': 0, 'S2': 0}  # library miRNA without reads keeps its row


@pytest.mark.parametrize('row, message', [
    (mapped_row('GATC', exact='mmu-miR-1a-3p', isomir='mmu-miR-2-5p', s1=1), 'both an exact'),
    (mapped_row('GATC', exact='mmu-miR-9-5p', s1=1), "'mmu-miR-9-5p' .* is not a row"),
])
def test_unexpected_reads_stop(tmp_path, row, message):
    with pytest.raises(ValueError, match=message):
        counts_module.count(io.StringIO(MAPPED + row), merges(tmp_path), NAMES, ['S1', 'S2'])


def test_samples_must_match_the_mirge3_table(tmp_path):
    with pytest.raises(ValueError, match="missing \\['S3'\\]"):
        counts_module.count(io.StringIO(MAPPED), merges(tmp_path), NAMES, ['S1', 'S2', 'S3'])


def test_command_line_keeps_the_mirge3_table_format(tmp_path):
    (tmp_path / 'mapped.csv').write_text(MAPPED)
    (tmp_path / 'merges.csv').write_text(MERGES)
    (tmp_path / 'table.csv').write_text('miRNA,S1,S2\n' + ''.join(name + ',0.0,0.0\n' for name in NAMES))
    result = subprocess.run([sys.executable, str(SCRIPT), '--mapped', str(tmp_path / 'mapped.csv'), '--merges',
                             str(tmp_path / 'merges.csv'), '--mirge3-table', str(tmp_path / 'table.csv')],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ['miRNA,S1,S2', 'mmu-miR-1a-3p,15,5', 'mmu-miR-2-5p,0,0',
                                          'mmu-miR-3a-5p/3b-5p,3,1', 'mmu-miR-4-3p,1,6', 'mmu-miR-5-5p,0,0',
                                          'mmu-miR-6-5p,1,1', 'mmu-miR-7-3p,0,0']


def test_spike_in_column_is_an_annotation_not_a_sample(tmp_path):
    """With --spikeIn, miRge3 adds a `spike-in` column after `isomiR miRNA`; spike-in reads are not miRNA reads."""
    columns = HEADER.strip().replace('isomiR miRNA,S1,S2', 'isomiR miRNA,spike-in,S1,S2').split(',')

    def row(**values):
        return ','.join(str(values.get(c.replace(' ', '_').replace('-', '_'), '')) for c in columns) + '\n'

    mapped = ','.join(columns) + '\n' + ''.join([
        row(Sequence='AAAA', annotFlag=1, exact_miRNA='mmu-miR-1a-3p', S1=10, S2=1),
        row(Sequence='AAAAC', annotFlag=1, isomiR_miRNA='mmu-miR-1a-3p', S1=5, S2=4),
        row(Sequence='TTTTGG', annotFlag=1, spike_in='ERCC-00002', S1=100, S2=200),  # spike-in only
    ])
    counts = counts_module.count(io.StringIO(mapped), merges(tmp_path), NAMES, ['S1', 'S2'])
    assert counts['mmu-miR-1a-3p'] == {'S1': 15, 'S2': 5}
    assert all(sum(c.values()) == 0 for name, c in counts.items() if name != 'mmu-miR-1a-3p')


def test_missing_and_extra_samples_are_still_detected(tmp_path):
    header = HEADER.replace('isomiR miRNA,S1,S2', 'isomiR miRNA,spike-in,S1,S3')
    with pytest.raises(ValueError, match=r"missing \['S2'\], extra \['S3'\]"):
        counts_module.count(io.StringIO(header), merges(tmp_path), NAMES, ['S1', 'S2'])
