"""unitas_modifications.py: every declared table is written, also when no sample has a modification type (#339)."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/unitas_modifications.py'
TABLES = ('modification_per_position', 'internal_modification', 'internal_modification_position', 'NT_3p_extension')

# UNITAS 1.8.0 layout: a title line per block, a blank line between blocks, empty blocks as title only.
TAILS_ONLY = ("3'-tailings (non-template nucleotides)\nA\t5\nAA\t0\nC\t0\n\n"
              "Internal modifications (?->[ATGCN] = mapped sequence exceeds precursor sequence)\n\n"
              "Positions were internal modifications occur\n\n"
              "Modifications per position\n")
WITH_INTERNAL = ("3'-tailings (non-template nucleotides)\nA\t5\nU\t2\n\n"
                 "Internal modifications (?->[ATGCN] = mapped sequence exceeds precursor sequence)\nA->G\t3\nC->T\t1\n\n"
                 "Positions were internal modifications occur\n5\t3\n12\t1\n\n"
                 "Modifications per position\n\tA->G\tC->T\n5\t3\t0\n12\t0\t1\n")


def run(tmp_path, contents, script=SCRIPT, out='out'):
    files = []
    for i, text in enumerate(contents, 1):
        path = tmp_path / 'unitas' / 'S{:02d}'.format(i) / 'unitas.miR-modifications_Mouse.txt'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        files.append(str(path))
    result = subprocess.run([sys.executable, str(script), '-o', str(tmp_path / out)] + files,
                            capture_output=True, text=True, timeout=60)
    return result, {name: (tmp_path / out / (name + '.tsv')) for name in TABLES}


def test_tables_without_rows_are_written_with_their_header(tmp_path):
    result, tables = run(tmp_path, [TAILS_ONLY, TAILS_ONLY])
    assert result.returncode == 0, result.stderr
    assert tables['internal_modification'].read_text() == 'Modification\tCount\tSample_ID\n'
    assert tables['internal_modification_position'].read_text() == 'Position\tCount\tSample_ID\n'
    assert tables['modification_per_position'].read_text() == 'Sample_ID\n'
    lines = tables['NT_3p_extension'].read_text().splitlines()
    assert lines[0] == 'NT_Tail\tCount\tSample_ID' and 'A\t5\tS01' in lines and 'A\t5\tS02' in lines


# Output of the script before #339 for [WITH_INTERNAL, TAILS_ONLY]; tables with rows must not change.
EXPECTED_WITH_ROWS = {
    'modification_per_position': 'A->G\tC->T\tSample_ID\n3\t0\tS01\n0\t1\tS01\n',
    'internal_modification': 'Modification\tCount\tSample_ID\nA->G\t3\tS01\nC->T\t1\tS01\n',
    'internal_modification_position': 'Position\tCount\tSample_ID\n5\t3\tS01\n12\t1\tS01\n',
    'NT_3p_extension': 'NT_Tail\tCount\tSample_ID\nA\t5\tS01\nU\t2\tS01\nA\t5\tS02\nAA\t0\tS02\nC\t0\tS02\n'}


def test_tables_with_rows_are_unchanged(tmp_path):
    result, tables = run(tmp_path, [WITH_INTERNAL, TAILS_ONLY])
    assert result.returncode == 0, result.stderr
    for name in TABLES:
        assert tables[name].read_text() == EXPECTED_WITH_ROWS[name], name
