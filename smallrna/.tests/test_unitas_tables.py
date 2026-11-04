"""Delivered UNITAS tables have one orientation and a run-independent row order (#265)."""
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'smallrna/rules/quant/scripts'
SAMPLES = {'S1': 1, 'S2': 2}

MIRS = [('mmu-miR-1a-3p', 10), ('mmu-let-7a-5p', 30), ('mmu-miR-21a-5p', 20)]
TRFS = [('tRNA-Gly-GCC', 5), ('tRNA-Glu-CTC', 7)]
ISOMIRS = [('TGAGGTAGTAGGTTGTATAGTT', 'mmu-let-7a-5p', 12), ('TGAGGTAGTAGGTTGTATAGT', 'mmu-let-7a-5p(3p-1)', 4),
           ('TAGCTTATCAGACTGATGTTGA', 'mmu-miR-21a-5p', 9)]
ISOTRFS = [('tRNA-Gly-GCC', "5'-tR-half", 3), ('tRNA-Glu-CTC', "5'-tRF", 2)]
ANNOTATION = [('miRNA', 0, 60), ('mmu-let-7a-5p', 1, 30), ('mmu-miR-21a-5p', 1, 20),
              ('tRNA', 0, 12), ("5'-tR-half", 1, 7), ('tRNA-Gly-GCC', 2, 5)]


def write_unitas(base, reverse):
    """Per-sample UNITAS-style files; `reverse` lists rows in the opposite order."""
    order = (lambda rows: rows[::-1]) if reverse else (lambda rows: rows)
    for sample, scale in SAMPLES.items():
        out = base / sample
        out.mkdir(parents=True)
        (out / 'miR-table_Mouse.simplified.txt').write_text(
            'miRNA\treads\n' + ''.join(f'{name}\t{count * scale}\n' for name, count in order(MIRS)))
        (out / 'unitas.tRF-table.simplified.txt').write_text(
            'tRNA\treads\tpre-tRNA\n' + ''.join(f'{name}\t{count * scale}\t0\n' for name, count in order(TRFS)))
        (out / 'unitas.miR-table_Mouse.txt').write_text(
            'sequence\tmiR-name\ttotal_reads\n' + ''.join(f'{seq}\t{name}\t{count * scale}\n' for seq, name, count in order(ISOMIRS)))
        (out / 'unitas.tRF-table.txt').write_text(
            'header 1\nheader 2\nheader 3\ntRNA\ttRF-type\treads\n'
            + ''.join(f'{trna}\t{kind}\t{count * scale}\n' for trna, kind, count in order(ISOTRFS)))
        # Annotation summary: three spaces per level. Children follow their parent; parents may be reordered.
        blocks = [[r for r in ANNOTATION[:3]], [r for r in ANNOTATION[3:]]]
        lines = [f"{'   ' * level}{name}\t{count * scale}\n" for block in order(blocks) for name, level, count in block]
        (out / 'unitas.annotation_summary.txt').write_text(''.join(lines))


def build(tmp_path, reverse):
    base = tmp_path / ('reversed' if reverse else 'forward')
    write_unitas(base / 'unitas', reverse)
    tables = base / 'tables'
    tables.mkdir()
    files = lambda name: [str(base / 'unitas' / s / name) for s in SAMPLES]
    def run(script, args, stdout=None):
        with open(stdout, 'w') if stdout else open('/dev/null', 'w') as out:
            subprocess.run([sys.executable, str(SCRIPTS / script), *args], check=True, stdout=out)
    run('unitas_mirtable.py', files('miR-table_Mouse.simplified.txt'), tables / 'mir_counts.tsv')
    run('unitas_trftable.py', files('unitas.tRF-table.simplified.txt'), tables / 'trf_counts.tsv')
    run('unitas_isomirtable.py', [*files('unitas.miR-table_Mouse.txt'), '--output', str(tables / 'isomir_counts.tsv')])
    run('unitas_isotrftable.py', files('unitas.tRF-table.txt'), tables / 'isotrf_counts.tsv')
    run('unitas_annotations.py', files('unitas.annotation_summary.txt'), tables / 'annotations.tsv')
    return tables


@pytest.fixture
def tables(tmp_path):
    return build(tmp_path, reverse=False), build(tmp_path, reverse=True)


@pytest.mark.parametrize('name', ['mir_counts.tsv', 'trf_counts.tsv', 'isomir_counts.tsv', 'isomir_counts_anno.tsv',
                                  'isotrf_counts.tsv', 'annotations.tsv'])
def test_tables_do_not_depend_on_unitas_row_order(tables, name):
    forward, reverse = tables
    assert (forward / name).read_bytes() == (reverse / name).read_bytes()


def test_count_tables_have_features_as_rows_and_samples_as_columns(tables):
    forward, _ = tables
    for name in ['mir_counts.tsv', 'trf_counts.tsv', 'isomir_counts.tsv']:
        table = pd.read_csv(forward / name, sep='\t', index_col=0)
        assert list(table.columns) == list(SAMPLES), name
        assert list(table.index) == sorted(table.index), name
    isomirs = pd.read_csv(forward / 'isomir_counts.tsv', sep='\t', index_col=0)
    assert isomirs.loc['TGAGGTAGTAGGTTGTATAGTT|mmu-let-7a-5p', 'S2'] == 24
    anno = pd.read_csv(forward / 'isomir_counts_anno.tsv', sep='\t', index_col=0)
    assert list(anno.index) == list(isomirs.index)


def test_isotrf_table_is_tab_separated(tables):
    forward, _ = tables
    header = (forward / 'isotrf_counts.tsv').read_text().splitlines()[0]
    assert header.split('\t') == ['tRNA', 'tRF-type', 'reads', 'Sample_ID']


def test_annotation_rows_carry_their_parent_levels(tables):
    forward, _ = tables
    table = pd.read_csv(forward / 'annotations.tsv', sep='\t', keep_default_na=False)
    s1 = table[table.Sample_ID == 'S1']
    assert list(zip(s1.Level1, s1.Level2, s1.Level3, s1.Count)) == [
        ('miRNA', '', '', 60), ('miRNA', 'mmu-let-7a-5p', '', 30), ('miRNA', 'mmu-miR-21a-5p', '', 20),
        ('tRNA', '', '', 12), ('tRNA', "5'-tR-half", '', 7), ('tRNA', "5'-tR-half", 'tRNA-Gly-GCC', 5)]
