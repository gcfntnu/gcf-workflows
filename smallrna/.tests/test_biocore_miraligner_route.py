"""biocore_miraligner route (#426): selection, the database from gcfdb's miRBase files, the commands of the earlier pipeline's
miraligner counting, and the counting script."""
import importlib.util
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from test_quant_methods import MIRBASE_GFF, MIRBASE_GFF_ID, dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/miraligner_counts.py'
DB = f'data/ext/biocore_miraligner/mirbase-23.1/mus_musculus/gff-{MIRBASE_GFF_ID}/db.v1'  # named by the corrected GFF3 (#440)
OUT = 'data/tmp/smallrna/quant/biocore_miraligner'

spec = importlib.util.spec_from_file_location('biocore_miraligner_counts', SCRIPT)
biocore_miraligner_counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(biocore_miraligner_counts)


def test_route_plans_database_and_counting_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,biocore_miraligner', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule biocore_miraligner_count:') == 4
    db = rule_block(output, 'biocore_miraligner_db')
    assert ('--species mmu --hairpin data/ext/mirbase/mus_musculus/release-23.1/fasta/hairpin.fa '
            '--mature data/ext/mirbase/mus_musculus/release-23.1/fasta/mature.fa '
            f'--gff {MIRBASE_GFF} --output ' + DB) in db
    collapse = rule_block(output, 'biocore_miraligner_collapse')
    assert 'fastx_collapser -Q 33 -i data/tmp/smallrna/filter/cutadapt/SR0' in collapse and f'> {OUT}/SR0' in collapse
    assert 'if [ -s data/tmp/smallrna/filter/cutadapt/SR0' in collapse and f'else : > {OUT}/SR0' in collapse  # no reads: empty table
    align = rule_block(output, 'biocore_miraligner_align')
    assert f'miraligner -Xmx3g -sub 1 -trim 3 -add 3 -s mmu -freq -i {OUT}/SR0' in align and f'-db {DB} -o {OUT}/SR0' in align
    assert f'if [ -s {OUT}/SR0' in align  # miraligner stops on an empty table; such a sample gets the header only
    assert '-minl' not in align  # miraligner's default minimum length (16 nt), as in the pipeline
    assert 'mem_mb=4000' in align.replace(' ', '')
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'biocore_miraligner_counts')
    returncode, output = dry_run(tmp_path, 'unitas,biocore_miraligner', target='bfq_all')
    assert returncode == 0, output
    assert 'biocore_miraligner' not in output


def planned_shell(output, rule):
    """The shell command of the first planned job of a rule (dry run with -p)."""
    return rule_block(output, rule).split('Shell command: ', 1)[1].split('\n', 1)[0]


def test_collapsed_reads_become_sequence_and_count(tmp_path):
    """The collapse rule's awk command, its arguments exactly as planned (OFS escape included), on fastx_collapser output."""
    returncode, output = dry_run(tmp_path, 'biocore_miraligner', printshell=True)
    assert returncode == 0, output
    awk = shlex.split('awk' + planned_shell(output, 'biocore_miraligner_collapse').split('| awk', 1)[1].split(' > ', 1)[0])
    result = subprocess.run(awk, input='>1-120\nACGTACGTACGTACGTACG\n>2-3\nGGGGCCCC\n', capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'ACGTACGTACGTACGTACG\t120\nGGGGCCCC\t3\n'


def test_a_sample_without_sequences_gets_a_header_the_counts_read(tmp_path):
    """The alignment's branch for an empty table writes miraligner's header, which miraligner_counts.py reads as no rows."""
    returncode, output = dry_run(tmp_path, 'biocore_miraligner', printshell=True)
    assert returncode == 0, output
    header = re.search(r"printf '([^']+)'", planned_shell(output, 'biocore_miraligner_align')).group(1)
    mirna = tmp_path / 'S1.mirna'
    subprocess.run(['bash', '-c', 'printf {} > {}'.format(shlex.quote(header), shlex.quote(str(mirna)))], check=True)
    assert mirna.read_text() == HEADER
    assert biocore_miraligner_counts.count(mirna) == ({}, {}, 0.0, 0.0)


HEADER = 'seq\tname\tfreq\tmir\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\n'


def row(seq, freq, mir, mism='0', add='0', t5='0', t3='0', ambiguity=1, n=1):
    return '\t'.join([seq, 'seq_{}_x{}'.format(n, freq), str(freq), mir, '1', '22', mism, add, t5, t3, 'x', 'y', 'miRNA',
                      mir.replace('miR', 'mir'), str(ambiguity)]) + '\n'


def test_counts_divide_by_ambiguity_and_keep_substitution_free_rows(tmp_path):
    mirna = tmp_path / 'S1.mirna'
    mirna.write_text(HEADER + row('AAAA', 30, 'mmu-miR-1a-3p', ambiguity=2, n=1) + row('AAAA', 30, 'mmu-miR-1b-3p', ambiguity=2, n=1)
                     + row('CCCC', 10, 'mmu-miR-1a-3p', add='A', n=2) + row('GGGG', 7, 'mmu-miR-1a-3p', mism='8GT', n=3)
                     + row('TTTT', 3, 'mmu-miR-2-5p', ambiguity=3, n=4))
    mature, isomir, annotated, perfect = biocore_miraligner_counts.count(mirna)
    assert mature == {'mmu-miR-1a-3p': 25.0, 'mmu-miR-1b-3p': 15.0, 'mmu-miR-2-5p': 1.0}
    assert isomir['mmu-miR-1a-3p_GGGG_8GT;0;0;0'] == 7.0 and isomir['mmu-miR-1a-3p_CCCC_0;A;0;0'] == 10.0
    assert (annotated, perfect) == (50.0, 43.0)


def test_ties_best_keeps_a_sequences_best_annotations(tmp_path):
    """quant.biocore_miraligner.ties (#524): `split` divides a sequence between all its annotations; `best` keeps only
    its best: exact before a 3' addition alone before a substitution, the divisor being the ambiguity less the dropped."""
    mirna = tmp_path / 'S1.mirna'
    # miR-151a-5p's canonical read: exact on 151a, on 151b its first 20 bases plus an added T
    mirna.write_text(HEADER + row('TCGAGG', 40, 'hsa-miR-151a-5p', ambiguity=2, n=1) + row('TCGAGG', 40, 'hsa-miR-151b', add='T', ambiguity=2, n=1)
                     # exact on 27b, a substitution on 27a
                     + row('TTCACA', 30, 'hsa-miR-27b-3p', ambiguity=2, n=2) + row('TTCACA', 30, 'hsa-miR-27a-3p', mism='19CT', ambiguity=2, n=2)
                     # an addition on one, a substitution on two others: the addition wins
                     + row('AGCAGC', 12, 'hsa-miR-103a-3p', add='A', ambiguity=3, n=3) + row('AGCAGC', 12, 'hsa-miR-107', mism='5GA', ambiguity=3, n=3)
                     + row('AGCAGC', 12, 'hsa-miR-107b', mism='5GA', add='A', ambiguity=3, n=3)
                     # two exact annotations (identical mature sequences) and a substitution: the two exact share it
                     + row('GGGGGG', 10, 'hsa-miR-199a-3p', ambiguity=3, n=4) + row('GGGGGG', 10, 'hsa-miR-199b-3p', ambiguity=3, n=4)
                     + row('GGGGGG', 10, 'hsa-miR-199c', mism='3AG', ambiguity=3, n=4)
                     # substitutions only, one with an addition: the same rank, still divided (mature table: none)
                     + row('ACACAC', 12, 'hsa-miR-9a', mism='4TC', ambiguity=2, n=6) + row('ACACAC', 12, 'hsa-miR-9b', mism='4TC', add='T', ambiguity=2, n=6)
                     # end shifts do not rank: still divided
                     + row('CCCCCC', 8, 'hsa-miR-1a', t3='g', ambiguity=2, n=5) + row('CCCCCC', 8, 'hsa-miR-1b', t5='a', ambiguity=2, n=5))
    split = biocore_miraligner_counts.count(mirna)
    assert split == biocore_miraligner_counts.count(mirna, 'split')
    assert split[0] == {'hsa-miR-151a-5p': 20.0, 'hsa-miR-151b': 20.0, 'hsa-miR-27b-3p': 15.0, 'hsa-miR-103a-3p': 4.0,
                        'hsa-miR-199a-3p': 10 / 3, 'hsa-miR-199b-3p': 10 / 3, 'hsa-miR-1a': 4.0, 'hsa-miR-1b': 4.0}
    mature, isomir, annotated, perfect = biocore_miraligner_counts.count(mirna, 'best')
    assert mature == {'hsa-miR-151a-5p': 40.0, 'hsa-miR-27b-3p': 30.0, 'hsa-miR-103a-3p': 12.0,
                      'hsa-miR-199a-3p': 5.0, 'hsa-miR-199b-3p': 5.0, 'hsa-miR-1a': 4.0, 'hsa-miR-1b': 4.0}
    assert isomir == {'hsa-miR-151a-5p_TCGAGG_0;0;0;0': 40.0, 'hsa-miR-27b-3p_TTCACA_0;0;0;0': 30.0, 'hsa-miR-103a-3p_AGCAGC_0;A;0;0': 12.0,
                      'hsa-miR-199a-3p_GGGGGG_0;0;0;0': 5.0, 'hsa-miR-199b-3p_GGGGGG_0;0;0;0': 5.0,
                      'hsa-miR-1a_CCCCCC_0;0;0;g': 4.0, 'hsa-miR-1b_CCCCCC_0;0;a;0': 4.0,
                      'hsa-miR-9a_ACACAC_4TC;0;0;0': 6.0, 'hsa-miR-9b_ACACAC_4TC;T;0;0': 6.0}
    assert {k: v for k, v in split[1].items() if 'ACACAC' in k} == {'hsa-miR-9a_ACACAC_4TC;0;0;0': 6.0, 'hsa-miR-9b_ACACAC_4TC;T;0;0': 6.0}
    assert (annotated, perfect) == split[2:] == (112.0, 100.0)
    # the isomiR table's substitution-free rows still add up to the mature table
    assert sum(v for k, v in isomir.items() if k.rsplit('_', 1)[1].startswith('0;')) == sum(mature.values())
    with pytest.raises(SystemExit, match="--ties must be one of split, best, got 'all'"):
        biocore_miraligner_counts.count(mirna, 'all')


@pytest.mark.parametrize('rows', [
    [dict(ambiguity=1), dict(add='T', ambiguity=1)],     # fewer than its annotations
    [dict(ambiguity=2), dict(add='T', ambiguity=3)],     # two ambiguities
])
def test_ties_best_stops_on_an_ambiguity_it_cannot_reduce(tmp_path, rows):
    mirna = tmp_path / 'S1.mirna'
    mirna.write_text(HEADER + ''.join(row('AAAA', 6, 'hsa-miR-{}'.format(i), **r) for i, r in enumerate(rows)))
    with pytest.raises(SystemExit, match='--ties best needs one ambiguity per sequence'):
        biocore_miraligner_counts.count(mirna, 'best')
    biocore_miraligner_counts.count(mirna, 'split')   # split divides by miraligner's ambiguity as it is


def test_ties_setting_reaches_the_count_rule_and_settings(tmp_path):
    """The default keeps the count command as before; `best` adds --ties best; the choice is recorded with the table;
    other values stop at load."""
    returncode, output = dry_run(tmp_path, 'biocore_miraligner', printshell=True)
    assert returncode == 0, output
    assert '--ties' not in rule_block(output, 'biocore_miraligner_count')
    best = lambda value: (lambda s: s['quant'].setdefault('biocore_miraligner', {}).update(ties=value))
    returncode, output = dry_run(tmp_path / 'best', 'biocore_miraligner', best('best'), printshell=True)
    assert returncode == 0, output
    assert re.search(r'--stats \S+ --ties best', rule_block(output, 'biocore_miraligner_count'))
    returncode, output = dry_run(tmp_path / 'bad', 'biocore_miraligner', best('all'))
    assert returncode != 0 and "quant.biocore_miraligner.ties must be split or best, got 'all'" in output
    for value in (['best'], []):   # an empty list too, not taken for the defaults
        returncode, output = dry_run(tmp_path / 'list{}'.format(len(value)), 'biocore_miraligner', lambda s: s['quant'].update(biocore_miraligner=value))
        assert returncode != 0 and 'quant.biocore_miraligner must be a mapping of settings' in output, output


def test_count_command_takes_ties(tmp_path):
    reads, mirna = tmp_path / 'reads.tab', tmp_path / 'S1.mirna'
    reads.write_text('AAAA\t30\n')
    mirna.write_text(HEADER + row('AAAA', 30, 'mmu-miR-1a-3p', ambiguity=2) + row('AAAA', 30, 'mmu-miR-1b-3p', add='T', ambiguity=2))
    out = {name: tmp_path / '{}.tsv'.format(name) for name in ('mature', 'isomir', 'stats')}
    command = [sys.executable, str(SCRIPT), 'count', '--reads', str(reads), '--mirna', str(mirna)] + ['--{}={}'.format(k, v) for k, v in out.items()]
    for ties, expected in (([], 'mmu-miR-1a-3p\t15\nmmu-miR-1b-3p\t15\n'), (['--ties', 'best'], 'mmu-miR-1a-3p\t30\n')):
        result = subprocess.run(command + ties, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert out['mature'].read_text() == 'miRNA\tcount\n' + expected


def test_count_command_writes_tables_and_stats(tmp_path):
    reads, mirna = tmp_path / 'reads.tab', tmp_path / 'S1.mirna'
    reads.write_text('AAAA\t30\nCCCC\t10\nTTTTTTTTTTTTTTTT\t1\n')
    mirna.write_text(HEADER + row('AAAA', 30, 'mmu-miR-1a-3p', ambiguity=3) + row('CCCC', 10, 'mmu-miR-1a-3p', n=2))
    out = {name: tmp_path / '{}.tsv'.format(name) for name in ('mature', 'isomir', 'stats')}
    result = subprocess.run([sys.executable, str(SCRIPT), 'count', '--reads', str(reads), '--mirna', str(mirna)]
                            + ['--{}={}'.format(k, v) for k, v in out.items()], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert out['mature'].read_text() == 'miRNA\tcount\nmmu-miR-1a-3p\t20\n'
    assert out['stats'].read_text() == ('statistic\tvalue\nreads\t41\nsequences\t3\nannotated_reads\t40\n'
                                        'substitution_free_reads\t40\nmature_counts\t20\n')


def test_fractional_counts_are_written_exactly(tmp_path):
    assert biocore_miraligner_counts.number(10.0) == '10' and biocore_miraligner_counts.number(100 / 3) == repr(100 / 3)


@pytest.mark.parametrize('text, message', [
    ('seq\tname\tmir\n', 'missing columns freq, mism, add, t5, t3, DB, ambiguity'),
    (HEADER + 'AAAA\tseq_1_x3\t3\n', '3 fields, expected 15'),
    (HEADER + row('AAAA', 3, 'mmu-miR-1a-3p', ambiguity=0), 'ambiguity 0'),
    (HEADER + row('AAAA', 3, 'mmu-miR-1a-3p').replace('\t1\n', '\tNA\n'), "ambiguity 'NA' is not a number"),
])
def test_count_stops_on_output_it_cannot_read(tmp_path, text, message):
    mirna = tmp_path / 'S1.mirna'
    mirna.write_text(text)
    with pytest.raises(SystemExit, match=message):
        biocore_miraligner_counts.count(mirna)


def test_join_fills_missing_mirnas_with_zero(tmp_path):
    first, second = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    first.write_text('miRNA\tcount\nmmu-miR-1a-3p\t2.5\n')
    second.write_text('miRNA\tcount\nmmu-miR-2-5p\t3\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'miRNA\tS1\tS2\nmmu-miR-1a-3p\t2.5\t0\nmmu-miR-2-5p\t0\t3\n'
    second.write_text('miRNA\tcount\nmmu-miR-2-5p\t3\nmmu-miR-2-5p\t1\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and 'mmu-miR-2-5p twice' in result.stderr
