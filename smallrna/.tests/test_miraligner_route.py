"""miraligner route (#426): selection, the database from gcfdb's miRBase files, the commands of the earlier pipeline's
miraligner counting, and the counting script."""
import importlib.util
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/miraligner_counts.py'
DB = 'data/ext/miraligner/mirbase-23.1/mus_musculus/db.v1'
OUT = 'data/tmp/smallrna/quant/miraligner'

spec = importlib.util.spec_from_file_location('miraligner_counts', SCRIPT)
miraligner_counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(miraligner_counts)


def test_route_plans_database_and_counting_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,miraligner', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule miraligner_count:') == 4
    db = rule_block(output, 'miraligner_db')
    assert ('--species mmu --hairpin data/ext/mirbase/mus_musculus/release-23.1/fasta/hairpin.fa '
            '--mature data/ext/mirbase/mus_musculus/release-23.1/fasta/mature.fa '
            '--gff data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3 --output ' + DB) in db
    collapse = rule_block(output, 'miraligner_collapse')
    assert 'fastx_collapser -Q 33 -i data/tmp/smallrna/filter/cutadapt/SR0' in collapse and f'> {OUT}/SR0' in collapse
    assert 'if [ -s data/tmp/smallrna/filter/cutadapt/SR0' in collapse and f'else : > {OUT}/SR0' in collapse  # no reads: empty table
    align = rule_block(output, 'miraligner_align')
    assert f'miraligner -Xmx3g -sub 1 -trim 3 -add 3 -s mmu -freq -i {OUT}/SR0' in align and f'-db {DB} -o {OUT}/SR0' in align
    assert f'if [ -s {OUT}/SR0' in align  # miraligner stops on an empty table; such a sample gets the header only
    assert '-minl' not in align  # miraligner's default minimum length (16 nt), as in the pipeline
    assert 'mem_mb=4000' in align.replace(' ', '')
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'miraligner_counts')
    returncode, output = dry_run(tmp_path, 'unitas,miraligner', target='bfq_all')
    assert returncode == 0, output
    assert 'miraligner' not in output


def planned_shell(output, rule):
    """The shell command of the first planned job of a rule (dry run with -p)."""
    return rule_block(output, rule).split('Shell command: ', 1)[1].split('\n', 1)[0]


def test_collapsed_reads_become_sequence_and_count(tmp_path):
    """The collapse rule's awk command, its arguments exactly as planned (OFS escape included), on fastx_collapser output."""
    returncode, output = dry_run(tmp_path, 'miraligner', printshell=True)
    assert returncode == 0, output
    awk = shlex.split('awk' + planned_shell(output, 'miraligner_collapse').split('| awk', 1)[1].split(' > ', 1)[0])
    result = subprocess.run(awk, input='>1-120\nACGTACGTACGTACGTACG\n>2-3\nGGGGCCCC\n', capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'ACGTACGTACGTACGTACG\t120\nGGGGCCCC\t3\n'


def test_a_sample_without_sequences_gets_a_header_the_counts_read(tmp_path):
    """The alignment's branch for an empty table writes miraligner's header, which miraligner_counts.py reads as no rows."""
    returncode, output = dry_run(tmp_path, 'miraligner', printshell=True)
    assert returncode == 0, output
    header = re.search(r"printf '([^']+)'", planned_shell(output, 'miraligner_align')).group(1)
    mirna = tmp_path / 'S1.mirna'
    subprocess.run(['bash', '-c', 'printf {} > {}'.format(shlex.quote(header), shlex.quote(str(mirna)))], check=True)
    assert mirna.read_text() == HEADER
    assert miraligner_counts.count(mirna) == ({}, {}, 0.0, 0.0)


HEADER = 'seq\tname\tfreq\tmir\tstart\tend\tmism\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\tambiguity\n'


def row(seq, freq, mir, mism='0', add='0', t5='0', t3='0', ambiguity=1, n=1):
    return '\t'.join([seq, 'seq_{}_x{}'.format(n, freq), str(freq), mir, '1', '22', mism, add, t5, t3, 'x', 'y', 'miRNA',
                      mir.replace('miR', 'mir'), str(ambiguity)]) + '\n'


def test_counts_divide_by_ambiguity_and_keep_substitution_free_rows(tmp_path):
    mirna = tmp_path / 'S1.mirna'
    mirna.write_text(HEADER + row('AAAA', 30, 'mmu-miR-1a-3p', ambiguity=2, n=1) + row('AAAA', 30, 'mmu-miR-1b-3p', ambiguity=2, n=1)
                     + row('CCCC', 10, 'mmu-miR-1a-3p', add='A', n=2) + row('GGGG', 7, 'mmu-miR-1a-3p', mism='8GT', n=3)
                     + row('TTTT', 3, 'mmu-miR-2-5p', ambiguity=3, n=4))
    mature, isomir, annotated, perfect = miraligner_counts.count(mirna)
    assert mature == {'mmu-miR-1a-3p': 25.0, 'mmu-miR-1b-3p': 15.0, 'mmu-miR-2-5p': 1.0}
    assert isomir['mmu-miR-1a-3p_GGGG_8GT;0;0;0'] == 7.0 and isomir['mmu-miR-1a-3p_CCCC_0;A;0;0'] == 10.0
    assert (annotated, perfect) == (50.0, 43.0)


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
    assert miraligner_counts.number(10.0) == '10' and miraligner_counts.number(100 / 3) == repr(100 / 3)


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
        miraligner_counts.count(mirna)


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
