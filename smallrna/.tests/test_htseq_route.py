"""htseq route (#425): selection, the reference (a Bowtie 2 index of TEsmall's shared genome, tiny-count's renamed
miRBase GFF), the commands of the earlier pipeline's genome counter, guards, and the counting script."""
import gzip
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/htseq_counts.py'
REF = 'data/ext/htseq/mirbase-23.1_ensembl-116/mus_musculus/genome/bowtie2'
OUT = 'data/tmp/smallrna/quant/htseq'

spec = importlib.util.spec_from_file_location('htseq_counts', SCRIPT)
htseq_counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(htseq_counts)


def test_route_plans_reference_and_counting_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,htseq', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule htseq_count:') == 4
    index = rule_block(output, 'htseq_bowtie2_index')
    assert 'input: data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus/genome/genome.fa' in index
    assert f'bowtie2-build --threads 1 data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus/genome/genome.fa {REF}/genome >' in index
    assert 'rule tesmall_db:' not in output  # only the shared genome, not the TEsmall database
    assert f'fastx_collapser -Q 33 -i data/tmp/smallrna/filter/cutadapt/SR01_R1.fastq -o {OUT}/SR01/collapsed.fa' in output
    align = rule_block(output, 'htseq_align')
    assert "length($0) >= 15 {print name; print}" in align
    assert f'bowtie2 -p 1 -k 10 -f -x {REF}/genome -U - 2>' in align and '| gzip > ' in align
    assert 'mem_mb=5000' in align.replace(' ', '')  # the index in memory, so --resources caps parallel jobs
    annotate = rule_block(output, 'htseq_annotate')
    # secondary alignments scored, as HTSeq 0.9.1 did by default (HTSeq 2 ignores them unless told)
    assert ('htseq-count -f sam -a 0 -s yes -t miRNA -i Name --secondary-alignments score --supplementary-alignments score '
            '--samout {0}/SR0'.format(OUT)) in annotate
    assert ' - data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/mirna.v2.gff3 > ' in annotate
    assert '--species mmu --collapsed' in rule_block(output, 'htseq_count')
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'htseq_counts')
    returncode, output = dry_run(tmp_path, 'unitas,htseq', target='bfq_all')
    assert returncode == 0, output
    assert 'htseq' not in output


def test_guards_stop_at_load(tmp_path):
    returncode, output = dry_run(tmp_path, 'htseq', lambda s: s.update(organism='danio_rerio'))
    assert returncode != 0 and 'organism' in output  # the mirflowz reference checks apply to this route


def test_length_filter_keeps_sequences_of_at_least_15_nt(tmp_path):
    """The alignment's awk program, as planned, on fastx_collapser records: 15 nt and longer go to Bowtie 2."""
    returncode, output = dry_run(tmp_path, 'htseq', printshell=True)
    assert returncode == 0, output
    program = re.search(r"awk '([^']+)'", rule_block(output, 'htseq_align')).group(1)
    reads = tmp_path / 'collapsed.fa'
    reads.write_text('>1-10\n' + 'A' * 15 + '\n>2-5\n' + 'C' * 14 + '\n>3-1\n' + 'G' * 22 + '\n')
    result = subprocess.run(['awk', program, str(reads)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == '>1-10\n' + 'A' * 15 + '\n>3-1\n' + 'G' * 22 + '\n'


def record(read, flag, sequence, md, feature, chrom='1', pos=1000):
    """An htseq-count --samout record of a Bowtie 2 alignment."""
    quality = 'I' * len(sequence)
    pos = 0 if flag & 4 else pos
    tags = ['YT:Z:UU'] if flag & 4 else ['AS:i:0', 'XN:i:0', 'XM:i:0', 'XO:i:0', 'XG:i:0', 'NM:i:0', 'MD:Z:' + md, 'YT:Z:UU']
    return '\t'.join([read, str(flag), '*' if flag & 4 else chrom, str(pos), '0' if flag & 4 else '1',
                      '*' if flag & 4 else '{}M'.format(len(sequence)), '*', '0', '0', sequence, quality] + tags +
                     ['XF:Z:' + feature]) + '\n'


MIR21 = 'TAGCTTATCAGACTGATGTTGA'
MIR194 = 'TGTAACAGCAACTCCATGTGGA'


def annotated(tmp_path, records):
    path = tmp_path / 'annotated.sam.gz'
    with gzip.open(path, 'wt') as fh:
        fh.write('@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:1\tLN:2000\n')
        fh.writelines(records)
    return path


def test_counts_perfect_alignments_and_every_name_a_read_overlaps(tmp_path):
    path = annotated(tmp_path, [
        record('1-100', 0, MIR21, '22', 'mmu-miR-21a-5p'),
        record('1-100', 256, MIR21, '22', 'mmu-miR-21a-5p', chrom='2'),          # the same name at a second locus: once
        record('2-40', 0, MIR194, '22', '__ambiguous[mmu-miR-194-5p+mmu-miR-x-5p]'),
        record('2-40', 256, MIR194, '22', 'mmu-miR-y-5p', chrom='3'),             # another name at another locus: in full too
        record('3-7', 0, MIR21[:-1] + 'C', '21A0', 'mmu-miR-21a-5p'),            # a mismatch: not counted
        record('4-5', 0, 'ACGTACGTACGTACGTAC', '18', '__no_feature'),
        record('5-3', 4, 'GGGGGGGGGGGGGGGGGG', '', '__not_aligned'),
        record('6-2', 16, MIR21, '22', 'mmu-miR-21a-5p'),                         # minus strand: the read as sequenced
    ])
    mature, isomir, aligned, counted = htseq_counts.count('mmu', path)
    assert mature == {'mmu-miR-21a-5p': 102, 'mmu-miR-194-5p': 40, 'mmu-miR-x-5p': 40, 'mmu-miR-y-5p': 40}
    reverse = MIR21[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    assert isomir['mmu-miR-21a-5p_' + MIR21] == 100 and isomir['mmu-miR-21a-5p_' + reverse] == 2
    assert aligned == 100 + 40 + 7 + 5 + 2 and counted == 100 + 40 + 2


def test_perfect_alignment_test_is_the_pipelines_md_pattern():
    perfect = lambda md: bool(htseq_counts.PERFECT.search('x\tMD:Z:{}\tYT:Z:UU'.format(md)))
    assert perfect('22') and perfect('9') and perfect('15')
    assert not perfect('0A21') and not perfect('10A11') and not perfect('100')
    assert perfect('21^A1')  # the pattern reads only the first run of matches: a deletion after it passes, as in the pipeline


def test_count_command_writes_tables_and_stats(tmp_path):
    collapsed = tmp_path / 'collapsed.fa'
    collapsed.write_text('>1-100\n{}\n>2-7\n{}\n>3-1\nACGTAC\n'.format(MIR21, MIR194))
    path = annotated(tmp_path, [record('1-100', 0, MIR21, '22', '__ambiguous[mmu-miR-21a-5p+mmu-miR-21b]'),
                                record('2-7', 0, MIR194, '20A1', 'mmu-miR-194-5p')])
    out = {name: tmp_path / '{}.tsv'.format(name) for name in ('mature', 'isomir', 'stats')}
    result = subprocess.run([sys.executable, str(SCRIPT), 'count', '--species', 'mmu', '--collapsed', str(collapsed),
                             '--annotated', str(path)] + ['--{}={}'.format(k, v) for k, v in out.items()],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert out['mature'].read_text() == 'miRNA\tcount\nmmu-miR-21a-5p\t100\nmmu-miR-21b\t100\n'
    assert out['isomir'].read_text().splitlines()[1] == 'mmu-miR-21a-5p_{}\t100'.format(MIR21)
    assert out['stats'].read_text() == ('statistic\tvalue\nreads\t108\nsequences\t3\naligned_reads\t107\n'
                                        'counted_reads\t100\nmature_counts\t200\n')


@pytest.mark.parametrize('line, message', [
    ('1-100\t0\t1\t1\t1\t22M\t*\t0\t0\tACGT\tIIII\tMD:Z:22\n', "without htseq-count's XF tag last"),
    ('read1\t0\t1\t1\t1\t22M\t*\t0\t0\tACGT\tIIII\tMD:Z:22\tXF:Z:mmu-miR-1\n', 'is not <rank>-<count>'),
])
def test_count_stops_on_records_it_cannot_read(tmp_path, line, message):
    with pytest.raises(SystemExit, match=message):
        htseq_counts.count('mmu', annotated(tmp_path, [line]))


def test_join_fills_missing_mirnas_with_zero(tmp_path):
    first, second = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    first.write_text('miRNA\tcount\nmmu-miR-1a-3p\t5\nmmu-miR-21a-5p\t7\n')
    second.write_text('miRNA\tcount\nmmu-miR-21a-5p\t3\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'miRNA\tS1\tS2\nmmu-miR-1a-3p\t5\t0\nmmu-miR-21a-5p\t7\t3\n'
    second.write_text('name\tcount\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and 'expected the header miRNA, count' in result.stderr
