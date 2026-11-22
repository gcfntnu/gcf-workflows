"""biocore_htseq route (#425): selection, the reference (a Bowtie 2 index of TEsmall's shared genome, tiny-count's renamed
miRBase GFF), the commands of the earlier pipeline's genome counter, guards, and the counting script."""
import gzip
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

from test_quant_methods import MIRBASE_GFF_ID, dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/htseq_counts.py'
REF = 'data/ext/biocore_htseq/mirbase-23.1_ensembl-116/mus_musculus/genome/bowtie2'
OUT = 'data/tmp/smallrna/quant/biocore_htseq'

spec = importlib.util.spec_from_file_location('biocore_htseq_counts', SCRIPT)
biocore_htseq_counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(biocore_htseq_counts)


def test_route_plans_reference_and_counting_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,biocore_htseq', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output and output.count('rule biocore_htseq_count:') == 4
    index = rule_block(output, 'biocore_htseq_bowtie2_index')
    assert 'input: data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus/genome/genome.fa' in index
    assert f'bowtie2-build --threads 1 data/ext/tesmall/mirbase-23.1_ensembl-116/mus_musculus/genome/genome.fa {REF}/genome >' in index
    assert 'rule tesmall_db:' not in output  # only the shared genome, not the TEsmall database
    collapse = rule_block(output, 'biocore_htseq_collapse')
    assert ('if [ -s data/tmp/smallrna/filter/cutadapt/SR0' in collapse and f'fastx_collapser -Q 33 -i data/tmp/smallrna/filter/cutadapt/SR0'
            in collapse and f'else : > {OUT}/SR0' in collapse)  # no filtered reads: no sequences, zero counts
    align = rule_block(output, 'biocore_htseq_align')
    assert "length($0) >= 15 {print name; print}" in align
    assert f'bowtie2 -p 1 -k 10 -f -x {REF}/genome -U - 2>' in align and '| gzip > ' in align
    assert 'mem_mb=5000' in align.replace(' ', '')  # the index in memory, so --resources caps parallel jobs
    annotate = rule_block(output, 'biocore_htseq_annotate')
    # secondary alignments scored, as HTSeq 0.9.1 did by default (HTSeq 2 ignores them unless told)
    assert ('htseq-count -f sam -a 0 -s yes -t miRNA -i Name --secondary-alignments score --supplementary-alignments score '
            '--samout {0}/SR0'.format(OUT)) in annotate
    assert f' - {OUT}/mirna.gff3 > ' in annotate  # tiny-count's renamed GFF without rows that have no strand
    assert '|| { rm -f ' in annotate  # a failed htseq-count leaves no uncompressed samout behind
    gff = f'gff --input data/ext/tinyrna/mirbase-23.1_ensembl-116/mus_musculus/gff-{MIRBASE_GFF_ID}/mirna.v2.gff3 > {OUT}/mirna.gff3'
    assert gff in rule_block(output, 'biocore_htseq_gff')
    assert '--species mmu --collapsed' in rule_block(output, 'biocore_htseq_count')
    assert '--samples SR01 SR02 SR03 SR04' in rule_block(output, 'biocore_htseq_counts')
    returncode, output = dry_run(tmp_path, 'unitas,biocore_htseq', target='bfq_all')
    assert returncode == 0, output
    assert 'biocore_htseq' not in output


def test_guards_stop_at_load(tmp_path):
    returncode, output = dry_run(tmp_path, 'biocore_htseq', lambda s: s.update(organism='danio_rerio'))
    assert returncode != 0 and 'organism' in output  # the mirflowz reference checks apply to this route


def test_length_filter_keeps_sequences_of_at_least_15_nt(tmp_path):
    """The alignment's awk program, as planned, on fastx_collapser records: 15 nt and longer go to Bowtie 2."""
    returncode, output = dry_run(tmp_path, 'biocore_htseq', printshell=True)
    assert returncode == 0, output
    program = re.search(r"awk '([^']+)'", rule_block(output, 'biocore_htseq_align')).group(1)
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
    mature, isomir, aligned, counted = biocore_htseq_counts.count('mmu', path)
    assert mature == {'mmu-miR-21a-5p': 102, 'mmu-miR-194-5p': 40, 'mmu-miR-x-5p': 40, 'mmu-miR-y-5p': 40}
    reverse = MIR21[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    assert isomir['mmu-miR-21a-5p_' + MIR21] == 100 and isomir['mmu-miR-21a-5p_' + reverse] == 2
    assert aligned == 100 + 40 + 7 + 5 + 2 and counted == 100 + 40 + 2


def test_perfect_alignment_test_is_the_pipelines_md_pattern():
    perfect = lambda md: bool(biocore_htseq_counts.PERFECT.search('x\tMD:Z:{}\tYT:Z:UU'.format(md)))
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
        biocore_htseq_counts.count('mmu', annotated(tmp_path, [line]))


def test_join_fills_missing_mirnas_with_zero(tmp_path):
    first, second = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    first.write_text('miRNA\tcount\nmmu-miR-1a-3p\t5\nmmu-miR-21a-5p\t7\n')
    second.write_text('miRNA\tcount\nmmu-miR-21a-5p\t3\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'miRNA\tS1\tS2\nmmu-miR-1a-3p\t5\t0\nmmu-miR-21a-5p\t7\t3\n'
    second.write_text('miRNA\tcount\nmmu-miR-21a-5p\t3\nmmu-miR-21a-5p\t4\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and 'mmu-miR-21a-5p twice' in result.stderr
    second.write_text('name\tcount\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and 'expected the header miRNA, count' in result.stderr


def test_gff_leaves_out_rows_without_a_strand(tmp_path):
    gff = tmp_path / 'mirna.gff3'
    gff.write_text('##gff-version 3\n# source x\n'
                   '1\t.\tmiRNA_primary_transcript\t10\t80\t.\t+\t.\tID=MI1;Name=hsa-mir-1\n'
                   '1\t.\tmiRNA_primary_transcript\t90\t150\t.\t.\t.\tID=MI2;Name=hsa-mir-13168\n'
                   '1\t.\tmiRNA\t100\t122\t.\t.\t.\tID=MIMAT2;Name=hsa-miR-13168-5p;Derives_from=MI2\n'
                   '1\t.\tmiRNA\t30\t51\t.\t-\t.\tID=MIMAT1;Name=hsa-miR-1-5p;Derives_from=MI1\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'gff', '--input', str(gff)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0] == '##gff-version 3'
    assert lines[1] == '# htseq_counts.py gff: rows without a strand left out (2): hsa-mir-13168,hsa-miR-13168-5p'
    assert [line.split('\t')[8].split(';')[0] for line in lines if not line.startswith('#')] == ['ID=MI1', 'ID=MIMAT1']


def test_local_alignment_mode(tmp_path):
    """quant.biocore_htseq.alignment: local (#417): Bowtie 2 --local with a minimum score of 2 per base for
    local_min_aligned of the read, and a filter that drops alignments clipping the read's 5' end or aligning too little
    of it; the default keeps the command as before; bad values stop at load."""
    local = lambda **kw: (lambda s: s['quant'].setdefault('biocore_htseq', {}).update(alignment='local', **kw))
    returncode, output = dry_run(tmp_path, 'biocore_htseq', local(), printshell=True)
    assert returncode == 0, output
    align = rule_block(output, 'biocore_htseq_align')
    assert re.search(r'bowtie2 -p \d+ -k 10 -f --local -L 16 --score-min L,0,1\.8 -x ', align)
    command = align.split('2> ', 1)[1].split(' | ', 1)[1].rsplit(' | gzip', 1)[0]
    assert command.startswith('awk ')
    header = '@SQ\tSN:1\tLN:1000\n'
    records = {'exact': ('r1', 0, '22M'), 'tail': ('r2', 0, '20M2S'), 'clipped5': ('r3', 0, '2S20M'), 'reverse_tail': ('r4', 16, '2S20M'),
               'reverse_clipped5': ('r5', 16, '20M2S'), 'short': ('r6', 0, '16M6S'), 'unmapped': ('r7', 4, '*'),
               'at_the_cutoff': ('r9', 0, '18M2S')}   # 18 of 20 aligned at share 0.9: kept (a 20-base read below)
    sam = header + ''.join('{}\t{}\t1\t100\t255\t{}\t*\t0\t0\t{}\t{}\tAS:i:40\n'.format(n, f, c, 'A' * (20 if n == 'r9' else 22),
                                                                                       'I' * (20 if n == 'r9' else 22)) for n, f, c in records.values())
    # a read aligned exactly at one locus and clipped at a relative's: only the best-scoring alignment is kept
    sam += 'r8\t0\t1\t500\t255\t22M\t*\t0\t0\t{0}\t{1}\tAS:i:44\nr8\t256\t2\t900\t255\t20M2S\t*\t0\t0\t{0}\t{1}\tAS:i:40\n'.format('A' * 22, 'I' * 22)
    kept = subprocess.run(['bash', '-c', command], input=sam, capture_output=True, text=True, check=True).stdout
    records_kept = [line.split('\t') for line in kept.splitlines() if not line.startswith('@')]
    assert kept.startswith('@SQ') and [f[0] for f in records_kept] == ['r1', 'r2', 'r4', 'r7', 'r9', 'r8']   # 5' clips and a 73% alignment dropped
    assert [f[2] for f in records_kept if f[0] == 'r8'] == ['1']
    returncode, output = dry_run(tmp_path / 'default', 'biocore_htseq', printshell=True)
    assert returncode == 0 and re.search(r'bowtie2 -p \d+ -k 10 -f -x ', output) and '--local' not in output
    # requested directly, without quant.method selecting the route: the same options
    returncode, output = dry_run(tmp_path / 'direct', 'unitas', local(), printshell=True, target=OUT + '/htseq_counts.tsv')
    assert returncode == 0 and '--local -L 16' in rule_block(output, 'biocore_htseq_align'), output
    returncode, output = dry_run(tmp_path / 'list', 'biocore_htseq', lambda s: s['quant'].update(biocore_htseq=['local']))
    assert returncode != 0 and 'quant.biocore_htseq must be a mapping of settings' in output
    for value, message in ((dict(local_min_aligned=0.3), 'quant.biocore_htseq.local_min_aligned must be a number from 0.5 to 1'),
                           (dict(local_max_5p_clip=-1), 'quant.biocore_htseq.local_max_5p_clip must be a whole number from 0 to 5')):
        returncode, output = dry_run(tmp_path / str(len(message)), 'biocore_htseq', local(**value))
        assert returncode != 0 and message in output, output
