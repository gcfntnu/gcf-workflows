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
    # mouse: no ncRNA tables (#503), RNAcentral 9.0 has mouse coordinates on GRCm38 only
    assert 'biocore_htseq_ncrna' not in output and 'biocore_htseq_rnacentral' not in output
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


def test_human_route_annotates_the_project_in_one_call_and_writes_ncrna_tables(tmp_path):
    """Human on the Ensembl reference (#503): the pipeline's annotation (miRBase GFF, RNAcentral 9.0 transcripts, the
    reference's Ensembl genes) and one htseq-count call for every sample of the project, the annotation loaded once,
    then the ncRNA tables per sample and their join."""
    returncode, output = dry_run(tmp_path, 'biocore_htseq', lambda s: s.update(organism='homo_sapiens'), printshell=True)
    assert returncode == 0, output
    out, ref = 'data/tmp/smallrna/quant/biocore_htseq', 'data/ext/biocore_htseq/mirbase-23.1_ensembl-116/homo_sapiens'
    assert output.count('rule biocore_htseq_annotate:') == 1 and 'rule biocore_htseq_gff:' not in output
    annotate = rule_block(output, 'biocore_htseq_annotate')
    annotation = re.search(r'({}/gff-[0-9a-f]+/rnacentral-9\.0/annotation\.v1\.gff3)'.format(ref), annotate).group(1)
    assert ('htseq-count -f sam -a 0 -s yes -t miRNA -i Name --secondary-alignments score --supplementary-alignments score '
            '-n 1 -o {0}/SR01/annotated.sam -o {0}/SR02/annotated.sam -o {0}/SR03/annotated.sam -o {0}/SR04/annotated.sam '
            '{0}/SR01/alignments.sam.gz'.format(out)) in annotate
    assert f'{out}/SR04/alignments.sam.gz {annotation} > {out}/htseq-count.log' in annotate
    assert 'gzip -f {0}/SR01/annotated.sam {0}/SR02/annotated.sam'.format(out) in annotate
    assert 'benchmark: data/tmp/smallrna/quant/benchmarks/biocore_htseq_annotate.tsv' in annotate
    build = rule_block(output, 'biocore_htseq_annotation')
    assert f'--rnacentral {ref}/sources/rnacentral-9.0/Homo_sapiens.GRCh38.gff.gz' in build
    assert '--fai data/ext/tesmall/mirbase-23.1_ensembl-116/homo_sapiens/genome/genome.fa.fai' in build and '--gtf ' in build
    assert '/rnacentral-9.0/benchmarks/biocore_htseq_annotation.tsv' in build
    download = rule_block(output, 'biocore_htseq_rnacentral')
    assert ('https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/9.0/genome_coordinates/Homo_sapiens.GRCh38.gff.gz'
            in download and 'gzip -t ' in download)
    assert output.count('rule biocore_htseq_ncrna:') == 4
    ncrna = rule_block(output, 'biocore_htseq_ncrna')
    assert 'ncrna --species hsa --annotated ' in ncrna and '/ncrna_types.tsv --features ' in ncrna
    assert f'join --key type --samples SR01 SR02 SR03 SR04 --tables {out}/SR01/ncrna_types.tsv' in rule_block(output, 'biocore_htseq_ncrna_types')


def test_custom_reference_keeps_the_mirna_annotation_and_writes_no_ncrna_tables(tmp_path):
    """A db.mirflowz.custom reference, also for human, annotates with the miRNA GFF alone, as before #503."""
    refs = tmp_path / 'refs'
    refs.mkdir()
    files = {key: refs / name for key, name in (('genome', 'genome.fa.gz'), ('gtf', 'genes.gtf.gz'), ('mirna_gff', 'mirna.gff3'))}
    for path in files.values():
        path.touch()
    custom = lambda s: (s.update(organism='homo_sapiens'),
                        s.setdefault('db', {}).setdefault('mirflowz', {}).update(custom=dict(name='reduced', **{k: str(v) for k, v in files.items()})))
    returncode, output = dry_run(tmp_path, 'biocore_htseq', custom, printshell=True)
    assert returncode == 0, output
    assert output.count('rule biocore_htseq_annotate:') == 4  # per sample
    assert f' - {OUT}/mirna.gff3 > ' in rule_block(output, 'biocore_htseq_annotate')
    assert 'biocore_htseq_ncrna' not in output and 'biocore_htseq_rnacentral' not in output and 'biocore_htseq_annotation' not in output


def test_ncrna_counts_classes_and_features_as_the_pipeline(tmp_path):
    """The pipeline's ncRNA path: records without the species in their feature, a feature, and a perfect alignment; a read
    divided evenly between its distinct types; each (name_sequence, read) pair once; stats by substring."""
    snord = 'ATGATGATGATGATGATGATGA'
    trf = 'TCCCTGGTGGTCTAGTGGTTAGGATTCG'
    path = annotated(tmp_path, [
        record('1-90', 0, snord, '22', 'URS000082BE64-snoRNA'),
        record('1-90', 256, snord, '22', '__ambiguous[URS00004219CE-snoRNA+ENSG00000207445-snoRNA]', chrom='2'),
        record('2-30', 0, trf, '28', '__ambiguous[URS0000A1B2C3-tRNA+ENSG00000276345-Mt_tRNA]'),
        record('2-30', 256, trf, '28', 'ENSG00000123456-lncRNA', chrom='3'),
        record('3-12', 0, MIR21, '22', '__ambiguous[mmu-miR-21a-5p+ENSG00000284190-miRNA]'),     # a miRNA: miRNA path only
        record('3-12', 256, MIR21, '22', 'URS0000D5E6F7-rRNA', chrom='4'),                        # its other locus: ncRNA too
        record('4-8', 0, snord, '20A1', 'URS000082BE64-snoRNA'),                                    # a mismatch: not counted
        record('5-4', 0, 'ACGTACGTACGTACGTAC', '18', '__no_feature'),
        record('6-2', 4, 'GGGGGGGGGGGGGGGGGG', '', '__not_aligned'),
        record('7-6', 16, snord, '22', 'URS000082BE64-snoRNA'),                                     # minus strand
    ])
    classes, features = biocore_htseq_counts.ncrna('mmu', path)
    # read 1 (90): snoRNA only; read 2 (30): tRNA, Mt_tRNA, lncRNA -> 10 each; read 3 (12): rRNA; read 7 (6): snoRNA
    assert classes == {'snoRNA': 96, 'tRNA': 10, 'Mt_tRNA': 10, 'lncRNA': 10, 'rRNA': 12}
    reverse = snord[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    assert features['URS000082BE64-snoRNA_' + snord] == 90 and features['URS000082BE64-snoRNA_' + reverse] == 6
    assert features['URS00004219CE-snoRNA_' + snord] == 90 and features['ENSG00000207445-snoRNA_' + snord] == 90
    assert features['URS0000D5E6F7-rRNA_' + MIR21] == 12 and not any(k.startswith('mmu') for k in features)
    out = {name: tmp_path / '{}.tsv'.format(name) for name in ('types', 'features', 'stats')}
    result = subprocess.run([sys.executable, str(SCRIPT), 'ncrna', '--species', 'mmu', '--annotated', str(path)] +
                            ['--{}={}'.format(k, v) for k, v in out.items()], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert out['types'].read_text() == 'type\tcount\nMt_tRNA\t10\nlncRNA\t10\nrRNA\t12\nsnoRNA\t96\ntRNA\t10\n'
    # by substring, as the pipeline's grep: tRNAs takes Mt_tRNA too
    assert out['stats'].read_text() == ('statistic\tvalue\nRNACentral\t138\ntRNAs\t20\nrRNAs\t12\nsnoRNAs\t96\n'
                                        'snRNAs\t0\nlncRNAs\t10\n')
    assert out['features'].read_text().splitlines()[0] == 'feature\tcount'


def test_ncrna_counts_can_be_fractions(tmp_path):
    path = annotated(tmp_path, [record('1-10', 0, 'ACGTTGCAACGTTGCAACGT', '20', '__ambiguous[URSA1-piRNA+URSB2-snoRNA+ENSG1-lncRNA]')])
    classes, _ = biocore_htseq_counts.ncrna('hsa', path)
    assert classes == pytest.approx({'piRNA': 10 / 3, 'snoRNA': 10 / 3, 'lncRNA': 10 / 3})
    assert biocore_htseq_counts.number(10 / 3) == '3.333333333' and biocore_htseq_counts.number(4.0) == '4'


def test_join_takes_class_tables(tmp_path):
    first, second = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    first.write_text('type\tcount\nsnoRNA\t2.5\ntRNA\t1\n')
    second.write_text('type\tcount\npiRNA\t3\n')
    result = subprocess.run([sys.executable, str(SCRIPT), 'join', '--key', 'type', '--samples', 'S1', 'S2', '--tables', str(first), str(second)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == 'type\tS1\tS2\npiRNA\t0\t3\nsnoRNA\t2.5\t0\ntRNA\t1\t0\n'


ANNOTATION = ROOT / 'gcfdb/scripts/biocore_htseq_annotation.py'


def test_annotation_joins_mirnas_rnacentral_transcripts_and_ensembl_genes(tmp_path):
    """The pipeline's annotation from current sources: RNAcentral transcripts named <URS>-<type> without the taxon suffix,
    Ensembl genes <gene id>-<biotype> (also unnamed ones), every row of type miRNA; other RNAcentral rows, rows without a
    strand and rows off the genome left out; chr prefixes stripped."""
    mirna = tmp_path / 'mirna.gff3'
    mirna.write_text('##gff-version 3\n'
                     '1\t.\tmiRNA\t30\t51\t.\t-\t.\tID=MIMAT1;Name=hsa-miR-1-5p\n'
                     '1\t.\tmiRNA\t100\t122\t.\t.\t.\tID=MIMAT2;Name=hsa-miR-13168-5p\n')
    rnacentral = tmp_path / 'rnacentral.gff3.gz'
    with gzip.open(rnacentral, 'wt') as fh:
        fh.write('##gff-version 3\n'
                 '1\tRNAcentral\ttranscript\t10\t40\t.\t+\t.\tName=URS000082BE64_9606;description=x;type=snoRNA;ID=URS000082BE64_9606.0\n'
                 '1\tRNAcentral\tnoncoding_exon\t10\t40\t.\t+\t.\tName=URS000082BE64_9606;type=snoRNA;Parent=URS000082BE64_9606.0\n'
                 'chr2\tRNAcentral\ttranscript\t5\t30\t.\t-\t.\tName=URS00004219CE_9606;type=piRNA;ID=URS00004219CE_9606.0\n'
                 '1\tRNAcentral\tpredicted_gene\t5\t30\t.\t-\t.\tName=URS0000AAAAAA_9606;type=lncRNA\n'
                 'KI270999.1\tRNAcentral\ttranscript\t5\t30\t.\t+\t.\tName=URS0000BBBBBB_9606;type=rRNA;ID=x\n'
                 '1\tRNAcentral\ttranscript\t200\t900\t.\t+\t.\tName=URS0000CCCCCC_9606;type=circRNA;ID=y\n'
                 # RNAcentral 9.0, the pipeline's release: attributes in quotes, IDs without a taxon suffix
                 '2\tRNAcentral\ttranscript\t600\t690\t.\t+\t.\tID "LN847466.1:1..91:ncRNA";Name "URS000082BE65";type "snoRNA"\n'
                 '2\tRNAcentral\tnoncoding_exon\t600\t690\t.\t+\t.\tID "LN847466.1:1..91:ncRNA_exon1";Name "URS000082BE65";type "snoRNA"\n')
    gtf = tmp_path / 'genes.gtf'
    gtf.write_text('#!genome-build GRCh38\n'
                   '1\thavana\tgene\t1\t500\t.\t+\t.\tgene_id "ENSG00000279493"; gene_version "1"; gene_source "havana"; gene_biotype "artifact";\n'
                   '1\thavana\ttranscript\t1\t500\t.\t+\t.\tgene_id "ENSG00000279493"; transcript_id "ENST1"; gene_biotype "artifact";\n'
                   '2\tensembl\tgene\t50\t90\t.\t-\t.\tgene_id "ENSG00000207445"; gene_name "SNORD1"; gene_biotype "snoRNA";\n')
    fai = tmp_path / 'genome.fa.fai'
    fai.write_text('1\t1000\t3\t60\t61\n2\t1000\t1100\t60\t61\n')
    result = subprocess.run([sys.executable, str(ANNOTATION), '--mirna-gff', str(mirna), '--rnacentral', str(rnacentral),
                             '--gtf', str(gtf), '--fai', str(fai)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    rows = [line.split('\t') for line in result.stdout.splitlines() if not line.startswith('#')]
    assert [(r[0], r[1], r[2], r[6], r[8]) for r in rows] == [
        ('1', '.', 'miRNA', '-', 'ID=MIMAT1;Name=hsa-miR-1-5p'),
        ('1', 'RNAcentral', 'miRNA', '+', 'Name=URS000082BE64-snoRNA'),
        ('2', 'RNAcentral', 'miRNA', '-', 'Name=URS00004219CE-piRNA'),
        ('2', 'RNAcentral', 'miRNA', '+', 'Name=URS000082BE65-snoRNA'),
        ('1', 'Ensembl', 'miRNA', '+', 'Name=ENSG00000279493-artifact'),
        ('2', 'Ensembl', 'miRNA', '-', 'Name=ENSG00000207445-snoRNA')]
    assert '# rows: 1 miRNA, 3 RNAcentral, 2 Ensembl' in result.stdout
    assert 'left out (1): hsa-miR-13168-5p' in result.stdout and '1 RNAcentral on a sequence the genome lacks' in result.stdout
    assert '1 RNAcentral circRNA' in result.stdout  # not in RNAcentral 9.0, the pipeline's release (#503)

