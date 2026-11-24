"""RNA classes and read lengths across routes (#504, smallrna/benchmark/classes.py): each route's parser on small files
in its own format, the common classes, the baseline and unassigned reads, the states and the report."""
import gzip
import struct
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmark'))
import classes  # noqa: E402

MIR21 = 'TAGCTTATCAGACTGATGTTGA'


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def bam(path, reads):
    """A minimal BAM (gzip) with one mapped primary record per (name, length), and an unmapped and a secondary one."""
    records = b''
    for name, length, flag in reads + [('u1', 30, 4), ('s1', 25, 256)]:
        rname = name.encode() + b'\0'
        body = struct.pack('<iiBBHHHiiii', 0, 100, len(rname), 60, 0, 0, flag, length, -1, -1, 0) + rname
        body += b'\0' * ((length + 1) // 2) + b'\xff' * length
        records += struct.pack('<i', len(body)) + body
    header = b'BAM\x01' + struct.pack('<i', 0) + struct.pack('<i', 1) + struct.pack('<i', 5) + b'rDNA\0' + struct.pack('<i', 45000)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wb') as fh:
        fh.write(header + records)
    return path


def test_native_classes_map_to_the_common_set():
    assert classes.common('exact miRNA') == classes.common('isomiR miRNA') == 'miRNA'
    assert classes.common('Mt_tRNA') == classes.common('mature tRNA') == 'tRNA'
    assert classes.common('rRNA (pre-filter)') == classes.common('rRNA_pseudogene') == 'rRNA'
    assert classes.common('scaRNA') == 'snoRNA' and classes.common('Y_RNA') == 'misc ncRNA'
    assert classes.common('ncrna others') == classes.LUMP_MIRGE3 and classes.common('structural_RNA') == classes.LUMP_TESMALL
    assert classes.common('structural_RNA:tRNA') == 'tRNA'  # TEsmall's structural RNA by type
    assert classes.common('transcribed_unprocessed_pseudogene') == 'pseudogene' and classes.common('IG_V_gene') == 'mRNA'
    assert classes.common('something new') == 'unmapped'  # never silently into `other`
    assert classes.common('exon') == classes.LUMP_EXONS  # TEsmall's gene exons include lncRNA and pseudogene exons
    assert classes.common('structural_RNA:7SK') == 'misc ncRNA'


@pytest.mark.parametrize('head, native', [
    ('miR-21-5p', 'miRNA'), ('let-7a-5p', 'miRNA'), ('miR-224-3p(11C->A)', 'miRNA'),
    ('mmu-mir-1983-5p(1G->A)', 'miRNA (other species)'), ('osa-MIR168a-5p', 'miRNA (other species)'),
    ('hsv1-mir-H15-3p(5C->G)', 'miRNA (other species)'), ('ebv-mir-BART19-5p', 'miRNA (other species)'),
    ('tRF-1 from tRNA-Ser-GCT-4-1 (1-17)', 'tRNA'), ('5p-tR-half from tRNA-Lys-CTT-1-2 (1-51)', 'tRNA'),
    ('3p-tR-half from MT-TE-ENSG00000210194.1 (20-70)', 'Mt_tRNA'), ('Mt_rRNA', 'Mt_rRNA'),
    ('piR-cluster:Hsap_0892_11:73927302..73932744', 'piRNA cluster'), ('lncRNA', 'lncRNA')])
def test_unitas_annotation_heads(head, native):
    assert classes.unitas_native(head, {'lncRNA', 'snoRNA', 'rRNA', 'tRNA', 'miRNA', 'no annotation'}) == native


def test_unitas_divides_a_sequence_between_its_annotations(tmp_path):
    folder = tmp_path / 'unitas/unitas/S1'
    write(folder / 'unitas.annotation_summary.txt', 'miRNA\t13\n   miRNA:homo_sapiens\t13\nlncRNA\t1.5\nsnoRNA\t1.5\n')
    write(folder / 'unitas.full_annotation_matrix.txt', 'sequence\treads\tannotation(s)\n'
          f'{MIR21}\t10\tmiR-21-5p\t\n'
          'ACGTACGTACGTACGTACGTAC\t3\tlncRNA|ENST1 (1-22)\tsnoRNA|SNORD1 (5-26)\n'
          'TTTTTTTTTTTTTTTTTTTT\t2\t\n')
    out = classes.unitas(tmp_path, ['S1'])
    assert out.classes['S1'] == {('miRNA', 22): 10, ('lncRNA', 22): 1.5, ('snoRNA', 22): 1.5}
    assert out.assigned['S1'] == {22: 13}  # the sequence without an annotation is unassigned


def mirge3_tables(tmp_path):
    header = 'Sequence,annotFlag,exact miRNA,hairpin miRNA,mature tRNA,primary tRNA,snoRNA,rRNA,ncrna others,mRNA,isomiR miRNA,S1,S2\n'
    write(tmp_path / 'mirge3/mapped.csv', header + f'{MIR21},1,hsa-miR-21-5p,,,,,,,,,10,4\n'
          'GCATTGGTGGTTCAGTGGTAGAATT,1,,,tRNA-Gly,,,,,,,3,0\n'
          'ACGTACGTACGTACGTACGTACGTACGT,1,,,,,,,ENST9,,,0,2\n')
    write(tmp_path / 'mirge3/unmapped.csv', header + 'GGGGGGGGGGGGGGGGGG,0,,,,,,,,,,5,1\n')


def test_mirge3_takes_each_sequences_one_class_and_its_tables_give_the_baseline(tmp_path):
    mirge3_tables(tmp_path)
    out = classes.mirge3(tmp_path, ['S1', 'S2'])
    assert out.classes['S1'] == {('exact miRNA', 22): 10, ('mature tRNA', 25): 3}
    assert out.classes['S2'] == {('exact miRNA', 22): 4, ('ncrna others', 28): 2}
    base, source = classes.baseline(tmp_path, ['S1', 'S2'])
    assert base['S1'] == {22: 10, 25: 3, 18: 5} and base['S2'] == {22: 4, 28: 2, 18: 1} and 'miRge3' in source


def test_tesmall_lengths_the_rrna_filters_reads_and_structural_rna_by_type(tmp_path):
    folder = tmp_path / 'tesmall/S1'
    write(folder / 'S1.anno.rlen.info', 'rlen\tmiRNA\tstructural_RNA\texon\n21\t5\t0\t1\n22\t20\t2\t0\n')
    bam(folder / 'S1.rRNA.bam', [('r1', 30, 0), ('r2', 30, 16), ('r1', 30, 0)])  # r1 aligned twice: once
    write(folder / 'count_summary.txt', 'fid\tftype\tS1\ntRNA:tRNA-Gly-1\tstructural_RNA\t1\nsnoRNA\tstructural_RNA\t1\nhsa-miR-21-5p\tmiRNA\t25\n')
    out = classes.tesmall(tmp_path, ['S1'])
    c = out.classes['S1']
    assert c[('miRNA', 22)] == 20 and c[('structural_RNA', 22)] == 2 and c[('exon', 21)] == 1
    assert c[('rRNA (pre-filter)', 30)] == 2  # the unmapped and the secondary record left out
    assert c[('structural_RNA:tRNA', None)] == 1 and c[('structural_RNA:snoRNA', None)] == 1  # totals only
    assert out.assigned['S1'] == {21: 6, 22: 22, 30: 2}


def test_tinyrna_assigned_reads_by_length(tmp_path):
    write(tmp_path / 'tinyrna/S1/assigned_nt_len_dist.csv', 'Length,A,T,G,C\n21,1.5,2,0,0\n22,0,10,0.5,0\n')
    out = classes.tinyrna(tmp_path, ['S1'])
    assert out.classes['S1'] == {('miRNA', 21): 3.5, ('miRNA', 22): 10.5}


MIRTOP = 'UID\tRead\tmiRNA\tVariant\tiso_5p\tiso_3p\tiso_add3p\tiso_snp\tiso_5p_nt\tiso_3p_nt\tiso_add3p_nt\tiso_snp_nt'


def test_smrnaseq_counts_a_read_once_and_falls_back_to_the_samples_tables(tmp_path):
    rows = f'u1\t{MIR21}\thsa-miR-21-5p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t7\t2\n' \
           f'u1\t{MIR21}\thsa-miR-21-3p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t7\t2\n'   # the same read again: its first row only
    write(tmp_path / 'joined/smrnaseq/joined_samples_mirtop.tsv', MIRTOP + '\tS1\tS2\n' + rows)
    assert classes.smrnaseq(tmp_path / 'joined', ['S1', 'S2']).classes['S1'] == {('miRNA', 22): 7}
    for s, n in (('S1', 7), ('S2', 2)):  # a cleanup removed the joined table: the per-sample tables
        write(tmp_path / f'samples/smrnaseq/mirtop/{s}/counts/mirtop.tsv',
              MIRTOP + f'\t{s}\nu1\t{MIR21}\thsa-miR-21-5p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t{n}\n')
    out = classes.smrnaseq(tmp_path / 'samples', ['S1', 'S2'])
    assert out.classes['S1'] == {('miRNA', 22): 7} and out.classes['S2'] == {('miRNA', 22): 2}


def test_biocore_divides_reads_between_types_and_assigns_a_read_once(tmp_path):
    folder = tmp_path / 'biocore_htseq/S1'
    write(folder / 'isomir.tsv', f'isomiR\tcount\nhsa-miR-21-5p_{MIR21}\t10\nhsa-miR-21b_{MIR21}\t10\n')
    write(folder / 'ncrna.tsv', 'feature\tcount\n'
          'URS000082BE64-snoRNA_ACGTACGTACGTACGTACGTAC\t6\nENSG00000207445-snoRNA_ACGTACGTACGTACGTACGTAC\t6\n'
          'URS0000A1B2C3-tRNA_TCCCTGGTGGTCTAGTGGTTAGGATTCG\t4\nENSG00000123456-lncRNA_TCCCTGGTGGTCTAGTGGTTAGGATTCG\t4\n'
          f'URS0000D5E6F7-rRNA_{MIR21}\t10\n')  # a read on a miRNA at one locus and on rRNA at another: in both
    out = classes.biocore_htseq(tmp_path, ['S1'])
    assert out.classes['S1'] == {('miRBase miRNA', 22): 10, ('snoRNA', 22): 6, ('tRNA', 28): 2, ('lncRNA', 28): 2, ('rRNA', 22): 10}
    assert out.assigned['S1'] == {22: 16, 28: 4}  # the union: the miRNA read assigned once


def test_summary_unassigned_states_and_a_route_that_cannot_be_read(tmp_path):
    mirge3_tables(tmp_path)
    write(tmp_path / 'tinyrna/S1/assigned_nt_len_dist.csv', 'Length,A,T,G,C\n22,0,10,0,0\n')
    write(tmp_path / 'tinyrna/S2/assigned_nt_len_dist.csv', 'Length,A,T,G,C\n22,0,4,0,0\n')
    write(tmp_path / 'biocore_htseq/S1/isomir.tsv', f'isomiR\tcount\nhsa-miR-21-5p_{MIR21}\t10\n')
    write(tmp_path / 'biocore_htseq/S2/isomir.tsv', f'isomiR\tcount\nhsa-miR-21-5p_{MIR21}\t4\n')
    rows, states, counts, base, notes = classes.summarize(tmp_path, ['mirge3', 'tinyrna', 'biocore_htseq', 'unitas'], ['S1', 'S2'])
    assert 'unitas' not in counts and any(n.startswith('unitas: not read') for n in notes)  # the evaluation goes on
    unassigned = {(m, s, length): r for m, s, cls, _, length, r in rows if cls == 'unassigned'}
    assert unassigned[('tinyrna', 'S1', 25)] == 3 and unassigned[('tinyrna', 'S1', 18)] == 5 and ('tinyrna', 'S1', 22) not in unassigned
    state = {(m, c): s for m, c, s in states}
    assert state[('mirge3', 'snRNA')] == 'lumped in ' + classes.LUMP_MIRGE3 and state[('mirge3', 'piRNA')] == 'not reported'
    assert state[('biocore_htseq', 'tRNA')] == 'not reported'  # no ncRNA tables (#503: human only): miRNA only
    assert state[('tinyrna', 'unassigned')] == 'reported'
    lines = classes.report_lines(rows, states, base, ['mirge3', 'tinyrna', 'biocore_htseq'], ['S1', 'S2'], notes)
    table = '\n'.join(lines)
    assert '| miRNA | 56.35 | 56.35 | 56.35 |' in table  # 10 of 18 and 4 of 7 filtered reads: the median of 55.56 and 57.14
    assert '| snRNA | lumped | n/r | n/r |' in table and '| tRNA | 8.33 | n/r | n/r |' in table  # 3 of 18 and 0 of 7


def test_a_broken_mirge3_table_falls_back_to_the_filtered_reads(tmp_path):
    mirge3_tables(tmp_path)
    write(tmp_path / 'mirge3/unmapped.csv', '')  # empty: the pair is incomplete
    fastq = {s: write(tmp_path / f'{s}.fastq', '@r1\nACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIII\n') for s in ('S1', 'S2')}
    rows, states, counts, base, notes = classes.summarize(tmp_path, ['mirge3'], ['S1', 'S2'], fastq)
    assert base['S1'] == {18: 1} and 'mirge3' in counts  # mapped.csv still gives miRge3's classes
    rows, states, counts, base, notes = classes.summarize(tmp_path, ['mirge3'], ['S1', 'S2'])  # no FASTQ either
    assert base is None and 'mirge3' in counts


def test_tesmall_without_its_rrna_bam_gives_neither_rrna_nor_unassigned(tmp_path):
    mirge3_tables(tmp_path)
    for s in ('S1', 'S2'):
        write(tmp_path / f'tesmall/{s}/{s}.anno.rlen.info', 'rlen\tmiRNA\texon\n22\t5\t1\n')
    rows, states, counts, base, notes = classes.summarize(tmp_path, ['tesmall'], ['S1', 'S2'])
    assert not [r for r in rows if r[0] == 'tesmall' and r[2] == 'unassigned']
    state = {(m, c): s for m, c, s in states}
    assert state[('tesmall', 'rRNA')].startswith('not available') and state[('tesmall', 'unassigned')].startswith('not available')
    assert state[('tesmall', 'mRNA')] == 'lumped in ' + classes.LUMP_EXONS
    table = '\n'.join(classes.report_lines(rows, states, base, ['tesmall'], ['S1', 'S2'], notes))
    assert '| rRNA | n/a |' in table and '| unassigned | n/a |' in table


def test_both_plots_are_drawn_with_totals_only_rows_and_routes_without_a_class(tmp_path):
    pytest.importorskip('matplotlib')
    rows = [('mirge3', 'S1', 'miRNA', 'exact miRNA', 22, 10.0), ('mirge3', 'S1', classes.LUMP_MIRGE3, 'ncrna others', 30, 2.0),
            ('mirge3', 'S1', 'unassigned', '', 18, 3.0), ('tesmall', 'S1', 'miRNA', 'miRNA', 22, 9.0),
            ('tesmall', 'S1', 'tRNA', 'structural_RNA:tRNA', '', 1.0), ('tesmall', 'S1', 'unassigned', '', 50, 6.0)]
    base = {'S1': {22: 10, 30: 2, 18: 3}}
    assert classes.plot(rows, base, ['mirge3', 'tesmall'], tmp_path / 'classes.png', tmp_path / 'class_lengths.png')
    assert (tmp_path / 'classes.png').stat().st_size > 1000 and (tmp_path / 'class_lengths.png').stat().st_size > 1000


def test_decoy_samples_score_where_a_route_puts_one_class_without_the_carriers(tmp_path):
    """P_control holds the carriers; a decoy sample adds reads of one class. A route with sequences is read without the
    carriers, so a carrier it places differently in the decoy sample does not count; tiny-count, without sequences,
    is scored by its classes minus P_control's, and leaves every decoy read unassigned."""
    carrier, trna = MIR21, 'GCATTGGTGGTTCAGTGGTAGAATTCTCGC'
    summary_text = 'miRNA\t10\nlncRNA\t1\ntRNA\t4\n'
    for sample, matrix in (('P_control', f'{carrier}\t10\tmiR-21-5p\t\n'),
                           ('P_decoy_tRNA-gtrnadb_1', f'{carrier}\t10\tlncRNA|ENST9 (1-22)\t\n{trna}\t4\ttRNA|tRNA-Gly (1-30)\t\n')):
        folder = tmp_path / 'unitas/unitas' / sample
        write(folder / 'unitas.annotation_summary.txt', summary_text)
        write(folder / 'unitas.full_annotation_matrix.txt', 'sequence\treads\tannotation(s)\n' + matrix)
        write(tmp_path / 'tinyrna' / sample / 'assigned_nt_len_dist.csv', 'Length,A,T,G,C\n22,0,10,0,0\n')
    samples = ['P_control', 'P_decoy_tRNA-gtrnadb_1']
    base = {'P_control': {22: 10}, 'P_decoy_tRNA-gtrnadb_1': {22: 10, 30: 5}}  # one decoy read no route annotates
    rows, states, counts, _, _ = classes.summarize(tmp_path, ['unitas', 'tinyrna'], samples)
    detail, summary = classes.decoy_scores(tmp_path, rows, states, base, counts, ['unitas', 'tinyrna'], samples, {carrier})
    score = {(m, o): v for m, s, o, v in summary if s == 'tRNA-gtrnadb'}
    assert score[('unitas', 'right class')] == 80.0 and score[('unitas', 'unassigned')] == 20.0
    assert score[('unitas', 'other class')] == 0 and score[('unitas', 'miRNA')] == 0  # the carrier's lncRNA left out
    assert score[('tinyrna', 'right class')] == 'n/r' and score[('tinyrna', 'unassigned')] == 100.0
    assert ('unitas', 'tRNA-gtrnadb', 'tRNA', 4.0, 80.0) in detail
    table = '\n'.join(classes.decoy_lines(summary, ['unitas', 'tinyrna']))
    assert '| tRNA-gtrnadb (tRNA) | 80.0 | n/r |' in table
    assert classes.decoy_scores(tmp_path, rows, states, base, counts, ['unitas'], ['P_control'], {carrier}) == ([], [])


def test_lumped_decoy_reads_are_shown_beside_the_right_class():
    summary = [('mirge3', 'Y_RNA-ensembl', o, v) for o, v in
               (('right class', 0.0), ('lumped', 100.0), ('miRNA', 0.0), ('other class', 0.0), ('unassigned', 0.0))]
    assert '| Y_RNA-ensembl (misc ncRNA) | 0.0 (+100.0 lumped) |' in '\n'.join(classes.decoy_lines(summary, ['mirge3']))


def test_a_decoy_source_is_also_scored_by_part_for_routes_with_sequences(tmp_path):
    """RefSeq's rRNA decoys mix mature rRNA and pre-rRNA spacers; each part is scored without the other's sequences."""
    carrier, mature, spacer = MIR21, 'TACCTGGTTGATCCTGCCAGTAGCAT', 'CGCGCCCCGCTTCTTCGGTCGCG'
    for sample, matrix in (('P_control', f'{carrier}\t10\tmiR-21-5p\t\n'),
                           ('P_decoy_rRNA-refseq_1', f'{carrier}\t10\tmiR-21-5p\t\n{mature}\t3\trRNA|18S (1-26)\t\n')):
        folder = tmp_path / 'unitas/unitas' / sample
        write(folder / 'unitas.annotation_summary.txt', 'miRNA\t10\nrRNA\t3\n')
        write(folder / 'unitas.full_annotation_matrix.txt', 'sequence\treads\tannotation(s)\n' + matrix)
    samples = ['P_control', 'P_decoy_rRNA-refseq_1']
    base = {'P_control': {22: 10}, 'P_decoy_rRNA-refseq_1': {22: 10, 26: 3, 23: 2}}
    parts = {'rRNA-refseq': {'18S': {'P_decoy_rRNA-refseq_1': {mature: 3}},
                             'spacers': {'P_decoy_rRNA-refseq_1': {spacer: 2}}}}
    rows, states, counts, _, _ = classes.summarize(tmp_path, ['unitas'], samples)
    _, summary = classes.decoy_scores(tmp_path, rows, states, base, counts, ['unitas'], samples, {carrier}, parts)
    score = {(s, o): v for _, s, o, v in summary}
    assert score[('rRNA-refseq', 'right class')] == 60.0
    assert score[('rRNA-refseq 18S', 'right class')] == 100.0 and score[('rRNA-refseq spacers', 'unassigned')] == 100.0
    table = '\n'.join(classes.decoy_lines(summary, ['unitas']))
    assert table.index('| rRNA-refseq (rRNA) | 60.0 |') < table.index('| rRNA-refseq 18S (rRNA) | 100.0 |')
    assert [classes.decoy_part('rRNA-refseq', o) for o in ('ITS1:fragment', '5.8S:fragment', '28S:x', 'other')] == [
        'spacers', '5.8S', '28S', None]
    assert [classes.decoy_part('rRNA-ensembl', o) for o in ('RNA5SP12:fragment', '5_8S_rRNA:fragment', 'MT-RNR2:fragment')] == [
        '5S', '5.8S', 'mitochondrial']


def test_tesmall_without_its_rrna_bam_is_not_scored_on_the_decoys(tmp_path):
    for sample in ('P_control', 'P_decoy_rRNA-refseq_1'):
        write(tmp_path / f'tesmall/{sample}/{sample}.anno.rlen.info', 'rlen\tmiRNA\n22\t10\n')
    samples = ['P_control', 'P_decoy_rRNA-refseq_1']
    base = {'P_control': {22: 10}, 'P_decoy_rRNA-refseq_1': {22: 10, 30: 5}}
    rows, states, counts, _, _ = classes.summarize(tmp_path, ['tesmall'], samples)
    _, summary = classes.decoy_scores(tmp_path, rows, states, base, counts, ['tesmall'], samples, set())
    assert {v for _, _, _, v in summary} == {'n/a'}  # not 100% unassigned
    table = '\n'.join(classes.decoy_lines(summary, ['tesmall']))
    assert '| rRNA-refseq (rRNA) | n/a |' in table and 'Decoy reads put in another class, %:' in table
