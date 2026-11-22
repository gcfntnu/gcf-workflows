"""Small-RNA benchmark: synthetic truth, method adapters, metrics, pools (truth groups, oracle, fill) and the prepare/evaluate/summary steps
on synthetic tables."""
import collections
import gzip
import io
import json
import math
import os
import statistics
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'smallrna/benchmark'
sys.path.insert(0, str(BENCH))
import adapters  # noqa: E402
import benchmark  # noqa: E402
import decoys  # noqa: E402
import fill  # noqa: E402
import metrics  # noqa: E402
import pool  # noqa: E402
import resources  # noqa: E402
import synthetic  # noqa: E402

MATURE = ('>mmu-miR-192-5p MIMAT0000517 Mus musculus miR-192-5p\nCUGACCUAUGAAUUGACAGCCA\n'
          '>mmu-miR-194-5p MIMAT0000224 Mus musculus miR-194-5p\nUGUAACAGCAACUCCAUGUGGA\n'
          '>mmu-miR-204-5p MIMAT0000237 Mus musculus miR-204-5p\nUUCCCUUUGUCAUCCUAUGCCU\n'
          '>hsa-miR-192-5p MIMAT0000222 Homo sapiens miR-192-5p\nCUGACCUAUGAAUUGACAGCCA\n')
GFF = ('##gff-version 3\n# genome-build-id: GRCm39\n'
       '19\t.\tmiRNA_primary_transcript\t1\t90\t.\t+\t.\tID=MI0000234;Alias=MI0000234;Name=mmu-mir-192\n'
       '19\t.\tmiRNA\t10\t31\t.\t+\t.\tID=MIMAT0000517;Alias=MIMAT0000517;Name=mmu-miR-192-5p;Derives_from=MI0000234\n'
       '19\t.\tmiRNA\t50\t71\t.\t+\t.\tID=MIMAT0000224;Alias=MIMAT0000224;Name=mmu-miR-194-5p;Derives_from=MI0000235\n')
NAME_TO_ID = {'mmu-miR-192-5p': 'MIMAT0000517', 'mmu-miR-194-5p': 'MIMAT0000224', 'mmu-miR-204-5p': 'MIMAT0000237'}


@pytest.fixture
def ref(tmp_path):
    ref = tmp_path / 'ext/benchmark/references/mmu-test'
    ref.mkdir(parents=True)
    (ref / 'mature.fa').write_text(MATURE)
    (ref / 'mirna.gff3').write_text(GFF)
    return ref


def test_mature_fasta_keeps_the_species_and_converts_to_dna(ref):
    mature = synthetic.read_mature(ref / 'mature.fa', 'mmu-')
    assert set(mature) == {'mmu-miR-192-5p', 'mmu-miR-194-5p', 'mmu-miR-204-5p'}
    assert mature['mmu-miR-192-5p'] == ('MIMAT0000517', 'CTGACCTATGAATTGACAGCCA')


def test_truth_is_deterministic_and_covers_every_read_class(ref):
    spec = {'seed': 3, 'mirnas': {'mmu-miR-192-5p': 4}, 'random_mirnas': {'n': 1, 'count': [2, 2]},
            'isomirs': [{'mirna': 'mmu-miR-192-5p', 'change': '3p_add_A', 'count': 2},
                        {'mirna': 'mmu-miR-192-5p', 'change': '5p_trim_1', 'count': 1}],
            'random_reads': {'n': 2, 'length': 20}}
    mature = synthetic.read_mature(ref / 'mature.fa', 'mmu-')
    rows = synthetic.truth_rows(spec, mature, synthetic.gff_mature_names(ref / 'mirna.gff3'))
    assert rows == synthetic.truth_rows(spec, mature, synthetic.gff_mature_names(ref / 'mirna.gff3'))
    assert rows[0] == ('canonical', 'mmu-miR-192-5p', 'MIMAT0000517', 'CTGACCTATGAATTGACAGCCA', 4)
    assert rows[1][1] == 'mmu-miR-194-5p'  # the only annotated miRNA left to draw; miR-204 is not on the reference
    assert rows[2][3] == 'CTGACCTATGAATTGACAGCCAA' and rows[3][3] == 'TGACCTATGAATTGACAGCCA'
    assert [r[0] for r in rows[4:]] == ['random', 'random'] and all(len(r[3]) == 20 for r in rows[4:])


def test_unknown_mirna_and_repeated_sequences_fail(ref):
    mature = synthetic.read_mature(ref / 'mature.fa', 'mmu-')
    with pytest.raises(ValueError, match='not in the mature miRNA FASTA'):
        synthetic.truth_rows({'mirnas': {'mmu-miR-1-3p': 2}}, mature, set())
    with pytest.raises(ValueError, match='repeat a sequence'):
        synthetic.truth_rows({'mirnas': {'mmu-miR-192-5p': 2}, 'isomirs': [{'mirna': 'mmu-miR-192-5p', 'change': '3p_trim_1', 'count': 1},
                                                                         {'mirna': 'mmu-miR-192-5p', 'change': '3p_trim_1', 'count': 1}]},
                             mature, set())


def test_reads_carry_umis_and_the_adapter(ref, tmp_path):
    rows = [('canonical', 'mmu-miR-192-5p', 'MIMAT0000517', 'CTGACC', 2)]
    import random
    path = tmp_path / 'S01.fastq.gz'
    assert synthetic.write_sample(path, 'S01', rows, 'TGGAATT', 4, random.Random(1)) == 2
    lines = gzip.open(path, 'rt').read().splitlines()
    assert len(lines[1]) == 4 + 6 + 4 + 7 and lines[1][4:10] == 'CTGACC' and lines[1].endswith('TGGAATT')


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def unitas_dump_dir(ext, organism):
    """Where the project's db.unitas (main.config) puts the organism's UNITAS reference dump."""
    release = yaml.safe_load((BENCH.parents[1] / 'main.config').read_text())['db']['unitas']['release']
    return Path(ext) / 'unitas' / release / 'UNITAS_refdump_{}'.format(organism)


def unitas_dump(ext, organism, species, hairpins):
    """A UNITAS reference dump with only the organism's miRBase hairpins, at unitas_dump_dir."""
    return write(unitas_dump_dir(ext, organism) / 'miRBase_hairpin_{}.fas'.format(species),
                 ''.join('>{}_stem-loop\n{}\n'.format(n, s) for n, s in hairpins.items())).parent


# mmu-mir-194-2: miR-194-5p on its 5' arm, nothing on its 3' arm
MMU_HAIRPINS = {'miR-194-2': 'GGTGTGTAACAGCAACTCCATGTGGACTGTGCTCGGATTCCAGTGGGGCTGCTGTTATCTGGGGCAACC'}


def write_records(quant, methods, database='mirbase', release='23.1', dump=None):
    """Each method's quant/<method>/reference.tsv, as unitas_reference, mirge3_reference and quant_reference write it
    (quant_reference with the SHA-256 of the count table as it is now; UNITAS's with the dump it read, when given)."""
    for method in methods:
        text = 'field\tvalue\tsource\nmirna_database\t{}\tx\nmirna_release\t{}\tx\n'.format(database, release)
        if method == 'unitas' and dump is not None:
            text += 'dump\t{}\tconfiguration\n'.format(dump)
        if method in adapters.TABLE_BOUND_RECORDS:
            text += 'table_sha256\t{}\tx\n'.format(adapters.table_sha256(method, quant))
        write(quant / method / 'reference.tsv', text)


UNITAS_MIR_TABLE = "miR-name\tsequence\ttotal_reads\t+N0-3'\t+N1-3'\t+N2-3'\tA-tailed_reads\tT-tailed_reads\n"


def isomir_outputs(quant, sequences=None):
    """isomiR outputs that agree with the UNITAS and miRge3 count tables written so far, as the routes write them
    (#454): each miRNA's reads on one read sequence, `sequences` ({name without the species prefix: sequence}) or one
    made from the name, all exact reads for miRge3."""
    def sequence(name):
        bare = name.strip().split('-', 1)[1] if name.strip()[:4] in ('hsa-', 'mmu-') else name.strip()
        return (sequences or {}).get(bare) or ''.join('ACGT'[ord(c) % 4] for c in (bare * 3))[:22]
    unitas = quant / 'unitas/mir_counts.tsv'
    if unitas.is_file():
        lines = [line.split('\t') for line in unitas.read_text().splitlines() if line]
        for i, sample in enumerate(lines[0][1:], 1):
            write(quant / 'unitas/unitas' / sample / 'unitas.miR-table_Human.txt', UNITAS_MIR_TABLE + ''.join(
                '{}\t{}\t{c}\t{c}\t0\t0\t0\t0\n'.format(row[0], sequence(row[0]), c=row[i]) for row in lines[1:] if float(row[i])))
    mirge3 = quant / 'mirge3/miR.Counts.csv'
    if mirge3.is_file():
        lines = mirge3.read_text().splitlines()
        annotation = 'Sequence,annotFlag,exact miRNA,hairpin miRNA,mature tRNA,primary tRNA,snoRNA,rRNA,ncrna others,mRNA,isomiR miRNA'
        write(quant / 'mirge3/mapped.csv', annotation + ',' + lines[0].split(',', 1)[1] + '\n' + ''.join(
            '{},1,{},,,,,,,,,{}\n'.format(sequence(line.split(',')[0]), line.split(',')[0], line.split(',', 1)[1]) for line in lines[1:] if line))
        write(quant / 'mirge3/lib/human/annotation.Libs/human_merges_miRBase.csv', '')


def native_outputs(quant):
    write(quant / 'unitas/mir_counts.tsv', '\tS01\tS02\nmiR-192-5p\t45\t45\nmiR-194-5p \t20\t25\nmiR-194-2-5p\t5\t0\nmiR-194-2-3p\t1\t0\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S01,S02\nmmu-miR-192-5p,45,45\nmmu-miR-194-5p/mmu-miR-204-5p,25,25\nmmu-miR-9999-5p,1,0\n')
    write(quant / 'tesmall/tesmall_counts.tsv', 'ftype\tfid\tS01\tS02\nhairpin\tmmu-mir-192\t3\t3\n'
          'miRNA\tmmu-miR-192-5p\t42\t44\nmiRNA\tmmu-miR-194-5p\t25\t0\nsense_TE\tL1:L1:LINE\t9\t9\n')
    # tiny-count: one row per locus (MIMAT0000224_1 is the second locus of miR-194-5p), fractional shares of multimappers
    write(quant / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS01\tS02\n'
          'MIMAT0000517\tmiRNA\tmmu-miR-192-5p\t40.0\t45.0\nMIMAT0000224\tmiRNA\tmmu-miR-194-5p\t12.5\t0.0\n'
          'MIMAT0000224_1\tmiRNA\tmmu-miR-194-5p\t12.5\t0.0\nMIMAT0000237\tmiRNA\tmmu-miR-204-5p\t0.0\t3.0\n')
    write(quant / 'smrnaseq/mirna.tsv', 'miRNA\tS02\tS01\nmmu-miR-192-5p\t44\t45\nmmu-miR-194-5p\t0\t25\nmmu-miR-9999-5p\t2\t0\n')
    # biocore_htseq: a read counts in full for every name it overlaps, so identical miR-194-5p and miR-204-5p reads count twice
    write(quant / 'biocore_htseq/htseq_counts.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\t45\nmmu-miR-194-5p\t25\t0\nmmu-miR-204-5p\t25\t0\n')
    # biocore_miraligner: a sequence annotated to two miRNAs is divided between them
    write(quant / 'biocore_miraligner/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\t44.5\nmmu-miR-194-5p\t12.5\t0\nmmu-miR-204-5p\t12.5\t0\n')
    write(quant / 'mirflowz/all_mirna_counts.tab',
          'ID\tS01\tS02\nmmu-miR-192-5p\t40\t40\nmmu-miR-192-5p|0|1|23M|23|CTGACCTATGAATTGACAGCCAA\t5\t5\nmmu-miR-194-2-5p\t25\t25\n'
          'mmu-miR-204-5p|0|0|22M|22|TTCCCTTTGTCATCCTATGCCT\tNA\t3\n')  # NA: no reads in that sample
    return write(quant.parent / 'mirflowz.gff3',
                 '19\t.\tmiRNA\t1\t30\t.\t+\t.\tID=MIMAT0000224_1;Alias=MIMAT0000224;Name=mmu-miR-194-2-5p;Derives_from=MI0000235\n')


def test_adapters_match_by_accession_whatever_the_method_names(tmp_path):
    gff = native_outputs(tmp_path / 'quant')
    rows = adapters.common_table(['unitas', 'mirge3', 'mirflowz', 'tesmall', 'tinyrna', 'smrnaseq', 'biocore_htseq', 'biocore_miraligner'], tmp_path / 'quant',
                                 NAME_TO_ID, 'mmu-', gff)
    got = {(m, s, n): (i, c) for m, s, n, i, c in rows}
    assert got[('unitas', 'S01', 'miR-194-5p')] == ('MIMAT0000224', 20)
    assert got[('unitas', 'S01', 'miR-194-2-5p')] == ('', 5)       # a hairpin name: placed only with the UNITAS hairpins
    assert got[('mirge3', 'S01', 'mmu-miR-194-5p/mmu-miR-204-5p')] == ('MIMAT0000224', 25)
    assert got[('mirge3', 'S01', 'mmu-miR-9999-5p')] == ('', 1)
    assert got[('mirflowz', 'S01', 'mmu-miR-192-5p')] == ('MIMAT0000517', 45)  # canonical and isomiR rows together
    assert got[('mirflowz', 'S01', 'mmu-miR-194-2-5p')] == ('MIMAT0000224', 25)
    assert got[('mirflowz', 'S02', 'mmu-miR-204-5p')] == ('MIMAT0000237', 3) and ('mirflowz', 'S01', 'mmu-miR-204-5p') not in got
    assert ('mirge3', 'S02', 'mmu-miR-9999-5p') not in got  # zero counts are dropped
    assert got[('tesmall', 'S01', 'mmu-miR-192-5p')] == ('MIMAT0000517', 42)  # mature miRNAs only, not hairpin or TE
    assert ('tesmall', 'S02', 'mmu-miR-194-5p') not in got and not any(n.startswith(('mmu-mir', 'L1')) for m, _, n in got if m == 'tesmall')
    assert got[('tinyrna', 'S01', 'mmu-miR-194-5p')] == ('MIMAT0000224', 25)  # loci summed by accession
    assert got[('tinyrna', 'S02', 'mmu-miR-192-5p')] == ('MIMAT0000517', 45) and got[('tinyrna', 'S02', 'mmu-miR-204-5p')] == ('MIMAT0000237', 3)
    assert ('tinyrna', 'S01', 'mmu-miR-204-5p') not in got
    assert got[('smrnaseq', 'S01', 'mmu-miR-192-5p')] == ('MIMAT0000517', 45)  # samples by column name, whatever their order
    assert got[('smrnaseq', 'S02', 'mmu-miR-9999-5p')] == ('', 2) and ('smrnaseq', 'S02', 'mmu-miR-194-5p') not in got
    assert got[('biocore_htseq', 'S01', 'mmu-miR-194-5p')] == ('MIMAT0000224', 25) and got[('biocore_htseq', 'S01', 'mmu-miR-204-5p')] == ('MIMAT0000237', 25)
    assert ('biocore_htseq', 'S02', 'mmu-miR-194-5p') not in got
    assert got[('biocore_miraligner', 'S02', 'mmu-miR-192-5p')] == ('MIMAT0000517', 44.5)
    assert got[('biocore_miraligner', 'S01', 'mmu-miR-204-5p')] == ('MIMAT0000237', 12.5)
    write(tmp_path / 'quant/smrnaseq/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\n')  # a truncated row
    with pytest.raises(ValueError, match='2 fields, expected 3'):
        adapters.smrnaseq(tmp_path / 'quant', NAME_TO_ID)
    write(tmp_path / 'quant/smrnaseq/mirna.tsv', 'miRNA\tS01\tS01\nmmu-miR-192-5p\t45\t44\n')
    with pytest.raises(ValueError, match='unique sample columns'):
        adapters.smrnaseq(tmp_path / 'quant', NAME_TO_ID)
    write(tmp_path / 'quant/biocore_htseq/htseq_counts.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\n')
    with pytest.raises(ValueError, match='htseq_counts.tsv:2: 2 fields, expected 3'):
        adapters.biocore_htseq(tmp_path / 'quant', NAME_TO_ID)
    write(tmp_path / 'quant/biocore_miraligner/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\n')
    with pytest.raises(ValueError, match='mirna.tsv:2: 2 fields, expected 3'):
        adapters.biocore_miraligner(tmp_path / 'quant', NAME_TO_ID)
    with pytest.raises(ValueError, match='no benchmark adapter'):
        adapters.common_table(['newtool'], tmp_path / 'quant', NAME_TO_ID, 'mmu-')


def test_unitas_hairpin_names_get_the_mature_miRNA_on_that_arm(tmp_path):
    mature = synthetic.read_mature(write(tmp_path / 'mature.fa', (
        '>hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p\nUAGCAGCACGUAAAUAUUGGCG\n'
        '>hsa-miR-16-1-3p MIMAT0004489 Homo sapiens miR-16-1-3p\nCCAGUAUUAACUGUGCUGCUGA\n'
        '>hsa-miR-325-5p MIMAT0000771 Homo sapiens miR-325-5p\nCCUAGUAGGUGUCCAGUAAGUGU\n'
        '>hsa-miR-y-5p MIMAT9999999 Homo sapiens miR-y-5p\nCCUAGCAGCACGUAAAUAUUGG\n')), 'hsa-')
    name_to_id = {n: a for n, (a, _) in mature.items()}
    hairpins = {'miR-16-1': 'GTCAGCAGTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAAGGTTGAC',  # hsa-mir-16-1
                'miR-325': 'ATACAGTGCTTGGTTCCTAGTAGGTGTCCAGTAAGTGTTTGTGACATAATTTGTTTATTGAGGACCTCCTATCAATCAAGCACTGTGCTAGGCTCTGG',
                'miR-x': 'TAGCAGCACGTAAATATTGGCGCCTAGCAGCACGTAAATATTGGCA' + 'A' * 50,
                'miR-y-2': 'TAGCAGCACGTAAATATTGGCGCCTAGCAGCACGTAAATATTGGCA' + 'A' * 50}
    place = lambda name: adapters.unitas_place(name, name_to_id, 'hsa-', hairpins, mature)
    assert place('miR-16-5p') == ('MIMAT0000069', 'name')
    assert place('miR-16-1-5p') == ('MIMAT0000069', 'hairpin')       # hsa-miR-16-5p on the 5' arm
    assert place('miR-16-1-3p') == ('MIMAT0004489', 'name')          # its own mature miRNA, not the 5' one
    assert place('miR-325-3p') == ('', '')                           # human miRBase has only miR-325-5p
    assert place('miR-x-5p') == ('', '')                             # two mature miRNAs of other families on the arm: ambiguous
    assert place('miR-y-2-5p') == ('MIMAT9999999', 'hairpin')        # the same two, but one is the hairpin family's (miR-y)
    assert place('miR-9999-1-5p') == ('', '')                        # no such hairpin
    dump = unitas_dump(tmp_path / 'ext', 'homo_sapiens', 'Human', hairpins)
    write(dump / 'miRBase_hairpin_non-Human.fas', '>cfa-mir-1_stem-loop\nACGT\n')
    assert adapters.unitas_hairpins(dump) == hairpins


def test_count_metrics():
    truth = {('S01', 'A'): 10, ('S01', 'B'): 5}
    common = [('m', 'S01', 'a', 'A', 10), ('m', 'S01', 'b', 'B', 3), ('m', 'S01', 'c', 'C', 2), ('m', 'S01', 'x', '', 4)]
    result = metrics.count_metrics(common, truth)['m']
    assert result == {'recovered': round(13 / 15, 4), 'exact': 0.5, 'misassigned': 2, 'unplaced': 4}


def test_exact_allows_floating_point_noise_from_divided_reads(tmp_path):
    """One read divided over ten loci sums to 0.9999999999999999 with plain addition; still exact."""
    assert sum([0.1] * 10) != 1
    common = [('m', 'S01', 'a', 'A', sum([0.1] * 10)), ('m', 'S01', 'b', 'B', 2.5)]
    assert metrics.count_metrics(common, {('S01', 'A'): 1, ('S01', 'B'): 3})['m']['exact'] == 0.5
    rows = ['MIMAT0000001{}\tmiRNA\tmmu-miR-1a-3p\t0.1'.format('_{}'.format(i) if i else '') for i in range(10)]
    write(tmp_path / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS01\n' + '\n'.join(rows) + '\n')
    assert adapters.tinyrna(tmp_path, {}) == [('tinyrna', 'S01', 'mmu-miR-1a-3p', 'MIMAT0000001', 1.0)]  # fsum


def test_prepare_then_evaluate(tmp_path, ref):
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({
        'name': 'syn', 'kind': 'synthetic', 'organism': 'mus_musculus', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
        'adapter': 'TGGAATTCTCGGGTGCCAAGG', 'methods': ['unitas', 'mirge3', 'mirflowz'],
        'reference': {'name': 'mmu-test', 'dir': 'benchmark/references/mmu-test'},
        'synthetic': {'seed': 1, 'samples': 2, 'mirnas': {'mmu-miR-192-5p': 40, 'mmu-miR-194-5p': 25},
                      'isomirs': [{'mirna': 'mmu-miR-192-5p', 'change': '3p_add_A', 'count': 5}], 'random_reads': {'n': 2}}}))
    project = tmp_path / 'project'
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(tmp_path / 'ext')],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert config['quant']['method'] == 'unitas,mirge3,mirflowz' and sorted(config['samples']) == ['S01', 'S02']
    assert config['db']['mirflowz']['custom']['genome'] == str(ref / 'genome.fa.gz')
    assert (project / 'data/raw/fastq/syn/S02_R1.fastq.gz').is_file()
    assert (project / 'pep/sample_table.csv').read_text().splitlines()[1] == 'S01,syn/S01_R1.fastq.gz,0'

    gff = native_outputs(project / 'data/tmp/smallrna/quant')
    write_records(project / 'data/tmp/smallrna/quant', ['unitas'], dump=unitas_dump_dir(tmp_path / 'ext', 'mus_musculus'))
    result = subprocess.run(cmd + ['evaluate', '--project', str(project), '--mirflowz-gff', str(gff)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'no UNITAS reference dump at' in result.stderr
    unitas_dump(tmp_path / 'ext', 'mus_musculus', 'Mouse', MMU_HAIRPINS)
    result = subprocess.run(cmd + ['evaluate', '--project', str(project), '--mirflowz-gff', str(gff)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'UNITAS names: 2 by name (135 reads), 1 by hairpin and arm (5 reads), 1 unplaced (1 reads)' in result.stdout
    assert (project / 'benchmark/unitas_names.tsv').read_text().splitlines()[1:3] == [
        'miR-192-5p\tname\tMIMAT0000517\t90', 'miR-194-2-3p\tunplaced\t\t1']
    report = (project / 'benchmark/report.md').read_text()
    assert '| recovered | 1.0 | 1.0 | 1.0 |' in report  # unitas's 5 reads of miR-194-5p in S01 on the hairpin name: 140/140
    assert '| mmu-miR-192-5p | 45 | 45 | 45 | 45 |' in report
    assert 'mirge3: `mmu-miR-9999-5p`' in report
    assert (project / 'benchmark/common.tsv').read_text().startswith('method\tsample\tmature_name\tmature_id\tcount\n')


@pytest.mark.parametrize('method', ['tesmall', 'tinyrna', 'biocore_htseq'])
def test_tesmall_or_tinyrna_only_project_keeps_the_reduced_reference(tmp_path, ref, method):
    """TEsmall, tiny-count and biocore_htseq use the mirflowz reference, so a project with any of them alone still gets the
    reduced reference, and its dry run builds the TEsmall genome (shared by all three) from it."""
    for name in ('genome.fa.gz', 'genes.gtf.gz'):
        (ref / name).touch()  # the reduced reference's genome and annotation (reduce_reference.py)
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({
        'name': 'syn', 'kind': 'synthetic', 'organism': 'mus_musculus', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
        'adapter': 'TGGAATTCTCGGGTGCCAAGG', 'methods': [method],
        'reference': {'name': 'mmu-test', 'dir': 'benchmark/references/mmu-test'},
        'synthetic': {'seed': 1, 'samples': 1, 'mirnas': {'mmu-miR-192-5p': 4}}}))
    project = tmp_path / 'project'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(tmp_path / 'ext')], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert config['quant']['method'] == method and config['db']['mirflowz']['custom']['genome'] == str(ref / 'genome.fa.gz')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '--cores', '1', 'quant_all'], cwd=project, env=env,
                            capture_output=True, text=True, timeout=120)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    genome = output.split('rule tesmall_genome:', 1)[1].split('\nrule ', 1)[0]
    assert 'input: {}'.format(ref / 'genome.fa.gz') in genome and 'tesmall/custom-mmu-test/mus_musculus/genome/genome.fa' in genome
    if method in ('tinyrna', 'biocore_htseq'):
        gff = output.split('rule tinyrna_gff:', 1)[1].split('\nrule ', 1)[0]
        assert 'input: {}'.format(ref / 'mirna.gff3') in gff and 'tinyrna/custom-mmu-test/mus_musculus/gff-custom/mirna.v2.gff3' in gff
    if method == 'biocore_htseq':
        index = output.split('rule biocore_htseq_bowtie2_index:', 1)[1].split('\nrule ', 1)[0]
        assert 'tesmall/custom-mmu-test/mus_musculus/genome/genome.fa' in index
        assert 'biocore_htseq/custom-mmu-test/mus_musculus/genome/bowtie2/genome.1.bt2' in index


def test_prepare_refuses_an_existing_project_and_a_missing_reference(tmp_path, ref):
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({'name': 'syn', 'kind': 'synthetic', 'organism': 'mus_musculus',
                                                         'libprepkit': 'x', 'adapter': 'TGG', 'methods': ['unitas'],
                                                         'reference': {'name': 'n', 'dir': 'benchmark/references/missing'},
                                                         'synthetic': {'mirnas': {}}}))
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--ext-dir', str(tmp_path / 'ext')]
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'p')], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'lacks mature.fa' in result.stderr
    (tmp_path / 'existing').mkdir()
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'existing')], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'exists' in result.stderr


# Pools: a human miRBase subset and a composition table in the study's two formats.
HSA = ('>hsa-miR-1-3p MIMAT0000416 Homo sapiens miR-1-3p\nUGGAAUGUAAAGAAGUAUGUAU\n'
       '>hsa-miR-199a-3p MIMAT0000232 Homo sapiens miR-199a-3p\nACAGUAGUCUGCACAUUGGUUA\n'
       '>hsa-miR-199b-3p MIMAT0004563 Homo sapiens miR-199b-3p\nACAGUAGUCUGCACAUUGGUUA\n'
       '>hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\nUAGCUUAUCAGACUGAUGUUGA\n'
       '>hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p\nUAGCAGCACGUAAAUAUUGGCG\n'
       '>hsa-miR-145-5p MIMAT0000437 Homo sapiens miR-145-5p\nGUCCAGUUUUCCCAGGAAUCCCU\n'
       '>mmu-miR-1a-3p MIMAT0000123 Mus musculus miR-1a-3p\nUGGAAUGUAAAGAAGUAUGUAU\n')
EQ = ('Sequence.ID\tSequence.type\tModification\tsequence\n'
      "MIR-1\tmiRNA\t5'-phosphorylation\tTGGAATGTAAAGAAGTATGTAT\n"            # exact
      "MIR-199A-3P\tmiRNA\t5'-phosphorylation\tACAGTAGTCTGCACATTGGTTA\n"      # exact, two paralog accessions
      "MIR-21\tmiRNA\t5'-phosphorylation\tTAGCTTATCAGACTGATGTTGAC\n"         # miR-21-5p + C: shifted
      "EBV-1\tmiRNA\t5'-phosphorylation\tTAACCTGATCAGCCCCGGAGTT\n"           # foreign
      "LONG-1\tother\t5'-phosphorylation\tACGTACGTACGTACGTACGTACGTACGTACG\n"  # 31 nt: excluded
      "MOD-1\tmiRNA\t5'-phosphorylation, 3'-biotin\tGGGGCCCCAAAATTTTGGGG\n")  # other modification: excluded
RATIO = ('ratio.seqID\tLength\tType\tsequence\tModification\tX.SAMPLE.A.\tX.SAMPLE.B.\n'
         "hsa-miR-1-3p\t22\tHuman miRNA\tTGGAATGTAAAGAAGTATGTAT\t5'-phosphorylation\t4x\t1x\n"
         "hsa-miR-199a-3p\t22\tHuman miRNA\tACAGTAGTCTGCACATTGGTTA\t5'-phosphorylation\t1x\t1x\n"
         "hsa-miR-21-5p\t22\tHuman miRNA\tTAGCTTATCAGACTGATGTTGA\t5'-phosphorylation\t1x\t8x\n")


def hsa_mature(tmp_path):
    return synthetic.read_mature(write(tmp_path / 'mature.fa', HSA), 'hsa-')


def test_pool_truth_groups(tmp_path):
    rows = {r[0]: r for r in pool.truth_groups(pool.read_composition(write(tmp_path / 'eq.csv', EQ)), hsa_mature(tmp_path))}
    assert {k: r[2] for k, r in rows.items()} == {'MIR-1': 'exact', 'MIR-199A-3P': 'exact', 'MIR-21': 'shifted',
                                                  'EBV-1': 'foreign', 'LONG-1': 'excluded', 'MOD-1': 'excluded'}
    assert rows['MIR-199A-3P'][3] == 'MIMAT0000232,MIMAT0004563' and rows['MIR-21'][3] == 'MIMAT0000076'
    ratio = pool.truth_groups(pool.read_composition(write(tmp_path / 'ratio.csv', RATIO)), hsa_mature(tmp_path))
    assert [(r[0], r[4], r[5]) for r in ratio] == [('hsa-miR-1-3p', '4', '1'), ('hsa-miR-199a-3p', '1', '1'), ('hsa-miR-21-5p', '1', '8')]
    assert pool.shifted('ATGGAATGTAAAGAAGTATG', 'TGGAATGTAAAGAAGTATGTAT')  # 5' end +1, 3' end -2 (older boundary)
    assert not pool.shifted('ATGGAATGTAAAGAAGTATC', 'TGGAATGTAAAGAAGTATGTAT')  # a mismatch in the overlap
    assert not pool.shifted('GGGGTGGAATGTAAAGAAGTATGTAT', 'TGGAATGTAAAGAAGTATGTAT')  # 4 nt more at the 5' end
    duplicated = EQ + "MIR-1-COPY\tmiRNA\t5'-phosphorylation\tTGGAATGTAAAGAAGTATGTAT\n"
    groups = [r[2] for r in pool.truth_groups(pool.read_composition(write(tmp_path / 'dup.csv', duplicated)), hsa_mature(tmp_path))]
    assert groups[0] == groups[-1] == 'excluded'  # a sequence repeated in the pool is not analysed


def test_oracle_counts_exact_matches_only(tmp_path):
    reads = ['TGGAATGTAAAGAAGTATGTAT'] * 3 + ['TGGAATGTAAAGAAGTATGTA', 'TGGAATGTAAAGAAGTATGTATA', 'ACAGTAGTCTGCACATTGGTTA']
    fastq = write(tmp_path / 'f.fastq', ''.join('@r{}\n{}\n+\n{}\n'.format(i, r, 'I' * len(r)) for i, r in enumerate(reads)))
    assert pool.oracle(fastq, ['TGGAATGTAAAGAAGTATGTAT', 'ACAGTAGTCTGCACATTGGTTA', 'TAACCTGATCAGCCCCGGAGTT']) == {
        'TGGAATGTAAAGAAGTATGTAT': 3, 'ACAGTAGTCTGCACATTGGTTA': 1, 'TAACCTGATCAGCCCCGGAGTT': 0}


def pool_entries(tmp_path, table):
    pool.write_pool(tmp_path / 'pool.tsv', pool.truth_groups(pool.read_composition(write(tmp_path / 'c.csv', table)), hsa_mature(tmp_path)))
    return pool.read_pool(tmp_path / 'pool.tsv')


def test_equimolar_metrics_against_the_oracle(tmp_path):
    entries = pool_entries(tmp_path, EQ)
    oracle = {'S1': {'TGGAATGTAAAGAAGTATGTAT': 100, 'ACAGTAGTCTGCACATTGGTTA': 10}}
    common = [('good', 'S1', 'miR-1-3p', 'MIMAT0000416', 100), ('good', 'S1', 'miR-199a-3p', 'MIMAT0000232', 6),
              ('good', 'S1', 'miR-199b-3p', 'MIMAT0004563', 4),  # paralog accessions of one pool sequence add up
              ('off', 'S1', 'miR-1-3p', 'MIMAT0000416', 2000), ('off', 'S1', 'miR-21-5p', 'MIMAT0000076', 50),
              ('off', 'S1', 'miR-155-5p', 'MIMAT0000646', 10),  # not in the pool: from a foreign sequence
              ('off', 'S1', 'miR-9999', '', 40)]  # a name without a miRBase accession: unplaced
    counts, reads = pool.method_counts(common, entries)
    assert counts[('good', 'S1', 'MIR-199A-3P')] == 10 and counts[('off', 'S1', 'MIR-21')] == 50
    assert reads[('off', 'S1')] == {'all': 2100, 'outside': 10, 'unplaced': 40}
    m = pool.equimolar_metrics(counts, reads, oracle, entries, ['good', 'off'], ['S1'])
    assert m['good']['median_log2_vs_oracle'] == 0 and m['good']['within_2fold'] == 1 and m['good']['outside'] == 0
    assert m['off']['beyond_10fold'] == 1 and m['off']['detected'] == 0.5
    assert m['off']['outside'] == round(10 / 2100, 4)  # among all reported reads: the shifted entry's and unplaced ones too
    assert m['off']['unplaced'] == round(40 / 2100, 4) and m['good']['unplaced'] == 0
    assert m['oracle']['detected'] == 1 and 'within_2fold' not in m['oracle']


def test_exact_entry_sharing_its_mirna_with_a_shifted_one_is_not_scored(tmp_path):
    """miR-1-3p one base shorter is a shifted entry with the exact entry's accession: a method's miR-1-3p count holds
    the reads of both sequences, the oracle of the exact entry one, so the exact entry is left out."""
    entries = pool_entries(tmp_path, EQ + "MIR-1-2\tmiRNA\t5'-phosphorylation\tTGGAATGTAAAGAAGTATGTA\n")
    assert {e['pool_id']: e['group'] for e in entries}['MIR-1-2'] == 'shifted'
    assert pool.shared(entries) == {'MIR-1': ['MIR-1-2'], 'MIR-1-2': ['MIR-1']}
    assert [e['pool_id'] for e in pool.scored(entries)] == ['MIR-199A-3P']
    oracle = {'S1': {'TGGAATGTAAAGAAGTATGTAT': 100, 'TGGAATGTAAAGAAGTATGTA': 100, 'ACAGTAGTCTGCACATTGGTTA': 10}}
    common = [('m', 'S1', 'miR-1-3p', 'MIMAT0000416', 200), ('m', 'S1', 'miR-199a-3p', 'MIMAT0000232', 10)]
    m = pool.equimolar_metrics(*pool.method_counts(common, entries), oracle, entries, ['m'], ['S1'])
    assert m['m']['median_log2_vs_oracle'] == 0 and m['m']['within_2fold'] == 1  # not 200 against 100 for MIR-1
    counts, _ = pool.method_counts(common, entries)
    kinds = pool.shift_kinds(entries, dict(hsa_mature(tmp_path).values()))
    rows = {r[1]: r for r in pool.per_sequence(counts, oracle, entries, ['m'], ['S1'], kinds)}  # per_sequence.tsv
    assert rows['MIR-1'] == ['S1', 'MIR-1', 'exact', '', 'MIR-1-2', 100, 'NA']
    assert rows['MIR-1-2'] == ['S1', 'MIR-1-2', 'shifted', '3p_shorter', 'MIR-1', 100, 'NA']
    assert rows['MIR-199A-3P'] == ['S1', 'MIR-199A-3P', 'exact', '', '', 10, '10']
    assert [e['pool_id'] for e in pool.scored(entries, 'shifted')] == ['MIR-21']  # MIR-1-2 shares miR-1-3p


SHIFTS = ("MIR-16-MOVED\tmiRNA\t5'-phosphorylation\tAGCAGCACGTAAATATTGGCG\n"     # miR-16-5p without its first base
          "MIR-145-SHORT\tmiRNA\t5'-phosphorylation\tGTCCAGTTTTCCCAGGAATCC\n")  # miR-145-5p 2 nt shorter at the 3' end


def test_shifted_entries_are_scored_by_kind_of_shift(tmp_path):
    """isomiRs of known amount: a method should count them on their miRNA; the kind of shift says whether the seed
    (5' end) is kept."""
    assert pool.shift('AGCAGCACGTAAATATTGGCG', 'TAGCAGCACGTAAATATTGGCG') == (1, 0)
    assert pool.shift('TAGCTTATCAGACTGATGTTGAC', 'TAGCTTATCAGACTGATGTTGA') == (0, 1)
    entries = pool_entries(tmp_path, EQ + SHIFTS)
    kinds = pool.shift_kinds(entries, dict(hsa_mature(tmp_path).values()))
    assert kinds == {'MIR-21': '3p_longer', 'MIR-16-MOVED': '5p_moved', 'MIR-145-SHORT': '3p_shorter'}
    two = [{'pool_id': 'X', 'group': 'shifted', 'sequence': 'ACGTACGTACGTACGTAC', 'mature_ids': ['A', 'B']}]
    assert pool.shift_kinds(two, {'A': 'ACGTACGTACGTACGTACT', 'B': 'GACGTACGTACGTACGTAC'}) == {'X': 'mixed'}
    oracle = {'S1': {'TAGCTTATCAGACTGATGTTGAC': 100, 'AGCAGCACGTAAATATTGGCG': 100, 'GTCCAGTTTTCCCAGGAATCC': 100}}
    common = [('m', 'S1', 'miR-21-5p', 'MIMAT0000076', 100), ('m', 'S1', 'miR-145-5p', 'MIMAT0000437', 40),
              ('m', 'S1', 'miR-16-2-5p', '', 100)]  # the moved 5' end counted on a name without an accession
    counts, _ = pool.method_counts(common, entries)
    shifted = pool.scored(entries, 'shifted')
    m = pool.agreement(counts, oracle, shifted, ['m'], ['S1'])
    assert m['oracle'] == {'detected': 1} and m['m']['detected'] == round(2 / 3, 4) and m['m']['within_2fold'] == round(1 / 3, 4)
    moved = pool.agreement(counts, oracle, [e for e in shifted if kinds[e['pool_id']] == '5p_moved'], ['m'], ['S1'])
    assert moved['m']['detected'] == 0 and moved['m']['beyond_10fold'] == 1


def test_ratiometric_metrics_against_the_design(tmp_path):
    entries = pool_entries(tmp_path, RATIO)
    design = {'A': (4000, 1000, 1000), 'B': (1000, 1000, 8000)}
    groups = {'A1': 'A', 'A2': 'A', 'B1': 'B'}
    seqs = [e['sequence'] for e in entries]
    oracle = {s: dict(zip(seqs, design[g])) for s, g in groups.items()}
    ids = [e['pool_id'] for e in entries]
    counts = {('flat', s, i): 100 for s in groups for i in ids}
    counts.update({('exact', s, i): v for s, g in groups.items() for i, v in zip(ids, design[g])})
    m, observed = pool.ratiometric_metrics(counts, oracle, entries, ['exact', 'flat'], list(groups.items()), pseudo=0)
    assert m['oracle'] == m['exact'] and m['exact']['median_abs_error'] == 0 and m['exact']['called'] == 1
    assert m['exact']['false_calls'] == 0 and m['flat']['called'] == 0
    assert observed['exact']['hsa-miR-21-5p'] < -2 and 'shifted_called' not in m['exact']  # no shifted entries

    # A shifted entry (miR-16-5p without its first base, 1:4) is scaled and centred as the exact ones: the exact
    # group's results do not change.
    shifted = pool_entries(tmp_path, RATIO + "MIR-16-MOVED\t21\tHuman miRNA\tAGCAGCACGTAAATATTGGCG\t5'-phosphorylation\t1x\t4x\n")
    for s, g in groups.items():
        oracle[s]['AGCAGCACGTAAATATTGGCG'] = counts[('exact', s, 'MIR-16-MOVED')] = {'A': 250, 'B': 1000}[g]
        counts[('flat', s, 'MIR-16-MOVED')] = 100
    with_shifted, observed = pool.ratiometric_metrics(counts, oracle, shifted, ['exact', 'flat'], list(groups.items()), pseudo=0)
    assert {k: str(v) for k, v in with_shifted['flat'].items() if not k.startswith('shifted_')} == {k: str(v) for k, v in m['flat'].items()}
    assert with_shifted['exact']['shifted_median_abs_error'] == 0 and with_shifted['exact']['shifted_called'] == 1
    assert with_shifted['flat']['shifted_called'] == 0 and observed['exact']['MIR-16-MOVED'] == -2


def test_fill_checks_and_places_files(tmp_path):
    source = tmp_path / 'source'
    data = {'SRR1.fastq.gz': b'reads', 'eq.csv': EQ.encode(), 'mature.fa': HSA.encode()}
    for name, body in data.items():
        write(source / name, '').write_bytes(body)
    import hashlib
    spec = {'name': 'p', 'kind': 'pool', 'data': 'benchmark/smallrna/p',
            'composition': {'file': 'eq.csv', 'url': (source / 'eq.csv').as_uri(), 'sha256': hashlib.sha256(data['eq.csv']).hexdigest()},
            'mirbase': {'release': '23.1', 'url': (source / 'mature.fa').as_uri(), 'sha256': hashlib.sha256(data['mature.fa']).hexdigest()},
            'runs': [{'run': 'SRR1', 'url': (source / 'SRR1.fastq.gz').as_uri(), 'md5': hashlib.md5(b'reads').hexdigest()}]}
    ext = tmp_path / 'ext'
    assert fill.fill(spec, ext, log=lambda _: None) == 3
    out = ext / 'benchmark/smallrna/p'
    assert (out / 'fastq/SRR1.fastq.gz').read_bytes() == b'reads' and (out / 'mirbase/23.1/mature.fa').read_text() == HSA
    assert fill.fill(spec, ext, log=lambda _: None) == 0  # in place: skipped
    spec['runs'].append({'run': 'SRR2', 'url': (source / 'SRR1.fastq.gz').as_uri(), 'md5': '0' * 32})
    with pytest.raises(ValueError, match='does not match'):
        fill.fill(spec, ext, log=lambda _: None)
    assert sorted(p.name for p in (out / 'fastq').iterdir()) == ['SRR1.fastq.gz']  # no partial or failed file left


def test_concurrent_downloads_of_one_file_keep_their_own_temporary_files(tmp_path, monkeypatch):
    """A second download of the same file starts while the first is still writing (fills of two declarations sharing a
    data folder): one fails its checksum, one succeeds; the first then completes. Only verified bytes end up in place."""
    import hashlib
    import io
    import urllib.request
    good = b'reads' * 1000
    path = tmp_path / 'fastq/SRR1.fastq.gz'
    opened = urllib.request.urlopen
    log = lambda _: None

    class Interleaved(io.BytesIO):
        first = True

        def read(self, size=-1):
            if Interleaved.first:
                Interleaved.first = False
                with pytest.raises(ValueError, match='does not match'):
                    fill.download('https://x/good', 'md5', '0' * 32, path, log)
                fill.download('https://x/good', 'md5', hashlib.md5(good).hexdigest(), path, log)
            return super().read(size)

    def urlopen(request, timeout):
        assert request.get_header('User-agent') == fill.USER_AGENT  # a client name of its own: GtRNAdb refuses Python's default
        url = request.full_url
        return Interleaved(good) if url == 'https://x/first' else io.BytesIO(good) if url == 'https://x/good' else opened(request, timeout=timeout)
    monkeypatch.setattr(fill.urllib.request, 'urlopen', urlopen)
    fill.download('https://x/first', 'md5', hashlib.md5(good).hexdigest(), path, log)
    assert path.read_bytes() == good and [p.name for p in path.parent.iterdir()] == ['SRR1.fastq.gz']


def test_a_download_the_server_cuts_is_resumed_from_where_it_stopped(tmp_path, monkeypatch):
    """A server that closes a long transfer early gives fewer bytes than it announced, which reads like the end of the
    file. The rest is asked for with a Range request on a new connection; the file in place is complete and checked.
    A server that keeps cutting ends the download after fill.CONNECTIONS connections, with nothing left behind."""
    import hashlib
    import io
    good = bytes(range(256)) * 40
    path = tmp_path / 'fastq/SRR1.fastq.gz'
    asked, messages, waits = [], [], []

    class Response(io.BytesIO):
        def __init__(self, data, start, announced):
            super().__init__(data)
            self.status = 206 if start else 200
            self.headers = {'Content-Length': str(announced)}

    def urlopen(request, timeout):
        assert request.get_header('User-agent') == fill.USER_AGENT  # on the first request and the resumed ones (#412)
        start = int(request.get_header('Range')[len('bytes='):-1]) if request.get_header('Range') else 0
        asked.append(start)
        return Response(good[start:start + cut], start, len(good) - start)  # at most `cut` bytes, then the connection closes
    monkeypatch.setattr(fill.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(fill.time, 'sleep', waits.append)
    cut = 4000
    fill.download('https://example.org/SRR1.fastq.gz', 'md5', hashlib.md5(good).hexdigest(), path, messages.append)
    assert path.read_bytes() == good and asked == [0, 4000, 8000] and waits == [30, 60]
    assert 'cut at 4000 of 10240 bytes; resuming' in messages[0] and messages[-1].endswith('md5 checked')
    path.unlink()
    cut, asked[:] = 100, []
    with pytest.raises(ValueError, match='incomplete after 20 connections \\(2000 of 10240 bytes\\)'):
        fill.download('https://example.org/SRR1.fastq.gz', 'md5', hashlib.md5(good).hexdigest(), path, messages.append)
    assert len(asked) == fill.CONNECTIONS and list(path.parent.iterdir()) == []

    class Full:  # a local write error is not a cut transfer: one connection, and the error ends the download
        def write(self, block):
            raise OSError(28, 'No space left on device')
    cut, asked[:] = 4000, []
    with pytest.raises(OSError, match='No space left'):
        fill.transfer('https://example.org/SRR1.fastq.gz', Full(), hashlib.md5(), messages.append)
    assert asked == [0]


def test_runs_filled_from_sra_are_checked_by_their_reads_not_their_file(tmp_path, monkeypatch):
    """fill --source sra: a run that declares its SRA file is downloaded in parts (range requests, a short part asked
    for again), checked against the file's md5, converted by the fasterq-dump command and written gzipped with a record
    of its source. The reads are checked by their number and by the MD5 over the sequence and quality lines, which
    does not depend on read names, block boundaries or compression. A run without `sra`, a wrong number of reads and a
    wrong content MD5 leave nothing behind."""
    import gzip
    import hashlib
    import json
    import stat
    reads = [('ACGT' * 5 + 'ACGT'[:i % 4 + 1], 'I' * (20 + i % 4 + 1)) for i in range(300)]
    fastq = ''.join('@SRR1.{0} {0} length={1}\n{2}\n+SRR1.{0} {0} length={1}\n{3}\n'.format(i + 1, len(s), s, q) for i, (s, q) in enumerate(reads)).encode()
    digest = hashlib.md5(''.join('{}\n{}\n'.format(s, q) for s, q in reads).encode()).hexdigest()
    for size in (1, 7, 4096, len(fastq)):  # blocks cut anywhere
        assert fill.content(fastq[i:i + size] for i in range(0, len(fastq), size)) == (300, digest)
    assert fill.content([fastq[:-1]]) == (300, digest)  # no newline after the last line
    with pytest.raises(ValueError, match='incomplete record after read 299'):
        fill.content([fastq[:-40]])
    with pytest.raises(ValueError, match='not FASTQ'):
        fill.content([b'>seq\nACGT\n>seq2\nACGT\n'])
    with pytest.raises(ValueError, match='one quality per base'):
        fill.content([b'@r\nACGT\n+\nIII\n'])

    sra = bytes(range(256)) * 100  # the "SRA file"; the fake fasterq-dump ignores it and prints the FASTQ
    (tmp_path / 'reads.fastq').write_bytes(fastq)
    dump = tmp_path / 'fasterq-dump'
    dump.write_text('#!/bin/sh\nif [ "$1" = --version ]; then echo; echo "fasterq-dump : 9.9.9"; exit 0; fi\n'
                    '[ -s "$1" ] || exit 3\ncat ' + str(tmp_path / 'reads.fastq') + '\n')
    dump.chmod(dump.stat().st_mode | stat.S_IXUSR)
    asked = []

    class Response(io.BytesIO):
        status = 206

    def urlopen(request, timeout):
        start, end = (int(v) for v in request.headers['Range'][len('bytes='):].split('-'))
        asked.append((request.full_url, start, end))
        short = start == 10000 and asked.count((request.full_url, start, end)) == 1
        return Response(sra[start:end + 1 - (5 if short else 0)])
    monkeypatch.setattr(fill.urllib.request, 'urlopen', urlopen)
    monkeypatch.setattr(fill.time, 'sleep', lambda seconds: None)
    monkeypatch.setattr(fill, 'PART', 10000)
    run = {'run': 'SRR1', 'reads': 300, 'md5': 'x', 'url': 'https://example.org/ena/SRR1.fastq.gz',
           'sra': {'md5': hashlib.md5(sra).hexdigest(), 'bytes': len(sra)}, 'content_md5': digest}
    spec = {'name': 'set', 'kind': 'biological', 'data': 'benchmark/set',
            'mirbase': {'release': '23.1', 'url': 'https://example.org/mature.fa', 'sha256': 'x'},
            'samples': [{'sample': 'A', 'group': 'g', 'runs': [run]}]}
    folder = tmp_path / 'ext/benchmark/set'
    (folder / 'mirbase/23.1').mkdir(parents=True)
    (folder / 'mirbase/23.1/mature.fa').write_text('>x\nACGT\n')
    messages = []
    monkeypatch.chdir(tmp_path)  # a relative EXT_DIR: the command, run in the temporary folder, still finds the file
    (tmp_path / 'link').symlink_to(tmp_path / 'ext')  # and a store reached through a link is named by the link
    dump.write_text(dump.read_text().replace('[ -s "$1" ] || exit 3', 'case "$1" in */link/*) ;; *) exit 4;; esac\n[ -s "$1" ] || exit 3'))
    assert fill.fill(spec, Path('link'), log=messages.append, source='sra', connections=3, fasterq_dump=str(dump)) == 1
    target = folder / 'fastq/SRR1.fastq.gz'
    assert sorted(p.name for p in target.parent.iterdir()) == ['SRR1.fastq.gz', 'SRR1.sra.json']
    assert fill.file_content(target) == (300, digest) and gzip.open(target, 'rt').read().split('\n')[2] == '+'
    record = json.loads((target.parent / 'SRR1.sra.json').read_text())
    assert record['source'] == fill.SRA_URL.format(run='SRR1') and record['sra_md5'] == run['sra']['md5']
    assert (record['reads'], record['content_md5'], record['fasterq_dump']) == (300, digest, '9.9.9')
    assert sorted(set(asked)) == [(fill.SRA_URL.format(run='SRR1'), start, min(start + 9999, len(sra) - 1)) for start in range(0, len(sra), 10000)]
    assert asked.count((fill.SRA_URL.format(run='SRR1'), 10000, 19999)) == 2  # the short part, asked for again
    assert fill.fill(spec, tmp_path / 'link', log=messages.append, source='sra', fasterq_dump=str(dump)) == 0  # in place

    for change, message in (({'reads': 299}, '300 reads from SRA, 299 declared'), ({'content_md5': '0' * 32}, 'content md5'),
                            ({'sra': {'md5': '0' * 32, 'bytes': len(sra)}}, 'md5 .* does not match the declared'), ({'sra': None}, 'declares `reads` and `sra')):
        target.unlink(missing_ok=True)
        (target.parent / 'SRR1.sra.json').unlink(missing_ok=True)
        spec['samples'][0]['runs'] = [dict(run, **change)]
        with pytest.raises(ValueError, match=message):
            fill.fill(spec, tmp_path / 'link', log=messages.append, source='sra', fasterq_dump=str(dump))
        assert list(target.parent.iterdir()) == []


def test_filtered_fastq_follows_the_project_filter_route(tmp_path):
    """The oracle reads the route's filtered reads: the trimming route's (cutadapt by default), or the spike-in
    prefilter's when the project lists spike-in sets, also when the trimming route's remain (quant_all --notemp)."""
    import benchmark
    filtered = tmp_path / 'data/tmp/smallrna/filter'
    config = write(tmp_path / 'config.yaml', yaml.safe_dump({'filter': {'spikein': {'ref': None}}}))
    with pytest.raises(SystemExit, match='quant_all --notemp'):
        benchmark.filtered_fastq(tmp_path, 'S1')
    write(filtered / 'cutadapt/S1_R1.fastq', '')
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'cutadapt/S1_R1.fastq'
    write(filtered / 'fastp/S1_R1.fastq', '')
    with pytest.raises(SystemExit, match='more than one trimming route'):
        benchmark.filtered_fastq(tmp_path, 'S1')
    (filtered / 'fastp/S1_R1.fastq').unlink()
    write(filtered / 'spikein/bowtie/S1_R1.fastq', '')  # left by an earlier run with spike-ins
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'cutadapt/S1_R1.fastq'
    config.write_text(yaml.safe_dump({'filter': {'spikein': {'ref': 'ercc'}}}))
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'spikein/bowtie/S1_R1.fastq'
    config.write_text(yaml.safe_dump({'filter': {'spikein': {'ref': 'ercc', 'quantifier': 'native'}}}))
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'cutadapt/S1_R1.fastq'  # UNITAS and miRge3 read these
    (filtered / 'spikein/bowtie/S1_R1.fastq').unlink()  # native with UNITAS and miRge3 only: no prefilter runs
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'cutadapt/S1_R1.fastq'
    config.write_text(yaml.safe_dump({'filter': {'spikein': {'ref': 'ercc'}}}))  # bowtie, the workflow default
    with pytest.raises(SystemExit, match='quant_all --notemp'):
        benchmark.filtered_fastq(tmp_path, 'S1')


def fastq(seqs):
    return ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)) for i, s in enumerate(seqs))


def write_stages(project, sample, adapter, filtered, final, cuts=(0, 0)):
    """The cutadapt route's reads after each pass and the filter pass's report (quant_all --notemp)."""
    stages = project / 'data/tmp/smallrna/filter'
    write(stages / 'cutadapt_trimmed/{}_R1.fastq'.format(sample), fastq(adapter))
    write(stages / 'cutadapt_filtered/{}_R1.fastq'.format(sample), fastq(filtered))
    write(stages / 'cutadapt/{}_R1.fastq'.format(sample), fastq(final))
    args = ['-j', '4'] + (['-u', str(cuts[0]), '-u', str(-cuts[1])] if any(cuts) else []) + ['--action=none', '-m', '15']
    write(stages / 'cutadapt_stages/{}.filter.json'.format(sample), json.dumps({'command_line_arguments': args}))


def test_end_cuts_from_the_filter_pass_report(tmp_path):
    """The reads after the adapter pass still carry the kit's random ends; the next pass's report says how many bases
    it cuts (cutadapt -u, fastp -f/-t)."""
    import benchmark
    stages = tmp_path / 'data/tmp/smallrna/filter'
    write(stages / 'cutadapt_stages/S1.filter.json', json.dumps({'command_line_arguments': ['-j', '4', '-u', '4', '-u', '-4', '-m', '15']}))
    assert benchmark.end_cuts(tmp_path, 'cutadapt', 'S1') == (4, 4)
    write(stages / 'cutadapt_stages/S1.filter.json', json.dumps({'command_line_arguments': ['--cut=3', '-m', '15']}))
    assert benchmark.end_cuts(tmp_path, 'cutadapt', 'S1') == (3, 0)
    write(stages / 'fastp_stages/S1.filter.json', json.dumps({'command': 'fastp -i a -o b -f 4 -t 4 -l 16 -u 40 -A'}))
    assert benchmark.end_cuts(tmp_path, 'fastp', 'S1') == (4, 4)
    write(stages / 'fastp_stages/S1.filter.json', json.dumps({'command': 'fastp -i a -o b --trim_front1=2 -l 16 -A'}))
    assert benchmark.end_cuts(tmp_path, 'fastp', 'S1') == (2, 0)
    with pytest.raises(SystemExit, match='quant_all --notemp'):
        benchmark.end_cuts(tmp_path, 'cutadapt', 'S2')


def test_filter_stages_and_what_each_pass_removes(tmp_path):
    """Pool sequences counted after each pass of the route: the adapter pass's reads with the random ends cut, and
    the reads after each later pass, the last being the oracle's."""
    import benchmark
    write(tmp_path / 'config.yaml', yaml.safe_dump({'filter': {'spikein': {'ref': None}}}))
    entries = pool_entries(tmp_path, EQ)
    mir1, mir21, ebv = 'TGGAATGTAAAGAAGTATGTAT', 'TAGCTTATCAGACTGATGTTGAC', 'TAACCTGATCAGCCCCGGAGTT'
    ends = lambda s: 'ACGT' + s + 'TTGA'  # NEXTflex-like random ends, cut by the filter pass
    adapter = [ends(mir1)] * 40 + [ends(mir21)] * 30 + [ends(ebv)] * 10 + ['ACGTTT']
    write_stages(tmp_path, 'S1', adapter, [mir1] * 40 + [mir21] * 20 + [ebv] * 10, [mir1] * 40 + [mir21] * 20 + [ebv] * 5, cuts=(4, 4))
    stages = benchmark.filter_stages(tmp_path, 'S1')
    assert [(st, path.parent.name, cuts) for st, path, cuts in stages] == [
        ('adapter', 'cutadapt_trimmed', (4, 4)), ('filter', 'cutadapt_filtered', (0, 0)), ('artefacts', 'cutadapt', (0, 0))]
    sequences = [e['sequence'] for e in entries if e['group'] != 'excluded']
    counts = {'S1': {st: pool.oracle(path, sequences, *cuts) for st, path, cuts in stages}}
    assert counts['S1']['adapter'][mir1] == 40 and counts['S1']['adapter'][mir21] == 30
    groups, lossy = pool.filter_losses(counts, entries, ['adapter', 'filter', 'artefacts'])
    assert groups['exact'] == {'sequences': 2, 'adapter': 40, 'filter': 40, 'artefacts': 40}
    assert groups['shifted'] == {'sequences': 1, 'adapter': 30, 'filter': 20, 'artefacts': 20}
    assert groups['foreign'] == {'sequences': 1, 'adapter': 10, 'filter': 10, 'artefacts': 5}
    assert [(e['pool_id'], worst) for e, _, worst in lossy] == [('MIR-21', 'filter')]  # EBV-1: fewer than 20 reads

    # With Bowtie spike-ins the prefilter's reads come last; the fastp route has two passes; skip filters nothing.
    write(tmp_path / 'config.yaml', yaml.safe_dump({'filter': {'spikein': {'ref': 'ercc'}}}))
    write(tmp_path / 'data/tmp/smallrna/filter/spikein/bowtie/S1_R1.fastq', fastq([mir1]))
    assert [st for st, _, _ in benchmark.filter_stages(tmp_path, 'S1')] == ['adapter', 'filter', 'artefacts', 'spikein']
    write(tmp_path / 'config.yaml', yaml.safe_dump({'filter': {'spikein': {'ref': None}}}))
    filtered = tmp_path / 'data/tmp/smallrna/filter'
    for folder in ('cutadapt_trimmed', 'cutadapt_filtered', 'cutadapt'):
        (filtered / folder / 'S1_R1.fastq').unlink()
    write(filtered / 'fastp/S1_R1.fastq', '')
    with pytest.raises(SystemExit, match='reads of the filter passes missing.*fastp_trimmed'):
        benchmark.filter_stages(tmp_path, 'S1')
    write(filtered / 'fastp_trimmed/S1_R1.fastq', '')
    write(filtered / 'fastp_stages/S1.filter.json', json.dumps({'command': 'fastp -f 4 -t 4 -A'}))
    assert [(st, cuts) for st, _, cuts in benchmark.filter_stages(tmp_path, 'S1')] == [('adapter', (4, 4)), ('filter', (0, 0))]
    for folder in ('fastp', 'fastp_trimmed'):
        (filtered / folder / 'S1_R1.fastq').unlink()
    write(filtered / 'fastq/S1_R1.fastq', '')
    assert [st for st, _, _ in benchmark.filter_stages(tmp_path, 'S1')] == ['input']  # skip: no pass to compare
    write(tmp_path / 'config.yaml', yaml.safe_dump({'filter': {'spikein': {'ref': 'ercc'}}}))
    assert [(st, cuts) for st, _, cuts in benchmark.filter_stages(tmp_path, 'S1')] == [('input', (0, 0)), ('spikein', (0, 0))]


def test_prepare_and_evaluate_a_pool(tmp_path):
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/p'
    write(data / 'composition/ratio.csv', RATIO)
    write(data / 'mirbase/23.1/mature.fa', HSA)
    reads = {'SRRA': ['TGGAATGTAAAGAAGTATGTAT'] * 8 + ['TAGCTTATCAGACTGATGTTGA'] * 2,
             'SRRB': ['TGGAATGTAAAGAAGTATGTAT'] * 2 + ['TAGCTTATCAGACTGATGTTGA'] * 8}
    for run, seqs in reads.items():
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/{}.fastq.gz'.format(run), 'wt') as fh:
            fh.write(''.join('@r{}\n{}TGGAATTCTC\n+\n{}\n'.format(i, s, 'I' * (len(s) + 10)) for i, s in enumerate(seqs)))
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({
        'name': 'p', 'kind': 'pool', 'pool': 'ratiometric', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
        'methods': ['unitas', 'mirge3', 'mirflowz'], 'data': 'benchmark/smallrna/p',
        'composition': {'file': 'ratio.csv', 'url': 'unused', 'sha256': 'unused'},
        'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'runs': [{'run': 'SRRA', 'group': 'A', 'md5': 'unused', 'url': 'unused'}, {'run': 'SRRB', 'group': 'B', 'md5': 'unused', 'url': 'unused'}]}))
    project = tmp_path / 'project'
    dump = unitas_dump(ext, 'homo_sapiens', 'Human', {})
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(ext),
                                   '--methods', 'unitas,mirge3', '--reads', '6'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert config['quant']['method'] == 'unitas,mirge3' and config['read_geometry'] == [32]
    assert config['db'] == {'mirbase': {'release': '23.1'}}  # the declaration's release, for every method
    with gzip.open(project / 'data/raw/fastq/p/SRRA_R1.fastq.gz', 'rt') as fh:
        assert len(fh.read().splitlines()) == 24  # first 6 reads
    assert '3 exact' in result.stdout

    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tSRRA\tSRRB\nmiR-1-3p\t8\t2\nmiR-21-5p\t2\t8\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,SRRA,SRRB\nhsa-miR-1-3p,8,2\nhsa-miR-21-5p,2,8\n')
    isomir_outputs(quant, {'miR-1-3p': 'TGGAATGTAAAGAAGTATGTAT', 'miR-21-5p': 'TAGCTTATCAGACTGATGTTGA'})
    write_records(quant, ['unitas', 'mirge3'], dump=dump)
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'quant_all --notemp' in result.stderr  # the oracle needs the filtered reads
    for run, seqs in reads.items():  # the adapter pass's reads hold one more read of each sequence, too short ones a filter removes
        write_stages(project, run, seqs + ['TGGAATGTAAAGAAGTATGTAT', 'TAGCTTATCAGACTGATGTTGA', 'ACGT'], seqs, seqs)
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    stages = (project / 'benchmark/filter_stages.tsv').read_text().splitlines()
    assert stages[0] == 'sample\tpool_id\tgroup\tadapter\tfilter\tartefacts' and 'SRRA\thsa-miR-1-3p\texact\t9\t8\t8' in stages
    report = (project / 'benchmark/report.md').read_text()
    isomir_lines = (project / 'benchmark/isomirs.tsv').read_text().splitlines()  # #454: the isomiR table and its metrics
    assert isomir_lines[0] == '\t'.join(adapters.ISOMIR_COLUMNS)
    assert 'mirge3\tSRRA\thsa-miR-1-3p\tMIMAT0000416\tTGGAATGTAAAGAAGTATGTAT\tcanonical\t8.0' in isomir_lines
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (project / 'benchmark/metrics.tsv').read_text().splitlines()}
    assert table['metric'] == ['oracle', 'unitas', 'mirge3'] and table['artefact_share'] == ['', '0.0', '0.0']
    assert table['exact_sequence_recovered'][0] == '' and all(0 < float(v) <= 1 for v in table['exact_sequence_recovered'][1:])
    assert '## isomiRs' in report and 'Without isomiR output' not in report
    ci = {(r[0], r[1]): r[2:] for r in (line.split('\t') for line in (project / 'benchmark/metrics_ci.tsv').read_text().splitlines())}
    for method in ('unitas', 'mirge3'):  # #459: intervals from the isomiR tallies, around the value
        value, low, high = ci[('exact_sequence_recovered', method)]
        assert low != '' and float(low) <= float(value) <= float(high)
    assert '| metric | oracle | unitas | mirge3 |' in report and '| called | 1.0 | 1.0 | 1.0 |' in report
    assert '| 1:8 | -3.00 | 1 |' in report
    per_sequence = (project / 'benchmark/per_sequence.tsv').read_text().splitlines()
    assert per_sequence[0] == 'sample\tpool_id\tgroup\tshift\tshared_with\toracle\tunitas\tmirge3' and 'SRRA\thsa-miR-1-3p\texact\t\t\t8\t8\t8' in per_sequence
    assert '## Filter route' in report and '| exact | 3 | 24 | 16.67% | 0.00% | 83.33% |' in report

    # With spike-ins the oracle reads the prefilter's reads (here one read fewer), while the trimmed reads remain.
    config['filter'] = {'spikein': {'ref': 'ercc'}}
    (project / 'config.yaml').write_text(yaml.safe_dump(config))
    for run, seqs in reads.items():
        write(project / 'data/tmp/smallrna/filter/spikein/bowtie/{}_R1.fastq'.format(run),
              ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)) for i, s in enumerate(seqs[1:])))
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'SRRA\thsa-miR-1-3p\texact\t\t\t7\t8\t8' in (project / 'benchmark/per_sequence.tsv').read_text().splitlines()

    # Native spike-ins with UNITAS and miRge3 only: the workflow runs no prefilter, and the oracle reads the trimmed reads.
    config['filter'] = {'spikein': {'ref': 'ERCC92', 'quantifier': 'native'}}
    (project / 'config.yaml').write_text(yaml.safe_dump(config))
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    def dry_run():
        result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '--cores', '1', '--notemp', '--forceall', 'quant_all'],
                                cwd=project, env=env, capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        return result
    dry = dry_run()
    assert 'rule mirge3_quant:' in dry.stdout and 'rule smallrna_cutadapt_artefacts:' in dry.stdout  # a real plan
    assert 'rule spikein_bowtie:' not in dry.stdout and 'data/tmp/smallrna/filter/cutadapt/SRRA_R1.fastq' in dry.stdout
    config['filter']['spikein'].pop('quantifier')  # bowtie, the default: the prefilter runs
    (project / 'config.yaml').write_text(yaml.safe_dump(config))
    assert 'rule spikein_bowtie:' in dry_run().stdout
    config['filter']['spikein']['quantifier'] = 'native'
    (project / 'config.yaml').write_text(yaml.safe_dump(config))
    for run in reads:
        (project / 'data/tmp/smallrna/filter/spikein/bowtie/{}_R1.fastq'.format(run)).unlink()
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'SRRA\thsa-miR-1-3p\texact\t\t\t8\t8\t8' in (project / 'benchmark/per_sequence.tsv').read_text().splitlines()

    # Bootstrap intervals (#413): published with checksums, the same for the same seed, and removed with --bootstrap 0.
    bench = project / 'benchmark'
    record = json.loads((bench / 'evaluation.json').read_text())
    assert record['bootstrap'] == {'resamples': 1000, 'seed': 1, 'skipped': record['bootstrap']['skipped']}
    assert {'metrics_ci.tsv', 'differences.tsv'} <= set(record) and '## Uncertainty' in (bench / 'report.md').read_text()
    ci = (bench / 'metrics_ci.tsv').read_text()
    assert ci.startswith('metric\trow\tvalue\tlow\thigh\n') and '\toracle\t' in ci
    assert (bench / 'differences.tsv').read_text().startswith('metric\trow_a\trow_b\tdifference\tlow\thigh\n')
    assert subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, timeout=60).returncode == 0
    assert (bench / 'metrics_ci.tsv').read_text() == ci
    assert subprocess.run(cmd + ['evaluate', '--project', str(project), '--bootstrap', '0'], capture_output=True, timeout=60).returncode == 0
    assert not (bench / 'metrics_ci.tsv').exists() and 'bootstrap' not in json.loads((bench / 'evaluation.json').read_text())


# Biological datasets: no truth; methods compared with each other, lane replicates, filtering and first bases.

import biological  # noqa: E402
import json  # noqa: E402
import random  # noqa: E402

BLOOD = ('>hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p\nUAGCAGCACGUAAAUAUUGGCG\n'
         '>hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\nUAGCUUAUCAGACUGAUGUUGA\n'
         '>hsa-miR-451a MIMAT0001631 Homo sapiens miR-451a\nAAACCGUUACCAUUACUGAGUU\n'
         '>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p\nUGAGGUAGUAGGUUGUAUAGUU\n')


def common_rows(method, sample, counts):
    return [(method, sample, mid, mid, c) for mid, c in counts.items()]


def test_method_summary_and_agreement():
    a = {'M1': 900, 'M2': 90, 'M3': 9, 'M4': 1}
    b = {'M1': 800, 'M2': 150, 'M3': 0, 'M5': 50}
    common = common_rows('x', 'S1', a) + common_rows('y', 'S1', b)
    summary = biological.method_summary(common, ['x', 'y'], ['S1'])
    assert summary['x']['detected'] == 2 and summary['y']['detected'] == 3       # at least 10 reads
    assert summary['x']['mirna_reads'] == 1000 and summary['x']['top10_share'] == 1.0
    pair = biological.agreement(common, ['x', 'y'], ['S1'])[('x', 'y')]
    assert pair['only_a'] == 0 and pair['only_b'] == 1                            # M5 only detected by y
    assert pair['spearman_either'] > 0 and round(pair['top100_overlap'], 3) == 0.667     # top sets {M1..M4} and {M1, M2, M5}: y's 0 for M3 is not in it
    blocked = biological.blocked(common, ['x', 'y'], ['S1'], {'miR-M2': 'M2'})
    assert blocked == {'x': {'miR-M2': 90000.0}, 'y': {'miR-M2': 150000.0}}
    # a method without reads in one sample: that sample's Spearman is left out, not NaN in the median
    empty = common + common_rows('x', 'S2', {'M1': 0}) + common_rows('y', 'S2', b)
    assert biological.agreement(empty, ['x', 'y'], ['S2', 'S1'])[('x', 'y')]['spearman_either'] == pair['spearman_either']


def test_identical_sequences_count_as_one_and_unplaced_reads_are_reported():
    common = [('a', 'S1', 'hsa-miR-199a-3p/hsa-miR-199b-3p', 'MIMAT0000232', 100), ('b', 'S1', 'miR-199b-3p', 'MIMAT0004563', 60),
              ('b', 'S1', 'miR-199a-2-3p', '', 40), ('b', 'S1', 'miR-16-5p', 'MIMAT0000069', 60)]
    seqs = {'MIMAT0000232': 'ACAGUAGUCUGCACAUUGGUUA', 'MIMAT0004563': 'ACAGUAGUCUGCACAUUGGUUA', 'MIMAT0000069': 'UAGCAGCACGUAAAUAUUGGCG'}
    collapsed, groups = biological.collapse(common, seqs)
    assert groups == [['MIMAT0000232', 'MIMAT0004563']] and {r[3] for r in collapsed} == {'MIMAT0000232', 'MIMAT0000069', ''}
    pair = biological.agreement(collapsed + [('a', 'S1', 'x', 'MIMAT0000069', 50)], ['a', 'b'], ['S1'])[('a', 'b')]
    assert pair['only_a'] == 0 and pair['only_b'] == 0
    assert biological.unplaced(common, ['a', 'b'], ['S1']) == {'a': {'share': 0.0, 'names': []},
                                                               'b': {'share': 0.25, 'names': [('miR-199a-2-3p', 40.0)]}}


def test_filtering_through_the_cutadapt_passes_and_their_chain(tmp_path):
    paths = write_cutadapt_reports(tmp_path, 'S1', 100)
    assert biological.filtering(paths) == {'reads_in': 100, 'with_3p_adapter': 95, 'discarded_5p_adapter': 1, 'too_short': 7,
                                           'too_many_n': 0, 'too_many_expected_errors': 2, 'discarded_3p_artefacts': 2,
                                           'discarded_reagent_oligos': 1, 'kept': 87}
    report = json.loads(paths['artefacts'].read_text())
    report['read_counts']['input'] = 80
    paths['artefacts'].write_text(json.dumps(report))
    with pytest.raises(ValueError, match='the artefacts pass read 80 reads, the filter pass wrote 90'):
        biological.filtering(paths)


def test_first_base_counts_substitutions_with_their_quality(tmp_path):
    reads = [('TAGCAGCACGTAAATATTGGCG', 'I'), ('AAGCAGCACGTAAATATTGGCG', '/'), ('TAGCAGCACGTAAATATTGGCGA', 'I'),
             ('CCCCCCCCCCCCCCCCCCCCCC', 'I')]
    fq = write(tmp_path / 'r.fastq', ''.join('@r{}\n{}\n+\n{}{}\n'.format(i, s, q, 'I' * (len(s) - 1)) for i, (s, q) in enumerate(reads)))
    counts, quality = biological.first_base(fq, {'hsa-miR-16-5p': 'TAGCAGCACGTAAATATTGGCG'})
    assert counts['hsa-miR-16-5p'] == {'T': 2, 'A': 1}
    assert quality == {'canonical': 40, 'substituted': 14}
    # miRNAs that differ only at the first base: each keeps its own reads, a third first base is not placed
    counts, quality = biological.first_base(fq, {'x': 'TAGCAGCACGTAAATATTGGCG', 'y': 'AAGCAGCACGTAAATATTGGCG', 'z': 'CCCCCCCCCCCCCCCCCCCCCC'})
    assert counts == {'x': {'T': 2}, 'y': {'A': 1}, 'z': {'C': 1}} and quality['substituted'] == ''


def test_lane_replicates_separate_poisson_from_overdispersion():
    rng = random.Random(5)
    def poisson(lam):  # Knuth for small lambdas, normal approximation above
        if lam > 30:
            return max(0, round(rng.gauss(lam, lam ** 0.5)))
        k, p, limit = 0, 1.0, pow(2.718281828, -lam)
        while True:
            p *= rng.random()
            if p <= limit:
                return k
            k += 1
    mirnas = {'M{}'.format(i): 20 + 30 * i for i in range(40)}
    common = []
    for lane in range(1, 5):
        common += common_rows('poisson', 'S1_L{}'.format(lane), {m: poisson(mu) for m, mu in mirnas.items()})
        common += common_rows('noisy', 'S1_L{}'.format(lane), {m: poisson(mu * rng.choice((0.5, 1.5))) for m, mu in mirnas.items()})
    rep = biological.lane_replicates(common, ['poisson', 'noisy'], {'S1': ['S1_L1', 'S1_L2', 'S1_L3', 'S1_L4']})
    assert 0.6 < rep['poisson']['dispersion_mean'] < 1.6 and rep['noisy']['dispersion_mean'] > 5
    assert rep['poisson']['expected_median'] == 0.794            # chi2(3)/3 for 4 lanes
    assert rep['poisson']['expected_overdispersed'] == round(40 * 0.0296)
    assert sum(n for _, n in rep['noisy']['overdispersed_by_mean'].values()) == 40
    median, above = biological.chi2_expectations([1, 3])         # 2 and 4 lanes: between the two medians
    assert biological.chi2_expectations([1])[0] < median < biological.chi2_expectations([3])[0] and 0.05 < above < 0.15
    assert rep['poisson']['overdispersed'] < rep['noisy']['overdispersed']


def test_fill_targets_for_a_biological_declaration(tmp_path):
    spec = {'name': 'b', 'kind': 'biological', 'data': 'benchmark/smallrna/b',
            'mirbase': {'release': '23.1', 'url': 'u', 'sha256': 's'},
            'samples': [{'sample': 'S1', 'runs': [{'run': 'R1', 'url': 'r1', 'md5': 'm1'}, {'run': 'R2', 'url': 'r2', 'md5': 'm2'}]}]}
    paths = [str(t[3].relative_to(tmp_path)) for t in fill.targets(spec, tmp_path)]
    assert paths == ['benchmark/smallrna/b/mirbase/23.1/mature.fa', 'benchmark/smallrna/b/fastq/R1.fastq.gz',
                     'benchmark/smallrna/b/fastq/R2.fastq.gz']


def release_declaration(tmp_path, **changes):
    """A pool declaration with miRBase 23.1 and, under mirbase_releases, 22 and an incomplete 21."""
    spec = {'name': 'p', 'kind': 'pool', 'pool': 'ratiometric', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
            'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/p', 'composition': {'file': 'ratio.csv', 'url': 'u', 'sha256': 's'},
            'mirbase': {'release': '23.1', 'url': 'u23', 'sha256': 's23'},
            'mirbase_releases': {'22': {'url': 'u22', 'sha256': 's22'}, 21: {'url': 'u21'}},
            'runs': [{'run': 'SRRA', 'group': 'A', 'md5': 'm', 'url': 'u'}, {'run': 'SRRB', 'group': 'B', 'md5': 'm', 'url': 'u'}]}
    spec.update(changes)
    return write(tmp_path / 'ds.yaml', yaml.safe_dump(spec)), spec


def test_load_dataset_selects_a_declared_mirbase_release(tmp_path):
    path, spec = release_declaration(tmp_path)
    assert benchmark.load_dataset(path)['mirbase'] == spec['mirbase'] == benchmark.load_dataset(path, '23.1')['mirbase']
    assert benchmark.load_dataset(path, '22')['mirbase'] == {'url': 'u22', 'sha256': 's22', 'release': '22'}
    with pytest.raises(SystemExit, match='miRBase 21 is not declared with url and sha256; declared: 23.1, 21, 22'):
        benchmark.load_dataset(path, '21')  # no sha256
    with pytest.raises(SystemExit, match='miRBase 20 is not declared'):
        benchmark.load_dataset(path, '20')
    genome = {'name': 'g', 'kind': 'synthetic', 'organism': 'homo_sapiens', 'libprepkit': 'k', 'methods': ['unitas'], 'adapter': 'A',
              'data': 'd', 'synthetic': {'mode': 'genome'},
              'mirbase': {'release': '23.1', 'url': 'u', 'sha256': 's', 'hairpin': {'url': 'h', 'sha256': 's'}, 'gff': {'url': 'g', 'sha256': 's'}},
              'mirbase_releases': {'22': {'url': 'u22', 'sha256': 's22'}}}
    path = write(tmp_path / 'genome.yaml', yaml.safe_dump(genome))
    with pytest.raises(SystemExit, match='a genome-wide synthetic dataset needs mirbase.hairpin and mirbase.gff'):
        benchmark.load_dataset(path, '22')


def test_prepare_with_another_mirbase_release_configures_it_for_every_method(tmp_path):
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/p'
    write(data / 'composition/ratio.csv', RATIO)
    write(data / 'mirbase/22/mature.fa', HSA)
    for run in ('SRRA', 'SRRB'):
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/{}.fastq.gz'.format(run), 'wt') as fh:
            fh.write('@r\nTGGAATGTAAAGAAGTATGTATTGGAATTCTC\n+\n{}\n'.format('I' * 32))
    unitas_dump(ext, 'homo_sapiens', 'Human', {})
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--project', str(tmp_path / 'project'), '--ext-dir', str(ext),
           '--mirbase-release', '22']
    dataset, _ = release_declaration(tmp_path)
    result = subprocess.run(cmd + ['--dataset', str(dataset)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((tmp_path / 'project/config.yaml').read_text())
    assert config['db']['mirbase'] == {'release': '22'}
    saved = yaml.safe_load((tmp_path / 'project/benchmark/dataset.yaml').read_text())
    assert saved['mirbase'] == {'url': 'u22', 'sha256': 's22', 'release': '22'} and 'mirbase_releases' not in saved
    dataset, _ = release_declaration(tmp_path, config={'db': {'mirbase': {'release': '23.1'}}})
    result = subprocess.run(cmd[:4] + [str(tmp_path / 'other')] + cmd[5:] + ['--dataset', str(dataset)], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and "the declaration's config sets db.mirbase.release 23.1, but the run uses miRBase 22" in result.stderr


def test_a_reduced_synthetic_dataset_keeps_its_reference_s_mirbase_release(tmp_path):
    reduced = {'name': 's', 'kind': 'synthetic', 'organism': 'mus_musculus', 'libprepkit': 'k', 'methods': ['unitas'], 'adapter': 'A',
               'synthetic': {'samples': 1}, 'reference': {'name': 'r', 'dir': 'd'},
               'mirbase': {'release': '23.1', 'url': 'u', 'sha256': 's'}, 'mirbase_releases': {'22': {'url': 'u22', 'sha256': 's22'}}}
    path = write(tmp_path / 'reduced.yaml', yaml.safe_dump(reduced))
    assert benchmark.load_dataset(path, '23.1')['mirbase']['release'] == '23.1'
    with pytest.raises(SystemExit, match="a reduced synthetic dataset's reference fixes its miRBase release"):
        benchmark.load_dataset(path, '22')


def test_evaluate_takes_the_unitas_dump_the_run_recorded_not_the_configured_one(tmp_path):
    project = tmp_path / 'project'
    quant = project / 'data/tmp/smallrna/quant'
    write(project / 'config.yaml', yaml.safe_dump({'db': {'unitas': {'release': 'another'}}}))
    with pytest.raises(SystemExit, match='no UNITAS reference record'):
        benchmark.unitas_refdump(project)
    write_records(quant, ['unitas'])
    with pytest.raises(SystemExit, match='has no dump row: run quant_all again'):
        benchmark.unitas_refdump(project)
    recorded = unitas_dump(tmp_path / 'ext', 'homo_sapiens', 'Human', {})
    write_records(quant, ['unitas'], dump=recorded)
    assert benchmark.unitas_refdump(project) == recorded  # not db.unitas.release 'another'
    write_records(quant, ['unitas'], dump=tmp_path / 'moved')
    with pytest.raises(SystemExit, match='no UNITAS reference dump at .*moved .*name one with --unitas-ref'):
        benchmark.unitas_refdump(project)
    write_records(quant, ['unitas'], database='mirgenedb', release='3.0', dump=recorded)
    with pytest.raises(SystemExit, match='the UNITAS dump was built from mirgenedb: the benchmark matches miRBase names'):
        benchmark.unitas_refdump(project)


def test_fill_targets_the_selected_mirbase_release(tmp_path):
    path, _ = release_declaration(tmp_path)
    paths = [str(t[3].relative_to(tmp_path)) for t in fill.targets(benchmark.load_dataset(path, '22'), tmp_path)]
    assert 'benchmark/smallrna/p/mirbase/22/mature.fa' in paths and not any('/23.1/' in p for p in paths)


def test_shipped_declarations_list_loadable_mirbase_releases():
    """Every release a declaration lists loads, with files from that release's miRBase folder (#431)."""
    listed = 0
    for path in sorted((BENCH / 'datasets').glob('*.yaml')):
        for release in (yaml.safe_load(path.read_text()).get('mirbase_releases') or {}):
            selected = benchmark.load_dataset(path, release)['mirbase']
            for files in [selected] + [selected[k] for k in ('hairpin', 'gff') if k in selected]:
                assert '/{}/'.format(release) in files['url'] and len(files['sha256']) == 64, (path.name, files)
            listed += 1
    assert listed


def write_cutadapt_reports(folder, sample, reads):
    """The cutadapt route's three JSON reports for a sample, chained: each pass reads what the one before wrote."""
    counts = {'adapter': {'input': reads, 'output': reads, 'read1_with_adapter': reads - 5, 'filtered': {}},
              'filter': {'input': reads, 'output': reads - 10,
                         'filtered': {'discard_trimmed': 1, 'too_short': 7, 'too_many_n': None, 'too_many_expected_errors': 2}},
              'artefacts': {'input': reads - 10, 'output': reads - 13, 'filtered': {'discard_trimmed': 3}}}
    paths = {}
    for name, c in counts.items():
        records = [{'name': 'three_prime_artefact_adapter_0_33', 'total_matches': 2}, {'name': 'reagent_oligo_1_22', 'total_matches': 1}]
        report = {'tag': 'Cutadapt report', 'command_line_arguments': [name], 'input': {'path1': 'in.fastq', 'path2': None, 'paired': False},
                  'read_counts': c, 'basepair_counts': {'input': c['input'] * 50, 'output': c['output'] * 20},
                  'adapters_read1': records if name == 'artefacts' else []}
        paths[name] = write(Path(folder) / '{}.{}.json'.format(sample, name), json.dumps(report))
    return paths


def write_benchmark(path, seconds, cpu, rss):
    """A Snakemake benchmark file (one run): wall time, peak RSS and CPU time; 'NA' as Snakemake writes it."""
    return write(path, 's\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n'
                 '{}\t0:00:{}\t{}\t1\t1\t1\t0\t0\t50\t{}\n'.format(seconds, seconds, rss, cpu))


def test_resources_per_method_from_the_benchmark_files(tmp_path):
    """Per method: quantification jobs summed, the largest job's memory, reference steps measured in this run apart;
    a method with a quantification step's file missing is not measured; a malformed file fails (#386)."""
    quant = tmp_path / 'quant'
    for sample in ('S1', 'S2'):
        for step, seconds in (('smrnaseq_collapse', 2), ('smrnaseq_align', 10), ('smrnaseq_mirtop', 20)):
            write_benchmark(quant / 'benchmarks' / step / '{}.tsv'.format(sample), seconds, seconds * 3, 100 * seconds)
        write_benchmark(quant / 'benchmarks/tesmall_quant/{}.tsv'.format(sample), 600, 590, 1500)
    write_benchmark(quant / 'benchmarks/smrnaseq_counts.tsv', 4, 'NA', 'NA')
    write_benchmark(quant / 'benchmarks/smrnaseq_hairpin_index.tsv', 3, 9, 50)
    write_benchmark(quant / 'benchmarks/tesmall_genome.tsv', 3600, 50000, 9000)   # a shared build's name: not a step here
    text, lines = resources.summarize(quant, ['smrnaseq', 'tesmall', 'tinyrna'], ['S1', 'S2'], 1000)
    table = '\n'.join(lines)
    assert '| smrnaseq | 7 | 68 | 192 | 2000 (smrnaseq_mirtop) | 14.7 |' in table
    assert '| tesmall | not measured |' in table and 'tesmall (1 missing, e.g. `tesmall_counts.tsv`)' in table
    assert 'tinyrna (7 missing, e.g. `tinyrna_collapse/S1.tsv`)' in table
    assert '| smrnaseq | smrnaseq_hairpin_index | 3 | 9 | 50 |' in table
    assert 'tesmall_genome' not in table and 'tesmall_genome' not in text
    assert text.splitlines()[0] == '\t'.join(resources.HEADER) and len(text.splitlines()) == 1 + 8 + 2  # header; smrnaseq and tesmall jobs
    _, lines = resources.summarize(quant, ['smrnaseq'], ['S1', 'S2'], None)
    assert '| smrnaseq | 7 | 68 | 192 | 2000 (smrnaseq_mirtop) |  |' in lines
    for bad in ('s\tcpu_time\n1\t2\n', 's\tmax_rss\tcpu_time\n1\t2\t3\n4\t5\t6\n', 's\tmax_rss\tcpu_time\nNA\t2\t3\n',
                's\tmax_rss\tcpu_time\n1\t-2\t3\n', 's\tmax_rss\tcpu_time\n1\tx\t3\n'):
        write(quant / 'benchmarks/smrnaseq_counts.tsv', bad)
        with pytest.raises(ValueError, match='smrnaseq_counts.tsv'):
            resources.summarize(quant, ['smrnaseq'], ['S1', 'S2'], None)


def biological_dataset(tmp_path, lanes=(1, 2)):
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/bio'
    write(data / 'mirbase/23.1/mature.fa', BLOOD)
    reads = ['TAGCAGCACGTAAATATTGGCG'] * 6 + ['AAGCAGCACGTAAATATTGGCG'] + ['TAGCTTATCAGACTGATGTTGA'] * 3
    runs = []
    for lane in lanes:
        run = 'SRR{}'.format(lane)
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/{}.fastq.gz'.format(run), 'wt') as fh:
            fh.write(''.join('@{}.{} NB1:1:FC:{}:11101:1:{}/1\n{}TGGAATTCTC\n+\n{}\n'.format(run, i, lane, i, s, 'I' * (len(s) + 10))
                             for i, s in enumerate(reads)))
        runs.append({'run': run, 'lane': lane, 'reads': len(reads), 'md5': 'unused', 'url': 'unused'})
    dataset = write(tmp_path / 'bio.yaml', yaml.safe_dump({
        'name': 'bio', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
        'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/bio', 'config': {'blockers': ['blood miRNA']},
        'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'samples': [{'sample': 'GSM1', 'group': 'control', 'runs': runs}]}))
    return ext, dataset, reads


def test_prepare_merges_runs_or_keeps_lanes_and_checks_lanes(tmp_path):
    ext, dataset, _ = biological_dataset(tmp_path)
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--ext-dir', str(ext)]
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'merged'), '--reads', '4'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((tmp_path / 'merged/config.yaml').read_text())
    assert config['samples']['GSM1']['R1'] == 'bio/SRR1_R1.fastq.gz,bio/SRR2_R1.fastq.gz' and config['blockers'] == ['blood miRNA']
    assert 'GSM1,"bio/SRR1_R1.fastq.gz,bio/SRR2_R1.fastq.gz",0' in (tmp_path / 'merged/pep/sample_table.csv').read_text()
    assert (tmp_path / 'merged/benchmark/samples.tsv').read_text().splitlines()[1] == 'GSM1\tGSM1\tcontrol\tmerged\tSRR1,SRR2\t8'
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'lanes'), '--lanes'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert sorted(yaml.safe_load((tmp_path / 'lanes/config.yaml').read_text())['samples']) == ['GSM1_L1', 'GSM1_L2']
    spec = yaml.safe_load(dataset.read_text())
    spec['samples'][0]['runs'][1]['lane'] = 3
    write(dataset, yaml.safe_dump(spec))
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'wrong')], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'SRR2: declared lane 3, but its reads are from lane 2' in result.stderr
    assert not (tmp_path / 'wrong').exists()                                          # checked before anything is written
    spec['samples'][0]['runs'][1]['lane'] = 1                                         # two runs in lane 1 (two flow cells)
    with gzip.open(ext / 'benchmark/smallrna/bio/fastq/SRR2.fastq.gz', 'wt') as fh:   # native Illumina read names
        fh.write('@NB1:1:FC2:1:11101:1:1 1:N:0:ACGT\nTAGCAGCACGTAAATATTGGCGTGGAATTCTC\n+\n{}\n'.format('I' * 32))
    write(dataset, yaml.safe_dump(spec))
    result = subprocess.run(cmd + ['--project', str(tmp_path / 'twice'), '--lanes'], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'two runs of GSM1 in lane 1' in result.stderr and not (tmp_path / 'twice').exists()
    for runs, message in (([{'run': 'SRR1', 'lane': '1'}], "lane must be a whole number"),
                          ([{'run': 'SRR1', 'lane': 1}, {'run': 'SRR1', 'lane': 2}], 'runs listed more than once: SRR1')):
        write(dataset, yaml.safe_dump(dict(spec, samples=[dict(spec['samples'][0], runs=runs)])))
        result = subprocess.run(cmd + ['--project', str(tmp_path / 'bad')], capture_output=True, text=True, timeout=60)
        assert result.returncode != 0 and message in result.stderr


def test_evaluate_a_biological_dataset_in_lanes(tmp_path):
    ext, dataset, reads = biological_dataset(tmp_path)
    project = tmp_path / 'project'
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(ext), '--lanes'],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tGSM1_L1\tGSM1_L2\nmiR-16-5p\t700\t690\nmiR-21-5p\t300\t310\nmiR-451a\t1\t0\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,GSM1_L1,GSM1_L2\nhsa-miR-16-5p,650,660\nhsa-miR-21-5p,320,300\n')
    write(quant / 'unitas/mir_counts.tsv', '\tGSM1_L1\tGSM1_L2\nmiR-16-5p\t700\t690\nmiR-21-5p\t300\t310\nmiR-451a\t1\t0\nmiR-16-1-5p\t30\t40\n'
          'miR-325-3p\t5\t0\n')  # miR-16-1-5p: hsa-miR-16-5p on the hairpin's 5' arm; miR-325-3p: no human mature on that arm
    isomir_outputs(quant)
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump_dir(ext, 'homo_sapiens'))
    unitas_dump(ext, 'homo_sapiens', 'Human', {
        'miR-16-1': 'GTCAGCAGTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAAGGTTGAC',
        'miR-325': 'ATACAGTGCTTGGTTCCTAGTAGGTGTCCAGTAAGTGTTTGTGACATAATTTGTTTATTGAGGACCTCCTATCAATCAAGCACTGTGCTAGGCTCTGG'})
    filtered = project / 'data/tmp/smallrna/filter'
    for lane in ('GSM1_L1', 'GSM1_L2'):
        write(filtered / 'cutadapt/{}_R1.fastq'.format(lane),
              ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, ('/' if s[0] == 'A' else 'I') + 'I' * (len(s) - 1)) for i, s in enumerate(reads)))
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'cutadapt reports missing' in result.stderr and not (project / 'benchmark/report.md').exists()
    for lane in ('GSM1_L1', 'GSM1_L2'):
        write_cutadapt_reports(filtered / 'cutadapt_stages', lane, 100)
    for lane, seconds in (('GSM1_L1', 30), ('GSM1_L2', 50)):   # UNITAS measured; miRge3 run before #386
        write_benchmark(quant / 'benchmarks/unitas/{}.tsv'.format(lane), seconds, 2 * seconds, 900 + seconds)
    write_benchmark(quant / 'benchmarks/unitas_mirtable.tsv', 1, 'NA', 'NA')
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    report = (project / 'benchmark/report.md').read_text()
    assert '2 lane samples' in report and '| unitas / mirge3 |' in report
    assert '| unitas | 3 | 81 | 160 | 950 (unitas) | 2.47 |' in report    # 200 input reads over 81 s
    assert '| mirge3 | not measured |' in report and 'mirge3 (2 missing, e.g. `mirge3_quant.tsv`)' in report
    record = json.loads((project / 'benchmark/evaluation.json').read_text())
    published = {name: (project / 'benchmark' / name).read_bytes() for name in ('common.tsv', 'report.md', 'evaluation.json')}
    # a malformed benchmark file stops a new evaluation before any scoring table is rewritten
    write(quant / 'unitas/mir_counts.tsv', '\tGSM1_L1\tGSM1_L2\nmiR-16-5p\t1\t2\n')
    write(quant / 'benchmarks/unitas_mirtable.tsv', 's\tmax_rss\tcpu_time\n1\tx\t1\n')
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'unitas_mirtable.tsv: max_rss is' in result.stderr
    assert {name: (project / 'benchmark' / name).read_bytes() for name in published} == published
    assert record['resources.tsv'] and (project / 'benchmark/resources.tsv').read_text().splitlines()[1:] == [
        'unitas\tquant\tunitas\tGSM1_L1\t30\t60\t930', 'unitas\tquant\tunitas\tGSM1_L2\t50\t100\t950',
        'unitas\tquant\tunitas_mirtable\t\t1\t\t']
    assert '| unplaced | 0.0024 | 0.0 |' in report and 'unitas: `miR-325-3p` (5)' in report      # median of 5/1036 and 0/1040
    assert '| isomir_share |' in report and '## isomiRs' in report and '| unitas / mirge3 |' in report.split('## isomiRs')[1]
    ci = {(r[0], r[1]): r[2:] for r in (line.split('\t') for line in (project / 'benchmark/metrics_ci.tsv').read_text().splitlines())}
    value, low, high = ci[('isomir_share', 'unitas')]  # #459: resampled samples; the detected count has no interval
    assert low != '' and float(low) <= float(value) <= float(high) and not any(k == 'isomirs_detected' for k, _ in ci)
    assert '| miR-451a | 485.0 | 0.0 |' in report     # median CPM of the blocked miRNA: 1 of 1031 placed reads in L1, none in L2
    assert '| GSM1_L1 | 100 | 95 | 1 | 7 | 0 | 2 | 2 | 1 | 87 |' in report                   # the filtering table
    assert '| GSM1_L1 | 10 | 10.00% | 40.0 | 14.0 |' in report                       # 1 of 10 top-miRNA reads has another first base
    assert '## Lane replicates' in report and '| mirna_sample_pairs | 2 | 2 |' in report and '| 100-999 | 0.0% of 2 | 0.0% of 2 |' in report
    assert 'resamples (seed 1) of the libraries, each with all its lanes' in report   # lanes of a library are not independent


def adapter_report(project, sample, count):
    """The first cutadapt pass's JSON report of a sample, as far as the summary reads it."""
    write(project / 'data/tmp/smallrna/filter/cutadapt_stages' / f'{sample}.adapter.json', json.dumps({'read_counts': {'input': count}}))


def summary_project(root, name, spec, metrics, samples=None, truth=None, published=True):
    """An evaluated project as the summary reads it: config.yaml, benchmark/dataset.yaml, metrics.tsv, samples.tsv
    (truth.tsv for synthetic data), report.md and each sample's first cutadapt report with its `reads` (the truth
    total for synthetic data; none when empty), published with evaluation.json as `evaluate` does (an older evaluation
    without it when not `published`)."""
    project = root / name
    names = [row[0] for row in samples[1:]] if samples else ['S01', 'S02']
    runs = {name: [name] for name in names}
    if samples and 'runs' in samples[0]:
        runs = {row[0]: row[samples[0].index('runs')].split(',') for row in samples[1:]}
    inputs = {s: {'R1': ','.join(f"{spec['name']}/{run}_R1.fastq.gz" for run in ids)} for s, ids in runs.items()}
    write(project / 'config.yaml', yaml.safe_dump({'workflow': 'smallrna', 'fastq_dir': 'data/raw/fastq', 'samples': inputs}))
    write(project / 'benchmark/dataset.yaml', yaml.safe_dump(dict(spec, organism='homo_sapiens')))
    write(project / 'benchmark/metrics.tsv', ''.join('\t'.join(row) + '\n' for row in metrics))
    report_rows = ['| ' + ' | '.join(row) + ' |' for row in metrics]
    report_rows.insert(1, '|---' * len(metrics[0]) + '|')
    write(project / 'benchmark/report.md', '# {}\n\n{}\n'.format(spec['name'], '\n'.join(report_rows)))
    if published:
        import evaluation
    if truth:
        write(project / 'benchmark/truth.tsv', 'sample\tclass\tmature_name\tmature_id\tsequence\tcount\n' + ''.join(
            '{}\tcanonical\tm\tMIMAT1\tACGT\t{}\n'.format(s, n) for s, n in truth))
        for sample, count in truth:
            adapter_report(project, sample, count)
    else:
        write(project / 'benchmark/samples.tsv', ''.join('\t'.join(row) + '\n' for row in samples))
        for row in samples[1:]:
            if row[samples[0].index('reads')]:
                adapter_report(project, row[0], int(row[samples[0].index('reads')]))
    if published:
        evaluation.publish(project / 'benchmark', metrics[0], metrics[1:], ['# {}'.format(spec['name'])],
                           evaluation.reads_processed(project))
    return project


POOL_SAMPLES = ['sample', 'group', 'lab', 'reads']


def summary_projects(root, published=True):
    kit = {'libprepkit': 'Kit A'}

    def project(*args, **kwargs):
        return summary_project(*args, published=published, **kwargs)
    return [
        project(root, 'syn', dict(kit, name='synthetic_set', kind='synthetic', methods=['unitas', 'tesmall']),
                        [['metric', 'unitas', 'tesmall'], ['recovered', '1.0', '0.9'], ['exact', '1.0', ''], ['misassigned', '0', '3']],
                        truth=[('S01', 60), ('S02', 40)]),
        project(root, 'eq3', dict(kit, name='eq_pool', kind='pool', pool='equimolar', methods=['unitas', 'mirge3', 'mirflowz']),
                        [['metric', 'oracle', 'unitas', 'mirge3', 'mirflowz'], ['detected', '0.9', '0.97', '0.95', '0.96'],
                         ['within_2fold', '', '0.77', '0.76', '0.74'], ['shifted_within_2fold', '', '0.79', '0.75', '0.79']],
                        samples=[POOL_SAMPLES, ['R1', 'EQ', 'Lab1', '500000'], ['R2', 'EQ', 'Lab2', '1500000']]),
        project(root, 'eq2', dict(libprepkit='Kit B', name='eq_pool_b', kind='pool', pool='equimolar', methods=['mirge3', 'unitas']),
                        [['metric', 'oracle', 'mirge3', 'unitas'], ['detected', '0.89', '0.89', '0.94'], ['within_2fold', '', '0.85', '0.86']],
                        samples=[POOL_SAMPLES, ['R1', 'EQ', 'Lab1', '1000'], ['R2', 'EQ', 'Lab1', '']]),
        project(root, 'ratio', dict(kit, name='ratio_pool', kind='pool', pool='ratiometric', methods=['unitas']),
                        [['metric', 'oracle', 'unitas'], ['called', '0.88', '0.89'], ['false_calls', '0.05', '0.06']],
                        samples=[POOL_SAMPLES, ['A1', 'A', 'Lab1', '10'], ['B1', 'B', 'Lab1', '10']]),
        project(root, 'bio_lanes', dict(kit, name='plasma', kind='biological', methods=['unitas', 'mirge3']),
                        [['metric', 'unitas', 'mirge3'], ['mirna_reads', '742128', '714978'], ['detected', '349.0', '373.5']],
                        samples=[['sample', 'library', 'group', 'lane', 'runs', 'reads'], ['G1_L1', 'G1', 'g', '1', 'SRR1', '100'],
                                 ['G1_L2', 'G1', 'g', '2', 'SRR2', '100'], ['G2_L1', 'G2', 'g', '1', 'SRR3', '100']]),
        project(root, 'bio_merged', dict(kit, name='plasma', kind='biological', methods=['unitas', 'mirge3', 'mirflowz']),
                        [['metric', 'unitas', 'mirge3', 'mirflowz'], ['mirna_reads', '2966861', '2858766', '2861471'], ['detected', '581.5', '640.5', '435.5']],
                        samples=[['sample', 'library', 'group', 'lane', 'runs', 'reads'], ['G1', 'G1', 'g', 'merged', 'SRR1,SRR2', '200'],
                                 ['G2', 'G2', 'g', 'merged', 'SRR3', '100']]),
    ]


def test_summary_across_kinds_and_method_sets(tmp_path):
    projects = summary_projects(tmp_path / 'projects')
    before = sorted((p, p.stat().st_mtime_ns) for project in projects for p in project.rglob('*'))
    output = tmp_path / 'out/summary.md'
    output.parent.mkdir()
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects'] + [str(p) for p in projects] +
                            ['--output', str(output)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert sorted((p, p.stat().st_mtime_ns) for project in projects for p in project.rglob('*')) == before   # nothing written into the projects
    text = output.read_text()
    assert [line for line in text.splitlines() if line.startswith('## ')] == [
        '## Synthetic data', '## Equimolar pools', '## Ratiometric pools', '## Biological samples']
    fixed = '| project | dataset | kit | samples | reads processed | methods |'
    # synthetic: a metric no project has (unplaced) is left out, one a method lacks is an empty cell
    assert fixed + ' recovered: unitas | recovered: tesmall | exact: unitas | misassigned: unitas | misassigned: tesmall |' in text
    assert '| [syn](../projects/syn/benchmark/report.md) | synthetic_set | Kit A | 2 | 100 | unitas, tesmall | 1.0 | 0.9 | 1.0 | 0 | 3 |' in text
    # equimolar: methods in order of appearance, the oracle only where it has values, missing methods and reads empty
    assert (fixed + ' detected: oracle | detected: unitas | detected: mirge3 | detected: mirflowz | within_2fold: unitas | within_2fold: mirge3 | '
            'within_2fold: mirflowz | shifted_within_2fold: unitas | shifted_within_2fold: mirge3 | shifted_within_2fold: mirflowz |') in text
    assert '| eq_pool | Kit A | 2 | 2,000,000 | unitas, mirge3, mirflowz | 0.9 | 0.97 | 0.95 | 0.96 | 0.77 | 0.76 | 0.74 | 0.79 | 0.75 | 0.79 |' in text
    assert '| eq_pool_b | Kit B | 2 |  | mirge3, unitas | 0.89 | 0.94 | 0.89 |  | 0.86 | 0.85 |  |  |  |  |' in text
    assert '| ratio_pool | Kit A | 2 | 20 | unitas | 0.88 | 0.89 | 0.05 | 0.06 |' in text
    assert '| plasma | Kit A | 3 lanes of 2 samples | 300 | unitas, mirge3 | 742128 | 714978 |  | 349.0 | 373.5 |  |' in text
    assert '| plasma | Kit A | 2 | 300 | unitas, mirge3, mirflowz | 2966861 | 2858766 | 2861471 | 581.5 | 640.5 | 435.5 |' in text


def test_summary_html_is_self_contained_and_bad_inputs_fail(tmp_path):
    import summary
    projects = summary_projects(tmp_path)
    output = tmp_path / 'summary.html'
    summary.write(projects, output)
    page = output.read_text()
    assert page.startswith('<!doctype html>') and '<th colspan="4">detected</th>' in page and '<a href="eq3/benchmark/report.md">eq3</a>' in page
    assert not any(token in page for token in ('http', 'src=', '<link', '@import'))
    assert '<td class="num">0.97</td><td class="num">0.95</td><td class="num">0.96</td>' in page
    (projects[1] / 'benchmark/metrics.tsv').unlink()
    for args, message in (((projects, tmp_path / 'summary.txt'), 'must be .md or .html'),
                          ((projects, projects[0] / 'summary.md'), 'inside the project'),
                          ((projects[:1] + projects[:1], output), 'projects named alike: syn'),
                          ((projects[1:2], output), 'run `benchmark.py evaluate')):
        with pytest.raises(ValueError, match=message):
            summary.write(*args)


def test_summary_requires_matching_report_and_metrics_for_legacy_and_new_evaluations(tmp_path):
    import evaluation
    import summary
    project = summary_projects(tmp_path, published=False)[0]
    bench = project / 'benchmark'
    output = tmp_path / 'overview.md'
    summary.write([project], output)  # Pre-marker report with matching metrics.
    (bench / 'report.md').write_text('unrelated evaluation')
    with pytest.raises(ValueError, match='incomplete or inconsistent'):
        summary.write([project], output)
    evaluation.publish(bench, ['metric', 'unitas'], [['exact', '1']], ['# Report'])
    summary.write([project], output)
    (bench / 'metrics.tsv').write_text('metric\tunitas\nexact\t0\n')
    with pytest.raises(ValueError, match='incomplete or inconsistent'):
        summary.write([project], output)
    (bench / 'report.md').unlink()
    with pytest.raises(ValueError, match='incomplete or inconsistent'):
        summary.write([project], output)


def test_summary_preserves_project_files_and_resolves_relative_report_link(tmp_path):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    renamed = project.rename(project.with_name('A|B]'))
    report = renamed / 'benchmark/report.md'
    output = tmp_path / 'summary.md'
    output.hardlink_to(report)
    original = report.read_bytes()
    through = tmp_path / 'link'
    through.symlink_to(renamed / 'benchmark', target_is_directory=True)
    summary.write([through / '..'], output)
    assert report.read_bytes() == original
    assert '[A&#124;B&#93;](projects/A%7CB%5D/benchmark/report.md)' in output.read_text()
    assert not os.path.samefile(output, report)


def test_prepare_actual_short_run_and_unknown_merged_counts(tmp_path):
    import benchmark
    ext, dataset, reads = biological_dataset(tmp_path)
    spec = yaml.safe_load(dataset.read_text())
    del spec['samples'][0]['runs'][0]['reads']
    write(dataset, yaml.safe_dump(spec))
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--ext-dir', str(ext)]
    for name, args, expected in [('capped', ['--reads', '100'], str(len(reads) * 2)), ('linked', [], '')]:
        project = tmp_path / name
        subprocess.run(cmd + ['--project', str(project), *args], check=True, capture_output=True, timeout=60)
        row, = benchmark.read_tsv(project / 'benchmark/samples.tsv')
        assert row['reads'] == expected


def test_summary_reads_processed_are_those_recorded_by_the_evaluation(tmp_path):
    """Not samples.tsv, and not the reports as they are now: a rerun without a new evaluation keeps the recorded total;
    a sample without a report gives no total; a report without a whole read count stops `evaluate`; an evaluation
    older than the record has no total."""
    import evaluation
    import summary
    project = summary_projects(tmp_path)[1]   # eq3: 500,000 and 1,500,000 reads
    table = project / 'benchmark/samples.tsv'
    table.write_text(table.read_text().replace('500000', '7'))
    assert summary.read_project(project, tmp_path)[1]['reads processed'] == '2,000,000'
    adapter_report(project, 'R1', 1000)                                   # rerun at another depth, not yet evaluated
    assert summary.read_project(project, tmp_path)[1]['reads processed'] == '2,000,000'
    assert evaluation.reads_processed(project) == 1501000
    (project / 'data/tmp/smallrna/filter/cutadapt_stages/R2.adapter.json').unlink()
    assert evaluation.reads_processed(project) is None
    for broken in ('{"read_counts": {}}', '{"read_counts": {"input": 2.5}}', '[]', 'not json'):
        write(project / 'data/tmp/smallrna/filter/cutadapt_stages/R2.adapter.json', broken)
        with pytest.raises(ValueError, match='R2.adapter.json'):
            evaluation.reads_processed(project)
    legacy = summary_project(tmp_path / 'old', 'eq', dict(libprepkit='Kit A', name='eq_pool', kind='pool', pool='equimolar',
                                                          methods=['unitas']), [['metric', 'unitas'], ['detected', '0.9']],
                             samples=[POOL_SAMPLES, ['R1', 'EQ', 'Lab1', '10']], published=False)
    assert summary.read_project(legacy, tmp_path)[1]['reads processed'] == ''
    record = json.loads((project / 'benchmark/evaluation.json').read_text())
    for bad in (-1, 2.5, '10'):
        (project / 'benchmark/evaluation.json').write_text(json.dumps(dict(record, reads_processed=bad)))
        with pytest.raises(ValueError, match='incomplete or inconsistent evaluation'):
            summary.read_project(project, tmp_path)


def test_short_and_incomplete_pool_fastqs(tmp_path):
    import benchmark
    records = '@read\nACGT\n+\nIIII\n'
    source = tmp_path / 'reads.fastq.gz'
    target = tmp_path / 'capped.fastq.gz'
    source.write_bytes(gzip.compress(records.encode() * 3))
    assert benchmark.head_fastq(source, target, 100) == 3
    assert gzip.decompress(target.read_bytes()) == records.encode() * 3
    source.write_bytes(gzip.compress((records + '@short\nACGT\n+\n').encode()))
    with pytest.raises(SystemExit, match='incomplete FASTQ record 2'):
        benchmark.head_fastq(source, target, 100)
    damaged = bytearray(gzip.compress(records.encode() * 3))
    damaged[-8] ^= 0xFF                                     # the CRC in the trailer, after exactly the 3 records asked for
    source.write_bytes(bytes(damaged))
    with pytest.raises(SystemExit, match='damaged gzip file'):
        benchmark.head_fastq(source, target, 3)


def test_summary_refuses_external_target_of_a_project_input(tmp_path):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    report = project / 'benchmark/report.md'
    external = report.rename(tmp_path / 'external.md')
    report.symlink_to(external)
    before = external.read_bytes()
    with pytest.raises(ValueError, match='aliases one of its inputs'):
        summary.write([project], external)
    assert external.read_bytes() == before


def test_interrupted_first_publication_needs_a_successful_retry(tmp_path, monkeypatch):
    import evaluation
    project = summary_projects(tmp_path, published=False)[0]
    bench = project / 'benchmark'
    rows = [('exact', 1)]
    lines = ['| metric | unitas |', '|---|---|', '| exact | 1 |']
    replace = Path.replace

    def interrupt(path, destination):
        if path.name == 'evaluation.json':  # Both data files have reached their destination.
            raise OSError('interrupted')
        return replace(path, destination)

    with monkeypatch.context() as patch:
        patch.setattr(Path, 'replace', interrupt)
        with pytest.raises(OSError, match='interrupted'):
            evaluation.publish(bench, ['metric', 'unitas'], rows, lines)
    with pytest.raises(ValueError, match='inconsistent evaluation'):
        evaluation.read_metrics(bench)
    (bench / 'evaluation.json').unlink()
    with pytest.raises(ValueError, match='inconsistent evaluation'):
        evaluation.read_metrics(bench)
    evaluation.publish(bench, ['metric', 'unitas'], rows, lines)
    assert evaluation.read_metrics(bench)[1] == ['exact', '1']


@pytest.mark.parametrize('problem', ['missing', 'duplicate', 'unexpected'])
def test_summary_requires_sample_table_to_cover_project(tmp_path, problem):
    """A biological project's samples.tsv (its lanes) must list the configured samples."""
    import summary
    project = summary_projects(tmp_path)[5]
    table = project / 'benchmark/samples.tsv'
    header, first, second = table.read_text().splitlines()
    rows = {'missing': [first], 'duplicate': [first, first], 'unexpected': [first, second.replace('G2', 'other')]}
    table.write_text('\n'.join([header] + rows[problem]) + '\n')
    with pytest.raises(ValueError, match='exactly one row'):
        summary.write([project], tmp_path / 'summary.md')


@pytest.mark.parametrize('bad', ['@bad\nACGT\n+\nI\n', 'bad\nA\n+\nI\n', '@bad\nA\n-\nI\n'])
def test_bad_capped_preparation_can_be_retried_at_the_same_path(tmp_path, bad):
    ext, dataset, _ = biological_dataset(tmp_path)
    spec = yaml.safe_load(dataset.read_text())
    for run in spec['samples'][0]['runs']:
        run.pop('lane')
    write(dataset, yaml.safe_dump(spec))
    source = ext / 'benchmark/smallrna/bio/fastq/SRR2.fastq.gz'
    original = source.read_bytes()
    source.write_bytes(gzip.compress(bad.encode()))
    project = tmp_path / 'prepared'
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--project', str(project),
           '--dataset', str(dataset), '--ext-dir', str(ext), '--reads', '100']
    failed = subprocess.run(cmd, capture_output=True, text=True)
    assert failed.returncode and 'malformed FASTQ record' in failed.stderr
    assert not project.exists()
    source.write_bytes(original)
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert str(project) in result.stdout and (project / 'benchmark/samples.tsv').is_file()


def test_summary_cli_prints_canonical_path_and_preserves_underscore_labels(tmp_path):
    from html import unescape
    project = summary_projects(tmp_path / 'projects')[0]
    renamed = project.rename(project.with_name('_set_'))
    target = tmp_path / 'depth/extra/output'
    target.mkdir(parents=True)
    alias = tmp_path / 'alias'
    alias.symlink_to(target, target_is_directory=True)
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(renamed), '--output', str(alias / 'page.md')]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    page = Path(result.stdout.strip())
    assert page == target / 'page.md'
    label, href = next(line.split('](', 1) for line in page.read_text().splitlines() if line.startswith('| ['))
    assert label == '| [&#95;set&#95;' and unescape(label) == '| [_set_'
    assert (page.parent / href.split(')', 1)[0]).resolve() == renamed / 'benchmark/report.md'


def test_summary_rejects_output_spelled_inside_project_through_external_directory(tmp_path):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    storage = tmp_path / 'external'
    storage.mkdir()
    (project / 'external-data').symlink_to(storage, target_is_directory=True)
    alias = tmp_path / 'project-alias'
    alias.symlink_to(project, target_is_directory=True)
    for entry in (project, alias):
        with pytest.raises(ValueError, match='inside the project'):
            summary.write([project], entry / 'external-data/overview.md')
    assert list(storage.iterdir()) == []


def test_evaluation_locks_are_regular_and_exclusive_descriptor_is_writable(tmp_path, monkeypatch):
    import evaluation
    import fcntl
    import stat
    bench = summary_projects(tmp_path, published=False)[0] / 'benchmark'   # an older evaluation, without a lock
    actual_flock = fcntl.flock
    modes = []

    def check_lock(fd, mode):
        assert stat.S_ISREG(os.fstat(fd).st_mode)
        access = fcntl.fcntl(fd, fcntl.F_GETFL) & os.O_ACCMODE
        assert access == (os.O_WRONLY if mode == fcntl.LOCK_EX else os.O_RDONLY)
        modes.append(mode)
        return actual_flock(fd, mode)

    monkeypatch.setattr(fcntl, 'flock', check_lock)
    assert evaluation.read_metrics(bench) and not modes
    assert not (bench / '.evaluation.lock').exists()
    evaluation.publish(bench, ['metric', 'unitas'], [['exact', '1']], ['# Report'])
    assert evaluation.read_metrics(bench)[1] == ['exact', '1']
    assert modes == [fcntl.LOCK_EX, fcntl.LOCK_SH]


def test_synthetic_summary_explains_metrics(tmp_path):
    import summary
    project = summary_projects(tmp_path)[0]
    output = tmp_path / 'overview.md'
    summary.write([project], output)
    assert 'smaller of observed and true counts' in output.read_text() and 'positive-truth pairs' in output.read_text()


def test_summary_follows_links_inside_link_targets_and_publishes_to_checked_parent(tmp_path, monkeypatch):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    storage = tmp_path / 'outside'
    storage.mkdir()
    (project / 'archive').symlink_to(storage)
    alias = tmp_path / 'report-alias'
    alias.symlink_to((project / 'archive').relative_to(tmp_path))
    with pytest.raises(ValueError, match='inside the project'):
        summary.write([project], alias / 'overview.md')
    assert list(storage.iterdir()) == []

    alias.unlink()
    alias.symlink_to(storage)
    render = summary.markdown

    def render_after_alias_change(sections):
        alias.unlink()
        alias.symlink_to(project)
        return render(sections)

    monkeypatch.setattr(summary, 'markdown', render_after_alias_change)
    published = summary.write([project], alias / 'overview.md')
    assert published == storage / 'overview.md' and published.is_file()
    assert not (project / 'overview.md').exists()


def test_first_legacy_reader_retries_if_evaluation_starts_during_snapshot(tmp_path, monkeypatch):
    import evaluation
    bench = summary_projects(tmp_path, published=False)[0] / 'benchmark'
    original = Path.read_bytes
    reads = 0

    def first_writer_after_metrics(path):
        nonlocal reads
        data = original(path)
        if path == bench / 'metrics.tsv':
            reads += 1
            if reads == 1:
                evaluation.publish(bench, ['metric', 'unitas'], [['exact', '0.25']], ['# New report'])
        return data

    monkeypatch.setattr(Path, 'read_bytes', first_writer_after_metrics)
    assert evaluation.read_metrics(bench) == [['metric', 'unitas'], ['exact', '0.25']]
    assert reads == 2


def test_legacy_biological_report_uses_its_primary_metrics_table(tmp_path):
    import evaluation
    import summary
    projects = summary_projects(tmp_path, published=False)
    bench = projects[4] / 'benchmark'
    secondary = [['metric', 'unitas', 'mirge3'], ['lane_agreement', '0.9', '0.88']]
    report = bench / 'report.md'
    with report.open('a') as fh:
        fh.write('\n## Lane replicates\n\n| metric | unitas | mirge3 |\n|---|---|---|\n| lane_agreement | 0.9 | 0.88 |\n')
    assert evaluation.read_metrics(bench)[1][0] == 'mirna_reads'
    (bench / 'metrics.tsv').write_text('\n'.join('\t'.join(row) for row in secondary) + '\n')
    with pytest.raises(ValueError, match='inconsistent'):
        evaluation.read_metrics(bench)
    page = tmp_path / 'overview.md'
    summary.write([projects[1]], page)
    assert 'unweighted means of per-sample shares' in page.read_text()


def test_summary_treats_double_slash_linux_paths_as_the_same_project(tmp_path):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    report = project / 'benchmark/report.md'
    original = report.read_bytes()
    alias = tmp_path / 'link-to-benchmark'
    alias.symlink_to('//' + str(report.parent).lstrip('/'))
    for destination in (Path('//' + str(report).lstrip('/')), alias / report.name):
        with pytest.raises(ValueError, match='inside the project'):
            summary.write([project], destination)
    assert report.read_bytes() == original


def test_summary_keeps_publication_bound_to_the_open_directory(tmp_path, monkeypatch):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    report = project / 'benchmark/report.md'
    original = report.read_bytes()
    out_dir = tmp_path / 'output'
    out_dir.mkdir()
    moved = tmp_path / 'moved'
    actual_replace = os.replace

    def move_parent_before_replace(src, dst, *, src_dir_fd, dst_dir_fd):
        assert src_dir_fd == dst_dir_fd
        out_dir.rename(moved)
        out_dir.symlink_to(report.parent, target_is_directory=True)
        return actual_replace(src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)

    monkeypatch.setattr(os, 'replace', move_parent_before_replace)
    with pytest.raises(ValueError, match='directory moved'):
        summary.write([project], out_dir / report.name)
    assert report.read_bytes() == original
    assert (moved / report.name).read_text().startswith('# Small-RNA benchmark summary')
    assert sorted(p.name for p in moved.iterdir()) == [report.name]


def test_project_input_link_entries_are_preserved(tmp_path):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    report = project / 'benchmark/report.md'
    archive = tmp_path / 'archive'
    archive.mkdir()
    entry = tmp_path / 'archive.md'
    entry.symlink_to(archive, target_is_directory=True)
    archived = report.rename(archive / 'saved.md')
    original = archived.read_bytes()
    report.symlink_to(entry / archived.name)
    with pytest.raises(ValueError, match='aliases one of its inputs'):
        summary.write([project], entry)
    assert entry.is_symlink() and report.read_bytes() == original


def test_write_rechecks_current_input_entries_before_publishing(tmp_path, monkeypatch):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    report = project / 'benchmark/report.md'
    outside = tmp_path / 'outside'
    outside.mkdir()
    destination = outside / report.name
    destination.write_bytes(report.read_bytes())
    original = destination.read_bytes()
    render = summary.markdown

    def render_after_report_relocation(sections):
        report.unlink()
        report.symlink_to(destination)
        return render(sections)

    monkeypatch.setattr(summary, 'markdown', render_after_report_relocation)
    with pytest.raises(ValueError, match='aliases one of its inputs'):
        summary.write([project], destination)
    assert destination.read_bytes() == original
    assert sorted(path.name for path in outside.iterdir()) == [report.name]


@pytest.mark.parametrize('kind', ['fifo', 'link'])
def test_evaluation_lock_requires_regular_file_without_following_links(tmp_path, kind):
    import evaluation
    bench = summary_projects(tmp_path, published=False)[0] / 'benchmark'
    lock = bench / '.evaluation.lock'
    if kind == 'fifo':
        os.mkfifo(lock)
    else:
        target = tmp_path / 'other'
        target.touch()
        lock.symlink_to(target)
    for writer in (False, True):
        with pytest.raises((OSError, ValueError)):
            with evaluation.evaluation_lock(bench, writer=writer):
                pytest.fail('accepted a nonregular lock')


def test_protected_paths_are_checked_after_the_output_directory_is_open(tmp_path, monkeypatch):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    old = tmp_path / 'archive'
    old.mkdir()
    output_dir = tmp_path / 'reports'
    output_dir.mkdir()
    report = project / 'benchmark/report.md'
    stored = report.rename(old / report.name)
    original = stored.read_bytes()
    report.symlink_to(stored)
    open_directory = summary.open_directory

    def open_after_relocation(path):
        output_dir.rmdir()
        old.rename(output_dir)
        old.symlink_to(output_dir, target_is_directory=True)
        return open_directory(path)

    monkeypatch.setattr(summary, 'open_directory', open_after_relocation)
    with pytest.raises(ValueError, match='aliases one of its inputs'):
        summary.write([project], output_dir / report.name)
    assert report.read_bytes() == original


@pytest.mark.parametrize('mounted_root', ['project', 'benchmark'])
def test_summary_refuses_output_below_an_alternate_project_mount(tmp_path, monkeypatch, mounted_root):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    root = project if mounted_root == 'project' else project / 'benchmark'
    alias = tmp_path / 'mount'
    output_dir = alias / 'results/reports'
    output_dir.mkdir(parents=True)
    # A bind mount exposes the same root inode under another name. Model that
    # identity locally without requiring mount privileges in the test runner.
    alias_stat = alias.stat()
    real_stat = Path.stat

    def stat_under_alternate_name(path, **kwargs):
        if path == root:
            return alias_stat
        return real_stat(path, **kwargs)

    monkeypatch.setattr(Path, 'stat', stat_under_alternate_name)
    with pytest.raises(ValueError, match='inside the project'):
        summary.write([project], output_dir / 'new-overview.html')
    assert not list(output_dir.iterdir())


def test_summary_waits_for_first_evaluation_before_checking_for_metrics(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from contextlib import contextmanager
    from threading import Event

    import evaluation
    import summary

    project = summary_projects(tmp_path)[0]
    bench = project / 'benchmark'
    for name in ('metrics.tsv', 'report.md'):
        (bench / name).unlink()
    writer_ready, reader_ready = Event(), Event()
    lock = evaluation.evaluation_lock

    @contextmanager
    def coordinate_first_publication(directory, *, writer):
        if not writer:
            reader_ready.set()
        with lock(directory, writer=writer) as acquired:
            if writer:
                writer_ready.set()
                assert reader_ready.wait(5)
            yield acquired

    monkeypatch.setattr(evaluation, 'evaluation_lock', coordinate_first_publication)
    output = tmp_path / 'summary.md'
    with ThreadPoolExecutor(max_workers=2) as workers:
        publication = workers.submit(evaluation.publish, bench, ['metric', 'unitas'], [['exact', '0.75']], ['# Report'])
        assert writer_ready.wait(5)
        summary.write([project], output)
        publication.result(timeout=5)
    assert '| 0.75 |' in output.read_text()


@pytest.mark.parametrize('prefix', ['missing', 'regular', 'link'])
def test_summary_requires_a_traversable_output_parent(tmp_path, prefix):
    import summary
    project = summary_projects(tmp_path / 'projects')[0]
    (tmp_path / 'regular').write_text('ordinary file')
    (tmp_path / 'link').symlink_to('missing/..', target_is_directory=True)
    output = tmp_path / 'overview.md'
    output.write_text('previous summary')
    with pytest.raises((FileNotFoundError, NotADirectoryError)):
        summary.write([project], tmp_path / prefix / '..' / output.name)
    assert output.read_text() == 'previous summary'


@pytest.mark.parametrize('trailing', ['/', '/.'])
def test_summary_cli_requires_a_filename_as_the_final_component(tmp_path, trailing):
    project = summary_projects(tmp_path / 'projects')[0]
    output = tmp_path / 'overview.md'
    output.write_text('previous summary')
    cmd = [sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project),
           '--output', str(output) + trailing]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    assert result.returncode and 'the extension must be .md or .html' in result.stderr
    assert output.read_text() == 'previous summary'


def test_summary_uses_one_selected_project_when_a_directory_link_is_retargeted(tmp_path, monkeypatch):
    import evaluation
    import summary

    first, second = summary_projects(tmp_path / 'projects')[:2]
    selection = tmp_path / 'selected'
    selection.symlink_to(first, target_is_directory=True)
    read_evaluation = evaluation.read_evaluation

    def retarget_after_metrics(bench):
        result = read_evaluation(bench)
        selection.unlink()
        selection.symlink_to(second, target_is_directory=True)
        return result

    expected = tmp_path / 'expected.md'
    summary.write([first], expected)
    monkeypatch.setattr(evaluation, 'read_evaluation', retarget_after_metrics)
    observed = tmp_path / 'observed.md'
    summary.write([selection], observed)
    assert selection.resolve() == second and observed.read_text() == expected.read_text()


def test_empty_tinyrna_counts_keep_synthetic_truth_rows(tmp_path):
    write(tmp_path / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS1\nMIMAT1\tmiRNA\tmir-1\t0\n')
    methods = ['tinyrna']
    common = adapters.common_table(methods, tmp_path, {}, 'mmu-')
    truth = {('S1', 'MIMAT1'): 12}
    assert metrics.per_mirna(common, truth, methods) == {('tinyrna', 'S1', 'MIMAT1'): (12, 0)}
    assert metrics.count_metrics(common, truth, methods)['tinyrna'] == {'recovered': 0, 'exact': 0, 'misassigned': 0, 'unplaced': 0}


def check_project(tmp_path, kind, samples, methods):
    """A project with config samples, the benchmark's sample file and a dataset declaration, for benchmark.check_methods."""
    project = tmp_path / 'project'
    spec = {'name': 'd', 'kind': kind, 'methods': methods}
    if kind == 'synthetic':
        write(project / 'benchmark/truth.tsv', 'sample\tmature_id\n' + ''.join('{}\tMIMAT1\n'.format(s) for s in samples))
    else:
        spec['mirbase'] = {'release': '23.1'}
        write(project / 'benchmark/samples.tsv', 'sample\tgroup\n' + ''.join('{}\tA\n'.format(s) for s in samples))
    return project, spec, {'samples': {s: {} for s in samples}}


def stopped_text(project, spec, config):
    """check_methods' stopped routes as `route: reason` text (#448: a failed route check no longer stops the evaluation)."""
    return '; '.join('{}: {}'.format(route, reason) for route, reason in benchmark.check_methods(project, spec, config).items())


@pytest.mark.parametrize('kind', ['pool', 'biological'])
def test_evaluate_stops_when_a_method_ran_on_other_samples_or_another_mirbase_release(tmp_path, kind):
    project, spec, config = check_project(tmp_path, kind, ['S1', 'S2'], ['unitas', 'mirge3', 'tinyrna'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tS2\tS1\nmiR-1-3p\t1\t2\n')  # column order free
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S1,S2\nhsa-miR-1-3p,1,2\n')
    write(quant / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS1\tS2\n')
    write_records(quant, ['unitas', 'mirge3', 'tinyrna'])
    benchmark.check_methods(project, spec, config)  # all consistent

    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S1\nhsa-miR-1-3p,1\n')  # a table from a run with fewer samples
    write(quant / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS1\tS2\tS3\n')
    write_records(quant, ['tinyrna'])
    write_records(quant, ['unitas'], release='22')
    err = stopped_text(project, spec, config)
    assert 'mirge3: table samples differ from the benchmark\'s; missing S2' in err
    assert 'tinyrna: table samples differ from the benchmark\'s; extra S3' in err
    assert 'unitas: annotated with mirbase 22, but the dataset declares miRBase 23.1' in err

    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S1,S2\nhsa-miR-1-3p,1,2\n')
    write(quant / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS1\tS2\n')
    write_records(quant, ['unitas'], database='mirgenedb', release='3.0')
    write_records(quant, ['tinyrna'], database='NA', release='NA')  # a db.mirflowz.custom reference
    err = stopped_text(project, spec, config)
    assert 'unitas: annotated with mirgenedb 3.0' in err and 'tinyrna: annotated with NA NA' in err and 'mirge3' not in err

    (quant / 'mirge3/reference.tsv').unlink()
    assert re.search('mirge3: reference record missing', stopped_text(project, spec, config))


def test_evaluate_of_a_project_prepared_before_the_biocore_rename_names_the_new_methods(tmp_path):
    project, spec, config = check_project(tmp_path, 'pool', ['S1'], ['unitas', 'htseq', 'miraligner'])
    with pytest.raises(SystemExit, match=r'names htseq, miraligner in benchmark/dataset.yaml, since renamed biocore_htseq, '
                                         r'biocore_miraligner \(#434\): rerun those routes under the new names'):
        benchmark.check_methods(project, spec, config)


@pytest.mark.parametrize('record, message', [
    ('field\tvalue\tsource\nmirna_database\tmirbase\tx\nmirna_release\t22\tx\nmirna_release\t23.1\tx\n', 'one non-empty mirna_release row, found 2'),
    ('field\tvalue\tsource\nmirna_database\tmirbase\tx\nmirna_release\t\tx\n', 'one non-empty mirna_release row, found 1'),
    ('name\tsetting\nmirna_database\tmirbase\nmirna_release\t23.1\n', 'expected unique field and value columns'),
    ('field\tvalue\tvalue\nmirna_database\tmirbase\tmirbase\nmirna_release\t22\t23.1\n', 'expected unique field and value columns'),
    ('field\tvalue\nmirna_database\tmirbase\nmirna_release\t22\t23.1\n', '3 fields, expected 2'),
])
def test_evaluate_stops_on_a_malformed_reference_record(tmp_path, record, message):
    project, spec, config = check_project(tmp_path, 'pool', ['S1'], ['mirge3'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S1\nhsa-miR-1-3p,1\n')
    write(quant / 'mirge3/reference.tsv', record)
    assert re.search(message, stopped_text(project, spec, config))


def test_evaluate_compares_configured_samples_and_literal_columns_with_the_scored_samples(tmp_path):
    project, spec, config = check_project(tmp_path, 'pool', ['S1', 'S2'], ['unitas'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tS1\tS2\nmiR-1-3p\t1\t2\n')
    write_records(quant, ['unitas'])
    benchmark.check_methods(project, spec, config)
    # requantified with one sample removed from the configuration, the prepared benchmark files kept
    with pytest.raises(SystemExit, match="configured samples differ from the benchmark's; not configured S2"):
        benchmark.check_methods(project, spec, {'samples': {'S1': {}}})
    write(quant / 'unitas/mir_counts.tsv', '\tS1 \tS2\nmiR-1-3p\t1\t2\n')  # a padded column name scores under 'S1 '
    assert re.search("unitas: table samples differ from the benchmark's; missing S1; extra S1 ", stopped_text(project, spec, config))


def test_synthetic_projects_check_samples_only(tmp_path):
    """Synthetic truth comes from the dataset's own reference, so no release is compared; samples still are."""
    project, spec, config = check_project(tmp_path, 'synthetic', ['S01'], ['tesmall', 'mirflowz'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'tesmall/tesmall_counts.tsv', 'ftype\tfid\tS01\n')
    write(quant / 'mirflowz/all_mirna_counts.tab', 'ID\tS01\n')
    benchmark.check_methods(project, spec, config)  # no reference.tsv needed
    write(quant / 'mirflowz/all_mirna_counts.tab', 'ID\tS01\tS01\n')
    assert re.search("mirflowz: table samples differ from the benchmark's; repeated columns", stopped_text(project, spec, config))
    with pytest.raises(SystemExit, match='add its count table to adapters.TABLES'):
        benchmark.check_methods(project, dict(spec, methods=['newtool']), config)


def test_a_table_rebuilt_after_its_record_was_written_stops_the_evaluation(tmp_path):
    """quant_reference records carry their count table's SHA-256: an older record cannot vouch for a newer table."""
    project, spec, config = check_project(tmp_path, 'pool', ['S1'], ['tesmall', 'smrnaseq'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'tesmall/tesmall_counts.tsv', 'ftype\tfid\tS1\nmiRNA\thsa-miR-1-3p\t5\n')
    write(quant / 'smrnaseq/mirna.tsv', 'miRNA\tS1\nhsa-miR-1-3p\t5\n')
    write_records(quant, ['tesmall', 'smrnaseq'])
    benchmark.check_methods(project, spec, config)
    write(quant / 'tesmall/tesmall_counts.tsv', 'ftype\tfid\tS1\nmiRNA\thsa-miR-1-3p\t7\n')  # rebuilt, record not yet rewritten
    assert re.search('tesmall: reference record was written for another count table', stopped_text(project, spec, config))
    write(quant / 'smrnaseq/reference.tsv', 'field\tvalue\tsource\nmirna_database\tmirbase\tx\nmirna_release\t23.1\tx\n')
    assert re.search('smrnaseq: .*one non-empty table_sha256 row, found 0', stopped_text(project, spec, config))


def test_every_method_the_shipped_datasets_declare_is_checked_before_scoring():
    """The declarations' default method lists (e.g. smrnaseq in the synthetic and Giraldez pools) must all be known to
    the preflight, or evaluation of a fully quantified project stops."""
    for path in sorted((BENCH / 'datasets').glob('*.yaml')):
        methods = benchmark.load_dataset(path)['methods']  # a scenario (#453) takes its base's
        assert set(methods) <= set(adapters.TABLES), (path.name, set(methods) - set(adapters.TABLES))
        assert set(methods) <= set(adapters.ADAPTERS), path.name


def test_bootstrap_intervals_are_paired_and_cover_the_value():
    """A constant offset between two methods gives a difference interval of exactly that offset (paired draws), equal
    methods one of 0; for the mean of normal units, the 95% interval covers the true mean in about 95% of datasets."""
    import random
    import uncertainty
    rng = random.Random(7)
    units = [rng.gauss(0, 1) for _ in range(40)]
    compute = lambda drawn: {'a': {'mean': statistics.mean(drawn)}, 'b': {'mean': statistics.mean(drawn) + 0.5},
                             'c': {'mean': statistics.mean(drawn)}, 'oracle': {'mean': 0}}
    point = compute(units)
    ci, differences, skipped = uncertainty.intervals(point, compute, lambda r: (uncertainty.resample(r, units),),
                                                     ['oracle', 'a', 'b', 'c'], ['mean'], 200, 3)
    assert skipped == 0 and [row[:2] for row in ci] == [['mean', 'oracle'], ['mean', 'a'], ['mean', 'b'], ['mean', 'c']]
    low, high = ci[1][3:]
    assert low < point['a']['mean'] < high
    by_pair = {(row[1], row[2]): row[3:] for row in differences}
    assert set(by_pair) == {('a', 'b'), ('a', 'c'), ('b', 'c')}                 # the oracle is not a method
    assert by_pair[('a', 'b')] == [-0.5, -0.5, -0.5] and by_pair[('a', 'c')] == [0.0, 0.0, 0.0]
    assert uncertainty.excludes_zero(['mean', 'a', 'b', -0.5, -0.5, -0.5]) and not uncertainty.excludes_zero(['mean', 'a', 'c', 0, 0, 0])
    covered = 0
    for seed in range(100):
        r = random.Random(seed)
        data = [r.gauss(0, 1) for _ in range(30)]
        f = lambda drawn: {'m': {'mean': statistics.mean(drawn)}}
        (_, _, _, low, high), = uncertainty.intervals(f(data), f, lambda g: (uncertainty.resample(g, data),), ['m'], ['mean'], 300, seed)[0]
        covered += low <= 0 <= high
    assert 85 <= covered <= 100


def test_failing_resamples_are_left_out_and_counted():
    import uncertainty
    def compute(drawn):
        if drawn.count(0) > 2:
            raise ValueError('not scorable')
        return {'m': {'x': sum(drawn)}}
    ci, _, skipped = uncertainty.intervals({'m': {'x': 3}}, compute, lambda r: (uncertainty.resample(r, [0, 1, 2]),), ['m'], ['x'], 100, 1)
    assert 0 < skipped < 100 and ci[0][:3] == ['x', 'm', 3]


def test_biological_bootstrap_summary_matches_method_summary():
    """summary_of, which each resample aggregates from per-sample values, gives method_summary's and unplaced's values."""
    common = [('u', 'S1', 'miR-1', 'MIMAT1', 20.0), ('u', 'S1', 'miR-2', 'MIMAT2', 5.0), ('u', 'S2', 'miR-1', 'MIMAT1', 12.0),
              ('u', 'S2', 'miR-2', 'MIMAT2', 15.0), ('m', 'S1', 'miR-1', 'MIMAT1', 30.0), ('m', 'S2', 'miR-3', 'MIMAT3', 11.0)]
    raw = common + [('u', 'S1', 'miR-x', '', 4.0), ('m', 'S2', 'miR-y', '', 2.0)]
    methods, samples = ['u', 'm'], ['S1', 'S2']
    expected = biological.method_summary(common, methods, samples)
    shares = biological.unplaced(raw, methods, samples)
    for m in methods:
        expected[m]['unplaced'] = shares[m]['share']
    stats = biological.sample_stats(common, raw, methods, samples)
    assert biological.summary_of(stats, methods, samples) == expected
    assert 'detected_all_samples' not in biological.summary_of(stats, methods, samples, all_samples=False)['u']
    assert biological.summary_of(stats, methods, ['S1', 'S1'])['u']['detected_all_samples'] == 1   # a resample repeating S1


def test_published_intervals_are_checked_and_shown_in_the_summary(tmp_path):
    import evaluation
    import summary
    project = summary_projects(tmp_path)[1]       # eq3
    bench = project / 'benchmark'
    metrics = [['metric', 'oracle', 'unitas', 'mirge3', 'mirflowz'], ['detected', '0.9', '0.97', '0.95', '0.96']]
    ci = 'metric\trow\tvalue\tlow\thigh\ndetected\tunitas\t0.97\t0.93\t0.99\ndetected\tmirge3\t0.95\t\t\n'
    evaluation.publish(bench, metrics[0], metrics[1:], ['# eq'], None,
                       {'metrics_ci.tsv': ci, 'differences.tsv': 'metric\trow_a\trow_b\tdifference\tlow\thigh\n'},
                       {'resamples': 10, 'seed': 1, 'skipped': 0})
    row = summary.read_project(project, tmp_path)[1]
    assert row['values'][('detected', 'unitas')] == '0.97 (0.93–0.99)' and row['values'][('detected', 'mirge3')] == '0.95'
    (bench / 'metrics_ci.tsv').write_text(ci.replace('0.99', '0.98'))
    with pytest.raises(ValueError, match='incomplete or inconsistent evaluation'):
        summary.read_project(project, tmp_path)
    with pytest.raises(ValueError, match='not an evaluation file'):
        evaluation.publish(bench, metrics[0], metrics[1:], ['# eq'], None, {'other.tsv': ''})
    for name, text in (('metrics_ci.tsv', 'metric\trow\n'), ('differences.tsv', 'metric\trow_a\n'),
                       ('differences.tsv', 'metric\trow_a\trow_b\tdifference\tlow\thigh\nx\ta\tb\t1\t0\n')):
        good = {'metrics_ci.tsv': ci, 'differences.tsv': 'metric\trow_a\trow_b\tdifference\tlow\thigh\n'}
        evaluation.publish(bench, metrics[0], metrics[1:], ['# eq'], None, dict(good, **{name: text}))
        with pytest.raises(ValueError, match='incomplete or inconsistent evaluation'):   # a wrong header, or a short row
            summary.read_project(project, tmp_path)


def test_synthetic_bootstrap_scores_and_options():
    """truth_scores over all truth pairs equals count_metrics; --bootstrap refuses a negative count."""
    import argparse
    common = [('m', 'S01', 'a', 'A', 10), ('m', 'S01', 'b', 'B', 1), ('m', 'S02', 'a', 'A', 10), ('m', 'S02', 'b', 'B', 3)]
    truth = {('S01', 'A'): 10, ('S01', 'B'): 3, ('S02', 'A'): 10, ('S02', 'B'): 3}
    table = metrics.per_mirna(common, truth, ['m'])
    full = metrics.count_metrics(common, truth, ['m'])['m']
    assert metrics.truth_scores(table, truth, ['m'], list(truth))['m'] == {k: full[k] for k in ('recovered', 'exact')}
    assert benchmark.nonnegative('0') == 0
    with pytest.raises(argparse.ArgumentTypeError, match='must be 0 or more'):
        benchmark.nonnegative('-5')


# Differential expression on the ratiometric pools (#414): inputs, scores, the de.smk rule and evaluate --de.

import de  # noqa: E402
import hashlib  # noqa: E402
import re  # noqa: E402
import math  # noqa: E402

DE_RUNS = {'A1': ('A', 'L1'), 'A2': ('A', 'L2'), 'B1': ('B', 'L1'), 'B2': ('B', 'L2')}
ONE, EQUAL, TWENTY_ONE = 'TGGAATGTAAAGAAGTATGTAT', 'ACAGTAGTCTGCACATTGGTTA', 'TAGCTTATCAGACTGATGTTGA'
DE_SEQUENCES = {'miR-1-3p': ONE, 'miR-199a-3p': EQUAL, 'miR-21-5p': TWENTY_ONE}


def de_entries(*designs):
    """Scored exact pool entries S0, S1, ... with designed amounts (a, b)."""
    return [{'pool_id': 'S{}'.format(i), 'sequence': 'SEQ{}'.format(i), 'group': 'exact', 'mature_ids': ['MIMAT{:07d}'.format(i)],
             'amount_a': a, 'amount_b': b} for i, (a, b) in enumerate(designs)]


def test_de_inputs_use_the_workflow_library_sizes_and_whole_counts():
    entries = de_entries((4, 1), (1, 1)) + [{'pool_id': 'SH', 'sequence': 'SEQSH', 'group': 'shifted', 'mature_ids': ['MIMAT9999999'],
                                             'amount_a': 1, 'amount_b': 1}]
    counts = {('m', 'A1', 'S0'): 0.5, ('m', 'A1', 'S1'): 1.5, ('m', 'B1', 'S0'): 2.5, ('m', 'B1', 'S1'): 3}
    reads = {('m', 'A1'): {'all': 1234.5}, ('m', 'B1'): {'all': 100}}  # every miRNA read, outside the pool too
    oracle = {'A1': {'SEQ0': 7, 'SEQ1': 2, 'SEQSH': 1, 'OTHER': 50}, 'B1': {'SEQ0': 1, 'SEQ1': 2, 'SEQSH': 0}}
    sizes = de.library_sizes(reads, oracle, entries, ['oracle', 'm'], ['A1', 'B1'])
    assert sizes == {('oracle', 'A1'): 10, ('oracle', 'B1'): 3, ('m', 'A1'): 1234.5, ('m', 'B1'): 100}  # oracle: exact + shifted reads
    count = de.pool_count(counts, oracle, entries)
    assert de.pool_features(entries[:2]) == [('S0', 2.0), ('S1', 0.0)]  # the scored sequences with their designed log2(A / B)
    texts = de.input_tables(count, ['oracle', 'm'], ['S0', 'S1'], [('A1', 'A', {'lab': 'L1'}), ('B1', 'B', {'lab': 'L2'})], sizes,
                            de.model_text('A', 'B'))
    assert list(texts) == list(de.INPUTS)
    assert texts['counts.tsv'] == 'row\tfeature\tA1\tB1\noracle\tS0\t7\t1\noracle\tS1\t2\t2\nm\tS0\t0\t2\nm\tS1\t2\t3\n'  # half to even
    assert texts['library_sizes.tsv'] == 'row\tA1\tB1\noracle\t10\t3\nm\t1234.5\t100\n'
    assert texts['design.tsv'] == 'sample\tgroup\tlab\nA1\tA\tL1\nB1\tB\tL2\n'
    assert texts['model.tsv'] == 'field\tvalue\nnumerator\tA\ndenominator\tB\nmin_mean\t0\n'  # pool A against pool B, every sequence with reads
    assert de.fractional_share(count, 'oracle', ['S0', 'S1'], ['A1', 'B1']) == 0.0
    assert de.fractional_share(count, 'm', ['S0', 'S1'], ['A1', 'B1']) == 1.0  # both sequences have a count that is rounded


def test_de_model_names_the_contrast_numeric_covariates_and_the_least_mean_count():
    assert de.model_text('case', 'Fertile control', ['nuisance'], 10) == (
        'field\tvalue\nnumerator\tcase\ndenominator\tFertile control\nmin_mean\t10\nnumeric\tnuisance\n')
    script = (BENCH / 'de.R').read_text()
    known = re.search(r"unknown <- setdiff\(settings\$field, c\(([^)]*)\)\)", script)
    assert known and known.group(1) == "'numerator', 'denominator', 'min_mean', 'numeric', 'normalisation'"  # the fields model_text writes
    assert de.model_text('A', 'B', (), 0, 'tmm').endswith('normalisation\ttmm\n') and 'calcNormFactors(y, method = \'TMM\')' in script
    assert 'lib.size = size * y$samples$norm.factors' in script  # voom on the TMM-normalised library sizes
    assert "levels = c(denominator, numerator)" in script and 'coefficient <- ncol(model)' in script  # log2(numerator / denominator)
    assert "deparse(formula, width.cutoff = 500L), collapse = ' ')" in script and script.count('deparse(') == 1  # one string in the record


def test_de_inputs_are_rewritten_only_when_changed(tmp_path):
    folder = tmp_path / 'de/input'
    de.write_inputs(folder, {'counts.tsv': 'a\n'})
    os.utime(folder / 'counts.tsv', (1, 1))
    de.write_inputs(folder, {'counts.tsv': 'a\n'})
    assert (folder / 'counts.tsv').stat().st_mtime == 1  # unchanged: de.smk does not run again
    de.write_inputs(folder, {'counts.tsv': 'b\n'})
    assert (folder / 'counts.tsv').read_text() == 'b\n'


def test_de_scores_call_in_the_designed_direction_and_false_calls_on_1to1_only():
    exact = de_entries((1, 1), (1, 1), (3, 2), (4, 1), (1, 8))  # designed log2(A/B): 0, 0, 0.58, 2, -3
    fits = {'S0': (0.2, 0.001, 0.01), 'S1': (-0.1, 0.5, 0.6), 'S2': (0.5, 0.001, 0.02), 'S3': (1.5, 0.001, 0.01), 'S4': (2.0, 0.001, 0.01)}
    results = {('m', test, pool_id): fit for test in de.TESTS for pool_id, fit in fits.items()}
    expected = {'called': 0.5,        # of S3 and S4 (at least 2-fold), S4 is called in the wrong direction
                'false_calls': 0.5,   # S0 of the 1:1 S0 and S1
                'fdr': 0.5,           # of the calls S0, S2, S3 and S4: S0 (1:1) and S4 (wrong direction); S2, 1.5-fold, is true
                'fc_bias': -0.5,      # median of -0.08 (S2), -0.5 (S3) and -5 (S4)
                'offset': 0.05,       # median of S0 and S1
                'false_calls_against': 0.0}  # 2 of the 3 changes go up (ONE_WAY); S0 is called up, not against
    assert de.scores(results, de.pool_features(exact), ['m']) == {'m': {'{}_{}'.format(t, k): v for t in de.TESTS for k, v in expected.items()}}  # no oracle row
    assert sorted(de.metric_keys()) == sorted('{}_{}'.format(t, k) for t in de.TESTS for k in list(expected) + list(de.ORACLE_KEYS))
    assert de.metric_keys(['', '.L1']) == de.metric_keys() + ['{}.L1'.format(k) for k in de.metric_keys()]


def test_de_scores_against_the_oracle_compare_calls_and_fold_changes_on_the_same_runs():
    exact = de_entries((1, 1), (1, 1), (4, 1), (1, 4))  # designed log2(A/B): 0, 0, 2, -2
    oracle = {'S0': (0.1, 0.5, 0.6), 'S1': (0.4, 0.001, 0.01), 'S2': (2.0, 0.001, 0.01), 'S3': (-1.5, 0.1, 0.2)}
    method = {'S0': (0.3, 0.001, 0.01), 'S1': (0.5, 0.001, 0.02), 'S2': (-1.0, 0.001, 0.01), 'S3': (None, None, None)}
    results = {(row, test, pool_id): fit for row, fits in (('oracle', oracle), ('m', method)) for test in de.TESTS for pool_id, fit in fits.items()}
    values = de.scores(results, de.pool_features(exact), ['oracle', 'm'], '.L1')
    assert values['m']['edger_own_false_calls.L1'] == 0.5       # S0 of the 1:1 S0 and S1; the oracle calls S1 too
    assert values['m']['voom_discordant_vs_oracle.L1'] == 0.5   # S0 (up, the oracle none) and S2 (down, the oracle up); S3 untested, none
    assert values['m']['edger_fc_shift_vs_oracle.L1'] == 0.1      # median of the differences 0.2, 0.1 and -3 over the sequences both test
    assert values['m']['edger_glm_fc_error_vs_oracle.L1'] == 0.1  # median of their distances from that shift: 0.1, 0 and 3.1
    assert not {'{}_{}.L1'.format(t, k) for t in de.TESTS for k in de.ORACLE_KEYS} & set(values['oracle'])
    # the synthetic dataset's reference row is the truth: the same metrics under its name
    renamed = {('truth' if row == 'oracle' else row, test, feature): fit for (row, test, feature), fit in results.items()}
    against_truth = de.scores(renamed, de.pool_features(exact), ['truth', 'm'], '.adjusted', 'truth')
    assert against_truth['m']['voom_discordant_vs_truth.adjusted'] == 0.5 and against_truth['m']['edger_fc_shift_vs_truth.adjusted'] == 0.1
    assert de.reference_keys('truth')[1:] == ('discordant_vs_truth', 'fc_shift_vs_truth', 'fc_error_vs_truth')
    assert 'voom_fc_error_vs_truth.adjusted' in de.metric_keys(['', '.adjusted'], 'truth')
    untested = {key: (None, None, None) if key[0] == 'm' else fit for key, fit in results.items()}
    values = de.scores(untested, de.pool_features(exact[2:]), ['oracle', 'm'])['m']  # no 1:1 sequence, nothing tested by both
    assert math.isnan(values['edger_own_false_calls']) and math.isnan(values['edger_fc_error_vs_oracle']) and math.isnan(values['edger_fc_shift_vs_oracle'])
    assert values['edger_discordant_vs_oracle'] == 0.5  # S2: the oracle calls it, the untested method does not


def test_de_fold_change_shift_against_the_oracle_is_not_counted_as_error():
    exact = de_entries((1, 1), (4, 1), (1, 4), (2, 1))
    oracle = {'S0': (0.1, 0.5, 0.6), 'S1': (2.1, 0.001, 0.01), 'S2': (-1.9, 0.001, 0.01), 'S3': (0.8, 0.01, 0.02)}
    method = {pool_id: (lfc - 0.25, p, padj) for pool_id, (lfc, p, padj) in oracle.items()}  # other library sizes: every log2fc 0.25 lower
    results = {(row, test, pool_id): fit for row, fits in (('oracle', oracle), ('m', method)) for test in de.TESTS for pool_id, fit in fits.items()}
    values = de.scores(results, de.pool_features(exact), ['oracle', 'm'])['m']
    assert values['voom_fc_shift_vs_oracle'] == -0.25 and values['voom_fc_error_vs_oracle'] == 0.0


def test_de_scores_of_a_resample_are_nan_not_errors():
    results = {('m', test, 'S0'): (None, None, None) for test in de.TESTS}  # untested: no reads
    values = de.scores(results, de.pool_features(de_entries((3, 2))), ['m'])['m']  # no 1:1 and no 2-fold sequence, no call
    assert len(values) == len(de.TESTS) * (len(de.KEYS) - len(de.ORACLE_KEYS)) and all(math.isnan(v) for v in values.values())


def test_de_results_must_cover_every_row_test_and_sequence(tmp_path):
    lines = ['\t'.join(de.RESULT_HEADER)] + ['m\t{}\t{}\t1.5\t0.01\t{}'.format(t, p, 'NA' if p == 'S1' else '0.02')
                                             for t in de.TESTS for p in ('S0', 'S1')]
    path = write(tmp_path / 'results.tsv', '\n'.join(lines) + '\n')
    read = lambda pool_ids: de.read_results(path.read_text(), ['m'], pool_ids, path)
    table = read(['S0', 'S1'])
    assert table[('m', 'edger', 'S0')] == (1.5, 0.01, 0.02) and table[('m', 'voom', 'S1')] == (1.5, 0.01, None)
    with pytest.raises(ValueError, match='and {} unexpected'.format(len(de.TESTS))):  # S1 in every test
        read(['S0'])
    with pytest.raises(ValueError, match='{} expected lines missing'.format(len(de.TESTS))):  # S2 in every test
        read(['S0', 'S1', 'S2'])
    write(path, '\n'.join(lines + [lines[1]]) + '\n')
    with pytest.raises(ValueError, match='repeated'):
        read(['S0', 'S1'])
    write(path, 'row\ttest\tpool_id\tlogFC\n')
    with pytest.raises(ValueError, match='columns'):
        read(['S0', 'S1'])


def de_declaration(tmp_path, change):
    spec = yaml.safe_load((BENCH / 'datasets/giraldez2018_ratiometric_truseq.yaml').read_text())
    change(spec)
    return write(tmp_path / 'ds.yaml', yaml.safe_dump(spec))


@pytest.mark.parametrize('change, message', [
    (lambda s: s.update(pool='equimolar'), 'is for ratiometric pools'),
    (lambda s: s.update(de={'covariates': ['batch']}), 'de.covariates batch: every run needs a value'),  # no such column of the runs
    (lambda s: s.update(de={'covariates': ['lab', 'lab']}), 'distinct column names'),
    (lambda s: s.update(de={'covariates': ['lab+group']}), 'distinct column names'),  # a term of de.R's formula
    (lambda s: s.update(de={'covariates': ['group']}), 'distinct column names'),
    (lambda s: s.update(de={'contrast': ['B', 'A']}), 'takes only `covariates` and `within`'),  # a pool is A against B, as its design
    (lambda s: s.update(de={'covariates': ['lab'], 'test': 'deseq2'}), 'takes only `covariates`'),
    (lambda s: [run.update(lab='Lab1') for run in s['runs']], 'at least two different ones'),
    (lambda s: s['runs'][0].pop('lab'), 'every run needs a value'),
    (lambda s: s.update(runs=[r for r in s['runs'] if r['group'] == 'A'] + [s['runs'][2]]), 'got 4 A and 1 B'),
    # two runs per pool over three labs, one lab shared: ~ lab + group has as many coefficients as runs
    (lambda s: s.update(runs=[dict(s['runs'][0], lab='L1'), dict(s['runs'][1], lab='L2'), dict(s['runs'][2], lab='L1'),
                              dict(s['runs'][3], lab='L3')]), '4 coefficients for 4 runs, so no residual degrees of freedom'),
])
def test_de_declaration_needs_lab_with_two_levels_and_replicates(tmp_path, change, message):
    with pytest.raises(SystemExit, match=message):
        benchmark.load_dataset(de_declaration(tmp_path, change))


def one_lab(spec):
    spec['de']['covariates'] = []
    for run in spec['runs']:
        run['lab'] = 'Lab1'


@pytest.mark.parametrize('change, message', [
    (lambda s: s['de'].update(within='instrument'), 'de.within instrument: every run needs a value that starts with a letter or digit'),  # HiSeq 4000
    (lambda s: s['de'].update(within='the lab'), "de.within must be a column name"),
    (one_lab, 'de.within lab: every run needs a value that starts with a letter or digit'),
    (lambda s: s['runs'][0].update(lab='Lab 1'), 'with at least two different ones'),
    (lambda s: s['runs'][0].update(lab='..'), 'de.within lab: every run needs a value that starts with a letter or digit'),  # not a folder of its own
    (lambda s: s['runs'][0].update(lab='lab1'), 'values that differ only in case would share a folder: Lab1, Lab5, lab1'),
    (lambda s: s['runs'].pop(0), 'de.within lab Lab1: differential expression needs at least two runs in each group, got 1 A and 2 B'),
])
def test_de_within_needs_safe_level_names_and_replicates_in_each_level(tmp_path, change, message):
    with pytest.raises(SystemExit, match=re.escape(message)):
        benchmark.load_dataset(de_declaration(tmp_path, change))


def test_shipped_ratiometric_declarations_declare_their_de_design():
    truseq = benchmark.load_dataset(BENCH / 'datasets/giraldez2018_ratiometric_truseq.yaml')
    nextflex = benchmark.load_dataset(BENCH / 'datasets/giraldez2018_ratiometric_nextflex.yaml')
    assert truseq['de'] == {'covariates': ['lab'], 'within': 'lab'} and nextflex['de'] == {'covariates': []}
    genome = benchmark.load_dataset(BENCH / 'datasets/synthetic_hsa_genome.yaml')
    assert genome['de'] == {} and benchmark.de_kind(genome) == 'synthetic' and benchmark.de_contrast(truseq) == ('A', 'B')
    for path in BENCH.glob('datasets/*.yaml'):
        spec = benchmark.load_dataset(path)  # a scenario (#453) takes its base's `de`
        assert ('de' in spec) == (spec.get('pool') == 'ratiometric' or spec['kind'] == 'mixture' or
                                  spec['name'] in ('synthetic_hsa_genome', 'synthetic_hsa_genome_null', 'synthetic_hsa_genome_unbalanced',
                                                   'gse317995_nextflex_v4', 'gse196633_truseq')), path.name


def test_de_sets_take_any_column_of_the_runs_as_a_covariate(tmp_path):
    """A covariate other than `lab` (#447): the instrument, from the declaration's runs, not from samples.tsv."""
    spec = benchmark.load_dataset(de_declaration(tmp_path, lambda s: s.update(de={'covariates': ['instrument']})))
    rows = [{'sample': run['run'], 'group': run['group'], 'lab': run['lab'], 'reads': ''} for run in spec['runs']]
    whole, fixed, tmm = benchmark.de_sets(spec, benchmark.de_units(spec, rows))
    assert whole.design[0] == ('SRR5234383', 'A', {'instrument': 'HiSeq 4000'}) and whole.design[-1][2] == {'instrument': 'HiSeq 2500'}
    assert (whole.suffix, whole.folder, whole.samples, whole.arm) == ('', '', '4 A and 4 B runs', None) and whole.model == de.model_text('A', 'B')
    # the normalisation arms (#453): the whole design again; tmm says so in its model, fixed gets other library sizes
    assert (fixed.suffix, fixed.folder, fixed.arm, tmm.suffix, tmm.folder, tmm.arm) == ('.fixed', 'arm/fixed', 'fixed', '.tmm', 'arm/tmm', 'tmm')
    assert fixed.design == tmm.design == whole.design and fixed.model == whole.model and tmm.model == de.model_text('A', 'B', (), 0, 'tmm')
    sets = benchmark.de_sets(benchmark.load_dataset(BENCH / 'datasets/giraldez2018_ratiometric_truseq.yaml'), benchmark.de_units(spec, rows))
    assert [(s.suffix, s.folder, s.samples, list(s.design[0][2])) for s in sets] == [
        ('', '', '4 A and 4 B runs', ['lab']), ('.Lab1', 'lab/Lab1', '2 A and 2 B runs', []), ('.Lab5', 'lab/Lab5', '2 A and 2 B runs', []),
        ('.fixed', 'arm/fixed', '4 A and 4 B runs', ['lab']), ('.tmm', 'arm/tmm', '4 A and 4 B runs', ['lab'])]


def biological_declaration(tmp_path, de_spec, groups=('case', 'case', 'control', 'control', 'other', 'other')):
    spec = {'name': 'bio', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
            'methods': ['unitas'], 'data': 'benchmark/smallrna/bio', 'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
            'samples': [{'sample': 'S{}'.format(i), 'group': group, 'patient': 'P{}'.format(i % 2), 'runs': [{'run': 'R{}'.format(i)}]}
                        for i, group in enumerate(groups)]}
    if de_spec is not None:
        spec['de'] = de_spec
    return write(tmp_path / 'bio.yaml', yaml.safe_dump(spec))


def test_a_biological_declaration_names_its_contrast_and_covariates(tmp_path):
    """`de:` for biological datasets (#447): a contrast of two of the samples' groups, with covariates such as the patient."""
    spec = benchmark.load_dataset(biological_declaration(tmp_path, {'contrast': ['case', 'control'], 'covariates': ['patient']}))
    assert benchmark.de_contrast(spec) == ('case', 'control')
    rows = [{'sample': s['sample'], 'library': s['sample'], 'group': s['group'], 'lane': 'merged'} for s in spec['samples']]
    (whole,) = benchmark.de_sets(spec, benchmark.de_units(spec, rows))
    assert whole.design == [('S0', 'case', {'patient': 'P0'}), ('S1', 'case', {'patient': 'P1'}), ('S2', 'control', {'patient': 'P0'}),
                            ('S3', 'control', {'patient': 'P1'})]  # the third group's samples are not in the set
    assert whole.samples == '2 case and 2 control samples' and whole.model == de.model_text('case', 'control', (), 10)  # detected miRNAs
    for de_spec, message in (({}, 'needs a `contrast`, `null` comparisons (#449) or both'),
                             ({'contrast': ['case', 'absent']}, 'de.contrast must name two different groups of the samples'),
                             ({'contrast': ['case', 'case']}, 'de.contrast must name two different groups'),
                             ({'contrast': ['case', 'control'], 'covariates': ['age']}, 'de.covariates age: every sample needs a value'),
                             ({'contrast': ['case', 'control'], 'test': 'deseq2'}, 'takes only `contrast`, `covariates`, `within` and `null`')):
        with pytest.raises(SystemExit, match=re.escape(message)):
            benchmark.load_dataset(biological_declaration(tmp_path, de_spec))
    with pytest.raises(SystemExit, match='at least two samples in each group, got 1 case and 2 control'):
        benchmark.load_dataset(biological_declaration(tmp_path, {'contrast': ['case', 'control']}, ('case', 'control', 'control')))
    # Group labels that YAML reads as numbers (a dose of 100 against 0): samples.tsv has them as text, and so has the contrast.
    spec = benchmark.load_dataset(biological_declaration(tmp_path, {'contrast': [100, 0]}, (100, 100, 0, 0, 10, 10)))
    assert benchmark.de_contrast(spec) == ('100', '0')
    rows = [{'sample': s['sample'], 'library': s['sample'], 'group': str(s['group']), 'lane': 'merged'} for s in spec['samples']]
    (whole,) = benchmark.de_sets(spec, benchmark.de_units(spec, rows))
    assert [(sample, group) for sample, group, _ in whole.design] == [('S0', '100'), ('S1', '100'), ('S2', '0'), ('S3', '0')]
    assert whole.model == de.model_text('100', '0', (), 10) and whole.samples == '2 100 and 2 0 samples'
    with pytest.raises(SystemExit, match=re.escape('benchmark/samples.tsv has 1 100 and 2 0 sample(s) for the DE set of the whole dataset')):
        benchmark.de_sets(spec, benchmark.de_units(spec, rows[1:]))  # a project prepared from another declaration
    for contrast in ([100, None], ['', 0], [[100], 0]):
        with pytest.raises(SystemExit, match='de.contrast must name two different groups'):
            benchmark.load_dataset(biological_declaration(tmp_path, {'contrast': contrast}, (100, 100, 0, 0)))


def de_snakemake(cwd, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(cwd), XDG_CACHE_HOME=str(cwd / '.cache'))
    return subprocess.run([str(Path(sys.executable).with_name('snakemake')), *args], cwd=cwd, env=env, capture_output=True, text=True, timeout=120)


def test_de_snakefile_runs_r_in_the_sartools_image_for_every_set(tmp_path):
    folders = ('', 'lab/L1/', 'lab/L2/')
    for folder in folders:
        for name in de.INPUTS:
            write(tmp_path / 'benchmark/de' / (folder + 'input') / name, 'x\n')
    result = de_snakemake(tmp_path, '-s', str(BENCH / 'de.smk'), '-n')
    assert result.returncode != 0 and 'benchmark/de/levels.tsv not found' in result.stdout + result.stderr  # evaluate writes it
    write(tmp_path / 'benchmark/de/levels.tsv', de.levels_text([('lab', 'L1'), ('lab', 'L2')]))
    write(tmp_path / 'Snakefile', 'include: {!r}\nimport json\nfrom pathlib import Path\n'
          'Path("rules.json").write_text(json.dumps({{r.name: r.container_img for r in workflow.rules}}))\n'.format(str(BENCH / 'de.smk')))
    result = de_snakemake(tmp_path, '--list-rules')
    assert result.returncode == 0, result.stdout + result.stderr
    sartools = 'docker://' + yaml.safe_load((ROOT / 'docker.config').read_text())['docker']['sartools']
    assert json.loads((tmp_path / 'rules.json').read_text()) == {'all': None, 'benchmark_de': sartools, 'benchmark_de_level': sartools}
    elsewhere = tmp_path / 'elsewhere'
    elsewhere.mkdir()
    result = de_snakemake(elsewhere, '-s', str(BENCH / 'de.smk'), '-d', str(tmp_path), '-n', '-p', '--cores', '1')
    assert result.returncode == 0, result.stdout + result.stderr
    for folder in folders:
        assert ('Rscript {0} --counts benchmark/de/{1}input/counts.tsv --library-sizes benchmark/de/{1}input/library_sizes.tsv '
                '--design benchmark/de/{1}input/design.tsv --model benchmark/de/{1}input/model.tsv --results benchmark/de/{1}results.tsv '
                '--record benchmark/de/{1}record.tsv > benchmark/de/{1}de.log 2>&1'.format(BENCH / 'de.R', folder)) in result.stdout
    nulls = tmp_path / 'nulls'  # null comparisons only (#449): no whole set
    for name in de.INPUTS:
        write(nulls / 'benchmark/de/null/case-1/input' / name, 'x\n')
    write(nulls / 'benchmark/de/levels.tsv', de.levels_text([]))
    result = de_snakemake(nulls, '-s', str(BENCH / 'de.smk'), '-n', '--cores', '1')
    assert result.returncode != 0 and 'has no DE set' in result.stdout + result.stderr
    write(nulls / 'benchmark/de/levels.tsv', de.levels_text([('null', 'case-1')]))
    result = de_snakemake(nulls, '-s', str(BENCH / 'de.smk'), '-n', '-p', '--cores', '1')
    assert result.returncode == 0, result.stdout + result.stderr
    assert '--results benchmark/de/null/case-1/results.tsv' in result.stdout and '--results benchmark/de/results.tsv' not in result.stdout


def test_de_snakefile_passes_the_options_the_r_script_defines():
    import re
    defined = set(re.findall(r"make_option\('(--[a-z-]+)'", (BENCH / 'de.R').read_text()))
    passed = set(re.findall(r"(--[a-z-]+) \{", (BENCH / 'de.smk').read_text()))
    assert defined == passed == {'--counts', '--library-sizes', '--design', '--model', '--results', '--record'}


def test_de_r_script_writes_the_tests_de_py_reads():
    import re
    tests = re.search(r"^tests <- c\(([^)]*)\)", (BENCH / 'de.R').read_text(), re.M)
    assert tests and tests.group(1) == ', '.join("'{}'".format(t) for t in de.TESTS) and len(de.TEST_NAMES) == len(de.TESTS)


DE_RUNS_BY_LAB = {'A1': ('A', 'L1'), 'A2': ('A', 'L1'), 'A3': ('A', 'L2'), 'A4': ('A', 'L2'),
                  'B1': ('B', 'L1'), 'B2': ('B', 'L1'), 'B3': ('B', 'L2'), 'B4': ('B', 'L2')}


def de_project(tmp_path, runs=DE_RUNS, declared=None):
    """A prepared and quantified ratiometric pool (RATIO: miR-1-3p 4:1, miR-199a-3p 1:1, miR-21-5p 1:8) with `runs`
    ({run: (pool, lab)}; the i-th run of a pool has i more reads), declaring `de: declared` (default: covariates lab);
    UNITAS and miRge3 outputs and the filtered reads, as quant_all --notemp leaves them."""
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/p'
    write(data / 'composition/ratio.csv', RATIO)
    write(data / 'mirbase/23.1/mature.fa', HSA)
    index = {run: [r for r, (g, _) in runs.items() if g == group].index(run) for run, (group, _) in runs.items()}
    amounts = {run: (8 + index[run], 4 + index[run], 1) if group == 'A' else (2, 4 + index[run], 8 + index[run])
               for run, (group, _) in runs.items()}  # reads of miR-1-3p, miR-199a-3p and miR-21-5p
    reads = {run: [ONE] * a + [EQUAL] * b + [TWENTY_ONE] * c for run, (a, b, c) in amounts.items()}
    for run, seqs in reads.items():
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/{}.fastq.gz'.format(run), 'wt') as fh:
            fh.write(''.join('@r{}\n{}TGGAATTCTC\n+\n{}\n'.format(i, s, 'I' * (len(s) + 10)) for i, s in enumerate(seqs)))
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({
        'name': 'p', 'kind': 'pool', 'pool': 'ratiometric', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
        'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/p', 'de': declared or {'covariates': ['lab']},
        'composition': {'file': 'ratio.csv', 'url': 'unused', 'sha256': 'unused'},
        'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'runs': [{'run': r, 'group': g, 'lab': lab, 'md5': 'unused', 'url': 'unused'} for r, (g, lab) in runs.items()]}))
    project = tmp_path / 'project'
    dump = unitas_dump(ext, 'homo_sapiens', 'Human', {})
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(ext)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    quant = project / 'data/tmp/smallrna/quant'
    names = list(runs)
    table = [(mirna, [amounts[run][i] for run in names]) for i, mirna in enumerate(['miR-1-3p', 'miR-199a-3p', 'miR-21-5p'])]
    write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(names) + '\n' + ''.join('{}\t{}\n'.format(m, '\t'.join(map(str, v))) for m, v in table))
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(names) + '\n' + ''.join('hsa-{},{}\n'.format(m, ','.join(map(str, v))) for m, v in table))
    isomir_outputs(quant, DE_SEQUENCES)
    write_records(quant, ['unitas', 'mirge3'], dump=dump)
    for run, seqs in reads.items():
        write_stages(project, run, seqs, seqs, seqs)
        write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', run, 100 + 10 * index[run])  # the `fixed` arm's sizes (#453)
    return project


def de_outputs(project, fits, change=None, drop_last=False, folder='', formula=None):
    """de.smk's outputs as de.R writes them for the current inputs of the DE set in benchmark/de/`folder`: every row and
    test gets `fits` ({pool_id: (log2fc, pvalue, padj)}); `change` edits the record; `drop_last` leaves the table's last
    line out (a table de.R recorded, but incomplete). The whole dataset's set (folder '') brings its normalisation arms
    (#453) with the same fits, written first, so that a check fails on the whole dataset's set as before."""
    if folder == '':
        for arm in de.ARMS:
            if (project / 'benchmark/de/arm' / arm / 'input').is_dir():
                de_outputs(project, fits, folder='arm/' + arm, formula='~lab + group')
    base = project / 'benchmark/de' / folder
    texts = {name: (base / 'input' / name).read_text() for name in de.INPUTS}
    rows = list(dict.fromkeys(line.split('\t')[0] for line in texts['counts.tsv'].splitlines()[1:]))
    lines = ['\t'.join(de.RESULT_HEADER)] + ['\t'.join([row, test, pool_id] + ['NA' if v is None else str(v) for v in fit])
                                             for row in rows for test in de.TESTS for pool_id, fit in fits.items()]
    results = '\n'.join(lines[:-1] if drop_last else lines) + '\n'
    write(base / 'results.tsv', results)
    record = {'md5_' + Path(name).stem: de.md5(text) for name, text in texts.items()}
    record.update(md5_script=de.md5(de.SCRIPT.read_bytes()), md5_results=de.md5(results), R='R version 4.5.3', edgeR='4.10.7',
                  limma='3.68.5', formula=formula or ('~group' if folder else '~lab + group'), edger='exact test', edger_glm='GLM test', voom='voom')
    if change:
        change(record)
    write(base / 'record.tsv', 'field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))


DE_FITS = {'hsa-miR-1-3p': (1.8, 0.001, 0.01), 'hsa-miR-199a-3p': (0.1, 0.4, 0.6), 'hsa-miR-21-5p': (-2.5, 0.001, 0.001)}


def evaluate_cmd(project, *args):
    return subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'evaluate', '--project', str(project), *args],
                          capture_output=True, text=True, timeout=120)


def test_evaluate_de_adds_metrics_intervals_and_report(tmp_path):
    project = de_project(tmp_path)
    bench = project / 'benchmark'
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    assert (bench / 'de/input/design.tsv').read_text() == 'sample\tgroup\tlab\nA1\tA\tL1\nA2\tA\tL2\nB1\tB\tL1\nB2\tB\tL2\n'
    assert (bench / 'de/input/library_sizes.tsv').read_text() == ('row\tA1\tA2\tB1\tB2\noracle\t13\t15\t14\t16\n'
                                                                  'unitas\t13.0\t15.0\t14.0\t16.0\nmirge3\t13.0\t15.0\t14.0\t16.0\n')
    assert 'edger_called' not in (bench / 'metrics.tsv').read_text()  # without --de, nothing is scored
    plain_ci = (bench / 'metrics_ci.tsv').read_text().splitlines()
    assert json.loads((bench / 'evaluation.json').read_text())['bootstrap']['skipped'] > 0  # draws without the 1:1 sequence
    de_outputs(project, DE_FITS)
    result = evaluate_cmd(project, '--de')
    assert result.returncode == 0, result.stderr
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['metric'] == ['oracle', 'unitas', 'mirge3']
    assert table['edger_called'] == table['edger_glm_called'] == table['voom_called'] == ['1.0'] * 3
    assert table['edger_false_calls'] == ['0.0'] * 3 and table['edger_offset'] == ['0.1'] * 3
    assert table['edger_fc_bias'] == ['-0.35'] * 3  # median of 1.8 - 2 and -(-2.5 + 3)
    ci = (bench / 'metrics_ci.tsv').read_text().splitlines()
    assert [line for line in ci if not line.startswith(('edger_', 'voom_'))] == plain_ci  # the existing intervals are unchanged
    assert json.loads((bench / 'evaluation.json').read_text())['bootstrap']['skipped'] == 0  # the DE scores use every draw
    assert any(line.startswith('edger_called\tunitas\t1.0\t') for line in ci)
    report = (bench / 'report.md').read_text()
    assert '## Differential expression' in report and '| 4× | 1 | 1.00 | 1.00 | 1.00 |' in report and '| 1× | 1 | 0.00 |' in report
    assert 'design `~lab + group`' in report and hashlib.sha256((bench / 'de/results.tsv').read_bytes()).hexdigest() in report
    assert "**edgeR's GLM likelihood-ratio test**" in report and 'GLM test: GLM test, design `~lab + group`' in report
    assert 'share of its sequences with reads that have a fractional count: unitas 0.0, mirge3 0.0.' in report
    summary_md = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'edger_called: unitas' in summary_md.read_text()


def test_evaluate_de_within_each_lab_scores_each_lab_and_against_the_oracle(tmp_path):
    project = de_project(tmp_path, DE_RUNS_BY_LAB, {'covariates': ['lab'], 'within': 'lab'})
    bench = project / 'benchmark'
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    assert (bench / 'de/levels.tsv').read_text() == 'column\tlevel\nlab\tL1\nlab\tL2\narm\tfixed\narm\ttmm\n'
    # #453: the fixed arm gives every row each run's filtered reads (the reports' kept reads); tmm keeps the workflow's sizes
    fixed = [line.split('\t') for line in (bench / 'de/arm/fixed/input/library_sizes.tsv').read_text().splitlines()]
    assert fixed[0] == ['row', 'A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4'] and fixed[1][1:] == fixed[2][1:] == fixed[3][1:]
    assert fixed[1][1:3] == ['87.0', '97.0']   # 100 and 110 reads in, 13 discarded by the passes
    assert (bench / 'de/arm/tmm/input/library_sizes.tsv').read_text() == (bench / 'de/input/library_sizes.tsv').read_text()
    assert (bench / 'de/arm/tmm/input/model.tsv').read_text() == de.model_text('A', 'B', (), 0, 'tmm')
    assert (bench / 'de/lab/L2/input/design.tsv').read_text() == 'sample\tgroup\nA3\tA\nA4\tA\nB3\tB\nB4\tB\n'  # ~ group
    assert (bench / 'de/lab/L2/input/library_sizes.tsv').read_text().splitlines()[1] == 'oracle\t17\t19\t18\t20'
    assert (bench / 'de/lab/L1/input/counts.tsv').read_text().splitlines()[0] == 'row\tfeature\tA1\tA2\tB1\tB2'
    assert (bench / 'de/lab/L1/input/model.tsv').read_text() == (bench / 'de/input/model.tsv').read_text() == de.model_text('A', 'B')
    assert (bench / 'de/input/design.tsv').read_text().splitlines()[0] == 'sample\tgroup\tlab'  # the whole pool's stays
    plain_ci = (bench / 'metrics_ci.tsv').read_text().splitlines()
    for folder in ('', 'lab/L1'):
        de_outputs(project, DE_FITS, folder=folder)
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and 'no lab/L2/results.tsv or lab/L2/record.tsv in' in result.stderr
    de_outputs(project, DE_FITS, folder='lab/L2')
    result = evaluate_cmd(project, '--de')
    assert result.returncode == 0, result.stderr
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['edger_called.L1'] == table['voom_called.L2'] == table['edger_glm_called'] == ['1.0'] * 3
    assert table['edger_fc_error_vs_oracle.L2'] == table['voom_own_false_calls'] == ['', '0.0', '0.0']  # nothing for the oracle
    ci = (bench / 'metrics_ci.tsv').read_text().splitlines()
    assert [line for line in ci if not line.startswith(('edger_', 'voom_'))] == plain_ci  # the existing intervals are unchanged
    assert any(line.startswith('edger_discordant_vs_oracle.L1\tunitas\t0.0\t') for line in ci)
    report = (bench / 'report.md').read_text()
    main = report.split('## Differential expression')[0]
    assert '| edger_fc_error_vs_oracle |  | 0.0 | 0.0 |' in main and '.L1 |' not in main  # each lab's metrics have their own table
    assert '### Within each lab' in report and '| L2 | called | 1.0 | 1.0 | 1.0 |' in report and '|  | fc_error_vs_oracle |  | 0.0 | 0.0 |' in report
    assert 'L1: 2 A and 2 B runs, design `~group` (edgeR 4.10.7, limma 3.68.5); `de/lab/L1/results.tsv` SHA-256 `{}`.'.format(
        hashlib.sha256((bench / 'de/lab/L1/results.tsv').read_bytes()).hexdigest()) in report
    assert '.L1' not in report.split('## Uncertainty')[1] and '\tunitas\tmirge3\t' in (bench / 'differences.tsv').read_text()
    arms = report.split('### Normalisation arms')[1]   # the whole design side by side in the three normalisations (#453)
    assert '| workflow | called | 1.0 | 1.0 | 1.0 |' in arms and '| fixed | called |' in arms and '| tmm | called |' in arms and '|  | offset |' in arms
    assert table['voom_offset.fixed'] and '### Within each lab' in report.split('### Normalisation arms')[0]
    assert any(line.startswith('edger_fc_shift_vs_oracle.L1\tunitas\tmirge3\t') for line in (bench / 'differences.tsv').read_text().splitlines())
    summary_md = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'edger_fc_error_vs_oracle: unitas' in summary_md.read_text() and 'edger_fc_shift_vs_oracle: unitas' in summary_md.read_text()


@pytest.mark.parametrize('case, message', [
    ('level missing', 'no lab/L1/results.tsv or lab/L1/record.tsv in'),
    ('level stale', 'lab/L2/results.tsv comes from other inputs or another de.R (counts)'),
])
def test_evaluate_de_within_each_lab_stops_before_writing(tmp_path, case, message):
    project = de_project(tmp_path, DE_RUNS_BY_LAB, {'covariates': ['lab'], 'within': 'lab'})
    bench = project / 'benchmark'
    assert evaluate_cmd(project).returncode == 0
    for folder in ('', 'lab/L2') if case == 'level missing' else ('', 'lab/L1', 'lab/L2'):
        de_outputs(project, DE_FITS, (lambda r: r.update(md5_counts='0' * 32)) if folder == 'lab/L2' and case == 'level stale' else None, folder=folder)
    names = ['metrics.tsv', 'report.md', 'evaluation.json', 'per_sequence.tsv', 'de/levels.tsv'] + \
        ['de/{}input/{}'.format(folder, name) for folder in ('', 'lab/L1/', 'lab/L2/') for name in de.INPUTS]
    before = {name: (bench / name).read_text() for name in names}
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and message in result.stderr, result.stderr
    assert {name: (bench / name).read_text() for name in names} == before


def no_de(project):
    path = project / 'benchmark/dataset.yaml'
    spec = yaml.safe_load(path.read_text())
    spec.pop('de')
    path.write_text(yaml.safe_dump(spec))


def equimolar(project):
    no_de(project)
    path = project / 'benchmark/dataset.yaml'
    path.write_text(path.read_text().replace('pool: ratiometric', 'pool: equimolar'))


@pytest.mark.parametrize('case, message', [
    ('equimolar', 'biological datasets, not this equimolar dataset'),
    ('undeclared', "declares no `de:`"),
    ('no results', 'no results.tsv or record.tsv'),
    ('other inputs', 'comes from other inputs or another de.R (counts)'),
    ('other script', 'comes from other inputs or another de.R (script)'),
    ('other results', 'is not the table de.R recorded'),
    ('new counts', 'comes from other inputs or another de.R (counts, library_sizes)'),
    ('incomplete', '1 expected lines missing'),
])
def test_evaluate_de_stops_early(tmp_path, case, message):
    project = de_project(tmp_path)
    bench = project / 'benchmark'
    assert evaluate_cmd(project).returncode == 0
    if case in ('other inputs', 'other script', 'other results', 'new counts', 'incomplete'):
        changes = {'other inputs': lambda r: r.update(md5_counts='0' * 32), 'other script': lambda r: r.update(md5_script='0' * 32)}
        de_outputs(project, DE_FITS, changes.get(case), drop_last=case == 'incomplete')
    if case == 'other results':  # an older table beside the current record, e.g. after a partial restore
        write(bench / 'de/results.tsv', (bench / 'de/results.tsv').read_text().replace('\t1.8\t', '\t2.8\t'))
    if case == 'new counts':  # the method's counts changed after de.smk ran
        write(project / 'data/tmp/smallrna/quant/mirge3/miR.Counts.csv',
              'miRNA,A1,A2,B1,B2\nhsa-miR-1-3p,9,9,2,2\nhsa-miR-199a-3p,4,5,4,5\nhsa-miR-21-5p,1,1,8,9\n')
        isomir_outputs(project / 'data/tmp/smallrna/quant', DE_SEQUENCES)  # the same run's isomiR table
    {'equimolar': equimolar, 'undeclared': no_de}.get(case, lambda p: None)(project)
    names = ['metrics.tsv', 'report.md', 'evaluation.json', 'common.tsv', 'per_sequence.tsv', 'filter_stages.tsv', 'unitas_names.tsv'] + \
        ['de/input/' + name for name in de.INPUTS]
    before = {name: (bench / name).read_text() for name in names}
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and message in result.stderr, result.stderr
    assert {name: (bench / name).read_text() for name in names} == before  # nothing published, no supporting table rewritten


def test_a_biological_dataset_goes_through_the_de_step_without_scores(tmp_path):
    """#447: a biological dataset that declares `de:` gets the inputs of de.smk for its contrast, with a covariate, and
    `evaluate --de` reports what each method tests and calls; nothing is scored, since no change is known."""
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/bio'
    write(data / 'mirbase/23.1/mature.fa', BLOOD)
    reads = ['TAGCAGCACGTAAATATTGGCG'] * 6 + ['TAGCTTATCAGACTGATGTTGA'] * 3
    names = ['GSM{}'.format(i) for i in range(1, 10)]
    groups = dict(zip(names, ['case'] * 4 + ['control'] * 4 + ['other']))
    labs = dict(zip(names, ['L1', 'L1', 'L2', 'L2'] * 2 + ['L1']))
    for i, name in enumerate(names):
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/SRR{}.fastq.gz'.format(i), 'wt') as fh:
            fh.write(''.join('@SRR{}.{} NB1:1:FC:1:11101:1:{}/1\n{}TGGAATTCTC\n+\n{}\n'.format(i, j, j, s, 'I' * (len(s) + 10))
                             for j, s in enumerate(reads)))
    dataset = write(tmp_path / 'bio.yaml', yaml.safe_dump({
        'name': 'bio', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
        'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/bio', 'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'de': {'contrast': ['case', 'control'], 'covariates': ['lab'], 'within': 'lab'},
        'samples': [{'sample': name, 'group': groups[name], 'lab': labs[name],
                     'runs': [{'run': 'SRR{}'.format(i), 'lane': 1, 'md5': 'unused', 'url': 'unused'}]} for i, name in enumerate(names)]}))
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]

    def prepared(folder, *options):
        project = tmp_path / folder
        result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(ext), *options],
                                capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stderr
        samples = list(yaml.safe_load((project / 'config.yaml').read_text())['samples'])
        quant = project / 'data/tmp/smallrna/quant'
        rows = [('miR-16-5p', [700 + 10 * k for k in range(len(samples))]), ('miR-21-5p', [300 - 10 * k for k in range(len(samples))]),
                ('miR-451a', [1] * len(samples))]
        write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(samples) + '\n' + ''.join('{}\t{}\n'.format(n, '\t'.join(map(str, v))) for n, v in rows))
        write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(samples) + '\n' + ''.join('hsa-{},{}\n'.format(n, ','.join(map(str, v)))
                                                                                         for n, v in rows[:2]))
        isomir_outputs(quant)
        write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(ext, 'homo_sapiens', 'Human', {}))
        for sample in samples:
            write(project / 'data/tmp/smallrna/filter/cutadapt/{}_R1.fastq'.format(sample),
                  ''.join('@r{}\n{}\n+\n{}\n'.format(j, s, 'I' * len(s)) for j, s in enumerate(reads)))
            write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', sample, 100)
        return project

    project = prepared('merged')
    bench = project / 'benchmark'
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    design = (bench / 'de/input/design.tsv').read_text().splitlines()
    assert design[:3] == ['sample\tgroup\tlab', 'GSM1\tcase\tL1', 'GSM2\tcase\tL1'] and len(design) == 9  # GSM9 is in neither group
    assert (bench / 'de/input/model.tsv').read_text() == de.model_text('case', 'control', (), 10)
    assert (bench / 'de/levels.tsv').read_text() == 'column\tlevel\nlab\tL1\nlab\tL2\n'
    assert (bench / 'de/lab/L2/input/design.tsv').read_text() == 'sample\tgroup\nGSM3\tcase\nGSM4\tcase\nGSM7\tcontrol\nGSM8\tcontrol\n'
    counts = (bench / 'de/input/counts.tsv').read_text().splitlines()
    assert counts[0] == 'row\tfeature\t' + '\t'.join(names[:8]) and 'unitas\tMIMAT0000069\t' + '\t'.join(str(700 + 10 * k) for k in range(8)) in counts
    assert 'mirge3\tMIMAT0001631' + '\t0' * 8 in counts  # a miRNA only UNITAS reports: every row lists every feature
    assert (bench / 'de/input/library_sizes.tsv').read_text() == ('row\t' + '\t'.join(names[:8]) + '\nunitas' + '\t1001.0' * 8 + '\n'
                                                                  'mirge3' + '\t1000.0' * 8 + '\n')
    assert '## Differential expression' not in (bench / 'report.md').read_text()
    metrics_before = (bench / 'metrics.tsv').read_text()
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and 'no results.tsv or record.tsv' in result.stderr
    fits = {'MIMAT0000069': (0.4, 0.001, 0.01), 'MIMAT0000076': (-0.6, 0.001, 0.02), 'MIMAT0001631': (None, None, None)}
    de_outputs(project, fits)
    de_outputs(project, fits, folder='lab/L1')
    de_outputs(project, dict(fits, MIMAT0000076=(-0.6, 0.2, 0.3)), folder='lab/L2')  # in L2, miR-21-5p is not called
    result = evaluate_cmd(project, '--de')
    assert result.returncode == 0, result.stderr
    report = (bench / 'report.md').read_text()
    assert '## Differential expression' in report and 'case against control over the 3 miRNAs any method reports' in report
    whole, lab1, lab2 = report.split('## Differential expression')[1].split('\n**')[1:4]  # a table per DE set, each lab's too
    assert whole.startswith('All samples of the contrast** (4 case and 4 control samples)') and lab2.startswith('L2** (2 case and 2 control samples)')
    assert "| edgeR's exact test | tested | 2 | 2 |" in whole and '|  | up | 1 | 1 |' in whole and '|  | down | 1 | 1 |' in whole
    assert '|  | down | 1 | 1 |' in lab1 and '|  | down | 0 | 0 |' in lab2 and '|  | up | 1 | 1 |' in lab2
    for folder in ('', 'lab/L1/', 'lab/L2/'):
        sha = hashlib.sha256((bench / 'de' / folder / 'results.tsv').read_bytes()).hexdigest()
        assert '`de/{}results.tsv` SHA-256 `{}`'.format(folder, sha) in report
    assert (bench / 'metrics.tsv').read_text() == metrics_before  # the calls are reported, not scored
    lanes = prepared('lanes', '--lanes')
    assert evaluate_cmd(lanes).returncode == 0 and not (lanes / 'benchmark/de').exists()  # lanes of one library are no replicates
    result = evaluate_cmd(lanes, '--de')
    assert result.returncode != 0 and 'was prepared with --lanes' in result.stderr


# isomiR level (#454): the routes' per-sequence counts, variants and scores.

import isomirs  # noqa: E402

MMU_192, MMU_194, MMU_204 = 'CTGACCTATGAATTGACAGCC', 'TGTAACAGCAACTCCATGTGGA', 'TTCCCTTTGTCATCCTATGCCT'
MMU_192_T = MMU_192[:4] + 'T' + MMU_192[5:]  # a substitution
MMU_MATURE = {'mmu-miR-192-5p': ('MIMAT0000517', MMU_192), 'mmu-miR-194-5p': ('MIMAT0000224', MMU_194), 'mmu-miR-204-5p': ('MIMAT0000237', MMU_204)}
MIRTOP_COLUMNS = '\t'.join(adapters.MIRTOP_ANNOTATION)


def native_isomir_outputs(quant):
    """Each isomiR route's count table and per-sequence output for samples S01 and S02, consistent as the routes write
    them; returns mirflowz's annotation."""
    write(quant / 'unitas/mir_counts.tsv', '\tS01\tS02\nmiR-192-5p\t45\t45\nmiR-194-2-5p\t5\t0\n')
    # the read without its non-templated tail, its reads split by tail; a substitution's name keeps its modification
    write(quant / 'unitas/unitas/S01/unitas.miR-table_Mouse.txt', UNITAS_MIR_TABLE + 'miR-192-5p\t{0}\t40\t35\t5\t0\t5\t0\n'
          'miR-192-5p(5C->T) \t{2}\t5\t5\t0\t0\t0\t0\nmiR-194-2-5p\t{1}\t5\t5\t0\t0\t0\t0\n'.format(MMU_192, MMU_194, MMU_192_T))
    write(quant / 'unitas/unitas/S02/unitas.miR-table_Mouse.txt', UNITAS_MIR_TABLE + 'miR-192-5p\t{}\t45\t45\t0\t0\t0\t0\n'.format(MMU_192))
    write(quant / 'unitas/unitas/S02/unitas.miR-table_non-Mouse.txt', UNITAS_MIR_TABLE + 'miR-192-5p\t{}\t9\t9\t0\t0\t0\t0\n'.format(MMU_192))
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S01,S02\nmmu-miR-192-5p,45,45\nmmu-miR-194-5p/204-5p,25,25\nmmu-miR-9999-5p,0,0\n')
    write(quant / 'mirge3/mapped.csv', 'Sequence,annotFlag,exact miRNA,hairpin miRNA,mature tRNA,primary tRNA,snoRNA,rRNA,ncrna others,'
          'mRNA,isomiR miRNA,S01,S02\n'
          '{0},1,mmu-miR-192-5p,,,,,,,,,40,45\n{0}AA,1,,,,,,,,,mmu-miR-192-5p,5,0\n'
          '{1},1,mmu-miR-204-5p,,,,,,,,,25,25\n'   # a merge group member: counted on the group
          'AAAACCCCGGGGTTTTACGTAC,1,mmu-miR-9999-5p,,,,,,,,,1,0\n'  # one exact read in the project: the route drops the miRNA
          'GCATTGGTGGTTCAGTGGTAGA,1,,,tRNA-Gly,,,,,,,9,9\n'.format(MMU_192, MMU_194))
    write(quant / 'mirge3/lib/mouse/annotation.Libs/mouse_merges_miRBase.csv', 'mmu-miR-194-5p/204-5p,mmu-miR-194-5p,mmu-miR-204-5p\n')
    write(quant / 'smrnaseq/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\t45\nmmu-miR-194-5p\t25\t0\n')
    write(quant / 'smrnaseq/joined_samples_mirtop.tsv', MIRTOP_COLUMNS + '\tS01\tS02\n'
          'iso-1\t{0}\tmmu-miR-192-5p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t40\t45\n'
          'iso-2\t{0}AA\tmmu-miR-192-5p\tiso_add3p:2\t0\t0\t2\t0\t0\t0\tAA\t0\t5\t0\n'
          'iso-3\t{1}\tmmu-miR-194-5p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t25\t0\n'
          'iso-3\t{1}\tmmu-miR-204-5p\tNA\t0\t0\t0\t0\t0\t0\t0\t0\t25\t0\n'.format(MMU_192, MMU_194))  # a UID's second row is not counted
    write(quant / 'biocore_miraligner/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t45\t44.5\nmmu-miR-194-5p\t12.5\t0\nmmu-miR-204-5p\t12.5\t0\n')
    write(quant / 'biocore_miraligner/S01/isomir.tsv', 'isomiR\tcount\nmmu-miR-192-5p_{0}_0;0;0;0\t40\nmmu-miR-192-5p_{0}AA_0;0;0;aa\t5\n'
          'mmu-miR-192-5p_{2}_5TA;0;0;0\t3\nmmu-miR-194-5p_{1}_0;0;0;0\t12.5\nmmu-miR-204-5p_{1}_0;0;0;0\t12.5\n'.format(
              MMU_192, MMU_194, MMU_192.replace('ACC', 'TCC', 1)))  # a substitution: out of the mature table
    write(quant / 'biocore_miraligner/S02/isomir.tsv', 'isomiR\tcount\nmmu-miR-192-5p_{}_0;0;0;0\t44.5\n'.format(MMU_192))
    # a plain row: the reads of the annotated mature; isomiR rows in genome orientation (miR-204-5p on the minus strand)
    write(quant / 'mirflowz/all_mirna_counts.tab', 'ID\tS01\tS02\nmmu-miR-192-5p\t40\t40\nmmu-miR-192-5p|0|2|23M|23|{}AA\t5\t5\n'
          'mmu-miR-204-5p|0|0|22M|22|{}\tNA\t3\n'.format(MMU_192, adapters.reverse_complement(MMU_204)))
    write_records(quant, ['unitas', 'mirge3', 'smrnaseq', 'biocore_miraligner', 'mirflowz'])
    return write(quant.parent / 'mirflowz.gff3', '19\t.\tmiRNA\t1\t30\t.\t+\t.\tID=MIMAT0000224_1;Alias=MIMAT0000224;Name=mmu-miR-194-2-5p\n'
                 '18\t.\tmiRNA\t1\t33\t.\t+\t.\tID=MIMAT0000517;Alias=MIMAT0000517;Name=mmu-miR-192-5p\n'
                 '2\t.\tmiRNA\t1\t34\t.\t-\t.\tID=MIMAT0000237;Alias=MIMAT0000237;Name=mmu-miR-204-5p\n')


def test_isomir_table_counts_as_each_route_counts_its_mirnas(tmp_path):
    quant = tmp_path / 'quant'
    gff = native_isomir_outputs(quant)
    methods = ['unitas', 'mirge3', 'smrnaseq', 'biocore_miraligner', 'mirflowz', 'tesmall']
    rows, without = adapters.isomir_table(methods, quant, NAME_TO_ID, 'mmu-', gff, mature=MMU_MATURE)
    assert without == ['tesmall']
    got = {(m, s, n, q): (i, label, c) for m, s, n, i, q, label, c in rows}
    assert got[('unitas', 'S01', 'miR-192-5p', MMU_192 + 'A')] == ('MIMAT0000517', '', 5.0)  # the read rebuilt with its tail; no labels
    assert got[('unitas', 'S01', 'miR-192-5p', MMU_192)][2] == 35 and got[('unitas', 'S01', 'miR-192-5p', MMU_192_T)][2] == 5  # modification dropped
    assert got[('mirge3', 'S01', 'mmu-miR-194-5p/204-5p', MMU_194)] == ('MIMAT0000224', 'canonical', 25.0)  # folded on the group
    assert got[('mirge3', 'S01', 'mmu-miR-192-5p', MMU_192 + 'AA')] == ('MIMAT0000517', 'isomiR', 5.0)
    assert not any(m == 'mirge3' and n == 'mmu-miR-9999-5p' for m, _, n, _ in got)  # fewer than 2 exact reads
    assert ('smrnaseq', 'S01', 'mmu-miR-204-5p', MMU_194) not in got and got[('smrnaseq', 'S01', 'mmu-miR-194-5p', MMU_194)][2] == 25
    assert got[('smrnaseq', 'S01', 'mmu-miR-192-5p', MMU_192 + 'AA')][1] == 'isomiR' and got[('smrnaseq', 'S01', 'mmu-miR-192-5p', MMU_192)][1] == 'canonical'
    assert got[('biocore_miraligner', 'S01', 'mmu-miR-192-5p', MMU_192 + 'AA')] == ('MIMAT0000517', 'isomiR', 5.0)
    assert not any(m == 'biocore_miraligner' and q.startswith('CTGTCC') for m, _, _, q in got)  # the substitution
    assert got[('mirflowz', 'S01', 'mmu-miR-192-5p', MMU_192)] == ('MIMAT0000517', 'canonical', 40.0)  # a plain row: the mature
    assert got[('mirflowz', 'S02', 'mmu-miR-204-5p', MMU_204)] == ('MIMAT0000237', 'isomiR', 3.0)  # minus strand: reverse complemented
    assert ('mirflowz', 'S01', 'mmu-miR-204-5p', MMU_204) not in got
    common = adapters.common_table(methods[:-1], quant, NAME_TO_ID, 'mmu-', gff)
    adapters.check_isomirs(rows, common, methods[:-1])
    with pytest.raises(ValueError, match='need its annotation'):
        adapters.isomir_table(['mirflowz'], quant, NAME_TO_ID, 'mmu-')
    with pytest.raises(ValueError, match='mirflowz miRNA mmu-miR-192-5p is not on a strand'):
        adapters.isomir_table(['mirflowz'], quant, NAME_TO_ID, 'mmu-', write(tmp_path / 'other.gff3', gff.read_text().splitlines()[0] + '\n'))
    write(quant / 'unitas/unitas/S02/unitas.miR-table_Mouse.txt', UNITAS_MIR_TABLE + 'miR-192-5p\t{}\t45\t40\t0\t0\t0\t0\n'.format(MMU_192))
    with pytest.raises(ValueError, match="total_reads 45, but \\+N0-3' and the tails add up to 40"):
        adapters.isomir_table(['unitas'], quant, NAME_TO_ID, 'mmu-')
    write(quant / 'smrnaseq/mirna.tsv', 'miRNA\tS01\tS02\nmmu-miR-192-5p\t46\t45\nmmu-miR-194-5p\t25\t0\n')
    with pytest.raises(ValueError, match=r'smrnaseq isomiR rows do not add up .* smrnaseq S01 mmu-miR-192-5p: 45 in the isomiR rows, 46'):
        adapters.check_isomirs(rows, adapters.common_table(methods[:-1], quant, NAME_TO_ID, 'mmu-', gff), methods[:-1])


LET7 = 'TGAGGTAGTAGGTTGTATAGTT'
LET7_HAIRPIN = 'CCGGA' + LET7 + 'TTCAGG'


@pytest.mark.parametrize('read, hairpin, expected, kind', [
    (LET7, LET7_HAIRPIN, (0, 0, '', 0), 'canonical'),
    (LET7[:-2], None, (0, -2, '', 0), '3p_trim'),
    (LET7 + 'T', LET7_HAIRPIN, (0, 1, '', 0), '3p_templated'),       # the hairpin's next base
    (LET7 + 'A', LET7_HAIRPIN, (0, 0, 'A', 0), '3p_add'),            # not the hairpin's
    (LET7 + 'T', None, (0, 0, 'T', 0), '3p_longer'),                 # without the hairpin, templated or not
    ('A' + LET7, LET7_HAIRPIN, (-1, 0, '', 0), '5p'),
    (LET7[1:], None, (1, 0, '', 0), '5p'),
    (LET7[:8] + 'C' + LET7[9:], None, (0, 0, '', 1), 'substitution'),
    (LET7[:-5], None, None, 'other'),                                # beyond 3 nt
    ('A' + LET7, None, None, 'other'),                               # upstream of the mature, no hairpin
])
def test_isomir_variant_places_the_read_on_its_mature(read, hairpin, expected, kind):
    v = isomirs.variant(read, LET7, hairpin)
    assert (tuple(v) if v else None) == expected and isomirs.kind(v, templated_known=hairpin is not None) == kind


def test_isomir_synthetic_scores_against_the_truth():
    a, b, random_read = LET7, 'TAGCTTATCAGACTGATGTTGA', 'ACGTACGTACGTACGTACGT'
    truth = [{'sample': 'M1', 'class': 'canonical', 'mature_id': 'A', 'sequence': a, 'count': '10', 'errors': '0'},
             {'sample': 'M1', 'class': '3p_trim_1', 'mature_id': 'A', 'sequence': a[:-1], 'count': '4', 'errors': '0'},
             {'sample': 'M1', 'class': '3p_add_A', 'mature_id': 'A', 'sequence': a + 'A', 'count': '2', 'errors': '0'},
             {'sample': 'M1', 'class': 'canonical', 'mature_id': 'A', 'sequence': 'C' + a[1:], 'count': '1', 'errors': '1'},
             {'sample': 'M1', 'class': 'random', 'mature_id': '', 'sequence': random_read, 'count': '3', 'errors': '0'},
             {'sample': 'M1', 'class': 'canonical', 'mature_id': 'B', 'sequence': b, 'count': '5', 'errors': '0'},
             {'sample': 'P1', 'class': 'canonical', 'mature_id': 'A', 'sequence': a, 'count': '99', 'errors': '0'}]
    samples = [{'sample': 'M1', 'kind': 'mixed'}, {'sample': 'P1', 'kind': 'probe'}]
    rows = [('m', 'M1', 'a', 'A', a, 'canonical', 10), ('m', 'M1', 'a', 'A', a[:-1], 'canonical', 3), ('m', 'M1', 'a', 'A', a + 'A', 'isomiR', 2),
            ('m', 'M1', 'b', 'B', random_read, 'isomiR', 1), ('m', 'M1', 'b', 'B', b, 'isomiR', 5), ('m', 'M1', 'a', 'A', 'GGGGGGGGGGGGGGGGGG', 'isomiR', 2),
            ('m', 'P1', 'a', 'A', a, 'canonical', 50), ('u', 'M1', 'a', 'A', a[:-1], '', 4)]
    values, classes = isomirs.synthetic_scores(truth, rows, samples, ['m', 'u'], {})
    assert classes == {'canonical': 15, '3p_trim_1': 4, '3p_add_A': 2, 'substituted': 1}  # mixed samples only
    m = values['m']
    assert [m['isomir_recovered_' + c] for c in ('canonical', '3p_trim_1', '3p_add_A', 'substituted')] == [1.0, 0.75, 1.0, 0.0]
    assert m['isomir_recovered'] == 0.8333  # (3 + 2) of the 6 true isomiR reads
    assert m['isomir_misassigned'] == 0.0435 and m['isomir_not_in_sample'] == 0.087  # 1 and 2 of its 23 reads in M1
    assert m['isomir_as_canonical'] == 0.2308 and m['canonical_as_isomir'] == 0.7143  # 3 of 13 canonical labels, 5 of 7 isomiR labels
    assert values['u']['isomir_recovered'] == 0.6667 and 'isomir_as_canonical' not in values['u']  # no labels, no label metrics
    assert isomirs.synthetic_keys(classes)[-4:] == ['isomir_recovered_canonical', 'isomir_recovered_3p_add_A', 'isomir_recovered_3p_trim_1',
                                                   'isomir_recovered_substituted']


def test_isomir_pool_scores_against_the_oracle_and_the_artefact_floor():
    exact = [{'pool_id': 'E1', 'sequence': LET7, 'mature_ids': ['A']}]
    shifted = [{'pool_id': 'H1', 'sequence': 'TAGCTTATCAGACTGATGTTG', 'mature_ids': ['B']}]
    oracle = {'S1': {LET7: 10, 'TAGCTTATCAGACTGATGTTG': 4}, 'S2': {LET7: 6}}
    rows = [('m', 'S1', 'a', 'A', LET7, '', 8), ('m', 'S2', 'a', 'A', LET7, '', 6), ('m', 'S1', 'b', 'B', 'TAGCTTATCAGACTGATGTTG', '', 5),
            ('m', 'S1', 'a', 'A', LET7 + 'A', '', 2), ('m', 'S1', 'c', 'C', 'TTTTTTTTTTTTTTTTTTTT', '', 9), ('m', 'S1', 'a', 'A', '', '', 7)]
    values, kinds = isomirs.pool_scores(rows, exact, shifted, exact + shifted, oracle, ['m'], {'A': LET7, 'B': 'TAGCTTATCAGACTGATGTTGA'})
    assert values['m'] == {'exact_sequence_recovered': 0.875, 'shifted_isomir_recovered': 1.0, 'shifted_isomir_within_2fold': 1.0,
                           'artefact_share': 0.0952}  # 14 of 16; 4 of 4; 6/5; 2 of the 21 reads on pool miRNAs
    assert kinds['m'] == {'3p_longer': 2}


def test_isomir_biological_scores_describe_the_routes():
    rows = [('m', 'S1', 'a', 'A', LET7, '', 30), ('m', 'S1', 'a', 'A', LET7[:-1], '', 10), ('m', 'S2', 'a2', 'A2', LET7, '', 5),
            ('n', 'S1', 'a', 'A', LET7, '', 40), ('n', 'S1', 'a', 'A', LET7 + 'T', '', 10), ('m', 'S1', 'x', '', 'ACGTACGTACGTACGTACGT', '', 7)]
    id_to_seq = {'A': LET7, 'A2': LET7}
    values, kinds, overlap = isomirs.biological_scores(rows, ['m', 'n'], id_to_seq, isomirs.representatives(id_to_seq))
    assert values == {'m': {'isomir_share': 0.2222, 'isomirs_detected': 2}, 'n': {'isomir_share': 0.2, 'isomirs_detected': 2}}
    assert kinds['m'] == {'canonical': 35, '3p_trim': 10} and kinds['n'] == {'canonical': 40, '3p_longer': 10}
    assert overlap == {('m', 'n'): 0.3333}  # the mature sequence of the 3 detected


def test_isomir_tallies_add_up_a_bootstrap_draw():
    """#459: a resample of the units (repeats allowed) adds up their tallies; fixed units count in every draw."""
    tally = isomirs.Tally({'a': {'x': 1}, 'b': {'x': 2}, '': {'x': 10}}, ['x'], fixed=('',))
    assert tally.sum() == {'x': 13} and tally.sum(['a', 'a']) == {'x': 12} and tally.sum(['b', 'not-a-unit']) == {'x': 12}
    exact = [{'pool_id': 'E1', 'sequence': LET7, 'mature_ids': ['A']}]
    shifted = [{'pool_id': 'H1', 'sequence': 'TAGCTTATCAGACTGATGTTG', 'mature_ids': ['B']}]
    oracle = {'S1': {LET7: 10, 'TAGCTTATCAGACTGATGTTG': 4}, 'S2': {LET7: 6}}
    rows = [('m', 'S1', 'a', 'A', LET7, '', 8), ('m', 'S2', 'a', 'A', LET7, '', 6), ('m', 'S1', 'b', 'B', 'TAGCTTATCAGACTGATGTTG', '', 5),
            ('m', 'S1', 'a', 'A', LET7 + 'A', '', 2)]
    tallies, _ = isomirs.pool_tallies(rows, exact, shifted, exact + shifted, oracle, ['m'], {'A': LET7})
    assert isomirs.pool_values(tallies) == isomirs.pool_scores(rows, exact, shifted, exact + shifted, oracle, ['m'], {'A': LET7})[0]
    assert isomirs.pool_values(tallies, exact * 2, shifted, ['S1']) == {'m': {
        'exact_sequence_recovered': 0.8,          # 8 of the oracle's 10, twice
        'shifted_isomir_recovered': 1.0, 'shifted_isomir_within_2fold': 1.0,
        'artefact_share': 0.16}}                  # 2 + 2 of the 10 + 10 + 5 reads on pool miRNAs in S1


# Null comparisons (#449): the samples of one group split in two, so every call is false.

def test_null_splits_halve_each_group_balanced_on_the_design():
    units = [('S{}'.format(i), 'case' if i < 6 else 'ctrl', {'lab': 'L1' if i % 2 else 'L2'}) for i in range(10)]
    units += [('T0', 'tiny', {'lab': 'L1'}), ('T1', 'tiny', {'lab': 'L2'})]
    splits = de.null_splits(units, 50, 1, ['lab'])
    assert splits == de.null_splits(units, 50, 1, ['lab']) and de.null_splits(units, 3, 1, ['lab']) != de.null_splits(units, 3, 2, ['lab'])
    assert set(splits) == {'case', 'ctrl'}  # tiny has fewer than 4 samples
    assert len(splits['case']) == 9 and len(splits['ctrl']) == 2  # every distinct balanced split: 3 x 3, and 2 x 2 halved by mirroring
    lab = {s: f['lab'] for s, _, f in units}
    for x, y in splits['case']:  # two odd strata (3 L1, 3 L2): their odd samples go to the halves in turn
        assert sorted(x + y) == ['S{}'.format(i) for i in range(6)] and sorted(lab[s] for s in x) == ['L1', 'L2', 'L2']
    assert len({tuple(sorted([tuple(x), tuple(y)])) for x, y in splits['ctrl']}) == 2  # never a split and its mirror
    assert len(de.null_splits(units, 3, 1)['case']) == 3  # without columns: plain random halves


def test_null_scores_count_every_call_as_false():
    features = ['F1', 'F2', 'F3']

    def results(fits):
        return {('m', test, f): fit for test in de.TESTS for f, fit in fits.items()}

    splits = [(results({'F1': (1.0, 0.001, 0.01), 'F2': (0.1, 0.04, 0.3), 'F3': (None, None, None)}), features),  # 1 of 2 tested called
              (results({'F1': (0.2, 0.5, 0.9), 'F2': (0.1, 0.6, 0.9), 'F3': (0.0, 0.9, 0.9)}), features)]          # nothing called
    values, spread, called = de.null_scores(splits, ['m'])
    assert values['m']['edger_null_called'] == 0.25 and spread['m']['voom'] == (0.0, 0.5)  # the median of 1/2 and 0/3
    assert values['m']['voom_null_any'] == 0.5 and values['m']['edger_glm_null_p05'] == 0.4  # 2 of the 5 tested p-values
    assert called['edger'] == {('m', 'F1'): 1} and sorted(de.null_keys())[:3] == ['edger_glm_null_any', 'edger_glm_null_called', 'edger_glm_null_p05']
    lines = de.null_lines(values, spread, called, ['m'], [('case', 6, 2)], lambda f: 'name-' + f, lambda row, f: 12.0)
    assert '| null_called (range) | 0.25 (0.0–0.5) |' in lines and '| name-F1 | 1/1/1 (12) |' in lines and 'case (6 samples) 2' in lines[3]


def test_null_comparisons_in_a_declaration(tmp_path):
    groups = ('case',) * 4 + ('control',) * 4 + ('other',)
    spec = benchmark.load_dataset(biological_declaration(tmp_path, {'null': {'splits': 5, 'seed': 1}, 'covariates': ['patient']}, groups))
    assert benchmark.de_contrast(spec) is None  # null comparisons only
    rows = [{'sample': s['sample'], 'library': s['sample'], 'group': str(s['group']), 'lane': 'merged'} for s in spec['samples']]
    sets = benchmark.de_sets(spec, benchmark.de_units(spec, rows))
    assert [s.folder for s in sets] == ['null/case-1', 'null/case-2', 'null/control-1', 'null/control-2']  # other: one sample
    assert sorted(g for _, g, _ in sets[0].design) == ['x', 'x', 'y', 'y'] and sets[0].samples == '2 x and 2 y samples of case'
    assert sets[0].design[0][2] == {'patient': 'P0'} and sets[0].model == de.model_text('x', 'y', (), biological.MIN_READS)
    for declared, message in (({'null': {'splits': 0, 'seed': 1}}, 'de.null takes `splits`'), ({'null': {'splits': 5}}, 'de.null takes `splits`'),
                              ({'null': {'splits': 5, 'seed': 1}}, 'de.null needs a group with at least 4 samples')):
        with pytest.raises(SystemExit, match=re.escape(message)):
            benchmark.load_dataset(biological_declaration(tmp_path, declared))
    spec = benchmark.load_dataset(biological_declaration(tmp_path, {'null': {'splits': 5, 'seed': 1}, 'within': 'patient'}, groups))
    (first, *_) = benchmark.de_sets(spec, benchmark.de_units(spec, rows))
    assert first.design[0][2] == {} and sorted(g for s, g, _ in first.design if s in ('S0', 'S2')) == ['x', 'y']  # within balances only
    members = [('S0', 'case', {'site': 'A', 'patient': 'P0'}), ('S1', 'case', {'site': 'A', 'patient': 'P1'})]
    assert benchmark.null_covariates(members, ['site', 'patient']) == ['patient']  # the group's samples share a site


def biological_de_project(tmp_path, declared):
    """A prepared biological project whose `de:` is `declared`: GSM1-GSM9, 4 case, 4 control and 1 other, over labs L1 and
    L2, with UNITAS and miRge3 counts and the filtered reads."""
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/bio'
    write(data / 'mirbase/23.1/mature.fa', BLOOD)
    reads = ['TAGCAGCACGTAAATATTGGCG'] * 6 + ['TAGCTTATCAGACTGATGTTGA'] * 3
    names = ['GSM{}'.format(i) for i in range(1, 10)]
    groups = dict(zip(names, ['case'] * 4 + ['control'] * 4 + ['other']))
    labs = dict(zip(names, ['L1', 'L1', 'L2', 'L2'] * 2 + ['L1']))
    for i, name in enumerate(names):
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/SRR{}.fastq.gz'.format(i), 'wt') as fh:
            fh.write(''.join('@SRR{}.{} NB1:1:FC:1:11101:1:{}/1\n{}TGGAATTCTC\n+\n{}\n'.format(i, j, j, s, 'I' * (len(s) + 10))
                             for j, s in enumerate(reads)))
    dataset = write(tmp_path / 'bio.yaml', yaml.safe_dump({
        'name': 'bio', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'Bioo Scientific NEXTflex Small RNA-Seq Kit v4',
        'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/bio', 'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'de': declared, 'samples': [{'sample': name, 'group': groups[name], 'lab': labs[name],
                                     'runs': [{'run': 'SRR{}'.format(i), 'lane': 1, 'md5': 'unused', 'url': 'unused'}]} for i, name in enumerate(names)]}))
    project = tmp_path / 'project'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(ext)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    samples = list(yaml.safe_load((project / 'config.yaml').read_text())['samples'])
    quant = project / 'data/tmp/smallrna/quant'
    rows = [('miR-16-5p', [700 + 10 * k for k in range(len(samples))]), ('miR-21-5p', [300 - 10 * k for k in range(len(samples))]),
            ('miR-451a', [1] * len(samples))]
    write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(samples) + '\n' + ''.join('{}\t{}\n'.format(n, '\t'.join(map(str, v))) for n, v in rows))
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(samples) + '\n' + ''.join('hsa-{},{}\n'.format(n, ','.join(map(str, v)))
                                                                                     for n, v in rows[:2]))
    isomir_outputs(quant)
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(ext, 'homo_sapiens', 'Human', {}))
    for sample in samples:
        write(project / 'data/tmp/smallrna/filter/cutadapt/{}_R1.fastq'.format(sample),
              ''.join('@r{}\n{}\n+\n{}\n'.format(j, s, 'I' * len(s)) for j, s in enumerate(reads)))
        write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', sample, 100)
    return project


def test_evaluate_reports_null_comparisons(tmp_path):
    project = biological_de_project(tmp_path, {'contrast': ['case', 'control'], 'covariates': ['lab'], 'within': 'lab',
                                               'null': {'splits': 3, 'seed': 1}})
    bench = project / 'benchmark'
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    nulls = ['null/case-1', 'null/case-2', 'null/control-1', 'null/control-2']  # labs halved: 2 distinct splits of each group
    assert (bench / 'de/levels.tsv').read_text() == 'column\tlevel\nlab\tL1\nlab\tL2\n' + ''.join('{}\t{}\n'.format(*f.split('/')) for f in nulls)
    design = (bench / 'de/null/case-1/input/design.tsv').read_text().splitlines()
    assert design[0] == 'sample\tgroup\tlab' and sorted(line.split('\t')[1] for line in design[1:]) == ['x', 'x', 'y', 'y']
    assert {line.split('\t')[2] for line in design[1:] if line.split('\t')[1] == 'x'} == {'L1', 'L2'}  # one of each lab in a half
    assert (bench / 'de/null/case-1/input/model.tsv').read_text() == de.model_text('x', 'y', (), biological.MIN_READS)
    assert (bench / 'de/input/design.tsv').read_text().splitlines()[:2] == ['sample\tgroup\tlab', 'GSM1\tcase\tL1']  # the contrast's, as before
    fits = {'MIMAT0000069': (0.4, 0.001, 0.01), 'MIMAT0000076': (-0.6, 0.3, 0.6), 'MIMAT0001631': (None, None, None)}
    for folder in ['', 'lab/L1', 'lab/L2'] + nulls[:-1]:
        de_outputs(project, fits, folder=folder)
    before = (bench / 'metrics.tsv').read_text()
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and 'null/control-2/results.tsv' in result.stderr and (bench / 'metrics.tsv').read_text() == before
    de_outputs(project, fits, folder=nulls[-1])
    result = evaluate_cmd(project, '--de')
    assert result.returncode == 0, result.stderr
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['edger_null_called'] == ['0.5', '0.5'] and table['voom_null_any'] == ['1.0', '1.0']  # 1 of 2 tested, in every split
    assert table['edger_glm_null_p05'] == ['0.5', '0.5'] and not any(line.split('\t')[0].endswith(('_null_called', '_null_any'))
                                                                    for line in (bench / 'metrics_ci.tsv').read_text().splitlines())
    report = (bench / 'report.md').read_text()
    assert '## Null comparisons' in report and 'case (4 samples) 2; control (4 samples) 2' in report
    assert '| null_called (range) | 0.5 (0.5–0.5) | 0.5 (0.5–0.5) |' in report and '| hsa-miR-16-5p | 4/4/4 (' in report
    assert '## Differential expression' in report  # the contrast is still described
    summary_md = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'edger_null_called: unitas' in summary_md.read_text()


def test_evaluate_null_comparisons_without_a_contrast(tmp_path):
    project = biological_de_project(tmp_path, {'null': {'splits': 1, 'seed': 7}})
    bench = project / 'benchmark'
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    assert not (bench / 'de/input').exists() and (bench / 'de/levels.tsv').read_text() == 'column\tlevel\nnull\tcase-1\nnull\tcontrol-1\n'
    assert (bench / 'de/null/case-1/input/design.tsv').read_text().splitlines()[0] == 'sample\tgroup'  # no covariates declared
    for folder in ('null/case-1', 'null/control-1'):
        de_outputs(project, {'MIMAT0000069': (0.1, 0.5, 0.9), 'MIMAT0000076': (0.2, 0.04, 0.9), 'MIMAT0001631': (None, None, None)}, folder=folder)
    result = evaluate_cmd(project, '--de')
    assert result.returncode == 0, result.stderr
    report = (bench / 'report.md').read_text()
    assert '## Differential expression' not in report and 'No route called any miRNA in any split.' in report
    table = {line.split('\t')[0]: line.split('\t')[1:] for line in (bench / 'metrics.tsv').read_text().splitlines()}
    assert table['voom_null_called'] == ['0.0', '0.0'] and table['edger_null_p05'] == ['0.5', '0.5']


def test_mirflowz_annotation_is_the_one_the_project_run_recorded(tmp_path):
    """#440: the extended annotation is in the mirflowz reference's gff-<ID> folder that quant_reference records; a
    record written before #440 has no mirna_gff_id, and its run read the annotation from the reference's prepared/."""
    project, reference = tmp_path / 'project', tmp_path / 'ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus'
    record = project / 'data/tmp/smallrna/quant/mirflowz/reference.tsv'
    name = 'extended_mirna_annotation_6_nt.gff3'
    with pytest.raises(SystemExit, match='--mirflowz-gff'):
        benchmark.mirflowz_annotation(project)
    write(record, 'field\tvalue\tsource\nmirna_database\tmirbase\tx\nreference\t{}\tconfiguration\nmirna_gff_id\tabc123\tconfiguration\n'.format(reference))
    with pytest.raises(SystemExit, match=re.escape('no mirflowz annotation at {}'.format(reference / 'gff-abc123/annotation' / name))):
        benchmark.mirflowz_annotation(project)
    for folder in ('gff-abc123/annotation', 'gff-other/annotation', 'prepared'):  # other corrections' and the old layout's too
        write(reference / folder / name, '')
    assert benchmark.mirflowz_annotation(project) == reference / 'gff-abc123/annotation' / name
    write(record, 'field\tvalue\tsource\nmirna_database\tmirbase\tx\nreference\t{}\tconfiguration\n'.format(reference))
    assert benchmark.mirflowz_annotation(project) == reference / 'prepared' / name


# Every route on every dataset, and a state for each declared route (#448)

def test_every_declaration_lists_the_eight_routes():
    for path in sorted((BENCH / 'datasets').glob('*.yaml')):
        assert sorted(benchmark.load_dataset(path)['methods']) == sorted(adapters.ADAPTERS), path.name


def test_gse317995_is_declared_in_full_and_per_group_takes_the_quick_run(tmp_path):
    spec = benchmark.load_dataset(BENCH / 'datasets/gse317995_nextflex_v4.yaml')
    assert len(spec['samples']) == 30 and collections.Counter(s['group'] for s in spec['samples']) == {
        'Adenomyosis': 10, 'Fertile control': 10, 'Infertile control': 10}
    assert all(sorted(run['lane'] for run in s['runs']) == [1, 2, 3, 4] for s in spec['samples'])
    assert len(fill.targets(spec, tmp_path)) == 1 + 3 + 120  # mature.fa, the other_rna sources (#415), the runs
    quick = benchmark.select_per_group(spec, 2)  # the two lowest-numbered of each group, as declared before #448
    assert [s['sample'] for s in quick['samples']] == ['GSM9484173', 'GSM9484174', 'GSM9484193', 'GSM9484194', 'GSM9484183', 'GSM9484184']
    assert spec['de'] == {'null': {'splits': 50, 'seed': 1}} and 'de' not in quick   # #449: halves need 4 samples per group
    assert benchmark.select_per_group(spec, 4)['de'] == spec['de']
    with pytest.raises(SystemExit, match='--per-group is for biological datasets'):
        benchmark.select_per_group(benchmark.load_dataset(BENCH / 'datasets/giraldez2018_equimolar_truseq.yaml'), 2)


def test_routes_not_evaluated_get_a_state_in_the_report_and_the_summary(tmp_path):
    """A route left out at prepare is `not run`; a route whose checks fail is `stopped` and the others are scored
    without it; benchmark/not_run.tsv gives the reason of a route not run by hand; the summary shows the state."""
    ext, dataset, reads = biological_dataset(tmp_path)
    spec = yaml.safe_load(dataset.read_text())
    write(dataset, yaml.safe_dump(dict(spec, methods=['unitas', 'mirge3', 'tinyrna'])))
    project = tmp_path / 'project'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(ext), '--methods', 'unitas,mirge3'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    saved = yaml.safe_load((project / 'benchmark/dataset.yaml').read_text())
    assert saved['methods'] == ['unitas', 'mirge3'] and saved['not_selected'] == ['tinyrna']
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tGSM1\nmiR-16-5p\t700\nmiR-21-5p\t300\n')   # miRge3 did not finish: no count table
    write_records(quant, ['unitas'], dump=unitas_dump(ext, 'homo_sapiens', 'Human', {}))
    isomir_outputs(quant)
    write(project / 'data/tmp/smallrna/filter/cutadapt/GSM1_R1.fastq', ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)) for i, s in enumerate(reads)))
    write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', 'GSM1', 100)
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    routes = (bench / 'routes.tsv').read_text().splitlines()
    assert routes == ['route\tstate\treason', 'unitas\tevaluated\t',
                      'mirge3\tstopped\tcount table missing (mirge3/miR.Counts.csv); run `snakemake quant_all` there',
                      'tinyrna\tnot run\tnot selected at prepare']
    assert 'benchmark.py: mirge3 stopped, scored without it: count table missing' in result.stderr
    report = (bench / 'report.md').read_text()
    assert '## Routes' in report and '| tinyrna | not run | not selected at prepare |' in report
    assert (bench / 'metrics.tsv').read_text().splitlines()[0] == 'metric\tunitas'   # scored without the stopped route
    write(bench / 'not_run.tsv', 'route\treason\nmirge3\tneeds more memory than the run had\n')
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    assert 'mirge3\tnot run\tneeds more memory than the run had' in (bench / 'routes.tsv').read_text() and 'stopped' not in result.stderr
    summary_md = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    row = next(line for line in summary_md.read_text().splitlines() if line.startswith('| [project]'))
    import summary
    keys = next(keys for kind, _, keys, _ in summary.SECTIONS if kind == 'biological')
    assert row.count('not run') == 2 * len(keys) and 'stopped' not in row   # mirge3 and tinyrna in each key metric
    write(bench / 'routes.tsv', 'route\tstate\treason\nunitas\tstopped\tedited\n')   # not the published file: the summary refuses it
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'incomplete or inconsistent evaluation' in result.stderr
    write(bench / 'not_run.tsv', 'route\treason\nmirflowz\tnot declared here\n')
    result = evaluate_cmd(project)
    assert result.returncode != 0 and 'not_run.tsv: mirflowz not declared for this project' in result.stderr


def test_a_bare_null_key_in_a_declaration_is_named(tmp_path):
    """#449: YAML reads `null:` as YAML's null, not as the key `null`; the declaration says to quote it."""
    groups = ('case',) * 4 + ('control',) * 4
    quoted = biological_declaration(tmp_path, {'null': {'splits': 2, 'seed': 1}}, groups)
    assert "'null':" in quoted.read_text()  # as PyYAML writes it, and as a declaration must
    assert benchmark.load_dataset(quoted)['de']['null'] == {'splits': 2, 'seed': 1}
    bare = write(tmp_path / 'bare.yaml', quoted.read_text().replace("'null':", 'null:'))
    with pytest.raises(SystemExit, match=re.escape("`de` has a key that YAML reads as null; write the null comparisons' key quoted, `'null':`")):
        benchmark.load_dataset(bare)


# Knockout cells as a negative control for miRNA counts (#452)

def knockout_common(rows):
    """Common-table rows from {(method, sample): {accession: count}}."""
    return [(m, s, i, i, c) for (m, s), counts in rows.items() for i, c in counts.items()]


def test_knockout_metrics_on_hand_made_counts():
    """Reads per million filtered reads, each knockout sample against the controls of its level; per miRNA a ratio
    for those with at least 10 control reads on average, surviving above 25%; survivors per method, only one method's,
    and every method's."""
    a = {'P1': {'X': 100, 'Y': 40, 'Z': 4, 'W': 20, 'V': 10, 'U': 50}, 'P2': {'X': 100, 'Y': 40, 'Z': 5, 'W': 20, 'V': 10, 'U': 50},
         'K1': {'X': 2, 'Y': 40, 'W': 20, 'U': 100}, 'K2': {'X': 1, 'Y': 30, 'U': 50}, 'P3': {'X': 100}, 'D1': {'X': 3}}
    b = {s: dict(c) for s, c in a.items()}
    b['K1'] = {'X': 2, 'W': 20, 'V': 20, 'U': 100}       # Y gone; V kept
    b['K2'] = {'X': 1, 'V': 5, 'U': 50}
    common = knockout_common({**{('a', s): c for s, c in a.items()}, **{('b', s): c for s, c in b.items()}})
    kept = {'P1': 1e6, 'P2': 1e6, 'K1': 2e6, 'K2': 1e6, 'P3': 1e6, 'D1': 1e6}   # K1 sequenced twice as deep
    units = [('P1', 'Parental', '1'), ('P2', 'Parental', '1'), ('K1', 'DROSHA KO', '1'), ('K2', 'DROSHA KO', '1'),
             ('P3', 'Parental', '2'), ('D1', 'DICER KO', '2')]
    result = biological.knockout(common, ['a', 'b'], kept, units, 'Parental')
    level1 = (224 + 225) / 2                                                       # set 1 controls' miRNA RPM
    assert result['a']['share']['DROSHA KO'] == pytest.approx((162 / 2 / level1 + 81 / level1) / 2)
    assert result['a']['share']['DICER KO'] == pytest.approx(3 / 100)              # against the set 2 control only
    assert result['a']['control_rpm'] == pytest.approx((224 + 225 + 100) / 3)
    ratios = result['a']['ratios']['DROSHA KO']
    assert set(ratios) == {'X', 'Y', 'W', 'V', 'U'}                                # Z: 4.5 control reads on average
    assert ratios['X'] == pytest.approx((100.0, 1.5, 0.01))                         # K1's 2 reads are 1 per million
    assert ratios['Y'][2] == pytest.approx((0.5 + 0.75) / 2) and ratios['W'][2] == pytest.approx(0.25)
    assert set(result['a']['ratios']['DICER KO']) == {'X'}                         # Y and the rest: no set 2 control reads
    per, every = biological.survivors(result, ['a', 'b'], 'DROSHA KO')
    assert per['a']['surviving'] == {'Y', 'U'} and per['b']['surviving'] == {'V', 'U'}   # W at exactly 25% does not survive
    assert every == {'U'} and per['a']['only'] == {'Y'} and per['b']['only'] == {'V'} and per['a']['tested'] == 5
    assert per['a']['only_reads'] == pytest.approx(35) and per['b']['only_reads'] == pytest.approx(12.5)


def knockout_declaration(tmp_path, change=None):
    spec = yaml.safe_load((BENCH / 'datasets/gse77989_hct116_knockouts.yaml').read_text())
    if change:
        change(spec)
    write(tmp_path / 'knockout_expected_hsa.tsv', (BENCH / 'datasets/knockout_expected_hsa.tsv').read_text())
    return write(tmp_path / 'ko.yaml', yaml.safe_dump(spec))


def test_the_knockout_declaration_lists_its_runs_and_pairs_each_knockout_with_its_set(tmp_path):
    spec = benchmark.load_dataset(BENCH / 'datasets/gse77989_hct116_knockouts.yaml')
    runs = [run for s in spec['samples'] for run in s['runs']]
    assert len(runs) == 9 and [t[3].name for t in fill.targets(spec, tmp_path)][4:] == ['{}.fastq.gz'.format(r['run']) for r in runs]  # after 3 other_rna
    assert collections.Counter(s['group'] for s in spec['samples']) == {'Parental': 3, 'DROSHA KO': 2, 'XPO5 KO': 2, 'DICER KO': 2}
    assert {s['group'] for s in spec['samples'] if s['set'] == 1} == {'Parental', 'DROSHA KO'}
    table = [{'sample': s['sample'], 'library': s['sample'], 'group': s['group']} for s in spec['samples']]
    units = benchmark.knockout_units(spec, table)
    assert ('GSM2064272', 'DROSHA KO', '1') in units and ('GSM2064274', 'Parental', '2') in units
    for change, message in ((lambda s: s.update(knockout={'control': 'Parental'}), 'must name two different groups'),
                            (lambda s: s['knockout'].update(survival='Parental'), 'must name two different groups'),
                            (lambda s: s['knockout'].update(control='Wild type'), 'must name two different groups'),
                            (lambda s: s['knockout'].update(normalise='miR-320a'), 'takes only `control`, `survival`, `within` and `expected`'),
                            (lambda s: s['knockout'].update(expected='lists/expected.tsv'), 'knockout.expected names a file next to the declaration'),
                            (lambda s: s['knockout'].update(expected='absent.tsv'), 'absent.tsv not found'),
                            (lambda s: s['samples'][0].pop('set'), "knockout.within is 'set', but not every sample has it"),
                            (lambda s: s['samples'][7].update(set=3), 'no Parental sample with the same set as GSM2064277')):
        with pytest.raises(SystemExit) as error:
            benchmark.load_dataset(knockout_declaration(tmp_path, change))
        assert message in str(error.value), (message, str(error.value))
    assert benchmark.load_dataset(knockout_declaration(tmp_path))['knockout']['expected'] == 'knockout_expected_hsa.tsv'
    pool = yaml.safe_load((BENCH / 'datasets/giraldez2018_equimolar_truseq.yaml').read_text())
    with pytest.raises(SystemExit, match='is for biological datasets'):
        benchmark.load_dataset(write(tmp_path / 'pool.yaml', yaml.safe_dump(dict(pool, knockout={'control': 'A', 'survival': 'B'}))))


def knockout_project(tmp_path, other_rna=None, sequences=None, expected=None, name='project', prepare=()):
    """A prepared knockout project (#452): P1, P2 (Parental) and K1, K2 (DROSHA KO), UNITAS and miRge3 counts; with
    `other_rna`, the declaration's other_rna (#415) and {source: file text} filled in place; `sequences`, the isomiR
    outputs' read sequences; `expected`, the text of a knockout.expected list next to the declaration; `name`, the project
    directory; `prepare`, more prepare options (e.g. --set, #417)."""
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/ko'
    write(data / 'mirbase/23.1/mature.fa', BLOOD)
    reads = ['TAGCAGCACGTAAATATTGGCG'] * 6 + ['TAGCTTATCAGACTGATGTTGA'] * 3
    names, groups = ['P1', 'P2', 'K1', 'K2'], ['Parental', 'Parental', 'DROSHA KO', 'DROSHA KO']
    for i in range(len(names)):
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/SRR{}.fastq.gz'.format(i), 'wt') as fh:
            fh.write(''.join('@SRR{}.{}\n{}TGGAATTCTC\n+\n{}\n'.format(i, j, s, 'I' * (len(s) + 10)) for j, s in enumerate(reads)))
    spec = {'name': 'ko', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
            'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/ko', 'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
            'knockout': {'control': 'Parental', 'survival': 'DROSHA KO'},
            'samples': [{'sample': n, 'group': g, 'runs': [{'run': 'SRR{}'.format(i), 'md5': 'unused', 'url': 'unused'}]}
                        for i, (n, g) in enumerate(zip(names, groups))]}
    if other_rna:
        spec['other_rna'] = {'sources': {key: {'url': 'https://x/{}{}'.format(key, '.gz' if key.startswith('ensembl') else ''), 'sha256': 'unused'}
                                         for key in other_rna}}
        for key, text in other_rna.items():
            path = fill.decoy_path(spec, ext, key, spec['other_rna']['sources'][key]['url'], 'other_rna')
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.name.endswith('.gz'):
                with gzip.open(path, 'wt') as fh:
                    fh.write(text)
            else:
                path.write_text(text)
    if expected:
        write(tmp_path / 'expected.tsv', expected)
        spec['knockout']['expected'] = 'expected.tsv'
    dataset = write(tmp_path / 'ko.yaml', yaml.safe_dump(spec))
    project = tmp_path / name
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(ext), *prepare], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    quant = project / 'data/tmp/smallrna/quant'
    # miR-16-5p vanishes in the knockout; miR-21-5p keeps half its reads in UNITAS only (e.g. reads of another RNA)
    rows = [('miR-16-5p', [600, 600, 3, 3]), ('miR-21-5p', [40, 40, 20, 20]), ('miR-451a', [1, 1, 1, 1])]
    write(quant / 'unitas/mir_counts.tsv', '\t' + '\t'.join(names) + '\n' + ''.join('{}\t{}\n'.format(n, '\t'.join(map(str, v))) for n, v in rows))
    mirge = [('miR-16-5p', [600, 600, 3, 3]), ('miR-21-5p', [40, 40, 2, 2])]
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,' + ','.join(names) + '\n' + ''.join('hsa-{},{}\n'.format(n, ','.join(map(str, v))) for n, v in mirge))
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(ext, 'homo_sapiens', 'Human', {}))
    isomir_outputs(quant, sequences)
    for sample in names:
        write(project / 'data/tmp/smallrna/filter/cutadapt/{}_R1.fastq'.format(sample),
              ''.join('@r{}\n{}\n+\n{}\n'.format(j, s, 'I' * len(s)) for j, s in enumerate(reads)))
        write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', sample, 100)   # 87 kept each
    return project


def test_evaluate_a_knockout_dataset(tmp_path):
    """#452 end to end: the knockout section, its metrics.tsv rows and knockout.tsv, per million filtered reads."""
    project = knockout_project(tmp_path)
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    report = (bench / 'report.md').read_text()
    section = report.split('## Knockouts (#452)')[1].split('\n## ')[0]
    assert '| Parental | 2 | 7367816 | 7356322 |' in section                    # 641 and 640 miRNA reads per 87 kept
    assert '| DROSHA KO | 2 | 0.0374 | 0.0078 |' in section                      # 24/641 and 5/640
    assert '| surviving | 1 | 0 |' in section and '| only | 1 | 0 |' in section and '| only_reads | 20.0 | 0.0 |' in section
    assert '0 survive in every method' in section and '| hsa-miR-21-5p | 1 | **0.50** | 0.05 |' in section
    metrics = {row.split('\t')[0]: row.split('\t')[1:] for row in (bench / 'metrics.tsv').read_text().splitlines()[1:]}
    assert metrics['knockout_mirna_share'] == ['0.0374', '0.0078'] and metrics['knockout_surviving'] == ['1', '0']
    assert '| knockout_mirna_share |' not in report.split('## Knockouts')[0]       # not in the methods table
    tsv = (bench / 'knockout.tsv').read_text().splitlines()
    assert tsv[0] == 'method\tgroup\tmature_id\tmirna\tcontrol_reads\tknockout_reads\tratio\tsurviving'
    assert 'unitas\tDROSHA KO\tMIMAT0000076\thsa-miR-21-5p\t40.0\t20.0\t0.5\t1' in tsv and len(tsv) == 1 + 4   # miR-451a not tested


def test_the_gse_declarations_name_each_run_s_sra_file():
    """Every run of the GSE declarations has its SRA file (`sra`: md5, bytes, from NCBI's locator) for fill --source sra,
    and the runs filled from ENA before have the content MD5 of their reads (#471)."""
    import re
    for name in ('gse317995_nextflex_v4', 'gse77989_hct116_knockouts', 'gse168416_nextflex_v3', 'gse196633_truseq', 'gse271895_nextflex_v4'):
        spec = benchmark.load_dataset(BENCH / 'datasets/{}.yaml'.format(name))
        runs = [run for sample in spec['samples'] for run in sample['runs']]
        assert all(re.fullmatch('[0-9a-f]{32}', run['sra']['md5']) and run['sra']['bytes'] > 0 and run['reads'] > 0 for run in runs), name
        assert len({run['run'] for run in runs}) == len(runs)
    knockouts = benchmark.load_dataset(BENCH / 'datasets/gse77989_hct116_knockouts.yaml')
    assert all('content_md5' in run for sample in knockouts['samples'] for run in sample['runs'])
    tissue = benchmark.load_dataset(BENCH / 'datasets/gse196633_truseq.yaml')
    assert collections.Counter(s['group'] for s in tissue['samples']) == {'tumour': 10, 'normal': 10}
    assert {s['patient'] for s in tissue['samples'] if s['group'] == 'tumour'} == {s['patient'] for s in tissue['samples'] if s['group'] == 'normal'}



# Settings variants (#417)

def test_prepare_set_takes_route_settings_of_the_workflow():
    """`prepare --set`: a setting smallrna.config has for a route the project runs, read as YAML, not its default."""
    variant = benchmark.settings_variant(['quant.tinyrna.mismatch=1', 'quant.tinyrna.normalize_by_genomic_hits=false'], ['tinyrna', 'unitas'])
    assert variant == {'quant.tinyrna.mismatch': 1, 'quant.tinyrna.normalize_by_genomic_hits': False}
    assert benchmark.settings_variant(['quant.tesmall.order=structural_RNA,miRNA,hairpin,exon,TE,intron,piRNA_cluster'], ['tesmall'])
    for assignments, message in ((['quant.tinyrna.mismatch'], 'expected quant.ROUTE.SETTING=VALUE'),
                                 (['db.mirbase.release=22'], 'expected quant.ROUTE.SETTING=VALUE'),
                                 (['quant.tesmall.mismatch=1'], 'tesmall is not one of the project\'s routes'),
                                 (['quant.tinyrna.mismatches=1'], 'smallrna.config has no setting mismatches for tinyrna'),
                                 (['quant.tinyrna.mismatch=0'], 'that is the default'),
                                 (['quant.tinyrna.mismatch=1', 'quant.tinyrna.mismatch=2'], 'set twice')):
        with pytest.raises(SystemExit, match=re.escape(message)):
            benchmark.settings_variant(assignments, ['tinyrna', 'unitas'])
    lines, _ = benchmark.settings_lines({}, {'quant.tinyrna.mismatch': 1})
    assert '**A settings variant** (`prepare --set`): quant.tinyrna.mismatch `1`.' in lines[2]
    spec = {'variant': {'quant.tinyrna.normalize_by_genomic_hits': False}, 'route_settings': {'tinyrna': [
        ('mismatch', '0', 'quant.tinyrna.mismatch'), ('normalize_by_genomic_hits', 'False', 'quant.tinyrna.normalize_by_genomic_hits')]}}
    benchmark.check_variant(spec)                                                   # the record confirms it
    for settings, message in (({'tinyrna': [('normalize_by_genomic_hits', 'True', 'quant.tinyrna.normalize_by_genomic_hits')]},
                               'counted with True, not False'), ({'tinyrna': None}, 'has no settings record'), ({}, 'tinyrna is not scored')):
        with pytest.raises(SystemExit, match=re.escape(message)):
            benchmark.check_variant(dict(spec, route_settings=settings))


def test_a_settings_variant_is_compared_with_the_defaults_in_the_summary(tmp_path):
    """#417 end to end: a project prepared with --set records the variant and sets it in its configuration; the summary
    compares its route with the project of the same dataset run with the defaults, and names the routes it passes."""
    import summary
    default = knockout_project(tmp_path / 'd', name='defaults')
    # any setting smallrna.config has for a route serves here; miRge3's database is its only one
    variant = knockout_project(tmp_path / 'v', name='other-db', prepare=['--set', 'quant.mirge3.db=mirgenedb'])
    assert yaml.safe_load((variant / 'benchmark/dataset.yaml').read_text())['variant'] == {'quant.mirge3.db': 'mirgenedb'}
    assert yaml.safe_load((variant / 'config.yaml').read_text())['quant']['mirge3'] == {'db': 'mirgenedb'}
    quant = variant / 'data/tmp/smallrna/quant'   # with the setting, miRge3 keeps half of miR-21-5p's reads in the knockout
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,P1,P2,K1,K2\nhsa-miR-16-5p,600,600,3,3\nhsa-miR-21-5p,40,40,20,20\n')
    write_records(quant, ['unitas', 'mirge3'], dump=unitas_dump(tmp_path / 'v/ext', 'homo_sapiens', 'Human', {}))
    isomir_outputs(quant)
    record = lambda db: write(quant / 'mirge3/settings.tsv', 'field\tvalue\tsource\nmir_db\t{}\tquant.mirge3.db\ntable_sha256\t{}\tx\n'.format(
        db, adapters.table_sha256('mirge3', quant)))
    record('miRBase')
    result = evaluate_cmd(variant)
    assert result.returncode != 0 and 'settings variant quant.mirge3.db: mirge3 counted with miRBase, not mirgenedb' in result.stderr
    record('MirGeneDB')
    for project in (default, variant):
        result = evaluate_cmd(project)
        assert result.returncode == 0, result.stderr
    assert '**A settings variant** (`prepare --set`): quant.mirge3.db `mirgenedb`.' in (variant / 'benchmark/report.md').read_text()
    for name in ('summary.md', 'summary.html'):
        result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(default), str(variant),
                                 '--output', str(tmp_path / name)], capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stderr
    text = (tmp_path / 'summary.md').read_text()
    assert '### ko: settings against defaults' in text
    assert '| knockout_surviving | mirge3 | quant.mirge3.db mirgenedb (other-db) | 0 | 1 | unitas 1 |' in text   # 0 < 1, then 1 = 1
    assert '<h3>ko: settings against defaults</h3>' in (tmp_path / 'summary.html').read_text()
    alone = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(variant), '--output', str(tmp_path / 'alone.md')],
                           capture_output=True, text=True, timeout=60)
    assert alone.returncode == 0 and 'settings variant other-db (ko) has no project of the same dataset' in alone.stderr
    rows = [{'project': p, 'dataset': 'ko', 'variant': {}, 'identity': 'x', 'values': {}} for p in ('a', 'b')]
    with pytest.raises(ValueError, match='projects a, b run ko with the defaults'):
        summary.settings_variants(rows + [{'project': 'c', 'dataset': 'ko', 'variant': {'quant.mirge3.db': 'mirgenedb'}, 'identity': 'x',
                                           'values': {}}], [])


# Sequencing depth (#409)

def test_depth_series_thins_every_method_alike_and_finds_saturation():
    """Binomial thinning to depths every sample reaches, then the full depth; detection rises with depth, and the depth
    reaching 90% of full-depth detection is interpolated on log depth."""
    import depth
    import numpy
    assert depth.parse_depths('2M,500k,100') == [100, 500000, 2000000]
    for bad in ('0', '1.5', 'x', '1M,1M'):
        with pytest.raises(ValueError, match='depths: '):
            depth.parse_depths(bad)
    thinned = depth.thin([('a', 'S1', 'n', 'M1', 1000.0), ('a', 'S2', 'n', 'M1', 2.4)], {'S1': 0.5, 'S2': 1.0}, numpy.random.default_rng(3))
    assert ('a', 'S2', 'n', 'M1', 2) in thinned and 420 < [r for r in thinned if r[1] == 'S1'][0][-1] < 580   # rounded; about half
    samples = ['A1', 'A2', 'B1', 'B2']
    common = [(m, s, 'n{}'.format(i), 'M{}'.format(i), 10.0 * 2 ** i) for m in ('a', 'b') for s in samples for i in range(12)]
    points, skipped = depth.series(common, ['a', 'b'], samples, dict.fromkeys(samples, 1_000_000), {'A': ['A1', 'A2'], 'B': ['B1', 'B2']},
                                   [1000, 10_000, 100_000, 2_000_000], seed=1)
    assert skipped == [2_000_000] and [d for d, _, _ in points] == [1000, 10_000, 100_000, depth.FULL]
    detected = [values['a']['detected_10'] for _, values, _ in points]
    assert detected == sorted(detected) and detected[-1] == 12 and detected[0] < 6
    assert points[-1][1]['a']['replicate_cv'] == 0.0 and points[-1][2][('a', 'b')]['spearman_either'] == 1.0   # identical tables
    reached = depth.saturation(points, ['a', 'b'], 1_000_000)
    assert 100_000 < reached['a'][0.9] < 1_000_000 and reached['a'][0.9] <= reached['a'][0.95]
    rising = depth.rise(points, ['a', 'b'], 1_000_000)   # the last step, 100k to 1M: (12 - detected at 100k) / 12 per doubling
    assert rising['a'] == round((12 - points[-2][1]['a']['detected_10']) / math.log2(10) / 12, 4) > 0
    fixed = depth.cv_sets(common, ['a'], {'A': ['A1', 'A2'], '': ['B1', 'B2']})
    assert len(fixed[('a', 'A')]) == 12 and ('a', '') not in fixed and depth.rise(points[-1:], ['a'], 1_000_000) == {'a': ''}
    # a depth every sample has is the full depth itself; with unequal depths the full depth is no common measurement
    assert depth.series(common, ['a'], samples, dict.fromkeys(samples, 1000), {}, [1000], seed=1)[1] == [1000]
    unequal = dict(zip(samples, (1_000_000, 1_000_000, 4_000_000, 4_000_000)))
    points, _ = depth.series(common, ['a', 'b'], samples, unequal, {}, [1000, 10_000], seed=1)
    assert points[1][1]['a']['detected_10'] < 0.9 * points[-1][1]['a']['detected_10']
    assert depth.saturation(points, ['a'], 2_500_000, equal=False)['a'][0.9] == depth.BEYOND
    assert depth.saturation_cell(depth.BEYOND, 10_000) == 'beyond 10k'
    rows, pairs = depth.tables(points, ['a', 'b'])
    assert ['a', 'full', 'detected_10', 12] in rows and ['a', 'b', '1k', 'spearman_either', points[0][2][('a', 'b')]['spearman_either']] in pairs


def test_depth_plot(tmp_path):
    """--plot's figure, where matplotlib is installed."""
    pytest.importorskip('matplotlib')
    import depth
    samples = ['A1', 'A2']
    common = [('a', s, 'n{}'.format(i), 'M{}'.format(i), 10.0 * 2 ** i) for s in samples for i in range(8)]
    points, _ = depth.series(common, ['a'], samples, dict.fromkeys(samples, 100_000), {'A': samples}, [1000, 10_000], seed=1)
    depth.plot(points, ['a'], 100_000, tmp_path / 'depth.png')
    assert (tmp_path / 'depth.png').read_bytes()[:4] == b'\x89PNG'


def test_depth_of_an_evaluated_knockout_project(tmp_path):
    """benchmark.py depth end to end: it needs an evaluation; writes its tables, report and plot to --output and nothing
    into the project; depths above the shallowest sample are left out."""
    project = knockout_project(tmp_path)
    command = [sys.executable, str(BENCH / 'benchmark.py'), 'depth', '--project', str(project), '--depths', '20,50,1k', '--output', str(tmp_path / 'depth')]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'run `benchmark.py evaluate' in result.stderr
    assert evaluate_cmd(project).returncode == 0
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    out = tmp_path / 'depth'
    assert (out / 'depth.tsv').read_text().startswith('method\tdepth\tmetric\tvalue\n') and not (out / 'depth.png').exists()
    report = (out / 'depth.md').read_text()
    assert 'Left out, above the shallowest sample: 1k.' in report and '| full |' in report and '**sequences**' in report
    assert 'under independent thinning' in report and '| method | 90% | 95% | rise per doubling |' in report
    write_cutadapt_reports(project / 'data/tmp/smallrna/filter/cutadapt_stages', 'P1', 200)   # the filtering changed since
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'the cutadapt reports give 500 reads, the evaluation recorded 400' in result.stderr
    assert 'unitas\tmirge3\t20\tspearman_either' in (out / 'depth_agreement.tsv').read_text()
    assert not (project / 'benchmark/depth').exists()



# miRNAs inside other RNA classes, on real data (#415)

MIR21 = 'TAGCTTATCAGACTGATGTTGA'
OTHER_RNA = {
    'gtrnadb': '>Homo_sapiens_tRNA-Leu-TAA-1-1 (tRNAscan-SE ID: chr6.trna1) Leu (TAA) 76 bp\n{}\n'.format(
        ('GCCGGATTAGCTCAGTTGGTAGAGC' + MIR21 + 'TTAAGCCTGTGGATCC').replace('T', 'U')),
    'ensembl_ncrna': '>ENST1.1 ncrna transcript_biotype:misc_RNA gene_symbol:RNY1\nGGCTGGTCCGAAGGTAGTGAGTTATCTCAATTGATTGTTCACAGTCAGTTACAGATCG\n'}


def test_other_rna_flags_the_mirnas_inside_other_classes(tmp_path):
    """#415 end to end: miR-21-5p lies inside a tRNA, so it is flagged; its reads, which UNITAS keeps in the knockout,
    lie inside the tRNA too."""
    sequences = {'miR-16-5p': 'TAGCAGCACGTAAATATTGGCG', 'miR-21-5p': MIR21}
    project = knockout_project(tmp_path, OTHER_RNA, sequences)
    result = evaluate_cmd(project, '--bootstrap', '20')
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    assert (bench / 'other_rna.tsv').read_text() == 'mature_id\tmirna\tclass\tmolecules\nMIMAT0000076\thsa-miR-21-5p\ttRNA\ttRNA-Leu-TAA-1-1\n'
    metrics = {row.split('\t')[0]: row.split('\t')[1:] for row in (bench / 'metrics.tsv').read_text().splitlines()[1:]}
    # per sample, UNITAS: P 40 of 641 miRNA reads on miR-21-5p, K 20 of 24; detected P {16, 21}, K {21}
    assert metrics['flagged_reads'][0] == str(round(statistics.median([40 / 641, 40 / 641, 20 / 24, 20 / 24]), 4))
    assert metrics['flagged_detected'] == ['1.0', '0.5'] and metrics['detected_unflagged'] == ['0.5', '0.5']  # miRge3: 2 reads in K
    assert metrics['read_in_other_class'] == [str(round(120 / 1330, 4)), str(round(84 / 1290, 4))]   # miR-21-5p's reads of all
    assert metrics['knockout_surviving_flagged'] == ['1', '0']
    report = (bench / 'report.md').read_text()
    assert '## miRNAs inside other RNA classes (#415)' in report and '- tRNA: hsa-miR-21-5p (tRNA-Leu-TAA-1-1).' in report
    assert 'Surviving miRNAs that lie inside molecules of other RNA classes (#415), per method: unitas: hsa-miR-21-5p.' in report
    ci = {row.split('\t')[0] for row in (bench / 'metrics_ci.tsv').read_text().splitlines()[1:]}  # metric, row, value, low, high
    assert 'detected_unflagged' in ci and 'read_in_other_class' not in ci
    tsv = (bench / 'knockout.tsv').read_text().splitlines()
    assert tsv[0].endswith('\tsurviving\tother_class') and 'unitas\tDROSHA KO\tMIMAT0000076\thsa-miR-21-5p\t40.0\t20.0\t0.5\t1\ttRNA:tRNA-Leu-TAA-1-1' in tsv
    out = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(out)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'detected_unflagged' in out.read_text() and 'read_in_other_class' in out.read_text()


def test_other_rna_declaration_and_reads_inside_molecules(tmp_path):
    molecules = {'tRNA-gtrnadb': [('tRNA-1', 'AAAA' + MIR21 + 'CCCC', None)], 'mRNA-ensembl': [('GENE', 'GGGG' + MIR21, None)]}
    assert decoys.contained(molecules, {MIR21, MIR21[2:], 'ACGTACGTACGTACGTAC'}) == {MIR21: ['tRNA'], MIR21[2:]: ['tRNA']}  # not mRNA
    assert decoys.embedded(molecules, {'MIMAT0000076': MIR21, 'MIMAT0000069': 'TAGCAGCACGTAAATATTGGCG'}) == {'tRNA': {'MIMAT0000076': ['tRNA-1']}}
    spec = {'name': 'b', 'kind': 'biological', 'organism': 'homo_sapiens', 'libprepkit': 'k', 'methods': ['unitas'], 'data': 'd',
            'mirbase': {'release': '23.1', 'url': 'u', 'sha256': 's'},
            'samples': [{'sample': 'S1', 'group': 'g', 'runs': [{'run': 'R1', 'url': 'u', 'md5': 'm'}]}]}
    for other, message in (({'sources': {'ensembl_cdna': {'url': 'u', 'sha256': 's'}}}, 'other_rna takes `sources`'),
                           ({'sources': {'gtrnadb': {'url': 'u'}}}, 'other_rna takes `sources`'), ({'files': {}}, 'other_rna takes `sources`')):
        with pytest.raises(SystemExit, match=re.escape(message)):
            benchmark.load_dataset(write(tmp_path / 'b.yaml', yaml.safe_dump(dict(spec, other_rna=other))))
    spec['other_rna'] = {'sources': {'gtrnadb': {'url': 'https://x/t.fa', 'sha256': 's'}, 'ensembl_ncrna': {'url': 'https://x/n.fa.gz', 'sha256': 't'}}}
    loaded = benchmark.load_dataset(write(tmp_path / 'b.yaml', yaml.safe_dump(spec)))
    targets = {str(t[3]) for t in fill.targets(loaded, '/ext')}
    assert {'/ext/d/other_rna/gtrnadb.fa', '/ext/d/other_rna/ensembl_ncrna.fa.gz'} <= targets


# The miRNAs expected to survive a DROSHA knockout (#478)

import knockout_expected  # noqa: E402

EXPECTED_HEADER = '# a test list\nmature_id\tmirna\treason\tsource\n'


def test_the_expected_list_from_the_study_and_mirtrondb(tmp_path):
    """knockout_expected.py: the study's six by accession, and every human mature exactly inside a mirtronDB precursor
    with the PubMed IDs behind it; one row per mature (miR-877-5p from both); mouse and outside matures left out."""
    rng = random.Random(7)
    seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    precursor = seq(60)
    inside, outside = precursor[5:27], seq(22)
    paper = {a: seq(22) for a, _, _ in knockout_expected.PAPER}
    paper['MIMAT0004949'] = precursor[30:52]                                     # miR-877-5p, a mirtron in both sources
    mature = ''.join('>hsa-miR-{} {} Homo sapiens\n{}\n'.format(i, a, s.replace('T', 'U')) for i, (a, s) in enumerate(paper.items()))
    mature += ('>hsa-miR-9001-5p MIMAT9000001 x\n{}\n>hsa-miR-9002 MIMAT9000002 x\n{}\n>mmu-miR-9003 MIMAT9000003 x\n{}\n'.format(
        inside.replace('T', 'U'), outside, inside))
    columns = 'specie\tid\tname\tpaper\ttype\thairpin arm\thost gene\tchromosome\tstart\tend\tstrand\tsequence\n'
    row = 'H. sapiens \t{} \thsa-mirtron-{} \t{} \tprecursor \t- \tGENE(1) \t1 \t1 \t60 \t+ \t{}\n'
    export = write(tmp_path / 'precursors.tsv', '##mirtronDB tabular format \n' + columns + row.format(1, 1, 26325366, precursor.replace('T', 'U')) +
                   row.format(2, 2, 22955976, 'GG' + inside + 'CC'))
    fa = write(tmp_path / 'mature.fa', mature)
    knockout_expected.main(['--precursors', str(export), '--exported', '2026-11-21', '--mature', str(fa), '--output', str(tmp_path / 'list.tsv')])
    listed = knockout_expected.read_list(tmp_path / 'list.tsv')
    assert set(listed) == set(paper) | {'MIMAT9000001'}
    assert listed['MIMAT9000001'] == ('hsa-miR-9001-5p', 'mirtron', 'mirtronDB v2 PMID 22955976,26325366')
    assert listed['MIMAT0004949'][1:] == ('mirtron', 'Kim et al. 2016 Fig. 3B; mirtronDB v2 PMID 26325366')
    assert listed['MIMAT0000510'][1:] == ("5' capped", 'Kim et al. 2016 Fig. 3B')
    text = (tmp_path / 'list.tsv').read_text()
    assert hashlib.sha256(export.read_bytes()).hexdigest() in text and '2 human precursors, exported 2026-11-21' in text
    for broken, message in ((mature.replace('MIMAT0000510', 'MIMAT0000000'), 'MIMAT0000510 of Kim et al. 2016 not in the mature'),
                            (mature + '>hsa-miR-1254 MIMAT0005905 x\nAAAA\n', "add Kim et al. 2016's miR-1254-5p"),
                            (mature.replace(paper['MIMAT0000510'].replace('T', 'U'), precursor[1:23].replace('T', 'U')), "is 5' capped in Kim")):
        write(fa, broken)
        with pytest.raises(SystemExit, match=re.escape(message)):
            knockout_expected.main(['--precursors', str(export), '--exported', 'd', '--mature', str(fa), '--output', str(tmp_path / 'x.tsv')])
    write(export, columns + row.format(1, 1, 1, precursor).replace('H. sapiens', 'M. musculus'))
    with pytest.raises(SystemExit, match='the human precursors only'):
        knockout_expected.precursors(export)
    for bad in ('mature_id\tmirna\treason\n', EXPECTED_HEADER + 'MIMAT1\tm\t\tx\n', EXPECTED_HEADER + 'MIMAT1\tm\tmirtron\tx\n' * 2):
        with pytest.raises(SystemExit, match='knockout_expected.py: '):
            knockout_expected.read_list(write(tmp_path / 'bad.tsv', bad))


def test_the_knockout_declaration_s_expected_list_is_made_from_its_mirbase_release():
    """GSE77989's list: the study's six and mirtronDB's mirtrons, from the mature.fa the declaration names."""
    spec = benchmark.load_dataset(BENCH / 'datasets/gse77989_hct116_knockouts.yaml')
    path = BENCH / 'datasets' / spec['knockout']['expected']
    listed = knockout_expected.read_list(path)
    assert 'miRBase mature.fa sha256 {};'.format(spec['mirbase']['sha256']) in path.read_text()
    assert {a for a, (_, _, source) in listed.items() if 'Kim et al. 2016' in source} == {a for a, _, _ in knockout_expected.PAPER}
    assert collections.Counter(reason for _, reason, _ in listed.values()) == {'mirtron': 469, "5' capped": 5}
    assert listed['MIMAT0002174'][:2] == ('hsa-miR-484-5p', "5' capped") and listed['MIMAT0027552'][:2] == ('hsa-miR-6826-5p', 'mirtron')


def test_knockout_survivors_split_by_the_expected_list():
    """Of a method's only survivors: expected; not expected but tested by another method (which does not list it);
    not tested elsewhere. Expected miRNAs tested and surviving per method."""
    result = {'a': {'ratios': {'KO': {'E': (10.0, 8.0, 0.8), 'T': (20.0, 10.0, 0.5), 'N': (30.0, 12.0, 0.4), 'M': (50.0, 1.0, 0.02)}}},
              'b': {'ratios': {'KO': {'T': (20.0, 1.0, 0.05), 'M': (50.0, 30.0, 0.6)}}}}
    per, _ = biological.survivors(result, ['a', 'b'], 'KO')
    split = biological.expected_split(result, per, ['a', 'b'], 'KO', {'E', 'M', 'Q'})
    assert split['a']['tested_expected'] == {'E', 'M'} and split['a']['surviving_expected'] == {'E'} and split['a']['only_expected'] == {'E'}
    assert split['a']['only_tested_elsewhere'] == {'T'} and split['a']['only_untested_elsewhere'] == {'N'}
    assert split['a']['only_unexpected_reads'] == pytest.approx(22.0)
    assert split['b']['surviving_expected'] == {'M'} and split['b']['only_expected'] == {'M'} and split['b']['only_unexpected_reads'] == 0.0


def test_evaluate_a_knockout_dataset_with_an_expected_list(tmp_path):
    """prepare copies the list into the project; the survivors are split by it; with other_rna, the flagged miRNAs are
    expected too."""
    listed = EXPECTED_HEADER + "MIMAT0000076\thsa-miR-21-5p\tmirtron\tx\nMIMAT0001631\thsa-miR-451a\t5' capped\tx\nMIMAT9999999\tother\tmirtron\tx\n"
    project = knockout_project(tmp_path, expected=listed)
    assert (project / 'benchmark/expected.tsv').read_text() == listed
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    section = (bench / 'report.md').read_text().split('## Knockouts (#452)')[1].split('\n## ')[0]
    assert "Expected to survive (`knockout.expected`, #478): 2 miRNAs (1 5' capped, 1 mirtron); 2 of the 3 listed are in this miRBase release." in section
    for row in ('| expected tested | 1 | 1 |', '| expected surviving | 1 | 0 |', '| only: expected | 1 | 0 |', '| only: tested elsewhere | 0 | 0 |',
                '| only_unexpected_reads | 0.0 | 0.0 |', '| hsa-miR-21-5p | 1 | mirtron | **0.50** | 0.05 |'):
        assert row in section, row
    assert 'Expected miRNAs tested but not surviving, per method (most control reads first): mirge3: hsa-miR-21-5p.' in section
    metrics = {row.split('\t')[0]: row.split('\t')[1:] for row in (bench / 'metrics.tsv').read_text().splitlines()[1:]}
    assert metrics['knockout_surviving_expected'] == ['1', '0'] and metrics['knockout_surviving_unexpected'] == ['0', '0']
    tsv = (bench / 'knockout.tsv').read_text().splitlines()
    assert tsv[0].endswith('\tsurviving\texpected') and 'unitas\tDROSHA KO\tMIMAT0000076\thsa-miR-21-5p\t40.0\t20.0\t0.5\t1\tmirtron' in tsv
    assert 'unitas\tDROSHA KO\tMIMAT0000069\thsa-miR-16-5p\t600.0\t3.0\t0.005\t0\t' in tsv
    out = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(out)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'knockout_surviving_expected' in out.read_text() and 'knockout_surviving_unexpected' in out.read_text()
    flagged = knockout_project(tmp_path / 'flagged', OTHER_RNA, {'miR-16-5p': 'TAGCAGCACGTAAATATTGGCG', 'miR-21-5p': MIR21},
                               EXPECTED_HEADER + "MIMAT0001631\thsa-miR-451a\t5' capped\tx\n")
    result = evaluate_cmd(flagged)
    assert result.returncode == 0, result.stderr
    tsv = (flagged / 'benchmark/knockout.tsv').read_text().splitlines()
    assert 'unitas\tDROSHA KO\tMIMAT0000076\thsa-miR-21-5p\t40.0\t20.0\t0.5\t1\ttRNA:tRNA-Leu-TAA-1-1\tinside another RNA class (#415)' in tsv



def test_each_routes_counting_settings_are_reported_with_its_results(tmp_path):
    """#417: the settings a route recorded with its count table (quant/<route>/settings.tsv) are in the report's Routes
    section and settings.tsv, listed in evaluation.json; a record for another table stops that route; a route without a
    record is "not recorded"."""
    project = knockout_project(tmp_path)
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/settings.tsv', 'field\tvalue\tsource\nrefseq\tnone\tfilter.spikein\noptions\tUNITAS defaults\tfixed\n'
          'table_sha256\t{}\tx\n'.format(adapters.table_sha256('unitas', quant)))
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    bench = project / 'benchmark'
    report = (bench / 'report.md').read_text()
    assert '- unitas: refseq `none`; options (fixed) `UNITAS defaults`' in report and '- mirge3: not recorded (quantified before #417)' in report
    assert (bench / 'settings.tsv').read_text() == ('route\tfield\tvalue\tsource\nunitas\trefseq\tnone\tfilter.spikein\n'
                                                   'unitas\toptions\tUNITAS defaults\tfixed\n')
    assert 'settings.tsv' in json.loads((bench / 'evaluation.json').read_text())
    out = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(out)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr  # read_evaluation accepts the new file
    write(quant / 'mirge3/settings.tsv', 'field\tvalue\tsource\nmir_db\tmiRBase\tquant.mirge3.db\ntable_sha256\t{}\tx\n'.format('0' * 64))
    result = evaluate_cmd(project)
    assert result.returncode == 0, result.stderr
    states = (bench / 'routes.tsv').read_text()
    assert 'mirge3\tstopped\tsettings record was written for another count table; run `snakemake quant_all` there' in states


def test_the_summary_groups_the_seeds_of_one_synthetic_declaration(tmp_path):
    """#455: projects of one synthetic declaration that differ only by the seed of the draw (or by the routes selected)
    get a table over the seeds: mean and range per key metric and route, and whether the routes rank the same in each."""
    import summary
    root = tmp_path / 'projects'
    base = {'libprepkit': 'Kit A', 'name': 'synthetic_set', 'kind': 'synthetic', 'methods': ['unitas', 'tesmall']}
    values = {1: ('0.9', '0.8', '0.5', '0.6'), 2: ('0.95', '0.85', '0.7', '0.6'), 3: ('0.92', '0.7', '0.6', '0.65')}
    projects = [summary_project(root, 'seed{}'.format(n), dict(base, synthetic={'mode': 'reduced', 'seed': n}),
                                [['metric', 'unitas', 'tesmall'], ['recovered', a, b], ['exact', c, d]], truth=[('S01', 60)])
                for n, (a, b, c, d) in values.items()]
    projects.append(summary_project(root, 'other', dict(base, synthetic={'mode': 'reduced', 'seed': 4, 'mirnas': 10}),  # another set
                                    [['metric', 'unitas', 'tesmall'], ['recovered', '0.1', '0.2'], ['exact', '0.1', '0.2']], truth=[('S01', 60)]))
    out = tmp_path / 'summary.md'
    summary.write(projects, out)
    text = out.read_text()
    table = text.split('### synthetic_set over seeds 1, 2, 3')[1]
    assert '| recovered | 0.9233 (0.9–0.95) | 0.7833 (0.7–0.85) | 1 of 1 pairs |' in table  # unitas above tesmall in each
    assert '| exact | 0.6 (0.5–0.7) | 0.6167 (0.6–0.65) | 0 of 1 pairs |' in table          # tesmall first in seed 1 only
    assert 'over seeds 1, 2, 3, 4' not in text and text.count('over seeds') == 1               # the other declaration stays apart
    summary.write(projects, tmp_path / 'summary.html')
    assert 'synthetic_set over seeds 1, 2, 3' in (tmp_path / 'summary.html').read_text()
    twice = summary_project(root, 'seed2_again', dict(base, synthetic={'mode': 'reduced', 'seed': 2}),
                            [['metric', 'unitas', 'tesmall'], ['recovered', '0.9', '0.8']], truth=[('S01', 60)])
    with pytest.raises(ValueError, match='with the same seed 2'):
        summary.write(projects + [twice], tmp_path / 'again.md')


# Scenarios of the genome-wide synthetic set (#453)

def test_a_scenario_extends_its_base_declaration(tmp_path):
    base = benchmark.load_dataset(BENCH / 'datasets/synthetic_hsa_genome.yaml')
    null = benchmark.load_dataset(BENCH / 'datasets/synthetic_hsa_genome_null.yaml')
    unbalanced = benchmark.load_dataset(BENCH / 'datasets/synthetic_hsa_genome_unbalanced.yaml')
    for scenario in (null, unbalanced):
        assert {k for k in scenario if scenario[k] != base[k]} == {'name', 'synthetic'} and 'extends' not in scenario
        assert {k for k in scenario['synthetic'] if scenario['synthetic'][k] != base['synthetic'][k]} == {'structure'}
        assert scenario['synthetic']['structure']['nuisance'] == base['synthetic']['structure']['nuisance']  # merged key by key
    assert null['synthetic']['structure']['group']['fraction'] == 0
    assert all(f > 0 for f in unbalanced['synthetic']['structure']['group']['log2']) and unbalanced['synthetic']['structure']['group']['abundant'] == 5
    assert benchmark.merged({'a': {'b': 1, 'c': [1, 2]}, 'd': 1}, {'a': {'c': [3]}}) == {'a': {'b': 1, 'c': [3]}, 'd': 1}  # lists are replaced
    write(tmp_path / 'base.yaml', yaml.safe_dump(dict(base, data='x')))
    for change, message in (({'extends': 'base.yaml'}, 'needs a name of its own'),
                            ({'extends': 'base.yaml', 'name': base['name']}, 'needs a name of its own'),
                            ({'extends': 'absent.yaml', 'name': 'n'}, 'which is not a declaration beside it'),
                            ({'extends': 'chained.yaml', 'name': 'n'}, 'one level only'),
                            ({'extends': '../datasets/base.yaml', 'name': 'n'}, 'by its file name'),
                            ({'extends': 'base.yaml', 'name': 'n', 'mirbase': {'release': '22'}}, 'keeps its dataset and references; it changes mirbase'),
                            ({'extends': 'base.yaml', 'name': 'n', 'data': 'other'}, 'it changes data'),
                            ({'extends': 'base.yaml', 'name': 'n', 'synthetic': {'decoys': {}}}, 'it changes synthetic.decoys')):
        write(tmp_path / 'chained.yaml', yaml.safe_dump({'extends': 'base.yaml', 'name': 'c'}))
        with pytest.raises(SystemExit, match=message):
            benchmark.load_dataset(write(tmp_path / 'scenario.yaml', yaml.safe_dump(change)))


def test_false_calls_against_the_changes_need_them_mostly_one_way():
    """#453: unchanged features called opposite to the known changes, when at least two thirds go one way."""
    features = [('U1', 0.0), ('U2', 0.0), ('U3', 0.0), ('C1', 1.0), ('C2', 2.0), ('C3', 1.0)]
    fit = {'U1': (-0.4, 0.001, 0.01), 'U2': (-0.3, 0.2, 0.3), 'U3': (0.3, 0.001, 0.01), 'C1': (1.0, 0.001, 0.01), 'C2': (2.0, 0.001, 0.01),
           'C3': (1.0, 0.001, 0.01), 'C4': (-1.0, 0.001, 0.01)}
    results = {('m', test, f): v for test in de.TESTS for f, v in fit.items()}
    values = de.scores(results, features, ['m'])['m']
    assert values['voom_false_calls_against'] == 0.3333 and values['voom_false_calls'] == 0.6667  # U1 against, U3 with
    balanced = [('U1', 0.0), ('U2', 0.0), ('U3', 0.0), ('C1', 1.0), ('C2', 2.0), ('C3', 1.0), ('C4', -1.0)]
    assert de.scores(results, balanced, ['m'])['m']['edger_false_calls_against'] == 0.3333  # 3 of 4 up: one way
    balanced[4] = ('C2', -2.0)
    assert math.isnan(de.scores(results, balanced, ['m'])['m']['edger_false_calls_against'])  # 2 up, 2 down: no direction
    assert de.scores(results, balanced, ['m'], way=1)['m']['edger_false_calls_against'] == 0.3333  # a bootstrap draw: the full design's way
    assert de.direction(features) == 1 and de.direction(balanced) == 0
