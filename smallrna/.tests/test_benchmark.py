"""Small-RNA benchmark: synthetic truth, method adapters, metrics, pools (truth groups, oracle, fill) and the prepare/evaluate/summary steps
on synthetic tables."""
import gzip
import json
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
        assert 'input: {}'.format(ref / 'mirna.gff3') in gff and 'tinyrna/custom-mmu-test/mus_musculus/mirna.v2.gff3' in gff
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
                    fill.download('good', 'md5', '0' * 32, path, log)
                fill.download('good', 'md5', hashlib.md5(good).hexdigest(), path, log)
            return super().read(size)

    def urlopen(url, timeout):
        return Interleaved(good) if url == 'first' else io.BytesIO(good) if url == 'good' else opened(url, timeout=timeout)
    monkeypatch.setattr(fill.urllib.request, 'urlopen', urlopen)
    fill.download('first', 'md5', hashlib.md5(good).hexdigest(), path, log)
    assert path.read_bytes() == good and [p.name for p in path.parent.iterdir()] == ['SRR1.fastq.gz']


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
    with pytest.raises(SystemExit) as stop:
        benchmark.check_methods(project, spec, config)
    err = str(stop.value)
    assert 'mirge3: table samples differ from the benchmark\'s; missing S2' in err
    assert 'tinyrna: table samples differ from the benchmark\'s; extra S3' in err
    assert 'unitas: annotated with mirbase 22, but the dataset declares miRBase 23.1' in err

    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S1,S2\nhsa-miR-1-3p,1,2\n')
    write(quant / 'tinyrna/tinyrna_counts.tsv', 'Feature ID\tClassifier\tFeature Name\tS1\tS2\n')
    write_records(quant, ['unitas'], database='mirgenedb', release='3.0')
    write_records(quant, ['tinyrna'], database='NA', release='NA')  # a db.mirflowz.custom reference
    with pytest.raises(SystemExit) as stop:
        benchmark.check_methods(project, spec, config)
    err = str(stop.value)
    assert 'unitas: annotated with mirgenedb 3.0' in err and 'tinyrna: annotated with NA NA' in err and 'mirge3' not in err

    (quant / 'mirge3/reference.tsv').unlink()
    with pytest.raises(SystemExit, match='mirge3: reference record missing'):
        benchmark.check_methods(project, spec, config)


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
    with pytest.raises(SystemExit, match=message):
        benchmark.check_methods(project, spec, config)


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
    with pytest.raises(SystemExit, match="unitas: table samples differ from the benchmark's; missing S1; extra S1 "):
        benchmark.check_methods(project, spec, config)


def test_synthetic_projects_check_samples_only(tmp_path):
    """Synthetic truth comes from the dataset's own reference, so no release is compared; samples still are."""
    project, spec, config = check_project(tmp_path, 'synthetic', ['S01'], ['tesmall', 'mirflowz'])
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'tesmall/tesmall_counts.tsv', 'ftype\tfid\tS01\n')
    write(quant / 'mirflowz/all_mirna_counts.tab', 'ID\tS01\n')
    benchmark.check_methods(project, spec, config)  # no reference.tsv needed
    write(quant / 'mirflowz/all_mirna_counts.tab', 'ID\tS01\tS01\n')
    with pytest.raises(SystemExit, match="mirflowz: table samples differ from the benchmark's; repeated columns"):
        benchmark.check_methods(project, spec, config)
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
    with pytest.raises(SystemExit, match='tesmall: reference record was written for another count table'):
        benchmark.check_methods(project, spec, config)
    write(quant / 'smrnaseq/reference.tsv', 'field\tvalue\tsource\nmirna_database\tmirbase\tx\nmirna_release\t23.1\tx\n')
    with pytest.raises(SystemExit, match='smrnaseq: .*one non-empty table_sha256 row, found 0'):
        benchmark.check_methods(project, spec, config)


def test_every_method_the_shipped_datasets_declare_is_checked_before_scoring():
    """The declarations' default method lists (e.g. smrnaseq in the synthetic and Giraldez pools) must all be known to
    the preflight, or evaluation of a fully quantified project stops."""
    for path in sorted((BENCH / 'datasets').glob('*.yaml')):
        methods = yaml.safe_load(path.read_text())['methods']
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
import math  # noqa: E402

DE_RUNS = {'A1': ('A', 'L1'), 'A2': ('A', 'L2'), 'B1': ('B', 'L1'), 'B2': ('B', 'L2')}
ONE, EQUAL, TWENTY_ONE = 'TGGAATGTAAAGAAGTATGTAT', 'ACAGTAGTCTGCACATTGGTTA', 'TAGCTTATCAGACTGATGTTGA'


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
    texts = de.input_tables(counts, oracle, entries, ['oracle', 'm'], [('A1', 'A', {'lab': 'L1'}), ('B1', 'B', {'lab': 'L2'})], sizes)
    assert texts['counts.tsv'] == 'row\tpool_id\tA1\tB1\noracle\tS0\t7\t1\noracle\tS1\t2\t2\nm\tS0\t0\t2\nm\tS1\t2\t3\n'  # half to even
    assert texts['library_sizes.tsv'] == 'row\tA1\tB1\noracle\t10\t3\nm\t1234.5\t100\n'
    assert texts['design.tsv'] == 'sample\tgroup\tlab\nA1\tA\tL1\nB1\tB\tL2\n'


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
                'offset': 0.05}       # median of S0 and S1
    assert de.scores(results, exact, ['m']) == {'m': {'{}_{}'.format(t, k): v for t in de.TESTS for k, v in expected.items()}}
    assert sorted(de.metric_keys()) == sorted('{}_{}'.format(t, k) for t in de.TESTS for k in expected)


def test_de_scores_of_a_resample_are_nan_not_errors():
    results = {('m', test, 'S0'): (None, None, None) for test in de.TESTS}  # untested: no reads
    values = de.scores(results, de_entries((3, 2)), ['m'])['m']  # no 1:1 and no 2-fold sequence, no call
    assert len(values) == len(de.metric_keys()) and all(math.isnan(v) for v in values.values())


def test_de_results_must_cover_every_row_test_and_sequence(tmp_path):
    lines = ['\t'.join(de.RESULT_HEADER)] + ['m\t{}\t{}\t1.5\t0.01\t{}'.format(t, p, 'NA' if p == 'S1' else '0.02')
                                             for t in de.TESTS for p in ('S0', 'S1')]
    path = write(tmp_path / 'results.tsv', '\n'.join(lines) + '\n')
    read = lambda pool_ids: de.read_results(path.read_text(), ['m'], pool_ids, path)
    table = read(['S0', 'S1'])
    assert table[('m', 'edger', 'S0')] == (1.5, 0.01, 0.02) and table[('m', 'voom', 'S1')] == (1.5, 0.01, None)
    with pytest.raises(ValueError, match='and 2 unexpected'):  # S1 in both tests
        read(['S0'])
    with pytest.raises(ValueError, match='2 expected lines missing'):
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
    (lambda s: s.update(de={'covariates': ['batch']}), 'distinct names from lab'),
    (lambda s: s.update(de={'covariates': ['lab', 'lab']}), 'distinct names from lab'),
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


def test_shipped_ratiometric_declarations_declare_their_de_design():
    truseq = benchmark.load_dataset(BENCH / 'datasets/giraldez2018_ratiometric_truseq.yaml')
    nextflex = benchmark.load_dataset(BENCH / 'datasets/giraldez2018_ratiometric_nextflex.yaml')
    assert truseq['de'] == {'covariates': ['lab']} and nextflex['de'] == {'covariates': []}
    for path in BENCH.glob('datasets/*.yaml'):
        spec = yaml.safe_load(path.read_text())
        assert ('de' in spec) == (spec.get('pool') == 'ratiometric'), path.name


def de_snakemake(cwd, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(cwd), XDG_CACHE_HOME=str(cwd / '.cache'))
    return subprocess.run([str(Path(sys.executable).with_name('snakemake')), *args], cwd=cwd, env=env, capture_output=True, text=True, timeout=120)


def test_de_snakefile_runs_r_in_the_sartools_image(tmp_path):
    for name in de.INPUTS:
        write(tmp_path / 'benchmark/de/input' / name, 'x\n')
    write(tmp_path / 'Snakefile', 'include: {!r}\nimport json\nfrom pathlib import Path\n'
          'Path("rules.json").write_text(json.dumps({{r.name: r.container_img for r in workflow.rules}}))\n'.format(str(BENCH / 'de.smk')))
    result = de_snakemake(tmp_path, '--list-rules')
    assert result.returncode == 0, result.stdout + result.stderr
    sartools = yaml.safe_load((ROOT / 'docker.config').read_text())['docker']['sartools']
    assert json.loads((tmp_path / 'rules.json').read_text()) == {'benchmark_de': 'docker://' + sartools}
    result = de_snakemake(tmp_path, '-s', str(BENCH / 'de.smk'), '-n', '-p', '--cores', '1')
    assert result.returncode == 0, result.stdout + result.stderr
    assert ('Rscript {} --counts benchmark/de/input/counts.tsv --library-sizes benchmark/de/input/library_sizes.tsv '
            '--design benchmark/de/input/design.tsv --results benchmark/de/results.tsv --record benchmark/de/record.tsv > '
            'benchmark/de/de.log 2>&1'.format(BENCH / 'de.R')) in result.stdout


def test_de_snakefile_passes_the_options_the_r_script_defines():
    import re
    defined = set(re.findall(r"make_option\('(--[a-z-]+)'", (BENCH / 'de.R').read_text()))
    passed = set(re.findall(r"(--[a-z-]+) \{", (BENCH / 'de.smk').read_text()))
    assert defined == passed == {'--counts', '--library-sizes', '--design', '--results', '--record'}


def de_project(tmp_path):
    """A prepared and quantified ratiometric pool (RATIO: miR-1-3p 4:1, miR-199a-3p 1:1, miR-21-5p 1:8) with two runs
    per pool over two labs, declaring `de:`; UNITAS and miRge3 outputs and the filtered reads, as quant_all --notemp
    leaves them."""
    ext = tmp_path / 'ext'
    data = ext / 'benchmark/smallrna/p'
    write(data / 'composition/ratio.csv', RATIO)
    write(data / 'mirbase/23.1/mature.fa', HSA)
    reads = {'A1': [ONE] * 8 + [EQUAL] * 4 + [TWENTY_ONE], 'A2': [ONE] * 9 + [EQUAL] * 5 + [TWENTY_ONE],
             'B1': [ONE] * 2 + [EQUAL] * 4 + [TWENTY_ONE] * 8, 'B2': [ONE] * 2 + [EQUAL] * 5 + [TWENTY_ONE] * 9}
    for run, seqs in reads.items():
        (data / 'fastq').mkdir(parents=True, exist_ok=True)
        with gzip.open(data / 'fastq/{}.fastq.gz'.format(run), 'wt') as fh:
            fh.write(''.join('@r{}\n{}TGGAATTCTC\n+\n{}\n'.format(i, s, 'I' * (len(s) + 10)) for i, s in enumerate(seqs)))
    dataset = write(tmp_path / 'ds.yaml', yaml.safe_dump({
        'name': 'p', 'kind': 'pool', 'pool': 'ratiometric', 'organism': 'homo_sapiens', 'libprepkit': 'TruSeq Small RNA Library Prep Kit',
        'methods': ['unitas', 'mirge3'], 'data': 'benchmark/smallrna/p', 'de': {'covariates': ['lab']},
        'composition': {'file': 'ratio.csv', 'url': 'unused', 'sha256': 'unused'},
        'mirbase': {'release': '23.1', 'url': 'unused', 'sha256': 'unused'},
        'runs': [{'run': r, 'group': g, 'lab': lab, 'md5': 'unused', 'url': 'unused'} for r, (g, lab) in DE_RUNS.items()]}))
    project = tmp_path / 'project'
    dump = unitas_dump(ext, 'homo_sapiens', 'Human', {})
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'prepare', '--dataset', str(dataset), '--project', str(project),
                             '--ext-dir', str(ext)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tA1\tA2\tB1\tB2\nmiR-1-3p\t8\t9\t2\t2\nmiR-199a-3p\t4\t5\t4\t5\nmiR-21-5p\t1\t1\t8\t9\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,A1,A2,B1,B2\nhsa-miR-1-3p,8,9,2,2\nhsa-miR-199a-3p,4,5,4,5\nhsa-miR-21-5p,1,1,8,9\n')
    write_records(quant, ['unitas', 'mirge3'], dump=dump)
    for run, seqs in reads.items():
        write_stages(project, run, seqs, seqs, seqs)
    return project


def de_outputs(project, fits, change=None, drop_last=False):
    """de.smk's outputs as de.R writes them for the project's current inputs: every row and test gets `fits`
    ({pool_id: (log2fc, pvalue, padj)}); `change` edits the record; `drop_last` leaves the table's last line out (a
    table de.R recorded, but incomplete)."""
    bench = project / 'benchmark'
    texts = {name: (bench / 'de/input' / name).read_text() for name in de.INPUTS}
    rows = list(dict.fromkeys(line.split('\t')[0] for line in texts['counts.tsv'].splitlines()[1:]))
    lines = ['\t'.join(de.RESULT_HEADER)] + ['\t'.join([row, test, pool_id] + ['NA' if v is None else str(v) for v in fit])
                                             for row in rows for test in de.TESTS for pool_id, fit in fits.items()]
    results = '\n'.join(lines[:-1] if drop_last else lines) + '\n'
    write(bench / 'de/results.tsv', results)
    record = {'md5_' + Path(name).stem: de.md5(text) for name, text in texts.items()}
    record.update(md5_script=de.md5(de.SCRIPT.read_bytes()), md5_results=de.md5(results), R='R version 4.5.3', edgeR='4.10.7',
                  limma='3.68.5', formula='~lab + group', edger='exact test', voom='voom')
    if change:
        change(record)
    write(bench / 'de/record.tsv', 'field\tvalue\n' + ''.join('{}\t{}\n'.format(k, v) for k, v in record.items()))


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
    assert table['edger_called'] == table['voom_called'] == ['1.0'] * 3
    assert table['edger_false_calls'] == ['0.0'] * 3 and table['edger_offset'] == ['0.1'] * 3
    assert table['edger_fc_bias'] == ['-0.35'] * 3  # median of 1.8 - 2 and -(-2.5 + 3)
    ci = (bench / 'metrics_ci.tsv').read_text().splitlines()
    assert [line for line in ci if not line.startswith(('edger_', 'voom_'))] == plain_ci  # the existing intervals are unchanged
    assert json.loads((bench / 'evaluation.json').read_text())['bootstrap']['skipped'] == 0  # the DE scores use every draw
    assert any(line.startswith('edger_called\tunitas\t1.0\t') for line in ci)
    report = (bench / 'report.md').read_text()
    assert '## Differential expression' in report and '| 4× | 1 | 1.00 | 1.00 | 1.00 |' in report and '| 1× | 1 | 0.00 |' in report
    assert 'design `~lab + group`' in report and hashlib.sha256((bench / 'de/results.tsv').read_bytes()).hexdigest() in report
    summary_md = tmp_path / 'summary.md'
    result = subprocess.run([sys.executable, str(BENCH / 'benchmark.py'), 'summary', '--projects', str(project), '--output', str(summary_md)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert 'edger_called: unitas' in summary_md.read_text()


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
    ('equimolar', 'is for ratiometric pools, not this equimolar dataset'),
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
    {'equimolar': equimolar, 'undeclared': no_de}.get(case, lambda p: None)(project)
    names = ['metrics.tsv', 'report.md', 'evaluation.json', 'common.tsv', 'per_sequence.tsv', 'filter_stages.tsv', 'unitas_names.tsv'] + \
        ['de/input/' + name for name in de.INPUTS]
    before = {name: (bench / name).read_text() for name in names}
    result = evaluate_cmd(project, '--de')
    assert result.returncode != 0 and message in result.stderr, result.stderr
    assert {name: (bench / name).read_text() for name in names} == before  # nothing published, no supporting table rewritten
