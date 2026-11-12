"""Small-RNA benchmark: synthetic truth, method adapters, metrics and the prepare/evaluate steps on synthetic tables."""
import gzip
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / 'smallrna/benchmark'
sys.path.insert(0, str(BENCH))
import adapters  # noqa: E402
import metrics  # noqa: E402
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


def native_outputs(quant):
    write(quant / 'unitas/mir_counts.tsv', '\tS01\tS02\nmiR-192-5p\t45\t45\nmiR-194-5p \t20\t25\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,S01,S02\nmmu-miR-192-5p,45,45\nmmu-miR-194-5p/mmu-miR-204-5p,25,25\nmmu-miR-9999-5p,1,0\n')
    write(quant / 'mirflowz/all_mirna_counts.tab',
          'ID\tS01\tS02\nmmu-miR-192-5p\t40\t40\nmmu-miR-192-5p|0|1|23M|23|CTGACCTATGAATTGACAGCCAA\t5\t5\nmmu-miR-194-2-5p\t25\t25\n')
    return write(quant.parent / 'mirflowz.gff3',
                 '19\t.\tmiRNA\t1\t30\t.\t+\t.\tID=MIMAT0000224_1;Alias=MIMAT0000224;Name=mmu-miR-194-2-5p;Derives_from=MI0000235\n')


def test_adapters_match_by_accession_whatever_the_method_names(tmp_path):
    gff = native_outputs(tmp_path / 'quant')
    rows = adapters.common_table(['unitas', 'mirge3', 'mirflowz'], tmp_path / 'quant', NAME_TO_ID, 'mmu-', gff)
    got = {(m, s, n): (i, c) for m, s, n, i, c in rows}
    assert got[('unitas', 'S01', 'miR-194-5p')] == ('MIMAT0000224', 20)
    assert got[('mirge3', 'S01', 'mmu-miR-194-5p/mmu-miR-204-5p')] == ('MIMAT0000224', 25)
    assert got[('mirge3', 'S01', 'mmu-miR-9999-5p')] == ('', 1)
    assert got[('mirflowz', 'S01', 'mmu-miR-192-5p')] == ('MIMAT0000517', 45)  # canonical and isomiR rows together
    assert got[('mirflowz', 'S01', 'mmu-miR-194-2-5p')] == ('MIMAT0000224', 25)
    assert ('mirge3', 'S02', 'mmu-miR-9999-5p') not in got  # zero counts are dropped
    with pytest.raises(ValueError, match='no benchmark adapter'):
        adapters.common_table(['newtool'], tmp_path / 'quant', NAME_TO_ID, 'mmu-')


def test_count_metrics():
    truth = {('S01', 'A'): 10, ('S01', 'B'): 5}
    common = [('m', 'S01', 'a', 'A', 10), ('m', 'S01', 'b', 'B', 3), ('m', 'S01', 'c', 'C', 2), ('m', 'S01', 'x', '', 4)]
    result = metrics.count_metrics(common, truth)['m']
    assert result == {'recovered': round(13 / 15, 4), 'exact': 0.5, 'misassigned': 2, 'unplaced': 4}


def test_equimolar_metrics():
    common = [('m', 'S01', 'a', 'A', 10), ('m', 'S01', 'b', 'B', 10), ('m', 'S01', 'c', 'C', 40), ('m', 'S01', 'z', 'Z', 99)]
    result = metrics.equimolar_metrics(common, {'A', 'B', 'C', 'D'})['m']
    assert result['detected'] == 0.75
    assert result['within_2fold'] == 0.5  # median 10: A and B within two-fold, C (40) and D (0) not
    assert result['cv'] > 0


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
    result = subprocess.run(cmd + ['evaluate', '--project', str(project), '--mirflowz-gff', str(gff)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    report = (project / 'benchmark/report.md').read_text()
    assert '| recovered | 0.9643 | 1.0 | 1.0 |' in report  # unitas misses 5 of miR-194-5p's 25 reads in S01: 135/140
    assert '| mmu-miR-192-5p | 45 | 45 | 45 | 45 |' in report
    assert 'mirge3: `mmu-miR-9999-5p`' in report
    assert (project / 'benchmark/common.tsv').read_text().startswith('method\tsample\tmature_name\tmature_id\tcount\n')


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
