"""Small-RNA benchmark: synthetic truth, method adapters, metrics, pools (truth groups, oracle, fill) and the prepare/evaluate steps
on synthetic tables."""
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
import fill  # noqa: E402
import metrics  # noqa: E402
import pool  # noqa: E402
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
          'ID\tS01\tS02\nmmu-miR-192-5p\t40\t40\nmmu-miR-192-5p|0|1|23M|23|CTGACCTATGAATTGACAGCCAA\t5\t5\nmmu-miR-194-2-5p\t25\t25\n'
          'mmu-miR-204-5p|0|0|22M|22|TTCCCTTTGTCATCCTATGCCT\tNA\t3\n')  # NA: no reads in that sample
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
    assert got[('mirflowz', 'S02', 'mmu-miR-204-5p')] == ('MIMAT0000237', 3) and ('mirflowz', 'S01', 'mmu-miR-204-5p') not in got
    assert ('mirge3', 'S02', 'mmu-miR-9999-5p') not in got  # zero counts are dropped
    with pytest.raises(ValueError, match='no benchmark adapter'):
        adapters.common_table(['newtool'], tmp_path / 'quant', NAME_TO_ID, 'mmu-')


def test_count_metrics():
    truth = {('S01', 'A'): 10, ('S01', 'B'): 5}
    common = [('m', 'S01', 'a', 'A', 10), ('m', 'S01', 'b', 'B', 3), ('m', 'S01', 'c', 'C', 2), ('m', 'S01', 'x', '', 4)]
    result = metrics.count_metrics(common, truth)['m']
    assert result == {'recovered': round(13 / 15, 4), 'exact': 0.5, 'misassigned': 2, 'unplaced': 4}


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


# Pools: a human miRBase subset and a composition table in the study's two formats.
HSA = ('>hsa-miR-1-3p MIMAT0000416 Homo sapiens miR-1-3p\nUGGAAUGUAAAGAAGUAUGUAU\n'
       '>hsa-miR-199a-3p MIMAT0000232 Homo sapiens miR-199a-3p\nACAGUAGUCUGCACAUUGGUUA\n'
       '>hsa-miR-199b-3p MIMAT0004563 Homo sapiens miR-199b-3p\nACAGUAGUCUGCACAUUGGUUA\n'
       '>hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\nUAGCUUAUCAGACUGAUGUUGA\n'
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
              ('off', 'S1', 'miR-155-5p', 'MIMAT0000646', 10)]  # not in the pool: from a foreign sequence
    counts, outside = pool.method_counts(common, entries)
    assert counts[('good', 'S1', 'MIR-199A-3P')] == 10 and counts[('off', 'S1', 'MIR-21')] == 50
    m = pool.equimolar_metrics(counts, outside, oracle, entries, ['good', 'off'], ['S1'])
    assert m['good']['median_log2_vs_oracle'] == 0 and m['good']['within_2fold'] == 1 and m['good']['outside'] == 0
    assert m['off']['beyond_10fold'] == 1 and m['off']['detected'] == 0.5
    assert m['off']['outside'] == round(10 / 2010, 4)
    assert m['oracle']['detected'] == 1 and 'within_2fold' not in m['oracle']


def test_ratiometric_metrics_against_the_design(tmp_path):
    entries = pool_entries(tmp_path, RATIO)
    design = {'A': (4000, 1000, 1000), 'B': (1000, 1000, 8000)}
    groups = {'A1': 'A', 'A2': 'A', 'B1': 'B'}
    seqs = [e['sequence'] for e in entries]
    oracle = {s: dict(zip(seqs, design[g])) for s, g in groups.items()}
    ids = [e['pool_id'] for e in entries]
    counts = {('flat', s, i): 100 for s in groups for i in ids}
    counts.update({('exact', s, i): v for s, g in groups.items() for i, v in zip(ids, design[g])})
    m, observed = pool.ratiometric_metrics(counts, oracle, entries, ['exact', 'flat'], groups, pseudo=0)
    assert m['oracle'] == m['exact'] and m['exact']['median_abs_error'] == 0 and m['exact']['called'] == 1
    assert m['exact']['false_calls'] == 0 and m['flat']['called'] == 0
    assert observed['exact']['hsa-miR-21-5p'] < -2


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


def test_filtered_fastq_follows_the_one_filter_route(tmp_path):
    """The oracle reads the route's filtered reads: cutadapt by default, the spike-in prefilter's with spike-ins."""
    import benchmark
    filtered = tmp_path / 'data/tmp/smallrna/filter'
    with pytest.raises(SystemExit, match='quant_all --notemp'):
        benchmark.filtered_fastq(tmp_path, 'S1')
    write(filtered / 'cutadapt/S1_R1.fastq', '')
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'cutadapt/S1_R1.fastq'
    write(filtered / 'spikein/bowtie/S1_R1.fastq', '')
    with pytest.raises(SystemExit, match='more than one filter route'):
        benchmark.filtered_fastq(tmp_path, 'S1')
    (filtered / 'cutadapt/S1_R1.fastq').unlink()
    assert benchmark.filtered_fastq(tmp_path, 'S1') == filtered / 'spikein/bowtie/S1_R1.fastq'


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
    cmd = [sys.executable, str(BENCH / 'benchmark.py')]
    result = subprocess.run(cmd + ['prepare', '--dataset', str(dataset), '--project', str(project), '--ext-dir', str(ext),
                                   '--methods', 'unitas,mirge3', '--reads', '6'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load((project / 'config.yaml').read_text())
    assert config['quant']['method'] == 'unitas,mirge3' and config['read_geometry'] == [32] and 'db' not in config
    with gzip.open(project / 'data/raw/fastq/p/SRRA_R1.fastq.gz', 'rt') as fh:
        assert len(fh.read().splitlines()) == 24  # first 6 reads
    assert '3 exact' in result.stdout

    quant = project / 'data/tmp/smallrna/quant'
    write(quant / 'unitas/mir_counts.tsv', '\tSRRA\tSRRB\nmiR-1-3p\t8\t2\nmiR-21-5p\t2\t8\n')
    write(quant / 'mirge3/miR.Counts.csv', 'miRNA,SRRA,SRRB\nhsa-miR-1-3p,8,2\nhsa-miR-21-5p,2,8\n')
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode != 0 and 'quant_all --notemp' in result.stderr  # the oracle needs the filtered reads
    for run, seqs in reads.items():
        write(project / 'data/tmp/smallrna/filter/cutadapt/{}_R1.fastq'.format(run),
              ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)) for i, s in enumerate(seqs)))
    result = subprocess.run(cmd + ['evaluate', '--project', str(project)], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    report = (project / 'benchmark/report.md').read_text()
    assert '| metric | oracle | unitas | mirge3 |' in report and '| called | 1.0 | 1.0 | 1.0 |' in report
    assert '| 1:8 | -3.00 | 1 |' in report
    per_sequence = (project / 'benchmark/per_sequence.tsv').read_text().splitlines()
    assert per_sequence[0] == 'sample\tpool_id\tgroup\toracle\tunitas\tmirge3' and 'SRRA\thsa-miR-1-3p\texact\t8\t8\t8' in per_sequence
