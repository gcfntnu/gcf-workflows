"""nf-core/smrnaseq's mirtop counting route (#403): planning with the pipeline's arguments, delivery unchanged, and the
hairpin and table scripts reproducing the pipeline's steps."""
from pathlib import Path
import re
import subprocess
import sys

import pytest

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'smallrna/rules/quant/scripts'
OUT = 'data/tmp/smallrna/quant/smrnaseq'
INDEX = ('1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt')
ANNOTATION = 'UID\tRead\tmiRNA\tVariant\tiso_5p\tiso_3p\tiso_add3p\tiso_snp\tiso_5p_nt\tiso_3p_nt\tiso_add3p_nt\tiso_snp_nt'


def test_route_plans_the_pipelines_steps_on_the_filtered_reads_and_delivery_stays_on_unitas(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas,smrnaseq', printshell=True)
    assert returncode == 0, output
    assert 'rule unitas_mirtable:' in output
    for rule, n in (('smrnaseq_hairpin', 1), ('smrnaseq_hairpin_index', 1), ('smrnaseq_collapse', 4), ('smrnaseq_align', 4),
                    ('smrnaseq_mirtop', 4), ('smrnaseq_counts', 1)):
        assert output.count(f'rule {rule}:') == n, rule
    hairpin = rule_block(output, 'smrnaseq_hairpin')
    assert 'input: data/ext/mirbase/mus_musculus/release-23.1/fasta/hairpin.fa' in hairpin  # db.mirbase.release
    assert '--species mmu' in hairpin
    assert re.search(rf'bowtie-build --threads \d+ {OUT}/reference/hairpin.fa {OUT}/reference/bowtie/hairpin ', output)
    collapse = rule_block(output, 'smrnaseq_collapse')
    assert 'input: data/tmp/smallrna/filter/cutadapt/SR0' in collapse  # the shared filtered reads
    assert 'seqcluster collapse -m 1 --min_size 15 -f data/tmp/smallrna/filter/cutadapt/SR0' in collapse
    align = rule_block(output, 'smrnaseq_align')
    assert f'output: {OUT}/bam/SR0' in align  # named by sample: mirtop's sample column
    assert re.search(rf'bowtie --threads \d+ --sam -x {OUT}/reference/bowtie/hairpin -q -t -k 50 --best --strata -e 99999 '
                     rf'--chunkmbs 2048 {OUT}/collapsed/SR0\d.fastq 2> {OUT}/bam/SR0\d.log \| samtools view -bS -o {OUT}/bam/SR0\d.bam -', align)
    for suffix in INDEX:
        assert f'{OUT}/reference/bowtie/hairpin.{suffix}' in align.split('output:', 1)[0]
    mirtop = rule_block(output, 'smrnaseq_mirtop')
    gff = 'data/ext/mirbase/mus_musculus/release-23.1/anno/genes.gff3'
    assert f'mirtop gff --sps mmu --hairpin {OUT}/reference/hairpin.fa --gtf {gff} -o {OUT}/mirtop/SR0' in mirtop
    assert f'mirtop counts --add-extra --hairpin {OUT}/reference/hairpin.fa --gtf {gff} --sps mmu --gff {OUT}/mirtop/SR0' in mirtop
    assert f'--joined {OUT}/joined_samples_mirtop.tsv --mirna {OUT}/mirna.tsv' in rule_block(output, 'smrnaseq_counts')
    returncode, output = dry_run(tmp_path, 'unitas,smrnaseq', target='bfq_all')
    assert returncode == 0, output
    assert 'smrnaseq' not in output


def test_route_quotes_paths_with_spaces(tmp_path):
    """A reference root with a space reaches every command of the route as one quoted argument."""
    ext = tmp_path / 'ref root'
    returncode, output = dry_run(tmp_path / 'project', 'smrnaseq', lambda s: s.update(ext_dir=str(ext)), printshell=True)
    assert returncode == 0, output
    hairpin = f"'{ext}/mirbase/mus_musculus/release-23.1/fasta/hairpin.fa'"
    assert f'--species mmu {hairpin} > ' in rule_block(output, 'smrnaseq_hairpin')
    assert f"--gtf '{ext}/mirbase/mus_musculus/release-23.1/anno/genes.gff3' -o " in rule_block(output, 'smrnaseq_mirtop')


def run(script, *args, stdin=None):
    return subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True, input=stdin)


def test_hairpins_of_one_species_as_the_pipeline_prepares_them(tmp_path):
    fasta = tmp_path / 'hairpin.fa'
    fasta.write_text('<p>&gt;cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop<br>UACACUGUGGAUCCGG<br>'
                     '&gt;hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop<br>UGGGAUGAGGUAGUAGG<br>UUGUAUAGuuNUU</p>\n'
                     '>hsa-mir-21 MI0000077 Homo sapiens miR-21 stem-loop\nUGUCGGGUAGCUUAUCAG\nACUGAUGUUGACUGU\n'
                     '>mmu-mir-21a MI0000569 Mus musculus\nUGUACCACCUUGUCGG\n>xhsa-mir-1 MI9 made up\nACGU\n\n\n')
    result = run('smrnaseq_hairpin.py', '--species', 'hsa', fasta)
    assert result.returncode == 0, result.stderr
    # headers cut at the first space, sequences on one line, U to T (lower case kept), other letters unchanged;
    # the species pattern is .*hsa-.*, as seqkit grep -r applies it
    assert result.stdout == ('>hsa-let-7a-1\nTGGGATGAGGTAGTAGGTTGTATAGttNTT\n>hsa-mir-21\nTGTCGGGTAGCTTATCAGACTGATGTTGACTGT\n'
                             '>xhsa-mir-1\nACGT\n')
    result = run('smrnaseq_hairpin.py', '--species', 'rno', fasta)
    assert result.returncode != 0 and 'no hairpin of species rno' in result.stderr
    result = run('smrnaseq_hairpin.py', '--species', 'Homo sapiens', fasta)
    assert result.returncode != 0 and 'species code' in result.stderr


def counts_table(path, sample, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'{ANNOTATION}\t{sample}\n' + ''.join('\t'.join(r) + '\n' for r in rows))
    return path


def iso(uid, mirna, variant, count, read='TAGCTTATCAGACTGATGTTGA'):
    return [uid, read, mirna, variant, '0', '0', '0', '0', '0', '0', '0', '0', str(count)]


def test_tables_join_samples_and_count_a_read_on_several_mirnas_once(tmp_path):
    s1 = counts_table(tmp_path / 'S1/counts/mirtop.tsv', 'S1', [
        iso('iso-22-B', 'hsa-miR-21-5p', 'NA', 10),
        iso('iso-22-A', 'hsa-miR-590-5p', 'iso_snv', 4),  # one collapsed read on two miRNAs:
        iso('iso-22-A', 'hsa-miR-21-5p', 'iso_snv', 4),   # only the line that sorts first counts
        iso('iso-23-C', 'hsa-miR-16-5p', 'iso_3p:+1', 7)])
    s2 = counts_table(tmp_path / 'S2/counts/mirtop.tsv', 'S2', [
        iso('iso-22-B', 'hsa-miR-21-5p', 'NA', 3), iso('iso-24-D', 'hsa-miR-16-5p', 'iso_add3p:2', 1)])
    joined, mirna = tmp_path / 'joined.tsv', tmp_path / 'mirna.tsv'
    result = run('smrnaseq_counts.py', '--tables', s1, s2, '--joined', joined, '--mirna', mirna)
    assert result.returncode == 0, result.stderr
    lines = joined.read_text().splitlines()
    assert lines[0] == f'{ANNOTATION}\tS1\tS2'
    assert [(l.split('\t')[0], l.split('\t')[2], l.split('\t')[-2:]) for l in lines[1:]] == [
        ('iso-22-A', 'hsa-miR-21-5p', ['4', '0']), ('iso-22-A', 'hsa-miR-590-5p', ['4', '0']),  # by UID, then the whole line
        ('iso-22-B', 'hsa-miR-21-5p', ['10', '3']), ('iso-23-C', 'hsa-miR-16-5p', ['7', '0']), ('iso-24-D', 'hsa-miR-16-5p', ['0', '1'])]
    assert mirna.read_text() == 'miRNA\tS1\tS2\nhsa-miR-21-5p\t14\t3\nhsa-miR-16-5p\t7\t1\n'  # miR-590-5p gets none


@pytest.mark.parametrize('content, message', [
    ('UID\tRead\tmiRNA\tS1\nx\tA\tm\t1\n', 'expected the columns'),
    (f'{ANNOTATION}\tS1\n' + '\t'.join(iso('u', 'm', 'NA', 1)[:-1]) + '\t1.5\n', 'whole number'),
    (f'{ANNOTATION}\tS1\n' + '\t'.join(iso('u', 'm', 'NA', 1)) + '\n' + '\t'.join(iso('u', 'm', 'NA', 2)) + '\n', 'row repeated'),
])
def test_tables_stop_on_unexpected_mirtop_output(tmp_path, content, message):
    table = tmp_path / 'mirtop.tsv'
    table.write_text(content)
    result = run('smrnaseq_counts.py', '--tables', table, '--joined', tmp_path / 'j.tsv', '--mirna', tmp_path / 'm.tsv')
    assert result.returncode != 0 and message in result.stderr


def test_tables_stop_on_a_sample_in_two_tables(tmp_path):
    a = counts_table(tmp_path / 'a/mirtop.tsv', 'S1', [iso('u', 'm', 'NA', 1)])
    b = counts_table(tmp_path / 'b/mirtop.tsv', 'S1', [iso('v', 'm', 'NA', 1)])
    result = run('smrnaseq_counts.py', '--tables', a, b, '--joined', tmp_path / 'j.tsv', '--mirna', tmp_path / 'm.tsv')
    assert result.returncode != 0 and 'sample S1 in two tables' in result.stderr
