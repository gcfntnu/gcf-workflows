"""UNITAS misc_RNA by gene family (#379): the family rules and table, its checks against annotations.tsv, the MultiQC
section and the rules in the quant_all and delivery dry runs. Synthetic UNITAS-style files only."""
import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
QUANT = ROOT / 'smallrna/rules/quant/scripts'
SCRIPT = QUANT / 'unitas_misc_rna.py'
spec = importlib.util.spec_from_file_location('unitas_misc_rna', SCRIPT)
unitas_misc_rna = importlib.util.module_from_spec(spec)
spec.loader.exec_module(unitas_misc_rna)
with open(ROOT / 'smallrna/.tests/pep/sample_table.csv', newline='') as fh:
    FIXTURE_SAMPLES = [row['sample_name'] for row in csv.DictReader(fh)]

# One read split between a Y RNA gene and a 7SK pseudogene (0.5 each), one between two Y RNA genes; a transcript-only
# name; other classes have targets with family-like names that must not count.
HITS = [('misc_RNA', 'RNY4-ENSG00000000001.1', 2.0), ('misc_RNA', 'Y_RNA-ENSG00000000002.1', 0.5),
        ('misc_RNA', 'RNY1P3-ENSG00000000003.1', 0.5), ('misc_RNA', 'RN7SKP9-ENSG00000000004.1', 0.5),
        ('misc_RNA', 'RN7SL1-ENSG00000000005.1', 1.0), ('misc_RNA', 'Metazoa_SRP-ENSG00000000006.1', 0.25),
        ('misc_RNA', 'Vault-ENSG00000000007.1', 0.5), ('misc_RNA', '>ENST00000000010.1', 0.25),
        ('lncRNA', 'RNY4-AS1-ENSG00000000008.1', 3.0), ('snRNA', 'RNU6-1-ENSG00000000009.1', 4.0)]
FAMILIES = {'Y_RNA': 3.0, '7SK': 0.5, '7SL': 1.25, 'other': 0.75}
SUMMARY = 'miRNA\t20\n   miRNA:homo_sapiens\t18\n   miRNA:other\t2\nlncRNA\t3\nmisc_RNA\t5.5\nsnRNA\t4\nno annotation\t7\n'
HITS_HEADER = 'TRANSCRIPT_CLASS\tTRANSCRIPT_NAME\tREAD_COUNT\n'
ZEROS = ['Y_RNA\t0.0', '7SK\t0.0', '7SL\t0.0', 'other\t0.0']


def write_sample(base, sample, summary=SUMMARY, hits=HITS, hits_text=None):
    out = base / sample
    out.mkdir(parents=True)
    (out / 'unitas.annotation_summary.txt').write_text(summary)
    (out / 'unitas.hits_per_target.txt').write_text(
        hits_text if hits_text is not None else HITS_HEADER + ''.join(f'{c}\t{n}\t{v}\n' for c, n, v in hits))
    return out


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True)


def annotations(base, samples):
    """annotations.tsv as the workflow makes it (unitas_annotations.py) from the samples' summaries."""
    result = run(QUANT / 'unitas_annotations.py', *[base / s / 'unitas.annotation_summary.txt' for s in samples])
    assert result.returncode == 0, result.stderr
    path = base / 'annotations.tsv'
    path.write_text(result.stdout)
    return path


def family_table(base, samples, annotations_path=None):
    """Run the script; the table and the MultiQC section are attached to the result when it succeeds."""
    out = base / 'out'
    out.mkdir(exist_ok=True)
    result = run(SCRIPT, '--annotations', annotations_path or annotations(base, samples),
                 '--hits-per-target', *[base / s / 'unitas.hits_per_target.txt' for s in samples],
                 '--table', out / 'misc_rna_families.tsv', '--mqc', out / 'unitas_misc_rna_mqc.yaml')
    if result.returncode == 0:
        result.table = (out / 'misc_rna_families.tsv').read_text()
        result.section = yaml.safe_load((out / 'unitas_misc_rna_mqc.yaml').read_text())
    else:
        assert not (out / 'misc_rna_families.tsv').exists() and not (out / 'unitas_misc_rna_mqc.yaml').exists()
    return result


@pytest.mark.parametrize('name, family', [
    ('RNY4-ENSG00000252316.1', 'Y_RNA'), ('RNY4P25-ENSG00000238711.1', 'Y_RNA'), ('Y_RNA-ENSG00000200615.1', 'Y_RNA'),
    ('Rny1-ENSMUSG00000065701.3', 'Y_RNA'), ('RN7SK-ENSG00000283293.1', '7SK'), ('RN7SKP279-ENSG00000000001.1', '7SK'),
    ('7SK-ENSG00000202198.1', '7SK'), ('Rn7sk-ENSMUSG00000000001.1', '7SK'), ('RN7SL1-ENSG00000276168.1', '7SL'),
    ('RN7SL5P-ENSG00000000002.1', '7SL'), ('Metazoa_SRP-ENSG00000000003.1', '7SL'), ('Rn7sl1-ENSRNOG00000085881.1', '7SL'),
    ('Rn7s1-ENSMUSG00000118866.1', '7SL'), ('Rn7s6-ENSMUSG00000092746.5', '7SL'), ('RN7SKP80-ENSG00000000011.1', '7SK'),
    ('Vault-ENSG00000000004.1', 'other'), ('VTRNA1-1-ENSG00000000005.1', 'other'), ('Gm22199-ENSMUSG00000084609.4', 'other'),
    ('>ENST00000621902.1', 'other'), ('XRNY-ENSG00000000006.1', 'other')])
def test_family_rules(name, family):
    assert unitas_misc_rna.family(name) == family


def test_table_keeps_fractional_counts_and_lists_every_family(tmp_path):
    write_sample(tmp_path, 'S2')
    write_sample(tmp_path, 'S1', summary='miRNA\t20\nmisc_RNA\t1\n', hits=[('misc_RNA', 'RNY4-ENSG00000000001.1', 1.0)])
    result = family_table(tmp_path, ['S2', 'S1'])
    assert result.returncode == 0, result.stderr
    assert result.table == ('Sample_ID\tfamily\tcount\nS1\tY_RNA\t1.0\nS1\t7SK\t0.0\nS1\t7SL\t0.0\nS1\tother\t0.0\n'
                            'S2\tY_RNA\t3.0\nS2\t7SK\t0.5\nS2\t7SL\t1.25\nS2\tother\t0.75\n')


def test_counts_are_written_without_losing_precision(tmp_path):
    write_sample(tmp_path, 'S1', summary='misc_RNA\t1234567.123456789\n',
                 hits=[('misc_RNA', 'RNY4-ENSG00000000001.1', '1234567.0'), ('misc_RNA', 'RNY1-ENSG00000000012.1', '0.123456789')])
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1] == 'S1\tY_RNA\t1234567.123456789'


def test_families_are_summed_without_float_noise(tmp_path):
    """0.3 + 0.3 + 0.3 + 0.4 added in turn is 1.2999999999999998; math.fsum gives 1.3."""
    hits = [('misc_RNA', f'RNY{i}-ENSG0000000000{i}.1', v) for i, v in enumerate([0.3, 0.3, 0.3, 0.4], start=1)]
    write_sample(tmp_path, 'S1', summary='misc_RNA\t1.3\n', hits=hits)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1] == 'S1\tY_RNA\t1.3'


@pytest.mark.parametrize('value, text', [(5e-05, '0.00005'), (1.5e-07, '0.00000015'), (3.0, '3.0'), (0.0, '0.0'),
                                         (1e16, '10000000000000000.0'), (3598991.109777921, '3598991.109777921')])
def test_counts_are_plain_decimals(value, text):
    assert unitas_misc_rna.decimal_count(value) == text


def test_small_counts_are_written_without_scientific_notation(tmp_path):
    write_sample(tmp_path, 'S1', summary='misc_RNA\t0.00005\n', hits=[('misc_RNA', 'RN7SK-ENSG00000283293.1', '5e-05')])
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1:] == ['S1\tY_RNA\t0.0', 'S1\t7SK\t0.00005', 'S1\t7SL\t0.0', 'S1\tother\t0.0']


@pytest.mark.parametrize('hits_text', ['', HITS_HEADER, '\n' + HITS_HEADER + '  \n'])
def test_sample_without_reads_has_no_targets(tmp_path, hits_text):
    write_sample(tmp_path, 'S1', summary='miRNA\t0\n', hits_text=hits_text)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1:] == ['S1\t' + z for z in ZEROS]


def test_blank_and_space_only_lines_in_the_hits_table_are_skipped(tmp_path):
    rows = [f'{c}\t{n}\t{v}\n' for c, n, v in HITS]
    write_sample(tmp_path, 'S1', hits_text=HITS_HEADER + ''.join(rows[:3]) + '\n   \n' + ''.join(rows[3:]) + '\n \n')
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.section['data'] == {'S1': FAMILIES}


def test_quote_in_a_target_name_does_not_merge_fields(tmp_path):
    """A `"` in a name is part of the name; with CSV quoting, it would join the following rows into one field."""
    write_sample(tmp_path, 'S1', summary='misc_RNA\t3\n', hits=[('misc_RNA', '"RNY4-ENSG00000000001.1', 1.0),
                                                                 ('misc_RNA', 'RN7SL1-ENSG00000000005.1"', 2.0)])
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1:] == ['S1\tY_RNA\t0.0', 'S1\t7SK\t0.0', 'S1\t7SL\t2.0', 'S1\tother\t1.0']


@pytest.mark.parametrize('summary, hits_text, message', [
    (SUMMARY.replace('misc_RNA\t5.5', 'misc_RNA\t9'), None, 'misc_RNA targets add up to 5.5 reads, annotations.tsv has 9'),
    ('miRNA\t20\n', None, 'misc_RNA targets add up to 5.5 reads, annotations.tsv has 0'),
    (SUMMARY, 'class\tname\tcount\nmisc_RNA\tRNY4\t1\n', 'unitas.hits_per_target.txt: line 1: not a UNITAS hits_per_target header'),
    (SUMMARY, HITS_HEADER + '\nmisc_RNA\tRNY4\n', "unitas.hits_per_target.txt: line 3: expected a class, a target name and a count"),
    (SUMMARY, HITS_HEADER + '\t\t\n', "unitas.hits_per_target.txt: line 2: expected a class, a target name and a count, got ['', '', '']"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\t \t1\n', 'unitas.hits_per_target.txt: line 2: expected a class, a target name and a count'),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\tx\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got 'x'"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\tnan\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got 'nan'"),
    (SUMMARY, HITS_HEADER + 'lncRNA\tRNY4-AS1\tinf\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got 'inf'"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\t-1\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got '-1'"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\t1_000\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got '1_000'"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\t1e999\n', "unitas.hits_per_target.txt: line 2: expected a non-negative decimal number, got '1e999'"),
    (SUMMARY, HITS_HEADER + 'misc_RNA\t' + 'R' * 200000 + '\t1\n', 'unitas.hits_per_target.txt: line 2: field larger than field limit'),
    ('miRNA\t20\nmisc_RNA\tnan\n', None, "annotations.tsv: line 3: expected a non-negative decimal number, got 'nan'")])
def test_table_stops_on_inconsistent_input(tmp_path, summary, hits_text, message):
    write_sample(tmp_path, 'S1', summary=summary, hits_text=hits_text)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode != 0 and message in result.stderr and 'Traceback' not in result.stderr


def test_total_check_fails_on_nan():
    counts = dict(FAMILIES, Y_RNA=math.nan)
    with pytest.raises(ValueError, match='add up to nan reads'):
        unitas_misc_rna.check_total(counts, 5.5, 'x')
    with pytest.raises(ValueError, match='annotations.tsv has nan'):
        unitas_misc_rna.check_total(FAMILIES, math.nan, 'x')


def test_sample_without_annotation_rows_and_without_misc_rna_hits_is_zero(tmp_path):
    """An empty summary (a sample without reads) gives no rows in annotations.tsv; with no misc_RNA hits, its families are 0."""
    write_sample(tmp_path, 'S1')
    write_sample(tmp_path, 'S2', summary='', hits_text='')
    write_sample(tmp_path, 'S3', summary='', hits=[('lncRNA', 'RNY4-AS1-ENSG00000000008.1', 3.0)])
    table = annotations(tmp_path, ['S1', 'S2', 'S3'])
    assert {line.split('\t')[-1] for line in table.read_text().splitlines()[1:]} == {'S1'}
    result = family_table(tmp_path, ['S1', 'S2', 'S3'], annotations_path=table)
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[5:] == ['S2\t' + z for z in ZEROS] + ['S3\t' + z for z in ZEROS]


@pytest.mark.parametrize('summary', ['', 'miRNA\t20\n'])
def test_sample_with_misc_rna_below_the_tolerance_passes_with_or_without_annotation_rows(tmp_path, summary):
    """An absent sample (empty summary) and a sample without a misc_RNA row are both checked against 0 with the tolerance."""
    write_sample(tmp_path, 'S1', summary=summary, hits=[('misc_RNA', 'RNY4-ENSG00000000001.1', 0.005)])
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.table.splitlines()[1] == 'S1\tY_RNA\t0.005'


def test_samples_must_match_annotations_and_hits_tables(tmp_path):
    for sample in ('S1', 'S2'):
        write_sample(tmp_path, sample)
    result = family_table(tmp_path, ['S1', 'S2'], annotations_path=annotations(tmp_path, ['S1']))
    assert result.returncode != 0
    assert 'misc_RNA targets add up to 5.5 reads, annotations.tsv has no rows for the sample, so 0' in result.stderr
    result = family_table(tmp_path, ['S1'], annotations_path=annotations(tmp_path, ['S1', 'S2']))
    assert result.returncode != 0 and "annotations.tsv: samples ['S2'] are in it without hits tables" in result.stderr


def test_annotations_table_with_two_misc_rna_rows_for_a_sample_fails(tmp_path):
    write_sample(tmp_path, 'S1')
    table = annotations(tmp_path, ['S1'])
    table.write_text(table.read_text() + 'misc_RNA\t\t\t5.5\tS1\n')
    result = family_table(tmp_path, ['S1'], annotations_path=table)
    assert result.returncode != 0 and 'annotations.tsv: line 9: a second misc_RNA row for sample S1' in result.stderr


def test_annotations_table_must_have_its_columns(tmp_path):
    write_sample(tmp_path, 'S1')
    bad = tmp_path / 'bad.tsv'
    bad.write_text('Sample_ID\tclass\tcount\nS1\tmisc_RNA\t5.5\n')
    result = family_table(tmp_path, ['S1'], annotations_path=bad)
    assert result.returncode != 0 and 'bad.tsv: line 1: not a UNITAS annotations table' in result.stderr


def test_multiqc_section_with_the_table(tmp_path):
    write_sample(tmp_path, 'S1')
    write_sample(tmp_path, 'S2', summary='miRNA\t20\n', hits=[])
    section = family_table(tmp_path, ['S1', 'S2']).section
    assert section['id'] == 'unitas_misc_rna' and section['plot_type'] == 'bargraph'
    assert section['categories'] == {name: {'name': label} for name, label, _ in unitas_misc_rna.FAMILIES}
    assert unitas_misc_rna.RULES in section['description'] and 'Metazoa_SRP...' in section['description']
    assert section['data'] == {'S1': FAMILIES, 'S2': dict.fromkeys(FAMILIES, 0.0)}


@pytest.mark.parametrize('summary, hits', [('miRNA\t20\n', []), ('misc_RNA\t0.004\n', [('misc_RNA', 'RNY4-ENSG00000000001.1', 0.004)])])
def test_multiqc_section_is_a_note_without_misc_rna_to_plot(tmp_path, summary, hits):
    """No misc_RNA, or only counts that round to 0 at the 2 decimals plotted: a note, not an empty plot."""
    write_sample(tmp_path, 'S1', summary=summary, hits=hits)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.section['plot_type'] == 'html'
    assert result.section['data'] == '<p>No sample has misc_RNA reads to plot (all below 0.005).</p>'


def test_quant_all_builds_the_family_table(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas')
    assert returncode == 0, output
    assert 'data/tmp/smallrna/quant/unitas/misc_rna_families.tsv' in rule_block(output, 'unitas_misc_rna')


def test_delivery_plans_the_family_split(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas', printshell=True, target='multiqc_report')
    assert returncode == 0, output
    families = rule_block(output, 'unitas_misc_rna')
    command = next(line for line in families.splitlines() if 'unitas_misc_rna.py' in line)
    assert '--annotations data/tmp/smallrna/quant/unitas/annotations.tsv --hits-per-target ' in command
    assert len(FIXTURE_SAMPLES) > 1 and 'annotation_summary' not in command
    assert command.count('unitas.hits_per_target.txt') == len(FIXTURE_SAMPLES)
    assert all(f'data/tmp/smallrna/quant/unitas/unitas/{s}/unitas.hits_per_target.txt' in command for s in FIXTURE_SAMPLES)
    assert ('--table data/tmp/smallrna/quant/unitas/misc_rna_families.tsv '
            '--mqc data/tmp/smallrna/quant/unitas/unitas_misc_rna_mqc.yaml') in command
    assert 'hits_per_target' not in rule_block(output, 'unitas_annotations')
    links = rule_block(output, 'bfq_level2_unitas_misc_rna')
    assert ('ln -srf data/tmp/smallrna/quant/unitas/misc_rna_families.tsv data/tmp/smallrna/bfq/exprs/misc_rna_families.tsv && '
            'ln -srf data/tmp/smallrna/quant/unitas/unitas_misc_rna_mqc.yaml data/tmp/smallrna/bfq/figs/unitas_misc_rna_mqc.yaml') in links
    multiqc = rule_block(output, 'multiqc_report')
    assert 'figs/unitas_misc_rna_mqc.yaml' in multiqc and 'exprs/misc_rna_families.tsv' in multiqc
