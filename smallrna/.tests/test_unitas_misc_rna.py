"""UNITAS misc_RNA by gene family (#379): the family rules and table, its checks, the MultiQC section and the rules in
the delivery dry run. Synthetic UNITAS-style files only."""
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
QUANT = ROOT / 'smallrna/rules/quant/scripts'
SCRIPT = QUANT / 'unitas_misc_rna.py'
sys.path.insert(0, str(QUANT))
import unitas_misc_rna  # noqa: E402

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


def write_sample(base, sample, summary=SUMMARY, hits=HITS, hits_text=None):
    out = base / sample
    out.mkdir(parents=True)
    (out / 'unitas.annotation_summary.txt').write_text(summary)
    (out / 'unitas.hits_per_target.txt').write_text(
        hits_text if hits_text is not None else HITS_HEADER + ''.join(f'{c}\t{n}\t{v}\n' for c, n, v in hits))
    return out


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True)


def family_table(base, samples):
    return run(SCRIPT, 'table', '--summaries', *[base / s / 'unitas.annotation_summary.txt' for s in samples],
               '--hits-per-target', *[base / s / 'unitas.hits_per_target.txt' for s in samples])


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
    assert result.stdout == ('Sample_ID\tfamily\tcount\nS1\tY_RNA\t1\nS1\t7SK\t0\nS1\t7SL\t0\nS1\tother\t0\n'
                             'S2\tY_RNA\t3\nS2\t7SK\t0.5\nS2\t7SL\t1.25\nS2\tother\t0.75\n')


@pytest.mark.parametrize('hits_text', ['', HITS_HEADER])
def test_sample_without_reads_has_no_targets(tmp_path, hits_text):
    write_sample(tmp_path, 'S1', summary='', hits_text=hits_text)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[1:] == ['S1\tY_RNA\t0', 'S1\t7SK\t0', 'S1\t7SL\t0', 'S1\tother\t0']


@pytest.mark.parametrize('summary, hits_text, message', [
    (SUMMARY.replace('misc_RNA\t5.5', 'misc_RNA\t9'), None, 'misc_RNA targets add up to 5.5 reads, the annotation summary has 9'),
    ('miRNA\t20\n', None, 'misc_RNA targets add up to 5.5 reads, the annotation summary has 0'),
    (SUMMARY, 'class\tname\tcount\nmisc_RNA\tRNY4\t1\n', 'not a UNITAS hits_per_target table'),
    (SUMMARY, HITS_HEADER + 'misc_RNA\tRNY4\n', 'expected 3 columns')])
def test_table_stops_on_inconsistent_input(tmp_path, summary, hits_text, message):
    write_sample(tmp_path, 'S1', summary=summary, hits_text=hits_text)
    result = family_table(tmp_path, ['S1'])
    assert result.returncode != 0 and message in result.stderr


def test_table_needs_one_hits_table_per_summary(tmp_path):
    for sample in ('S1', 'S2'):
        write_sample(tmp_path, sample)
    result = run(SCRIPT, 'table', '--summaries', tmp_path / 'S1/unitas.annotation_summary.txt',
                 tmp_path / 'S2/unitas.annotation_summary.txt', '--hits-per-target', tmp_path / 'S1/unitas.hits_per_target.txt')
    assert result.returncode != 0 and 'one summary and one hits table per sample directory' in result.stderr


def test_quote_in_a_target_name_does_not_merge_fields(tmp_path):
    """A `"` in a name is part of the name; with CSV quoting, it would join the following rows into one field."""
    write_sample(tmp_path, 'S1', summary='misc_RNA\t3\n', hits=[('misc_RNA', '"RNY4-ENSG00000000001.1', 1.0),
                                                                 ('misc_RNA', 'RN7SL1-ENSG00000000005.1"', 2.0)])
    result = family_table(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[1:] == ['S1\tY_RNA\t0', 'S1\t7SK\t0', 'S1\t7SL\t2', 'S1\tother\t1']


def mqc_section(tmp_path, samples):
    table = tmp_path / 'misc_rna_families.tsv'
    table.write_text(family_table(tmp_path, samples).stdout)
    result = run(SCRIPT, 'mqc', table)
    assert result.returncode == 0, result.stderr
    return yaml.safe_load(result.stdout)


def test_multiqc_section_from_the_family_table(tmp_path):
    write_sample(tmp_path, 'S1')
    write_sample(tmp_path, 'S2', summary='miRNA\t20\n', hits=[])
    section = mqc_section(tmp_path, ['S1', 'S2'])
    assert section['id'] == 'unitas_misc_rna' and section['plot_type'] == 'bargraph'
    assert section['categories'] == {name: {'name': label} for name, label, _ in unitas_misc_rna.FAMILIES}
    assert unitas_misc_rna.RULES in section['description'] and 'Metazoa_SRP...' in section['description']
    assert section['data'] == {'S1': FAMILIES, 'S2': dict.fromkeys(FAMILIES, 0.0)}


def test_multiqc_section_is_a_note_without_misc_rna(tmp_path):
    write_sample(tmp_path, 'S1', summary='miRNA\t20\n', hits=[])
    section = mqc_section(tmp_path, ['S1'])
    assert section['plot_type'] == 'html' and section['data'] == '<p>No sample has misc_RNA reads.</p>'


def test_delivery_plans_the_family_split(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas', printshell=True, target='multiqc_report')
    assert returncode == 0, output
    families = rule_block(output, 'unitas_misc_rna')
    command = next(line for line in families.splitlines() if 'unitas_misc_rna.py table' in line)
    summaries, hits = command.split('--hits-per-target')
    assert summaries.count('unitas.annotation_summary.txt') == hits.count('unitas.hits_per_target.txt') > 0
    assert 'data/tmp/smallrna/quant/unitas/misc_rna_families.tsv' in families
    assert 'hits_per_target' not in rule_block(output, 'unitas_annotations')
    section = rule_block(output, 'bfq_level2_unitas_misc_rna')
    assert 'unitas_misc_rna.py mqc data/tmp/smallrna/quant/unitas/misc_rna_families.tsv' in section
    assert 'data/tmp/smallrna/bfq/exprs/misc_rna_families.tsv' in rule_block(output, 'bfq_level2_misc_rna_families')
    multiqc = rule_block(output, 'multiqc_report')
    assert 'figs/unitas_misc_rna_mqc.yaml' in multiqc and 'exprs/misc_rna_families.tsv' in multiqc
