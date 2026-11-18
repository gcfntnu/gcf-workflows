"""UNITAS misc_RNA by gene family (#379): the family rules, the misc_RNA rows of annotations.tsv and the annotation JSON,
the MultiQC section, and the rules in the delivery dry run. Synthetic UNITAS-style files only."""
import io
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
import yaml

from test_quant_methods import dry_run, rule_block

ROOT = Path(__file__).resolve().parents[2]
QUANT = ROOT / 'smallrna/rules/quant/scripts'
BFQ = ROOT / 'smallrna/rules/bfq/scripts'
sys.path.insert(0, str(QUANT))
import unitas_misc_rna  # noqa: E402

# One read split between a Y RNA gene and a 7SK pseudogene (0.5 each), one between two Y RNA genes; other classes
# have targets with family-like names that must not count.
HITS = [('misc_RNA', 'RNY4-ENSG00000000001.1', 2.0), ('misc_RNA', 'Y_RNA-ENSG00000000002.1', 0.5),
        ('misc_RNA', 'RNY1P3-ENSG00000000003.1', 0.5), ('misc_RNA', 'RN7SKP9-ENSG00000000004.1', 0.5),
        ('misc_RNA', 'RN7SL1-ENSG00000000005.1', 1.0), ('misc_RNA', 'Metazoa_SRP-ENSG00000000006.1', 0.25),
        ('misc_RNA', 'Vault-ENSG00000000007.1', 0.75), ('lncRNA', 'RNY4-AS1-ENSG00000000008.1', 3.0),
        ('snRNA', 'RNU6-1-ENSG00000000009.1', 4.0)]
MISC_RNA = 5.5
SUMMARY = f'miRNA\t20\n   miRNA:homo_sapiens\t18\n   miRNA:other\t2\nlncRNA\t3\nmisc_RNA\t{MISC_RNA}\nsnRNA\t4\nno annotation\t7\n'


def write_sample(base, sample, summary=SUMMARY, hits=HITS):
    out = base / sample
    out.mkdir(parents=True)
    (out / 'unitas.annotation_summary.txt').write_text(summary)
    (out / 'unitas.hits_per_target.txt').write_text(
        'TRANSCRIPT_CLASS\tTRANSCRIPT_NAME\tREAD_COUNT\n' + ''.join(f'{c}\t{n}\t{v}\n' for c, n, v in hits))
    return out


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True)


def annotations(base, samples):
    summaries = [base / s / 'unitas.annotation_summary.txt' for s in samples]
    hits = [base / s / 'unitas.hits_per_target.txt' for s in samples]
    return run(QUANT / 'unitas_annotations.py', *summaries, '--hits-per-target', *hits)


@pytest.mark.parametrize('name, family', [
    ('RNY4-ENSG00000252316.1', 'Y_RNA'), ('RNY4P25-ENSG00000238711.1', 'Y_RNA'), ('Y_RNA-ENSG00000200615.1', 'Y_RNA'),
    ('Rny1-ENSMUSG00000065701.3', 'Y_RNA'), ('RN7SK-ENSG00000283293.1', '7SK'), ('RN7SKP279-ENSG00000000001.1', '7SK'),
    ('7SK-ENSG00000202198.1', '7SK'), ('Rn7sk-ENSMUSG00000000001.1', '7SK'), ('RN7SL1-ENSG00000276168.1', '7SL'),
    ('RN7SL5P-ENSG00000000002.1', '7SL'), ('Metazoa_SRP-ENSG00000000003.1', '7SL'), ('Rn7sl1-ENSRNOG00000085881.1', '7SL'),
    ('Vault-ENSG00000000004.1', 'other'), ('VTRNA1-1-ENSG00000000005.1', 'other'), ('Gm22199-ENSMUSG00000084609.4', 'other'),
    ('>ENST00000621902.1', 'other'), ('XRNY-ENSG00000000006.1', 'other')])
def test_family_rules(name, family):
    assert unitas_misc_rna.family(name) == family


def test_families_keep_fractional_counts_and_add_up(tmp_path):
    hits = write_sample(tmp_path, 'S1') / 'unitas.hits_per_target.txt'
    counts = unitas_misc_rna.misc_rna_families(hits)
    assert counts == {'Y_RNA': 3.0, '7SK': 0.5, '7SL': 1.25, 'other': 0.75}
    unitas_misc_rna.check_total(counts, MISC_RNA, hits)
    with pytest.raises(ValueError, match='add up to 5.5 reads, the annotation summary has 6'):
        unitas_misc_rna.check_total(counts, 6.0, hits)


def test_annotations_add_family_rows_under_unchanged_misc_rna(tmp_path):
    write_sample(tmp_path, 'S1')
    result = annotations(tmp_path, ['S1'])
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(io.StringIO(result.stdout), sep='\t', keep_default_na=False)
    misc = table[table.Level1 == 'misc_RNA']
    assert list(zip(misc.Level2, misc.Count)) == [('', MISC_RNA), ('7SK', 0.5), ('7SL', 1.25), ('Y_RNA', 3.0), ('other', 0.75)]
    rest = table[table.Level1 != 'misc_RNA']
    assert list(zip(rest.Level1, rest.Level2, rest.Count)) == [
        ('lncRNA', '', 3), ('miRNA', '', 20), ('miRNA', 'miRNA:homo_sapiens', 18), ('miRNA', 'miRNA:other', 2),
        ('no annotation', '', 7), ('snRNA', '', 4)]
    top = table[table.Level2 == '']
    assert top.Count.sum() == 20 + 3 + MISC_RNA + 4 + 7


def test_annotations_stop_when_families_and_summary_disagree(tmp_path):
    write_sample(tmp_path, 'S1', summary=SUMMARY.replace(f'misc_RNA\t{MISC_RNA}', 'misc_RNA\t9'))
    result = annotations(tmp_path, ['S1'])
    assert result.returncode != 0
    assert 'misc_RNA targets add up to 5.5 reads, the annotation summary has 9' in result.stderr


def test_sample_without_misc_rna_gets_no_family_rows(tmp_path):
    write_sample(tmp_path, 'S1')
    write_sample(tmp_path, 'S2', summary='miRNA\t20\n', hits=[('lncRNA', 'RNY4-AS1-ENSG00000000008.1', 3.0)])
    result = annotations(tmp_path, ['S1', 'S2'])
    assert result.returncode == 0, result.stderr
    table = pd.read_csv(io.StringIO(result.stdout), sep='\t', keep_default_na=False)
    assert list(table[table.Sample_ID == 'S2'].Level1) == ['miRNA']
    write_sample(tmp_path, 'S3', summary='miRNA\t20\n')
    result = annotations(tmp_path, ['S3'])
    assert result.returncode != 0 and 'the annotation summary has 0' in result.stderr


def test_annotations_need_one_hits_table_per_summary(tmp_path):
    for sample in ('S1', 'S2'):
        write_sample(tmp_path, sample)
    result = run(QUANT / 'unitas_annotations.py', tmp_path / 'S1/unitas.annotation_summary.txt',
                 tmp_path / 'S2/unitas.annotation_summary.txt', '--hits-per-target', tmp_path / 'S1/unitas.hits_per_target.txt')
    assert result.returncode != 0 and 'one file per summary' in result.stderr


def test_annotation_json_nests_families_under_misc_rna(tmp_path):
    out = write_sample(tmp_path, 'S1')
    result = run(QUANT / 'unitas_annotation2json.py', out / 'unitas.annotation_summary.txt',
                 '--hits-per-target', out / 'unitas.hits_per_target.txt')
    assert result.returncode == 0, result.stderr
    tree = json.loads(result.stdout)
    assert tree['misc_RNA'] == {'Y_RNA': 3.0, '7SK': 0.5, '7SL': 1.25, 'other': 0.75}
    assert tree['lncRNA'] == 3.0 and tree['snRNA'] == 4.0


def test_multiqc_section_from_annotations(tmp_path):
    write_sample(tmp_path, 'S1')
    write_sample(tmp_path, 'S2', summary='miRNA\t20\n', hits=[])
    table = tmp_path / 'annotations.tsv'
    table.write_text(annotations(tmp_path, ['S1', 'S2']).stdout)
    result = run(BFQ / 'unitas_misc_rna_mqc.py', table)
    assert result.returncode == 0, result.stderr
    section = yaml.safe_load(result.stdout)
    assert section['id'] == 'unitas_misc_rna' and section['plot_type'] == 'bargraph'
    assert list(section['categories']) == ['Y_RNA', '7SK', '7SL', 'other']
    assert section['data'] == {'S1': {'Y_RNA': 3.0, '7SK': 0.5, '7SL': 1.25, 'other': 0.75},
                               'S2': {'Y_RNA': 0.0, '7SK': 0.0, '7SL': 0.0, 'other': 0.0}}


def test_delivery_plans_the_family_split(tmp_path):
    returncode, output = dry_run(tmp_path, 'unitas', printshell=True, target='multiqc_report')
    assert returncode == 0, output
    annotations = rule_block(output, 'unitas_annotations')
    command = next(line for line in annotations.splitlines() if 'unitas_annotations.py' in line)
    summaries, hits = command.split('--hits-per-target')
    assert summaries.count('unitas.annotation_summary.txt') == hits.count('unitas.hits_per_target.txt') > 0
    section = rule_block(output, 'bfq_level2_unitas_misc_rna')
    assert 'data/tmp/smallrna/bfq/figs/unitas_misc_rna_mqc.yaml' in section
    assert 'unitas_misc_rna_mqc.py data/tmp/smallrna/quant/unitas/annotations.tsv' in section
    assert 'unitas_misc_rna_mqc.yaml' in rule_block(output, 'multiqc_report')
