"""Unsupported small-RNA routes stop at workflow load instead of failing mid-run."""
import os
import shlex
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

from test_quant_methods import rule_block

ROOT = Path(__file__).resolve().parents[2]


def dry_run(tmp_path, update, printshell=False, target='multiqc_report', dry=True):
    shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    update(settings)
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                      f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    flags = (['-n'] if dry else []) + (['-p'] if printshell else [])
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), *flags, '--cores', '1', target],
                            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    return result.returncode, result.stdout + result.stderr


def set_aligner(name):
    return lambda settings: settings.setdefault('align', {}).setdefault('genome', {}).update(aligner=name)


def set_dedup(value):
    return lambda settings: settings.setdefault('filter', {}).update(dedup=value)


@pytest.mark.parametrize('aligner', ['bowtie2', 'bwa'])
def test_genome_aligner_other_than_star_fails(tmp_path, aligner):
    returncode, output = dry_run(tmp_path, set_aligner(aligner))
    assert returncode != 0
    assert f"align.genome.aligner: smallrna supports only star, got '{aligner}'" in output


def test_dedup_true_fails(tmp_path):
    returncode, output = dry_run(tmp_path, set_dedup(True))
    assert returncode != 0
    assert 'filter.dedup: read deduplication is not supported in smallrna' in output


@pytest.mark.parametrize('update', [set_aligner('star'), set_dedup(False)], ids=['star', 'dedup-false'])
def test_supported_settings_plan_delivery(tmp_path, update):
    returncode, output = dry_run(tmp_path, update)
    assert returncode == 0, output
    assert 'rule star_align:' in output


def set_spikein(ref):
    return lambda settings: settings.setdefault('filter', {}).setdefault('spikein', {}).update(quantifier='unitas', ref=ref)


@pytest.mark.parametrize('ref, fastas', [('ERCC92,smallrna_calibrators', ['ERCC92.fa', 'smallrna_calibrators.fa']),
                                         ('ERCC92', ['ERCC92.fa']), ('smallrna_calibrators', ['smallrna_calibrators.fa'])])
def test_spikein_sets_select_their_fasta_files(tmp_path, ref, fastas):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode == 0, output
    inputs = output.split('rule spikein_unitas:')[1].split('output:')[0]
    assert sorted(path.rsplit('/', 1)[1] for path in inputs.split('input:')[1].strip().split(', ')) == sorted(fastas)


@pytest.mark.parametrize('ref', ['ercc', 'ERCC92,mirbase'])
def test_unknown_spikein_set_fails(tmp_path, ref):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode != 0
    assert 'filter.spikein.ref: unsupported value' in output


@pytest.mark.parametrize('ref', ['', None])
def test_spikein_quantification_without_sets_fails(tmp_path, ref):
    returncode, output = dry_run(tmp_path, set_spikein(ref))
    assert returncode != 0
    assert 'filter.spikein.ref lists no spike-in sets' in output


def test_default_has_no_spikeins(tmp_path):
    returncode, output = dry_run(tmp_path, lambda settings: None)
    assert returncode == 0, output
    assert 'rule spikein_unitas:' not in output and 'spikein/fasta' not in output


def shell_lines(output, needle):
    return [line for line in output.splitlines() if line.startswith('Shell command:') and needle in line]


def trim_route(kit, quantifier):
    def update(settings):
        settings['libprepkit'] = kit
        trim = settings.setdefault('filter', {}).setdefault('trim', {})
        trim['quantifier'] = quantifier
        if quantifier == 'fastp':  # no kit default since #350: the project sets it
            trim['fastp'] = {'args': '-l 15 -q 15 -u 40 -n 5 -e 0 --low_complexity_filter --complexity_threshold 30'}
    return update


MIRTRACE_ROUTES = [('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', 'fastp', 'fastp'),
                   ('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', 'cutadapt', 'cutadapt'),
                   ('Bioo Scientific NEXTflex Small RNA-Seq Kit v4', 'fastp', 'fastp'),
                   ('TruSeq Small RNA Library Prep Kit', 'cutadapt', 'cutadapt'),
                   ('TruSeq Small RNA Library Prep Kit', 'skip', 'fastq')]


@pytest.mark.parametrize('kit, quantifier, stage', MIRTRACE_ROUTES)
def test_mirtrace_reads_the_filtered_reads_without_trimming(tmp_path, kit, quantifier, stage):
    """miRTrace reads the filtered reads the quantifiers read, with no adapter and protocol illumina (#352)."""
    returncode, output = dry_run(tmp_path, trim_route(kit, quantifier), printshell=True)
    assert returncode == 0, output
    mirtrace_input = rule_block(output, 'qc_mirtrace').split('output:', 1)[0]
    assert all(f'data/tmp/smallrna/filter/{stage}/{sample}_R1.fastq' in mirtrace_input for sample in ['SR01', 'SR02', 'SR03', 'SR04'])
    assert '_R1.fastq.gz' not in mirtrace_input  # not the merged reads
    assert f'input: data/tmp/smallrna/filter/{stage}/SR0' in rule_block(output, 'unitas')  # the quantifiers' reads
    mirtrace = shell_lines(output, 'mirtrace qc ')
    assert len(mirtrace) == 1 and ' -p illumina ' in mirtrace[0] and ' -a ' not in mirtrace[0] and '--force' in mirtrace[0]
    assert 'rule bfq_level2_mirtrace:' in output
    assert 'bfq/logs/mirtrace/mirtrace-stats-mirna-complexity.tsv' in output
    assert any('-m mirtrace' in line for line in shell_lines(output, 'multiqc '))


def test_mirtrace_filter_route_is_rejected(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: s.setdefault('filter', {}).setdefault('trim', {}).update(quantifier='mirtrace'))
    assert returncode != 0
    assert "filter.trim.quantifier: smallrna supports fastp, cutadapt or skip, got 'mirtrace'" in output


@pytest.mark.parametrize('quantifier, stage', [('fastp', 'fastp'), ('cutadapt', 'cutadapt')])
def test_mirtrace_config_lists_filtered_reads_without_adapter(tmp_path, quantifier, stage):
    update = trim_route('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', quantifier)
    returncode, output = dry_run(tmp_path, update, target='data/tmp/smallrna/qc/mirtrace/mirtrace.config', dry=False)
    assert returncode == 0, output
    rows = [line.split(',') for line in (tmp_path / 'data/tmp/smallrna/qc/mirtrace/mirtrace.config').read_text().splitlines()]
    assert [row[1] for row in rows] == ['SR01', 'SR02', 'SR03', 'SR04']
    assert all(len(row) == 2 and row[0].endswith(f'/filter/{stage}/{row[1]}_R1.fastq') for row in rows)  # no adapter column


def test_report_shows_mirtrace_contamination_and_complexity_only():
    template = yaml.safe_load((ROOT / 'misc/multiqc/multiqc_config-smallrna.yaml').read_text())
    removed = set(template['remove_sections'])
    assert {'mirtrace_qc', 'mirtrace_rna_categories', 'mirtrace_length'} <= removed
    assert not {'mirtrace_contamination_check', 'mirtrace_complexity'} & removed
    assert 'filtered reads used for quantification' in template['section_comments']['mirtrace']


def test_mirtrace_title_is_one_argument_for_several_projects(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: s.update(project_id=['GCF-0000-006', 'GCF-0000-007']), printshell=True)
    assert returncode == 0, output
    argv = shlex.split(shell_lines(output, 'mirtrace qc ')[0].split('Shell command:', 1)[1])
    title = argv.index('--title')
    assert argv[title + 1] == 'GCF-0000-006,GCF-0000-007'
    assert argv[title + 2] == '-t' and len(argv) == title + 4


@pytest.mark.parametrize('organism', ['danio_rerio', 'human'])
def test_organism_without_unitas_reference_fails(tmp_path, organism):
    returncode, output = dry_run(tmp_path, lambda s: s.update(organism=organism))
    assert returncode != 0
    assert f"UNITAS supports organism homo_sapiens, mus_musculus, rattus_norvegicus, got '{organism}'" in output


def test_unitas_ribosomal_filter_is_rejected(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: s.setdefault('filter', {}).setdefault('ribosomal', {}).update(quantifier='unitas'))
    assert returncode != 0
    assert "filter.ribosomal.quantifier: smallrna supports only skip, got 'unitas'" in output


def test_unitas_uses_the_container_seqmap_and_writes_nothing_in_the_project_root(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: None, printshell=True)
    assert returncode == 0, output
    assert 'seqmap' not in output and '.unitas_dummy' not in output
    unitas = shell_lines(output, 'run_unitas.sh')
    assert len(unitas) == 4 and all('--log data/tmp/smallrna/quant/unitas/logs/SR0' in line for line in unitas)


def set_unitas_refseq(spikein, contaminants):
    def update(settings):
        filters = settings.setdefault('filter', {})
        if spikein:
            filters.setdefault('spikein', {}).update(quantifier='unitas', ref='ERCC92')
        if contaminants:
            filters.setdefault('contaminants', {}).update(quantifier='unitas')
    return update


@pytest.mark.parametrize('spikein, contaminants', [(False, False), (True, False), (False, True), (True, True)],
                         ids=['none', 'spikein', 'contaminants', 'both'])
def test_extra_unitas_reference_only_when_configured(tmp_path, spikein, contaminants):
    returncode, output = dry_run(tmp_path, set_unitas_refseq(spikein, contaminants), printshell=True)
    assert returncode == 0, output
    unitas = shell_lines(output, 'run_unitas.sh')
    assert len(unitas) == 4
    if not (spikein or contaminants):
        assert 'rule unitas_extra_reference:' not in output
        assert all(line.rstrip().endswith(' --') for line in unitas)
        return
    refseq = str(tmp_path / 'data/tmp/smallrna/quant/unitas/refseq.fa')
    assert all(line.endswith(f' -- -refseq {refseq}') for line in unitas)
    concat = shell_lines(output, '> data/tmp/smallrna/quant/unitas/refseq.fa')
    assert len(concat) == 1
    assert ('_spikein_unitas_formatted.fa' in concat[0]) == spikein
    assert ('_contaminants_unitas_formatted.fa' in concat[0]) == contaminants


def test_report_template_orders_sections_and_tables_for_smallrna():
    """#358: section order by the IDs MultiQC reports, cutadapt columns placed as fastp's, its comment, and the
    General Statistics hide and rename."""
    template = yaml.safe_load((ROOT / 'misc/multiqc/multiqc_config-smallrna.yaml').read_text())
    order = {key: value['order'] for key, value in template['report_section_order'].items()}
    assert order['cutadapt'] == order['fastp'] > order['fastqc_rnaseq'] > order['mirtrace'] > order['star']
    assert order['star'] > order['unitas'] > order['gene_high'] > order['exploratory']
    assert 'pca' not in order and 'umap' not in order  # the PCA's parent section is `exploratory`
    assert template['table_columns_placement']['cutadapt'] == {k: v for k, v in template['table_columns_placement']['fastp'].items()
                                                               if k in ('before_filtering_total_reads', 'pct_surviving', 'pct_adapter')}
    assert '5′ adapter artefacts' in template['section_comments']['cutadapt']
    assert 'three passes' in template['section_comments']['cutadapt'] and 'Illumina primer' in template['section_comments']['cutadapt']
    assert template['table_columns_visible']['Custom content: gcf'] == {'paired_end': False}
    assert template['table_columns_name'] == {'unitas': {'miRNA_Detected': '# miRNAs'}}
