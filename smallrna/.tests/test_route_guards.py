"""Unsupported small-RNA routes stop at workflow load instead of failing mid-run."""
import os
import shlex
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

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


@pytest.mark.parametrize('kit, protocol', [('Bioo Scientific NEXTflex Small RNA-Seq Kit v3', 'nextflex'),
                                           ('Bioo Scientific NEXTflex Small RNA-Seq Kit v4', 'illumina'),
                                           ('TruSeq Small RNA Library Prep Kit', 'illumina')])
def test_mirtrace_qc_is_delivered_and_reported(tmp_path, kit, protocol):
    returncode, output = dry_run(tmp_path, lambda s: s.update(libprepkit=kit), printshell=True)
    assert returncode == 0, output
    mirtrace = shell_lines(output, 'mirtrace qc ')
    assert len(mirtrace) == 1 and f'-p {protocol} ' in mirtrace[0] and '--force' in mirtrace[0]
    assert 'rule bfq_level2_mirtrace:' in output
    assert 'bfq/logs/mirtrace/mirtrace-stats-mirna-complexity.tsv' in output
    assert any('-m mirtrace' in line for line in shell_lines(output, 'multiqc '))


def test_mirtrace_filter_route_is_rejected(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: s.setdefault('filter', {}).setdefault('trim', {}).update(quantifier='mirtrace'))
    assert returncode != 0
    assert "filter.trim.quantifier: smallrna supports fastp, cutadapt or skip, got 'mirtrace'" in output


def test_missing_mirtrace_protocol_fails(tmp_path):
    update = lambda s: s.setdefault('filter', {}).setdefault('trim', {}).setdefault('mirtrace', {}).update(protocol=None)
    returncode, output = dry_run(tmp_path, update)
    assert returncode != 0
    assert 'filter.trim.mirtrace.protocol: set the miRTrace protocol' in output


def test_mirtrace_config_lists_merged_fastq_per_sample(tmp_path):
    returncode, output = dry_run(tmp_path, lambda s: None, target='data/tmp/smallrna/qc/mirtrace/mirtrace.config', dry=False)
    assert returncode == 0, output
    rows = [line.split(',') for line in (tmp_path / 'data/tmp/smallrna/qc/mirtrace/mirtrace.config').read_text().splitlines()]
    assert [row[1] for row in rows] == ['SR01', 'SR02', 'SR03', 'SR04']
    assert all(row[0].endswith(f'/filter/fastq/{row[1]}_R1.fastq.gz') and row[2] == 'TGGAATTCTCGGGTGCCAAGG' for row in rows)


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
