"""Synthetic paired-read checks through both workflow CLIs and their real rules."""

import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys

from Bio.Seq import Seq
import pytest
import yaml


ROOT = Path(__file__).resolve().parents[2]
PRIMERS = {'first': ('ACGTCAGTACGATCGATGCA', 'GTACGTAGCTAGTCAGTCGA'),
           'second': ('TGCATCGTAGCTACGTCAGT', 'CATGACTGATCGTAGCATGC')}
PAYLOAD = 'GATCTAGGCTAACGTTACCGATGCTAGGATCA'


@pytest.fixture(autouse=True)
def installed_tools(monkeypatch):
    monkeypatch.setenv('PATH', str(Path(sys.executable).parent) + os.pathsep + os.environ['PATH'])


def fastq(records):
    return ''.join(f'@{name}\n{sequence}\n+\n{"I" * len(sequence)}\n' for name, sequence in records).encode()


def prepare(work, sample='sample', regions=None, empty=False):
    work.mkdir(parents=True, exist_ok=True)
    regions = PRIMERS if regions is None else regions
    kit = 'Synthetic kit; $(touch CONFIG_EXECUTED)'
    config = work / 'primer settings.yaml'
    config.write_text(yaml.safe_dump({kit: {'db': {'primers': {r: '-'.join(p) for r, p in regions.items()}}}}, sort_keys=False))
    # The mate deliberately has a different basename, exercising --in2.
    inputs = [work / f'{sample}_R1.fastq.gz', work / 'explicit mate.fastq.gz']
    for mate, path in enumerate(inputs):
        records = []
        if not empty:
            for region, primers in regions.items():
                tail = str(Seq(primers[1 - mate]).reverse_complement())
                records.append((region, primers[mate] + PAYLOAD + tail))
            records.append(('unknown', 'T' * 45))
            records.append(('too_short', next(iter(regions.values()))[mate] + 'A'))
        with gzip.open(path, 'wb') as handle:
            handle.write(fastq(records))
    return ['--in1', str(inputs[0]), '--in2', str(inputs[1]), '--libkit', kit, '--libprepconfig', str(config)]


def invoke(workflow, work, output, logdir, args, compress=False):
    command = [sys.executable, str(ROOT / workflow / 'rules/filter/scripts/demultiplex_16s_its.py'), *args,
               '--output-dir', str(output), '--log-dir', str(logdir)]
    if compress:
        command.append('--compress')
    return subprocess.run(command, cwd=work, capture_output=True, timeout=60)


@pytest.mark.parametrize('workflow', ['microbiome', 'metagenome'])
@pytest.mark.parametrize('compress', [False, True])
def test_real_cutadapt_assignments_and_literal_paths(tmp_path, workflow, compress):
    if compress and not shutil.which('pigz'):
        pytest.fail('pigz is required for the compression contract')
    sample = 'sample [x]; $(touch EXECUTED) `touch BACKTICK`'
    work = tmp_path / 'input_R1.fastq directory with spaces and \'quotes\''
    args = prepare(work, sample)
    output, logdir = work / 'output [literal]; space', work / 'separate logs $literal'
    result = invoke(workflow, work, output, logdir, args, compress)
    assert result.returncode == 0, result.stderr.decode()
    expected = fastq([('unknown:region=no_adapter', 'T' * 45),
                      ('first:region=first', PAYLOAD), ('second:region=second', PAYLOAD)])
    suffix = '.gz' if compress else ''
    for mate in (1, 2):
        path = output / f'{sample}_R{mate}.fastq{suffix}'
        actual = gzip.decompress(path.read_bytes()) if compress else path.read_bytes()
        assert actual == expected
    log = logdir / f'{sample}_qiaseq_demultiplex.log'
    assert log.read_text().count('Total read pairs processed:') == 2
    for region in PRIMERS:
        reverse_log = output / 'rev_comp_log' / f'{sample}_{region}_rev_comp.log'
        assert 'Total read pairs processed:' in reverse_log.read_text()
    assert not list(output.glob('.qiaseq-*'))
    assert not any((work / name).exists() for name in ('EXECUTED', 'BACKTICK', 'CONFIG_EXECUTED'))

    # A retry must not ingest old unmatched reads or remove another sample's
    # files. Main and reverse-complement logs describe this invocation only.
    stale = output / f'{sample}_unknown_R1.fastq'
    stale.write_bytes(b'old interrupted run')
    neighbor = output / 'sample x_unrelated_R1.fastq'
    neighbor.write_bytes(b'keep')
    result = invoke(workflow, work, output, logdir, args, compress)
    assert result.returncode == 0, result.stderr.decode()
    assert stale.read_bytes() == b'old interrupted run'
    assert neighbor.read_bytes() == b'keep'
    assert log.read_text().count('Total read pairs processed:') == 2
    for mate in (1, 2):
        path = output / f'{sample}_R{mate}.fastq{suffix}'
        assert (gzip.decompress(path.read_bytes()) if compress else path.read_bytes()) == expected


@pytest.mark.parametrize('workflow', ['microbiome', 'metagenome'])
@pytest.mark.parametrize('empty', [False, True])
def test_one_region_and_empty_pairs(tmp_path, workflow, empty):
    args = prepare(tmp_path, regions={'first': PRIMERS['first']}, empty=empty)
    output = tmp_path / 'output'
    result = invoke(workflow, tmp_path, output, output / 'log', args)
    assert result.returncode == 0, result.stderr.decode()
    expected = b'' if empty else fastq([('unknown:region=no_adapter', 'T' * 45), ('first:region=first', PAYLOAD)])
    for mate in (1, 2):
        assert (output / f'sample_R{mate}.fastq').read_bytes() == expected


@pytest.mark.parametrize('failure', ['cutadapt', 'pigz'])
def test_failure_keeps_existing_outputs_and_inputs(tmp_path, monkeypatch, failure):
    args = prepare(tmp_path)
    inputs_before = {p: p.read_bytes() for p in tmp_path.glob('*.gz')}
    output = tmp_path / 'output'
    output.mkdir()
    compress = failure == 'pigz'
    suffix = '.gz' if compress else ''
    finals = [output / f'sample_R{mate}.fastq{suffix}' for mate in (1, 2)]
    for path in finals:
        path.write_bytes(b'previous successful output')
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    executable = bindir / failure
    executable.write_text('#!/bin/sh\necho synthetic tool failure >&2\nexit 17\n')
    executable.chmod(0o755)
    monkeypatch.setenv('PATH', str(bindir) + os.pathsep + os.environ['PATH'])
    result = invoke('microbiome', tmp_path, output, output / 'log', args, compress)
    assert result.returncode != 0
    assert b'synthetic tool failure' in result.stderr
    assert all(path.read_bytes() == b'previous successful output' for path in finals)
    assert all(path.read_bytes() == before for path, before in inputs_before.items())
    assert not list(output.glob('.qiaseq-*'))


def test_configuration_cannot_execute_shell_code(tmp_path):
    args = prepare(tmp_path)
    config = Path(args[-1])
    value = yaml.safe_load(config.read_text())
    next(iter(value.values()))['db']['primers']['first'] = '$(touch EXECUTED)-ACTG'
    config.write_text(yaml.safe_dump(value))
    result = invoke('metagenome', tmp_path, tmp_path / 'output', tmp_path / 'log', args)
    assert result.returncode != 0
    assert not (tmp_path / 'EXECUTED').exists()


@pytest.mark.parametrize('region', ['../escape', 'nested/region', 'back\\slash', '{name}', 'a=b', 'unknown', ''])
def test_rejects_unsafe_region_paths_before_processing(tmp_path, region):
    args = prepare(tmp_path, regions={region: PRIMERS['first']})
    output = tmp_path / 'output'
    result = invoke('microbiome', tmp_path, output, tmp_path / 'log', args)
    assert result.returncode != 0
    assert not output.exists()


def test_first_matching_region_wins(tmp_path):
    args = prepare(tmp_path, regions={'first': PRIMERS['first'], 'second': PRIMERS['first']})
    output = tmp_path / 'output'
    result = invoke('metagenome', tmp_path, output, output / 'log', args)
    assert result.returncode == 0, result.stderr.decode()
    expected = fastq([('unknown:region=no_adapter', 'T' * 45),
                      ('first:region=first', PAYLOAD), ('second:region=first', PAYLOAD)])
    assert (output / 'sample_R1.fastq').read_bytes() == expected
    assert (output / 'sample_R2.fastq').read_bytes() == expected


@pytest.mark.parametrize('workflow', ['microbiome', 'metagenome'])
def test_real_rule_quotes_arguments(tmp_path, workflow):
    work = tmp_path / 'rule space; $(touch EXECUTED)'
    args = prepare(work, sample='sample [x]')
    output = work / 'filter'
    rule = ROOT / workflow / 'rules/filter/cutadapt.smk'
    config = {'libprepkit': 'Synthetic kit; $(touch CONFIG_EXECUTED)', 'read_geometry': [100, 100],
              'docker': {'cutadapt': 'example.invalid/cutadapt:test'}}
    # Keep the real rule's PE suffix and source resolution, with synthetic
    # configuration and inputs; no container, scientific data or network.
    primer_config = yaml.safe_load(Path(args[-1]).read_text())
    primer_config[config['libprepkit'] + ' PE'] = primer_config.pop(config['libprepkit'])
    Path(args[-1]).write_text(yaml.safe_dump(primer_config, sort_keys=False))
    snakefile = work / 'Snakefile'
    snakefile.write_text(
        'from os.path import join\n'
        f'config.update({config!r})\nFILTER_INTERIM = {str(output)!r}\nSAMPLES = ["sample [x]"]\n'
        f'def get_raw_fastq(wildcards):\n    return {{"R1": {args[1]!r}, "R2": {args[3]!r}}}\n'
        'get_filtered_fastq = get_raw_fastq\n'
        f'def src_gcf(path):\n    return {str(Path(args[-1]))!r} if path.endswith("libprep.config") else join({str(rule.parent)!r}, path)\n'
        f'include: {str(rule)!r}\n'
    )
    result = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile), '--directory', str(work),
                             '--cores', '1', '--printshellcmds', 'cutadapt_all'], capture_output=True, timeout=60)
    assert result.returncode == 0, (result.stdout + result.stderr).decode()
    assert not (work / 'EXECUTED').exists()
    assert not (work / 'CONFIG_EXECUTED').exists()
    log = output / 'cutadapt_demultiplex/log/sample [x]_qiaseq_demultiplex.log'
    assert log.read_text().count('Total read pairs processed:') == 2
