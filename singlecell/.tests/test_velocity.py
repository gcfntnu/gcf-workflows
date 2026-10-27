"""Velocity defaults and actual single-cell DAGs; no scientific commands are run."""
from copy import deepcopy
from pathlib import Path
import json
import os
import re
import shlex
import subprocess
import sys

import pytest
import yaml
from configmaker.libprep import LibprepConfig


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / 'singlecell/singlecell.smk'
CONVERTER = ROOT / 'singlecell/rules/quant/scripts/convert_scanpy.py'
LIBPREP = ROOT / 'libprep.config'
PARSE_KIT = 'Parse Biosciences Evercode WT v4'
TENX_KIT = '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3'
OMITTED = object()
SPLITPIPE_WT_KITS = [name.removesuffix(' PE') for name, config in yaml.safe_load(LIBPREP.read_text()).items()
    if name.startswith('Parse Biosciences') and config.get('quant', {}).get('method') == 'splitpipe'
    and config.get('quant', {}).get('assay', 'wt') == 'wt']


def protected_env(tmp_path):
    """Keep subprocess caches local and prevent network/SMTP access during dry runs."""
    support = tmp_path / 'support'
    support.mkdir(exist_ok=True)
    (support / 'sitecustomize.py').write_text(
        'import sys\n'
        'def deny_network(event, args):\n'
        '    if event in {"smtplib.connect", "socket.connect"}:\n'
        '        raise RuntimeError("Network and SMTP are forbidden in velocity tests")\n'
        'sys.addaudithook(deny_network)\n')
    env = os.environ.copy()
    env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
    for name in ['GCF_SECRET', 'GCF_EXT', 'GCF_INTERIM', 'GCF_FASTQ', 'GCF_REF']:
        env.pop(name, None)
    env['PYTHONPATH'] = str(support) + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'), NUMBA_CACHE_DIR=str(tmp_path / '.numba'))
    return env


def workflow_fixture(tmp_path, methods=None, use_velo=OMITTED, kit=PARSE_KIT):
    """Use the real libprep selector and workflow default merge for synthetic inputs."""
    parse = kit.startswith('Parse Biosciences')
    geometry = [100, 100] if parse else [28, 90]
    selection = LibprepConfig.load(LIBPREP).select(kit, geometry)
    source = deepcopy(selection.parameters)
    source.update(workflow='singlecell', organism='mus_musculus', project_id=['SYNTHETIC'], libprepkit=kit,
        machine='Synthetic sequencer', read_geometry=geometry, ext_dir='data/ext', interim_dir='data/tmp',
        fastq_dir='data/raw/fastq', preprocessing={'enabled': False},
        samples={sample: {'R1': f'{sample}_R1.fastq.gz', 'R2': f'{sample}_R2.fastq.gz'} for sample in ['libA', 'libB']})
    if methods is not None:
        source['quant']['method'] = methods
    if use_velo is not OMITTED:
        source['quant']['use_velo'] = use_velo
    if parse:
        source['wells'] = {'bioA': {'Wells': 'A1'}, 'bioB': {'Wells': 'B1'}}
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(source))
    (tmp_path / 'pep.yaml').write_text('pep_version: 2.0.0\nsample_table: samples.csv\n')
    (tmp_path / 'samples.csv').write_text('sample_name,R1,R2,paired_end\n'
        'libA,libA_R1.fastq.gz,libA_R2.fastq.gz,1\nlibB,libB_R1.fastq.gz,libB_R2.fastq.gz,1\n')
    fastqs = tmp_path / 'data/raw/fastq'
    fastqs.mkdir(parents=True)
    for sample in source['samples']:
        for read in ['R1', 'R2']:
            (fastqs / f'{sample}_{read}.fastq.gz').touch()

    # Inspect actual rule declarations/parameters after the full workflow merge.
    inspection = '''
from types import SimpleNamespace
import json
velocity_contract = {'attachment': {}, 'star_args': {}, 'producer_outputs': {}}
for method in QUANT_AVAILABLE_METHODS:
    wc = SimpleNamespace(method=method, aggr_id='all_samples')
    velocity_contract['attachment'][method] = scanpy_aggr_velo(wc)
if '10x_starsolo' in METHODS:
    velocity_contract['star_args']['10x_starsolo'] = STARSOLO_10X_ARGS
    velocity_contract['producer_outputs']['10x_starsolo'] = list(rules.starsolo_quant.output)
if any(method in PARSEBIO_STARSOLO_MODES for method in METHODS):
    velocity_contract['star_args']['parsebio_starsolo'] = STARSOLO_PARSEBIO_ARGS
    velocity_contract['producer_outputs']['parsebio_starsolo'] = list(rules.parsebio_starsolo_quant.output)
    velocity_contract['parse_sample_outputs'] = list(rules.parsebio_starsolo_sample_raw.output)
Path('velocity-contract.json').write_text(json.dumps(velocity_contract))
'''
    (tmp_path / 'Snakefile').write_text('from pathlib import Path\npepfile: "pep.yaml"\nconfigfile: "config.yaml"\n'
        + f'include: {str(WORKFLOW)!r}\n' + inspection)
    return source


def dry_run(tmp_path, methods=None, use_velo=OMITTED, kit=PARSE_KIT, target='expression', expect_success=True):
    source = workflow_fixture(tmp_path, methods, use_velo, kit)
    selected = [method.strip() for method in source['quant']['method'].split(',')]
    if target == 'expression':
        targets = [f'data/tmp/singlecell/quant/aggregate/{method}/scanpy/all_samples_filtered.h5ad' for method in selected]
    else:
        targets = [target]
    command = [str(Path(sys.executable).with_name('snakemake')), '--dry-run', '--printshellcmds', '--cores', '1', *targets]
    result = subprocess.run(command, cwd=tmp_path, env=protected_env(tmp_path), capture_output=True, text=True)
    transcript = result.stdout + result.stderr
    (tmp_path / 'dry-run.log').write_text(transcript)
    if expect_success:
        assert result.returncode == 0, transcript
    else:
        assert result.returncode != 0, transcript
    return transcript


def conversion_commands(transcript, method):
    commands = []
    for line in transcript.splitlines():
        if 'convert_scanpy.py ' not in line:
            continue
        args = shlex.split(line)
        if '-f' in args and args[args.index('-f') + 1] == method and '-F' in args and args[args.index('-F') + 1] == 'anndata':
            commands.append(args)
    return commands


def assert_velocity_contract(tmp_path, transcript, expected):
    observed = json.loads((tmp_path / 'velocity-contract.json').read_text())
    for method, enabled in expected.items():
        assert observed['attachment'][method] == ('--use-velo' if enabled else '')
        command_method = 'parsebio_starsolo' if method == 'parsebio_starsolo_rt' else method
        commands = conversion_commands(transcript, command_method)
        assert commands, transcript
        assert all(('--use-velo' in command) == enabled for command in commands)
        if method in {'10x_starsolo', 'parsebio_starsolo', 'parsebio_starsolo_rt'}:
            family = command_method
            args = shlex.split(observed['star_args'][family])
            feature_start = args.index('--soloFeatures') + 1
            features = []
            for arg in args[feature_start:]:
                if arg.startswith('--'):
                    break
                features.append(arg)
            assert ('Velocyto' in features) == enabled
            outputs = observed['producer_outputs'][family]
            for matrix in ['spliced.mtx', 'unspliced.mtx', 'ambiguous.mtx']:
                assert any('/Velocyto/' in output and output.endswith(matrix) for output in outputs) == enabled
            if family == 'parsebio_starsolo':
                for matrix in ['spliced.mtx', 'unspliced.mtx', 'ambiguous.mtx']:
                    assert any('/Velocyto/' in output and output.endswith(matrix)
                        for output in observed['parse_sample_outputs']) == enabled
    if 'splitpipe' in expected:
        assert 'rule splitpipe_splice_from_tscp:' not in transcript
        assert 'splitpipe_velo.py' not in transcript
    if 'parsebio_starsolo' in expected:
        # The automatically available RT representation must not regain a disabled flag.
        assert observed['attachment']['parsebio_starsolo_rt'] == observed['attachment']['parsebio_starsolo']


@pytest.mark.parametrize('kit', SPLITPIPE_WT_KITS)
def test_splitpipe_wt_presets_default_to_expression_only(tmp_path, kit):
    transcript = dry_run(tmp_path, kit=kit, target='multiqc_report')
    assert_velocity_contract(tmp_path, transcript, {'splitpipe': False})
    for rule in ['splitpipe_quant', 'splitpipe_aggr', 'scanpy_aggr_filtered', 'multiqc_report']:
        assert re.search(rf'^(?:local)?rule {rule}:$', transcript, re.MULTILINE), transcript


@pytest.mark.parametrize('value', [False, {}, {'splitpipe': False}])
def test_splitpipe_explicit_disabled_configuration(tmp_path, value):
    transcript = dry_run(tmp_path, use_velo=value)
    assert_velocity_contract(tmp_path, transcript, {'splitpipe': False})


@pytest.mark.parametrize('method,kit', [('10x_starsolo', TENX_KIT), ('parsebio_starsolo', PARSE_KIT)])
@pytest.mark.parametrize('setting', ['omitted', 'true', 'false', 'map-false'])
def test_starsolo_default_and_explicit_velocity(tmp_path, method, kit, setting):
    value = {'omitted': OMITTED, 'true': True, 'false': False, 'map-false': {method: False}}[setting]
    transcript = dry_run(tmp_path, methods=method, use_velo=value, kit=kit)
    assert_velocity_contract(tmp_path, transcript, {method: setting in {'omitted', 'true'}})


@pytest.mark.parametrize('methods', ['splitpipe,parsebio_starsolo', 'parsebio_starsolo,splitpipe'])
@pytest.mark.parametrize('setting', ['omitted', 'explicit-map', 'partial-map', 'false', 'alias-false'])
def test_mixed_parse_quantifiers_resolve_velocity_independently(tmp_path, methods, setting):
    value = {'omitted': OMITTED, 'explicit-map': {'splitpipe': False, 'parsebio_starsolo': True},
        'partial-map': {'splitpipe': False}, 'false': False, 'alias-false': {'parsebio_starsolo_rt': False}}[setting]
    transcript = dry_run(tmp_path, methods=methods, use_velo=value)
    assert_velocity_contract(tmp_path, transcript,
        {'splitpipe': False, 'parsebio_starsolo': setting in {'omitted', 'explicit-map', 'partial-map'}})


@pytest.mark.parametrize('methods,value,error_terms', [
    ('splitpipe', True, ['Split-Pipe', 'parsebio_starsolo']),
    ('splitpipe', {'splitpipe': True}, ['Split-Pipe', 'parsebio_starsolo']),
    ('splitpipe,parsebio_starsolo', True, ['Split-Pipe', 'parsebio_starsolo']),
    ('splitpipe', 'false', ['use_velo', 'bool']),
    ('splitpipe', 0, ['use_velo', 'bool']),
    ('splitpipe', None, ['use_velo', 'bool']),
    ('splitpipe', {'splitpipe': 'false'}, ['use_velo', 'bool']),
    ('splitpipe', {'splitpipe': 0}, ['use_velo', 'bool']),
    ('splitpipe', {'typo': False}, ['use_velo', 'typo']),
    ('splitpipe', {'parsebio_starsolo': True}, ['use_velo', 'parsebio_starsolo']),
    ('parsebio_starsolo', {'parsebio_starsolo': True, 'parsebio_starsolo_rt': False}, ['use_velo', 'parsebio_starsolo']),
])
def test_invalid_velocity_requests_fail_before_workflow_dag(tmp_path, methods, value, error_terms):
    transcript = dry_run(tmp_path, methods=methods, use_velo=value, expect_success=False)
    for term in error_terms:
        assert term.lower() in transcript.lower(), transcript
    assert 'Building DAG of jobs' not in transcript
    assert not (tmp_path / 'velocity-contract.json').exists()
    assert not (tmp_path / 'data/tmp/singlecell/quant').exists()


@pytest.mark.parametrize('input_format', ['splitpipe', 'splitpipe_aggr'])
def test_direct_splitpipe_velocity_request_fails_before_loading_input(tmp_path, input_format):
    output = tmp_path / 'result.h5ad'
    result = subprocess.run([sys.executable, str(CONVERTER), str(tmp_path / 'missing.mtx'),
        '-f', input_format, '--use-velo', '-F', 'anndata', '-o', str(output)],
        cwd=tmp_path, env=protected_env(tmp_path), capture_output=True, text=True)
    transcript = result.stdout + result.stderr
    (tmp_path / 'converter.log').write_text(transcript)
    assert result.returncode != 0
    assert 'Split-Pipe' in transcript and 'parsebio_starsolo' in transcript, transcript
    assert 'No such file or directory' not in transcript
    assert not output.exists()
