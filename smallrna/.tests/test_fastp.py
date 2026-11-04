"""Synthetic small-RNA stage contracts and read outcomes with fastp 0.20.1."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'smallrna/rules/filter/scripts'
KITS = yaml.safe_load((ROOT / 'libprep.config').read_text())
SMALLRNA = {name: kit for name, kit in KITS.items() if kit['workflow'] == 'smallrna'}
ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'
INSERT = 'ACGTACGTACGTACGTA'


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


policy = load('fastp_config')
reports = load('merge_fastp')


@pytest.mark.parametrize('kit', SMALLRNA.values(), ids=SMALLRNA.keys())
def test_presets_have_one_explicit_final_policy(kit):
    settings = kit['filter']['trim']['fastp']
    assert set(settings) == {'args'}
    args = shlex.split(policy.final_fastp_args(settings))
    assert args[args.index('-l') + 1] == ('16' if 'v3' in kit['libprep_name'] else '15')
    assert all(flag in args for flag in ('-q', '-u', '-n', '-e', '--low_complexity_filter', '--complexity_threshold'))
    assert '-Q' not in args and '-L' not in args
    assert policy.ADAPTER_ONLY_ARGS == '-Q -L -G'


@pytest.mark.parametrize('settings', [{}, {'args': ''}, {'args': {}}, {'args': '-A'},
    {'args': '-l 16 -L'}, {'args': '-l 16', 'trim_args': ''}, {'args': '-l 16', 'params': '-l 0'}])
def test_invalid_or_legacy_overrides_fail_with_migration_error(settings):
    with pytest.raises(ValueError, match='Small-RNA|filter.trim.fastp.args'):
        policy.final_fastp_args(settings)


@pytest.mark.parametrize('args', ['-l 16', '--length_required 16', '-L', '--disable_length_filtering', '-l 0'])
def test_explicit_final_length_overrides(args):
    assert policy.final_fastp_args({'args': args, 'params': ''}) == args


@pytest.mark.parametrize('name', SMALLRNA)
def test_resolved_workflow_commands(tmp_path, name):
    shutil.copytree(ROOT / 'smallrna/.tests', tmp_path, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
    settings = yaml.safe_load((tmp_path / 'config.yaml').read_text())
    settings['libprepkit'] = name
    (tmp_path / 'config.yaml').write_text(yaml.safe_dump(settings))
    (tmp_path / 'Snakefile').write_text('configfile: "config.yaml"\npepfile: "pep/pep_config.yaml"\n'
                                      f'include: {str(ROOT / "smallrna/smallrna.smk")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    env.update(TMPDIR=str(tmp_path), XDG_CACHE_HOME=str(tmp_path / '.cache'))
    result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '4', 'multiqc_report'],
                            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    commands = [shlex.split(line.split('fastp ', 1)[1]) for line in output.splitlines() if 'Shell command: fastp ' in line]
    assert len(commands) == 2 * len(settings['samples'])
    paths = []
    for args in commands:
        paths.extend([args[args.index('--json') + 1], args[args.index('--html') + 1]])
        if '--adapter_sequence' in args:
            assert all(flag in args for flag in ('-Q', '-L', '-G'))
            assert '-f' not in args and '-l' not in args and '-A' not in args
        else:
            assert '-A' in args and '-Q' not in args and '-L' not in args
            assert args[args.index('-l') + 1] == ('16' if 'v3' in name else '15')
    assert len(paths) == len(set(paths))
    assert 'provenance/fastp/' in output


@pytest.fixture
def fastp():
    binary = os.environ.get('FASTP_0201')
    if not binary:
        pytest.skip('Set FASTP_0201 to the pinned fastp 0.20.1 binary for synthetic execution')
    version = subprocess.run([binary, '--version'], capture_output=True, text=True, check=True)
    assert 'fastp 0.20.1' in version.stdout + version.stderr
    return binary


def write_reads(path, records, nextseq=False):
    text = ''
    for name, seq, quality in records:
        prefix = 'NS500000:1:SYNTHETIC:1:1:1:' if nextseq else 'SYNTHETIC:'
        text += f'@{prefix}{name}\n{seq}\n+\n{quality or "I" * len(seq)}\n'
    path.write_text(text)


def run_fastp(binary, work, name, source, args):
    out = work / (name + '.fq')
    report = work / (name + '.json')
    command = [binary, '-i', str(source), '-o', str(out), '-j', str(report), '-h', str(work / (name + '.html')),
               '-w', '1', *shlex.split(args)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    lines = out.read_text().splitlines()
    reads = {lines[i].split(':')[-1]: (lines[i+1], lines[i+3]) for i in range(0, len(lines), 4)}
    return out, reads, json.loads(report.read_text())


def run_stages(binary, work, records, args, nextseq=False):
    raw = work / 'input.fq'
    write_reads(raw, records, nextseq=nextseq)
    trimmed, intermediate, first = run_fastp(binary, work, 'adapter', raw,
                                            policy.ADAPTER_ONLY_ARGS + ' --adapter_sequence ' + ADAPTER)
    _, retained, final = run_fastp(binary, work, 'filter', trimmed, args + ' -A')
    return intermediate, retained, first, final


@pytest.mark.parametrize('name', SMALLRNA)
def test_read_outcomes_and_reconciled_reports(fastp, tmp_path, name):
    v3 = 'v3' in name
    def raw(insert):
        return ('TGCA' + insert + 'CATG' if v3 else insert) + ADAPTER
    records = [(str(n), raw(INSERT[:n]), None) for n in (14, 15, 16, 17)]
    records += [('empty', ADAPTER, None), ('lowcomplexity', raw('A' * 20), None),
                ('N6', raw('ACGTNNNNNNACGTACGT'), None), ('N5', raw('ACGTNNNNNACGTACGT'), None),
                ('lowquality', raw(INSERT), '!' * len(raw(INSERT))),
                ('short', raw('ACGT'), None)]
    if v3:
        # N bases and poor quality in the random ends must not discard a valid insert in pass one.
        records += [('tags', 'NNNN' + INSERT + 'NNNN' + ADAPTER, '!!!!' + 'I' * 17 + '!!!!' + 'I' * len(ADAPTER))]
    intermediate, retained, first, final = run_stages(fastp, tmp_path, records, SMALLRNA[name]['filter']['trim']['fastp']['args'])
    expected = {'16', '17', 'N5'} | ({'tags'} if v3 else {'15'})
    assert set(retained) == expected
    for n in ('16', '17'):
        assert retained[n] == (INSERT[:int(n)], 'I' * int(n))
    if v3:
        assert retained['tags'] == (INSERT, 'I' * 17)
    assert set(intermediate) == {r[0] for r in records} - {'empty'}
    assert first['filtering_result']['too_short_reads'] == 1
    assert final['filtering_result']['low_quality_reads'] == 1
    assert final['filtering_result']['too_many_N_reads'] == 1
    assert final['filtering_result']['low_complexity_reads'] == 1
    merged = reports.merge_reports(first, final)
    assert merged['summary']['before_filtering']['total_reads'] == len(records)
    assert merged['summary']['after_filtering']['total_reads'] == len(expected)
    assert sum(merged['filtering_result'].values()) == len(records)
    assert merged['command'] == first['command'].strip() + '\n' + final['command'].strip()
    assert [s['filtering_result'] for s in merged['processing_stages']] == [first['filtering_result'], final['filtering_result']]
    cli = subprocess.run([sys.executable, str(SCRIPTS / 'merge_fastp.py'), str(tmp_path / 'adapter.json'),
                          str(tmp_path / 'filter.json'), '--html', str(tmp_path / 'combined.html')],
                         capture_output=True, text=True, check=True)
    assert json.loads(cli.stdout) == merged
    assert 'Original reads:' in (tmp_path / 'combined.html').read_text()
    for field in ('total_reads', 'total_bases'):
        wrong = copy.deepcopy(final)
        wrong['summary']['before_filtering'][field] += 1
        with pytest.raises(ValueError):
            reports.merge_reports(first, wrong)


@pytest.mark.parametrize('length', ['-L', '-l 0'])
def test_disabled_length_retains_nonempty_reads_only(fastp, tmp_path, length):
    records = [('empty', ADAPTER, None), ('one', 'C' + ADAPTER, None), ('short', 'CGTA' + ADAPTER, None),
               ('Ns', 'N' * 20 + ADAPTER, None), ('quality', INSERT + ADAPTER, '!' * (17 + len(ADAPTER)))]
    _, retained, first, final = run_stages(fastp, tmp_path, records, '-Q -G ' + length)
    assert set(retained) == {'one', 'short', 'Ns', 'quality'}
    assert first['filtering_result']['too_short_reads'] == 1
    assert final['filtering_result']['too_short_reads'] == 0


@pytest.mark.parametrize('name', [name for name in SMALLRNA if 'v3' not in name])
def test_other_kits_match_previous_route_and_one_pass_changes_polyg(fastp, tmp_path, name):
    records = [('polyg', 'ACGTACGTACGTACGT' + 'G' * 12 + ADAPTER, None),
               ('normal', INSERT + ADAPTER, None), ('short', INSERT[:14] + ADAPTER, None), ('empty', ADAPTER, None)]
    args = SMALLRNA[name]['filter']['trim']['fastp']['args']
    _, retained, first, final = run_stages(fastp, tmp_path, records, args, nextseq=True)
    old_trim, _, old_first = run_fastp(fastp, tmp_path, 'old_adapter', tmp_path / 'input.fq',
                                      '-Q -n 0 -l 0 -G --adapter_sequence ' + ADAPTER)
    _, old_retained, old_final = run_fastp(fastp, tmp_path, 'old_filter', old_trim, '--low_complexity_filter -A')
    assert retained == old_retained
    assert first['filtering_result'] == old_first['filtering_result']
    assert final['filtering_result'] == old_final['filtering_result']
    _, one_pass, _ = run_fastp(fastp, tmp_path, 'one_pass', tmp_path / 'input.fq', args + ' --adapter_sequence ' + ADAPTER)
    assert 'polyg' in one_pass and 'polyg' not in retained


def test_report_accounting_without_tool_execution():
    def report(before, after, command, **rejected):
        return {'summary': {'before_filtering': {'total_reads': before, 'total_bases': before * 20},
                            'after_filtering': {'total_reads': after, 'total_bases': after * 20}},
                'filtering_result': dict(passed_filter_reads=after, **rejected),
                'adapter_cutting': {'adapter_trimmed_reads': before}, 'read1_before_filtering': {'total_reads': before},
                'duplication': {'rate': 0.25}, 'command': command}
    first = report(10, 8, 'fastp adapter', too_short_reads=2)
    final = report(8, 5, 'fastp filter', too_short_reads=1, low_quality_reads=2)
    final['duplication']['rate'] = 0.5
    original = copy.deepcopy((first, final))
    combined = reports.merge_reports(first, final)
    assert combined['filtering_result'] == {'passed_filter_reads': 5, 'too_short_reads': 3, 'low_quality_reads': 2}
    assert combined['duplication']['rate'] == 0.25
    assert (first, final) == original
    assert len(combined['processing_stages']) == 2
    final['filtering_result']['passed_filter_reads'] = 6
    with pytest.raises(ValueError, match='counts do not reconcile'):
        reports.merge_reports(first, final)
