"""run_unitas.sh keeps every UNITAS run in its own work directory (#263), using a stand-in unitas.pl."""
import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'smallrna/rules/quant/scripts/run_unitas.sh'

# Behaves like unitas.pl where it matters here: needs seqmap.exe and UNITAS_refdump_<species> in the
# current directory, leaves scratch files there and writes UNITAS_<date>_<input name>_#<n>, taking the
# next free <n> when an earlier result directory is still present.
STUB = r'''#!/usr/bin/env bash
set -eu
while (( $# )); do case $1 in -input) input=$2; shift 2;; -s) species=$2; shift 2;; -threads) shift 2;; *) extra+=("$1"); shift;; esac; done
[[ -x seqmap.exe ]] || { echo "no seqmap.exe" >&2; exit 3; }
[[ -d UNITAS_refdump_$species ]] || { echo "no reference dump" >&2; exit 3; }
name=$(basename "$input"); n=1
while [[ -e UNITAS_04-11-2026_${name}_#$n ]]; do n=$((n + 1)); done
for copy in $(seq 1 "${STUB_RESULTS:-1}"); do
    dir=UNITAS_04-11-2026_${name}_#$((n + copy - 1)); mkdir "$dir"
    printf '%s\n' "$(cat "$input")" > "$dir/results.html"
    printf '%s\n' "${extra[@]:-}" > "$dir/args.txt"
done
sleep "${STUB_SLEEP:-0}"
touch "$name.fas.collapsed"
'''


@pytest.fixture
def env(tmp_path):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'unitas.pl').write_text(STUB)
    (bin_dir / 'unitas.pl').chmod(0o755)
    (tmp_path / 'ext/UNITAS_refdump_mus_musculus').mkdir(parents=True)
    (tmp_path / 'ext/seqmap.exe').write_text('seqmap\n')
    return dict(os.environ, PATH=f'{bin_dir}:{os.environ["PATH"]}')


def start(tmp_path, env, sample, extra=(), **settings):
    fastq = tmp_path / 'filter' / f'{sample}_R1.fastq'
    fastq.parent.mkdir(exist_ok=True)
    fastq.write_text(f'reads of {sample}\n')
    command = ['bash', str(SCRIPT), '--input', str(fastq.relative_to(tmp_path)), '--species', 'mus_musculus',
               '--threads', '2', '--seqmap', 'ext/seqmap.exe', '--refdump', 'ext/UNITAS_refdump_mus_musculus',
               '--workdir', f'quant/unitas/work/{sample}', '--outdir', f'quant/unitas/unitas/{sample}', '--', *extra]
    return subprocess.Popen(command, cwd=tmp_path, env=dict(env, **settings), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)


def run(tmp_path, env, sample, extra=(), **settings):
    process = start(tmp_path, env, sample, extra, **settings)
    output, _ = process.communicate(timeout=60)
    return process.returncode, output


def test_results_move_to_outdir_and_scratch_stays_in_workdir(tmp_path, env):
    returncode, output = run(tmp_path, env, 'S1')
    assert returncode == 0, output
    assert (tmp_path / 'quant/unitas/unitas/S1/results.html').read_text() == 'reads of S1\n'
    assert not (tmp_path / 'quant/unitas/work/S1').exists()
    assert not list(tmp_path.glob('UNITAS_*')) and not list(tmp_path.glob('*.fas.collapsed'))
    assert not (tmp_path / 'seqmap.exe').exists()


def test_concurrent_samples_with_overlapping_names_keep_their_own_results(tmp_path, env):
    samples = ['S1', 'XS1', 'S1X']
    processes = {sample: start(tmp_path, env, sample, STUB_SLEEP='1') for sample in samples}
    for sample, process in processes.items():
        output, _ = process.communicate(timeout=60)
        assert process.returncode == 0, output
    for sample in samples:
        assert (tmp_path / f'quant/unitas/unitas/{sample}/results.html').read_text() == f'reads of {sample}\n'


def test_leftover_from_interrupted_run_is_not_delivered(tmp_path, env):
    stale = tmp_path / 'quant/unitas/work/S1/UNITAS_03-11-2026_S1_R1.fastq_#1'
    stale.mkdir(parents=True)
    (stale / 'results.html').write_text('stale\n')
    returncode, output = run(tmp_path, env, 'S1')
    assert returncode == 0, output
    assert (tmp_path / 'quant/unitas/unitas/S1/results.html').read_text() == 'reads of S1\n'


def test_more_than_one_result_directory_fails_without_touching_outdir(tmp_path, env):
    previous = tmp_path / 'quant/unitas/unitas/S1'
    previous.mkdir(parents=True)
    (previous / 'results.html').write_text('previous\n')
    returncode, output = run(tmp_path, env, 'S1', STUB_RESULTS='2')
    assert returncode != 0
    assert 'expected one UNITAS result directory' in output and 'found 2' in output
    assert (previous / 'results.html').read_text() == 'previous\n'


def test_extra_arguments_reach_unitas(tmp_path, env):
    returncode, output = run(tmp_path, env, 'S1', extra=('--riborase', '-refseq', '/abs/refseq.fa'))
    assert returncode == 0, output
    assert (tmp_path / 'quant/unitas/unitas/S1/args.txt').read_text().split() == ['--riborase', '-refseq', '/abs/refseq.fa']


def test_reference_dump_keeps_its_name_when_it_is_a_symlink(tmp_path, env):
    target = tmp_path / 'ext/release/dump-content'
    target.mkdir(parents=True)
    (tmp_path / 'ext/UNITAS_refdump_mus_musculus').rmdir()
    (tmp_path / 'ext/UNITAS_refdump_mus_musculus').symlink_to(target)
    returncode, output = run(tmp_path, env, 'S1')
    assert returncode == 0, output
