"""mirge3_collect.py (#404): miRge3's run moved into the project's miRge3 folder, also over an earlier run's entries,
with miRge3's own miR.Counts.csv never at the path mirge3_counts owns."""
from pathlib import Path
import subprocess
import sys

from test_quant_methods import dry_run, rule_block

SCRIPT = Path(__file__).resolve().parents[1] / 'rules/quant/scripts/mirge3_collect.py'
RESERVED = ['fastq', 'lib', 'mirbase', 'spikein.tsv', 'reference.tsv', 'miR.Counts.mirge3.csv', 'unmapped.fasta']


def mirge3_run(out, stamp='2026-01-01_00-00-00', tag='new', extra=()):
    """A miRge3 run folder as miRge3 writes it, files holding tag."""
    run = out / 'mirbase' / 'miRge.{}'.format(stamp)
    (run / 'tRFs.samples.tmp').mkdir(parents=True)
    (run / 'tRFs.samples.tmp' / 'S1.aligned_tRFs.report').write_text(tag)
    for name in ('miR.Counts.csv', 'mapped.csv', 'isomirs.csv', 'run.log') + tuple(extra):
        (run / name).write_text(tag)
    return run


def project(tmp_path):
    out = tmp_path / 'quant' / 'mirge3'
    (out / 'lib' / 'human').mkdir(parents=True)
    (out / 'reference.tsv').write_text('reference')
    return out


def collect(out):
    return subprocess.run([sys.executable, str(SCRIPT), '--run-root', str(out / 'mirbase'), '--output', str(out),
                           '--mirtable', str(out / 'miR.Counts.mirge3.csv'), '--reserved', *RESERVED], capture_output=True, text=True)


def test_run_is_moved_and_mirge3_table_kept_apart(tmp_path):
    out = project(tmp_path)
    mirge3_run(out)
    result = collect(out)
    assert result.returncode == 0, result.stderr
    assert (out / 'miR.Counts.mirge3.csv').read_text() == 'new'
    assert not (out / 'miR.Counts.csv').exists()  # mirge3_counts' path
    assert (out / 'mapped.csv').read_text() == 'new' and (out / 'tRFs.samples.tmp' / 'S1.aligned_tRFs.report').read_text() == 'new'
    assert not (out / 'mirbase').exists()
    assert (out / 'reference.tsv').read_text() == 'reference' and (out / 'lib' / 'human').is_dir()


def test_rerun_replaces_an_earlier_runs_entries(tmp_path):
    """The earlier run's tRFs.samples.tmp (undeclared, not removed by Snakemake) used to stop mv; the workflow's
    miR.Counts.csv stays as it is until mirge3_counts remakes it."""
    out = project(tmp_path)
    mirge3_run(out, tag='old')
    assert collect(out).returncode == 0
    (out / 'tRFs.samples.tmp' / 'S9.only_in_the_old_run').write_text('old')
    (out / 'miR.Counts.csv').write_text('workflow table')
    mirge3_run(out, stamp='2026-01-02_00-00-00')
    result = collect(out)
    assert result.returncode == 0, result.stderr
    assert sorted(p.name for p in (out / 'tRFs.samples.tmp').iterdir()) == ['S1.aligned_tRFs.report']
    assert (out / 'tRFs.samples.tmp' / 'S1.aligned_tRFs.report').read_text() == 'new'
    assert (out / 'mapped.csv').read_text() == 'new' and (out / 'miR.Counts.mirge3.csv').read_text() == 'new'
    assert (out / 'miR.Counts.csv').read_text() == 'workflow table'


def test_a_file_replaces_a_directory_of_the_same_name(tmp_path):
    out = project(tmp_path)
    (out / 'isomirs.csv').mkdir()
    mirge3_run(out)
    assert collect(out).returncode == 0 and (out / 'isomirs.csv').read_text() == 'new'


def check_stops_before_moving(out, message):
    result = collect(out)
    assert result.returncode != 0 and message in result.stderr and 'Traceback' not in result.stderr
    assert not (out / 'mapped.csv').exists() and not (out / 'miR.Counts.mirge3.csv').exists()


def test_two_runs_or_none_stop(tmp_path):
    out = project(tmp_path)
    (out / 'mirbase').mkdir()
    check_stops_before_moving(out, 'expected one miRge3 run in')
    mirge3_run(out)
    mirge3_run(out, stamp='2026-01-02_00-00-00')
    check_stops_before_moving(out, "found ['miRge.2026-01-01_00-00-00', 'miRge.2026-01-02_00-00-00']")


def test_a_run_without_mirge3_table_stops(tmp_path):
    out = project(tmp_path)
    (mirge3_run(out) / 'miR.Counts.csv').unlink()
    check_stops_before_moving(out, 'miRge3 wrote no miR.Counts.csv')


def test_a_reserved_name_in_the_run_stops(tmp_path):
    out = project(tmp_path)
    mirge3_run(out, extra=('reference.tsv',))
    check_stops_before_moving(out, 'miRge3 wrote reference.tsv in')
    assert (out / 'reference.tsv').read_text() == 'reference'


def test_rule_clears_a_failed_run_and_collects(tmp_path):
    returncode, output = dry_run(tmp_path, 'mirge3', printshell=True)
    assert returncode == 0, output
    quant = ' '.join(rule_block(output, 'mirge3_quant').split())
    assert 'rm -rf data/tmp/smallrna/quant/mirge3/fastq data/tmp/smallrna/quant/mirge3/mirbase && mkdir -p' in quant
    assert ('mirge3_collect.py --run-root data/tmp/smallrna/quant/mirge3/mirbase --output data/tmp/smallrna/quant/mirge3 '
            '--mirtable data/tmp/smallrna/quant/mirge3/miR.Counts.mirge3.csv --reserved ' + ' '.join(RESERVED)) in quant
    assert 'mv ' not in quant.split('miRge3.0', 1)[1]
