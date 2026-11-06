"""Shared Kraken2 screening contracts on synthetic FASTQs; no classifier or database execution."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'common/scripts/k2_screen.sh'
FASTQ = b'@synthetic_read\nACGT\n+\nIIII\n'
ZERO_REPORT = '0.00\t0\t0\tU\t0\tunclassified\n'


class KrakenScreenTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='kraken-screen-')
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.fastq = self.work / "sample 'quoted'.fastq.gz"
        self.classifications = self.work / 'sample.kraken.out'
        self.report = self.work / 'sample.kraken.kreport'
        self.calls = self.work / 'calls.json'
        self.database = self.work / 'database with spaces'
        binary = self.work / 'bin'
        binary.mkdir()
        stub = binary / 'kraken2'
        stub.write_text(f'#!{sys.executable}\n' + '''
import json, os
from pathlib import Path
import sys
Path(os.environ['CALLS']).write_text(json.dumps(sys.argv[1:]))
if os.environ.get('FAIL_CLASSIFIER'):
    print('synthetic classifier failure', file=sys.stderr)
    sys.exit(23)
Path(sys.argv[sys.argv.index('--output') + 1]).write_text('C\\tsynthetic_read\\t2\\t4\\t2:1\\n')
Path(sys.argv[sys.argv.index('--report') + 1]).write_text('100.00\\t1\\t1\\tD\\t2\\tBacteria\\n')
''')
        stub.chmod(0o755)
        self.env = dict(os.environ, PATH=str(binary)+os.pathsep+os.environ['PATH'], CALLS=str(self.calls))

    def write_fastq(self, content):
        self.fastq.write_bytes(gzip.compress(content, mtime=0))

    def run_screen(self, success=True, threads='6'):
        result = subprocess.run(['bash', str(SCRIPT), '--input', str(self.fastq), '--db', str(self.database),
                                 '--output', str(self.classifications), '--report', str(self.report), '--threads', threads],
                                env=self.env, cwd=self.work, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode == 0, success, result.stdout+result.stderr)
        return result

    def assert_rejected(self):
        self.run_screen(success=False)
        self.assertFalse(self.calls.exists(), 'invalid input must not invoke Kraken2')
        self.assertFalse(self.classifications.exists())
        self.assertFalse(self.report.exists())

    def test_empty_gzip_writes_zero_results_without_running_classifier(self):
        self.write_fastq(b'')
        self.assertGreater(self.fastq.stat().st_size, 0)
        result = self.run_screen()
        self.assertFalse(self.calls.exists())
        self.assertEqual(self.classifications.read_bytes(), b'')
        self.assertEqual(self.report.read_text(), ZERO_REPORT)
        self.assertIn('zero reads', result.stderr)
        self.assertIn(str(self.fastq), result.stderr)

    def test_nonempty_gzip_keeps_arguments_and_tool_outputs(self):
        self.write_fastq(FASTQ)
        original = self.fastq.read_bytes()
        self.run_screen()
        self.assertEqual(json.loads(self.calls.read_text()), [
            '--db', str(self.database), '--output', str(self.classifications), '--report', str(self.report),
            '--threads', '6', '--gzip-compressed', '--memory-mapping', str(self.fastq)])
        self.assertEqual(self.classifications.read_text(), 'C\tsynthetic_read\t2\t4\t2:1\n')
        self.assertEqual(self.report.read_text(), '100.00\t1\t1\tD\t2\tBacteria\n')
        self.assertEqual(self.fastq.read_bytes(), original)

    def test_concatenated_gzip_with_empty_first_member_is_not_empty(self):
        self.fastq.write_bytes(gzip.compress(b'')+gzip.compress(FASTQ))
        self.run_screen()
        self.assertTrue(self.calls.exists())

    def test_multiple_empty_gzip_members_are_empty(self):
        self.fastq.write_bytes(gzip.compress(b'')+gzip.compress(b''))
        self.run_screen()
        self.assertEqual(self.report.read_text(), ZERO_REPORT)
        self.assertFalse(self.calls.exists())

    def test_empty_gzip_with_optional_header_is_empty(self):
        with self.fastq.open('wb') as raw:
            with gzip.GzipFile(filename='synthetic.fastq', mode='wb', fileobj=raw, mtime=0):
                pass
        self.run_screen()
        self.assertEqual(self.report.read_text(), ZERO_REPORT)

    def test_corrupt_empty_gzip_does_not_publish_results(self):
        self.fastq.write_bytes(gzip.compress(b'')[:-4])
        self.assert_rejected()

    def test_corrupt_later_gzip_member_is_rejected(self):
        self.fastq.write_bytes(gzip.compress(FASTQ)+gzip.compress(FASTQ)[:-4])
        self.assert_rejected()

    def test_missing_or_uncompressed_input_is_rejected(self):
        self.assert_rejected()
        self.fastq.write_bytes(FASTQ)
        self.assert_rejected()

    def test_malformed_fastq_is_rejected(self):
        for content in (b'not a FASTQ\n', b'\n', b'@read\nACGT\n+\n', b'@read\nACGT\n+\nIII\n',
                        b'@read\nACGT\nwrong\nIIII\n', b'@read\n\n+\n\n', FASTQ+b'@truncated\n'):
            with self.subTest(content=content):
                self.write_fastq(content)
                self.assert_rejected()

    def test_crlf_and_quality_characters_are_valid(self):
        self.write_fastq(b'@read comment\r\nACGT\r\n+read comment\r\n@+!I\r\n')
        self.run_screen()
        self.assertTrue(self.calls.exists())

    def test_classifier_failure_is_propagated(self):
        self.write_fastq(FASTQ)
        self.env['FAIL_CLASSIFIER'] = '1'
        result = self.run_screen(success=False)
        self.assertEqual(result.returncode, 23)
        self.assertIn('synthetic classifier failure', result.stderr)
        self.assertFalse(self.classifications.exists())
        self.assertFalse(self.report.exists())

    def test_invalid_threads_or_output_path_fail(self):
        self.write_fastq(b'')
        self.run_screen(success=False, threads='0')
        self.classifications = self.work/'absent/output'
        self.run_screen(success=False)
        self.assertFalse(self.report.exists())

    def test_mixed_and_all_empty_reports_keep_zero_counts(self):
        for contents in ((b'', FASTQ), (b'', b'')):
            reports = []
            for content in contents:
                self.write_fastq(content)
                self.run_screen()
                # KrakenTools consumes this six-column format, including U/0 rows.
                columns = self.report.read_text().strip().split('\t')
                self.assertEqual(len(columns), 6)
                self.assertEqual(int(columns[1]), int(columns[2]))
                reports.append((int(columns[1]), columns[3], int(columns[4])))
            self.assertEqual(reports[0], (0, 'U', 0))
            self.assertEqual(sum(row[0] for row in reports), int(bool(contents[1])))


class KrakenScreenRuleTests(unittest.TestCase):
    def test_default_dry_run_wires_helper_and_keeps_krona_outputs(self):
        with tempfile.TemporaryDirectory(prefix='kraken-plan-') as directory:
            work = Path(directory)
            shutil.copytree(ROOT/'default/.tests', work, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
            (work/'Snakefile').write_text('pepfile: "pep/pep_config.yaml"\nconfigfile: "config.yaml"\n'
                                        f'include: {str(ROOT/"default/default.smk")!r}\n')
            env = {k:v for k,v in os.environ.items() if not k.startswith('GCF_')}
            env.update(TMPDIR=str(work), XDG_CACHE_HOME=str(work/'.cache'))
            result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '1', 'multiqc_report'],
                                    cwd=work, env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn(str(SCRIPT), result.stdout)
        self.assertIn('kreport2krona.py', result.stdout)
        self.assertIn('ktImportText', result.stdout)
        self.assertIn('kraken.kreport', result.stdout)
        self.assertIn('kraken.out', result.stdout)


if __name__ == '__main__':
    unittest.main()
