"""Exercise the screening subsampler with real gzip and a tiny stand-in for seqkit."""
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'common/scripts/k2_subsample.sh'
READ = b'@example\nACGT\n+\nIIII\n'


class SubsampleTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.input = self.work / 'filtered reads.fastq.gz'
        self.output = self.work / 'screening reads.fastq.gz'
        seqkit = self.work / 'seqkit'
        seqkit.write_text(f'#!{sys.executable}\n' + """
import gzip, os, shutil, sys
args = sys.argv[1:]
failure = os.environ.get('FAIL_STAGE')
if args[0] == 'sample':
    assert args[1:] == ['--proportion', '0.4', '--rand-seed', '1234']
    shutil.copyfileobj(sys.stdin.buffer, sys.stdout.buffer)
    if failure == 'sample':
        sys.exit(23)
else:
    assert args[:4] == ['head', '--number', '1000000', '--out-file']
    with gzip.open(args[4], 'wb') as output:
        # Stand in for reaching the cap after a single record.
        for _ in range(4):
            output.write(sys.stdin.buffer.readline())
    if failure == 'head':
        sys.exit(29)
""")
        seqkit.chmod(0o755)
        self.env = dict(os.environ, PATH=str(self.work) + os.pathsep + os.environ['PATH'])

    def run_helper(self, fail_stage='', expect_success=True):
        result = subprocess.run(['bash', str(SCRIPT), str(self.input), str(self.output)], cwd=self.work,
                                env=dict(self.env, FAIL_STAGE=fail_stage), capture_output=True, timeout=15)
        self.assertEqual(result.returncode == 0, expect_success, result.stderr)
        self.assertEqual(list(self.work.glob('.k2-subsample.*')), [])
        return result.returncode

    def test_reaching_cap_and_empty_gzip(self):
        for data in [READ * 100000, b'']:
            with self.subTest(size=len(data)):
                self.input.write_bytes(gzip.compress(data))
                self.run_helper()
                self.assertEqual(gzip.decompress(self.output.read_bytes()), READ if data else b'')

    def test_uncompressed_filtered_fastq(self):
        self.input = self.input.with_suffix('')
        self.input.write_bytes(READ)
        self.run_helper()
        self.assertEqual(gzip.decompress(self.output.read_bytes()), READ)

    def test_missing_and_corrupt_input_even_after_cap(self):
        self.run_helper(expect_success=False)
        self.assertFalse(self.output.exists())
        self.input.write_bytes(gzip.compress(READ * 100000)[:-8])
        self.run_helper(expect_success=False)
        self.assertFalse(self.output.exists())

    def test_sampling_and_head_failures_preserve_prior_output(self):
        self.input.write_bytes(gzip.compress(READ))
        self.output.write_bytes(b'previous result')
        for stage, code in [('sample', 23), ('head', 29)]:
            with self.subTest(stage=stage):
                self.assertEqual(self.run_helper(stage, expect_success=False), code)
                self.assertEqual(self.output.read_bytes(), b'previous result')
