"""Managed raw genome and annotation files become visible only after their producer succeeds."""
import gzip
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time
import unittest


PUBLISH = Path(__file__).resolve().parents[1] / 'scripts/publish_reference.sh'


class ManagedReferenceDownload(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name)
        self.wget = self.work / 'wget'
        # The synthetic source arguments are local files, including concatenated compressed members.
        self.wget.write_text('#!/bin/bash\nset -eu\n[[ $1 = -O- && $2 = -- ]]\nshift 2\ncat -- "$@"\nexit "${WGET_STATUS:-0}"\n')
        self.wget.chmod(0o755)
        self.env = dict(os.environ, PATH=str(self.work) + os.pathsep + os.environ['PATH'])
        self.output = self.work / 'shared reference.fa'

    def fetch(self, mode, *sources, **env):
        return subprocess.run(['bash', PUBLISH, self.output, mode, '--', *sources], env=dict(self.env, **env),
                              capture_output=True, text=True, timeout=15)

    def test_plain_and_multiple_compressed_sources_preserve_source_order(self):
        first, second = self.work / 'part 1', self.work / 'part 2'
        first.write_bytes(b'>synthetic\nACGT\n')
        self.assertEqual(self.fetch('plain', first).returncode, 0)
        self.assertEqual(self.output.read_bytes(), first.read_bytes())
        self.assertEqual(self.output.stat().st_mode & 0o777, first.stat().st_mode & 0o777)
        first.write_bytes(gzip.compress(b'#header\nfirst\n'))
        second.write_bytes(gzip.compress(b'#second header\nsecond\n'))
        self.assertEqual(self.fetch('gzip', first, second).returncode, 0)
        self.assertEqual(self.output.read_bytes(), b'#header\nfirst\n#second header\nsecond\n')
        self.assertEqual(self.fetch('refseq-gtf', first, second).returncode, 0)
        self.assertEqual(self.output.read_bytes(), b'first\nsecond\n')

    def test_network_error_or_damaged_gzip_preserves_complete_destination(self):
        self.output.write_bytes(b'published data\n')
        source = self.work / 'download'
        source.write_bytes(gzip.compress(b'replacement\n'))
        self.assertNotEqual(self.fetch('gzip', source, WGET_STATUS='1').returncode, 0)
        self.assertEqual(self.output.read_bytes(), b'published data\n')
        source.write_bytes(source.read_bytes()[:-8])
        self.assertNotEqual(self.fetch('gzip', source).returncode, 0)
        self.assertEqual(self.output.read_bytes(), b'published data\n')
        self.assertEqual(list(self.work.glob('.gcf-reference.*')), [])

    def test_concurrent_writers_stage_separately_and_interrupted_writer_cleans_up(self):
        complete = self.work / 'complete download'
        complete.write_bytes(b'complete\n')
        ready = self.work / 'ready'
        # The command-mode publisher supplies its private output path after producer arguments.
        blocked = subprocess.Popen(['bash', PUBLISH, self.output, 'command', 'bash', '-c',
                                    'printf partial > "$2"; touch "$1"; sleep 60', 'bash', ready], start_new_session=True)
        try:
            deadline = time.monotonic() + 10
            while not ready.exists() and blocked.poll() is None and time.monotonic() < deadline:
                time.sleep(.02)
            self.assertTrue(ready.exists())
            self.assertFalse(self.output.exists())
            staged, = self.work.glob('.gcf-reference.*/content')
            self.assertEqual(staged.read_bytes(), b'partial')
            self.assertEqual(self.fetch('plain', complete).returncode, 0)
            self.assertEqual(self.output.read_bytes(), complete.read_bytes())
        finally:
            os.killpg(blocked.pid, signal.SIGTERM)
            blocked.wait(timeout=10)
        self.assertNotEqual(blocked.returncode, 0)
        self.assertEqual(self.output.read_bytes(), complete.read_bytes())
        self.assertEqual(list(self.work.glob('.gcf-reference.*')), [])
