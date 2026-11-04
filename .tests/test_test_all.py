"""Regression tests for test-all.sh exit status and invocations (issue #227).

Uses a fake `snakemake` executable and a temporary checkout copy; no workflow, container or scientific tool runs.
Run from the repository root: python -m unittest discover -s .tests -p 'test_*.py' -v
"""

import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "test-all.sh"
WORKFLOWS = ["smallrna", "singlecell", "metagenome", "rnaseq", "microbiome", "default", "wgs"]
SYSTEM_PATH = "/usr/bin:/bin"
BASH = shutil.which("bash")

FAKE_SNAKEMAKE = """#!/bin/sh
# Record the invocation, simulate Snakemake metadata, and fail for workflows listed in FAKE_FAIL.
wf=$(basename "$(dirname "$PWD")")
printf '%s|%s\\n' "$PWD" "$*" >> "$FAKE_LOG"
mkdir -p .snakemake
case " $FAKE_FAIL " in
    *" $wf "*) echo "fake snakemake: simulated failure in $wf" >&2; exit 23 ;;
esac
exit 0
"""


class TestAllScript(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="test-all "))
        self.addCleanup(shutil.rmtree, self.tmp)
        self.root = self.tmp / "checkout with space"
        self.root.mkdir()
        shutil.copy2(SCRIPT, self.root / "test-all.sh")
        for wf in WORKFLOWS:
            (self.root / wf / ".tests").mkdir(parents=True)
        self.bin = self.tmp / "bin"
        self.bin.mkdir()
        fake = self.bin / "snakemake"
        fake.write_text(FAKE_SNAKEMAKE)
        fake.chmod(fake.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        self.log = self.tmp / "invocations.log"

    def run_script(self, fail=(), fake=True):
        env = {
            "PATH": f"{self.bin}:{SYSTEM_PATH}" if fake else SYSTEM_PATH,
            "FAKE_LOG": str(self.log),
            "FAKE_FAIL": " ".join(fail),
        }
        # Invoke from an unrelated directory: the script must resolve its own checkout root.
        return subprocess.run([BASH, str(self.root / "test-all.sh")], cwd=self.tmp, env=env,
                              capture_output=True, text=True)

    def invocations(self):
        if not self.log.exists():
            return []
        return [line.split("|", 1) for line in self.log.read_text().splitlines()]

    def assert_all_invoked(self):
        expected = [[str(self.root / wf / ".tests"), "-n -p --cores 1 multiqc_report"] for wf in WORKFLOWS]
        self.assertEqual(self.invocations(), expected)

    def test_all_pass(self):
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_all_invoked()

    def test_early_failure_is_nonzero_and_named(self):
        result = self.run_script(fail=["smallrna"])
        self.assertNotEqual(result.returncode, 0)
        self.assert_all_invoked()
        self.assertIn("simulated failure in smallrna", result.stderr)
        self.assertIn("failed for: smallrna", result.stderr)

    def test_later_failures_are_nonzero_and_named(self):
        result = self.run_script(fail=["default", "wgs"])
        self.assertNotEqual(result.returncode, 0)
        self.assert_all_invoked()
        self.assertIn("failed for: default wgs", result.stderr)

    @unittest.skipIf(shutil.which("snakemake", path=SYSTEM_PATH), "snakemake installed in system PATH")
    def test_missing_snakemake_is_nonzero(self):
        result = self.run_script(fake=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("snakemake not found", result.stderr)

    def test_missing_fixture_directory_is_nonzero_before_running(self):
        shutil.rmtree(self.root / "metagenome")
        result = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("metagenome/.tests", result.stderr)
        self.assertEqual(self.invocations(), [])

    def test_cleanup_only_removes_metadata_created_by_invocation(self):
        existing = self.root / "rnaseq" / ".tests" / ".snakemake"
        existing.mkdir()
        (existing / "keep").write_text("pre-existing")
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((existing / "keep").exists())
        for wf in WORKFLOWS:
            if wf != "rnaseq":
                self.assertFalse((self.root / wf / ".tests" / ".snakemake").exists(), wf)


if __name__ == "__main__":
    unittest.main()
