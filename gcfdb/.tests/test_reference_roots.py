"""Reference getters resolve their root through the shared EXT_DIR contract (issue #228).

The contract in utils.py: a nonempty project `ext_dir` wins over `GCF_EXT`, which wins over the default `data/ext`.
These tests dry-run a real getter (gcfdb/homd.db) included after utils.py, so nothing is downloaded or executed.

Run from the repository root: python -m unittest discover -s gcfdb/.tests -p 'test_*.py' -v
Requires Snakemake and PyYAML.
"""

import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RELEASE = "15.2"


class ReferenceRootTests(unittest.TestCase):
    def dry_run(self, work, ext_dir_config=None, ext_dir_env=None, expected_root=None):
        settings = {"workflow": "default", "read_geometry": [50]}
        if ext_dir_config is not None:
            settings["ext_dir"] = ext_dir_config
        snakefile = work / "Snakefile"
        snakefile.write_text(
            f"config.update({settings!r})\n"
            f"include: {str(ROOT / 'utils.py')!r}\n"
            f"include: {str(ROOT / 'gcfdb' / 'homd.db')!r}\n"
        )
        env = {k: v for k, v in os.environ.items() if k not in ("GCF_EXT", "GCF_SECRET", "GCF_INTERIM", "GCF_FASTQ")}
        if ext_dir_env is not None:
            env["GCF_EXT"] = ext_dir_env
        target = os.path.join(expected_root, "homdb", RELEASE, "fasta", "repseq.fa")
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile), "--directory", str(work),
             "--dry-run", "--printshellcmds", "--cores", "1", target],
            cwd=work, env=env, capture_output=True, text=True, timeout=120,
        )
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output)
        self.assertIn("rule homdb_db:", output)
        self.assertIn(target, output)
        return output

    def test_default_root(self):
        with tempfile.TemporaryDirectory() as d:
            self.dry_run(Path(d), expected_root="data/ext")

    def test_environment_root(self):
        with tempfile.TemporaryDirectory() as d:
            env_root = str(Path(d) / "env-ext")
            self.dry_run(Path(d), ext_dir_env=env_root, expected_root=env_root)

    def test_config_root(self):
        with tempfile.TemporaryDirectory() as d:
            config_root = str(Path(d) / "config-ext")
            self.dry_run(Path(d), ext_dir_config=config_root, expected_root=config_root)

    def test_config_wins_over_environment(self):
        with tempfile.TemporaryDirectory() as d:
            config_root, env_root = str(Path(d) / "config-ext"), str(Path(d) / "env-ext")
            output = self.dry_run(Path(d), ext_dir_config=config_root, ext_dir_env=env_root, expected_root=config_root)
            self.assertNotIn(env_root + os.sep, output)


class NoHardcodedReferenceRootTests(unittest.TestCase):
    """Executable workflow sources must not pin facility-specific absolute reference roots."""

    def test_no_absolute_ext_roots(self):
        files = subprocess.run(["git", "ls-files", "*.smk", "*.db", "*/Snakefile", "utils.py", "common.smk"],
                               cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
        files = [f for f in files if not f.startswith("sars-cov2/external_pipelines/")]
        self.assertTrue(files)
        assignment = re.compile(r"^\s*EXT_DIR\s*=", re.M)
        absolute_mnt = re.compile(r"['\"]/mnt/")
        offenders = []
        for f in files:
            text = (ROOT / f).read_text(errors="replace")
            if f != "utils.py" and assignment.search(text):
                offenders.append(f"{f}: assigns EXT_DIR")
            if absolute_mnt.search(text):
                offenders.append(f"{f}: absolute /mnt/ path literal")
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
