"""Use the metagenome entry point to check project database paths at the consuming commands."""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "metagenome/.tests"
DATABASE = Path("references/panel v1")


class ProjectDatabaseTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name)
        shutil.copytree(FIXTURE / "data/raw", self.work / "data/raw")
        self.settings = yaml.safe_load((FIXTURE / "config.yaml").read_text())
        self.settings.update(ref_dir=str(DATABASE), skip_peppy=True)
        for provider in ("ensembl", "langmead", "ncbi_16s"):
            self.settings["db"][provider] = {"release": None}
        (self.work / "Snakefile").write_text(
            'configfile: "config.yaml"\n'
            f'include: {str(ROOT / "metagenome/metagenome.smk")!r}\n'
        )
        db = self.work / DATABASE
        db.mkdir(parents=True)
        for name in ("hash.k2d", "opts.k2d", "taxo.k2d", "database150mers.kmer_distrib"):
            (db / name).touch()

    def dry_run(self, target):
        (self.work / "config.yaml").write_text(yaml.safe_dump(self.settings))
        return subprocess.run(
            [sys.executable, "-m", "snakemake", "-n", "-p", "--cores", "1", target],
            cwd=self.work, env={k: v for k, v in os.environ.items() if not k.startswith("GCF_")},
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60,
        )

    def test_reporting_uses_project_database(self):
        for provider in ("langmead", "ncbi_16s"):
            with self.subTest(provider=provider):
                self.settings["db"]["reference_db"] = provider
                self.check_reporting_commands()

    def check_reporting_commands(self):
        result = self.dry_run("multiqc_report")
        self.assertEqual(result.returncode, 0, result.stdout)
        commands = [shlex.split(line.removeprefix("Shell command: ")) for line in result.stdout.splitlines()
                    if line.startswith("Shell command: kraken2 ") or "run_bracken.py -d " in line]
        kraken = [cmd for cmd in commands if cmd[0] == "kraken2"]
        bracken = [cmd for cmd in commands if len(cmd) > 1 and cmd[1].endswith("/run_bracken.py")]
        self.assertEqual(len(kraken), len(self.settings["samples"]))
        self.assertEqual(len(bracken), len(self.settings["samples"]))
        for cmd in kraken:
            self.assertEqual(cmd[cmd.index("--db") + 1], str(DATABASE))
        for cmd in bracken:
            self.assertEqual(cmd[cmd.index("-d") + 1], str(DATABASE))
        self.assertNotIn("rule langmead_kraken_prebuild:", result.stdout)
        self.assertNotIn("rule langmead_shmem:", result.stdout)
        self.assertNotIn("rule ncbi_16s_kraken_build:", result.stdout)
        self.assertNotIn("rule ncbi_16s_shmem:", result.stdout)

    def test_single_end_kraken_only_needs_the_index(self):
        (self.work / DATABASE / "database150mers.kmer_distrib").unlink()
        self.settings["read_geometry"] = [151]
        self.settings["delta_readlen"] = [0]
        for sample in self.settings["samples"].values():
            del sample["R2"]
        result = self.dry_run("kraken_classify_all")
        self.assertEqual(result.returncode, 0, result.stdout)
        command = next(line for line in result.stdout.splitlines() if line.startswith("Shell command: kraken2 "))
        args = shlex.split(command)
        self.assertEqual(args[args.index("--db") + 1], str(DATABASE))
        self.assertNotIn("--paired", args)

    def test_missing_project_inputs_name_the_path(self):
        for name, target in (("taxo.k2d", "kraken_classify_all"), ("database150mers.kmer_distrib", "multiqc_report")):
            with self.subTest(name=name):
                path = self.work / DATABASE / name
                path.unlink()
                result = self.dry_run(target)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(f"Missing required database file in project ref_dir: {DATABASE / name}", result.stdout)
                path.touch()


if __name__ == "__main__":
    unittest.main()
