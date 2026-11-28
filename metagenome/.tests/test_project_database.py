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

    def dry_run(self, target, *options, work=None):
        work = work or self.work
        (work / "config.yaml").write_text(yaml.safe_dump(self.settings))
        return subprocess.run(
            [sys.executable, "-m", "snakemake", "-n", "-p", "--cores", "1", *options, target],
            cwd=work, env={k: v for k, v in os.environ.items() if not k.startswith("GCF_")},
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60,
        )

    def test_reporting_uses_project_database(self):
        for provider in ("langmead", "ncbi_16s"):
            with self.subTest(provider=provider):
                self.settings["db"]["reference_db"] = provider
                self.check_reporting_commands()

    def check_reporting_commands(self, *options):
        result = self.dry_run("multiqc_report", *options)
        self.assertEqual(result.returncode, 0, result.stdout)
        commands = [shlex.split(line.removeprefix("Shell command: ")) for line in result.stdout.splitlines()
                    if line.startswith("Shell command: kraken2 ") or "run_bracken.py -d " in line]
        kraken = [cmd for cmd in commands if cmd[0] == "kraken2"]
        bracken = [cmd for cmd in commands if len(cmd) > 1 and cmd[1].endswith("/run_bracken.py")]
        self.assertEqual(len(kraken), len(self.settings["samples"]))
        self.assertEqual(len(bracken), len(self.settings["samples"]))
        expected = str((self.work / self.settings["ref_dir"]).resolve())
        for cmd in kraken:
            self.assertEqual(cmd[cmd.index("--db") + 1], expected)
        for cmd in bracken:
            self.assertEqual(cmd[cmd.index("-d") + 1], expected)
        self.assertNotIn("rule langmead_kraken_prebuild:", result.stdout)
        self.assertNotIn("rule langmead_shmem:", result.stdout)
        self.assertNotIn("rule ncbi_16s_kraken_build:", result.stdout)
        self.assertNotIn("rule ncbi_16s_bracken:", result.stdout)
        self.assertNotIn("rule ncbi_16s_shmem:", result.stdout)

    def test_explicit_databases_at_managed_paths_are_inputs_only(self):
        for provider, release, assembly in (("langmead", "20260626", "k2_pluspf"), ("ncbi_16s", "222", "ncbi_16s")):
            with self.subTest(provider=provider):
                path = Path("data/ext") / provider / f"release-{release}" / "metagenome" / assembly
                db = self.work / path
                shutil.copytree(self.work / DATABASE, db)
                self.settings["ref_dir"] = str(path)
                self.settings["db"]["reference_db"] = provider
                self.settings["db"][provider]["release"] = release
                if provider == "ncbi_16s":
                    # A library marker newer than the supplied index triggers a rebuild if the getter is registered.
                    marker = db / "library/added/libraries.added"
                    marker.parent.mkdir(parents=True)
                    marker.touch()
                    for name in ("hash.k2d", "taxo.k2d", "opts.k2d"):
                        os.utime(db / name, (1, 1))
                    self.check_reporting_commands()
                else:
                    # The archive getter has no inputs, so force all jobs to expose any producer for this path.
                    self.check_reporting_commands("--forceall")

    def test_independent_screening_respects_its_configured_database(self):
        self.settings["db"]["langmead"]["release"] = "20260626"
        ram_db = "/dev/shm/langmead/release-20260626/metagenome/k2_pluspf"
        disk_db = "data/ext/langmead/release-20260626/metagenome/k2_pluspf"
        shutil.copytree(self.work / DATABASE, self.work / disk_db)
        alias = self.work / "current"
        alias.symlink_to(disk_db, target_is_directory=True)
        for ref_dir, expected_db in ((str(DATABASE), ram_db), (disk_db, str(self.work / disk_db)),
                                     (alias.name, str(self.work / disk_db))):
            with self.subTest(ref_dir=ref_dir):
                self.settings["ref_dir"] = ref_dir
                result = self.dry_run("data/tmp/metagenome/qc/kraken2/Bact.kraken.kreport", "--forceall")
                self.assertEqual(result.returncode, 0, result.stdout)
                command = next(line for line in result.stdout.splitlines() if "k2_screen.sh --input" in line)
                args = shlex.split(command)
                self.assertEqual(args[args.index("--db") + 1], expected_db)
                self.assertEqual("rule langmead_kraken_prebuild:" in result.stdout, expected_db == ram_db)

    def test_database_symlink_retarget_between_kraken_and_bracken(self):
        first, second = self.work / DATABASE, self.work / "another database"
        shutil.copytree(first, second)
        for db, label in ((first, "first"), (second, "second")):
            (db / "hash.k2d").write_text(label)
            (db / "database150mers.kmer_distrib").write_text(label)
        alias = self.work / "current"
        alias.symlink_to(first, target_is_directory=True)
        self.settings["ref_dir"] = alias.name
        target = Path("data/tmp/metagenome/quant/kraken2/k2_pluspf/Bact/Bact.bracken_out")

        def plan(work, expected):
            result = self.dry_run(str(target), work=work)
            self.assertEqual(result.returncode, 0, result.stdout)
            lines = result.stdout.splitlines()
            kraken = shlex.split(next(line for line in lines if line.startswith("Shell command: kraken2 ")))
            bracken = shlex.split(next(line for line in lines if "run_bracken.py -d " in line))
            self.assertEqual(kraken[kraken.index("--db") + 1], str(expected))
            self.assertEqual(bracken[bracken.index("-d") + 1], str(expected))
            self.assertEqual(bracken[bracken.index("-r") + 1], "150")
            for rule, files in (("kraken_classify", ("hash.k2d", "opts.k2d", "taxo.k2d")),
                                ("bracken", ("database150mers.kmer_distrib",))):
                inputs = lines[lines.index(f"rule {rule}:") + 1].strip().removeprefix("input: ").split(", ")
                for name in files:
                    self.assertIn(str(expected / name), inputs)
            return [line for line in lines if line.startswith("rule ")]

        # Supply the filtered reads so the only jobs run are the two command stand-ins.
        reads = self.work / "data/tmp/metagenome/filter/fastp"
        reads.mkdir(parents=True)
        for mate in (1, 2):
            (reads / f"Bact_R{mate}.fastq.gz").touch()
        self.assertEqual(plan(self.work, first), ["rule kraken_classify:", "rule bracken:"])
        tools = self.work / "bin"
        tools.mkdir()
        stub = ('#!/usr/bin/env python3\nfrom pathlib import Path\nimport os, sys\n'
                'def arg(flag): return sys.argv[sys.argv.index(flag) + 1]\n')
        kraken = tools / "kraken2"
        kraken.write_text(stub +
            'db = Path(arg("--db"))\n'
            'label = (db / "hash.k2d").read_text()\n'
            'Path(arg("--output")).write_text(label)\n'
            'Path(arg("--report")).write_text(f"100\\t1\\t1\\tS\\t1\\t{label}\\n")\n'
            'alias = Path(os.environ["PROJECT_ALIAS"])\n'
            'alias.unlink()\n'
            'alias.symlink_to(os.environ["NEXT_DB"], target_is_directory=True)\n')
        bracken = tools / "bracken"
        bracken.write_text(stub +
            'assert arg("-l") == "S"\n'
            'label = (Path(arg("-d")) / ("database" + arg("-r") + "mers.kmer_distrib")).read_text()\n'
            'classified = Path(arg("-i")).read_text().split()[-1]\n'
            'Path(arg("-o")).write_text(classified + " " + label)\n'
            'Path(arg("-w")).write_text(label)\n')
        for tool in (kraken, bracken):
            tool.chmod(0o755)
        env = {k: v for k, v in os.environ.items() if not k.startswith("GCF_")}
        env.update(PATH=str(tools) + os.pathsep + env["PATH"], PROJECT_ALIAS=str(alias), NEXT_DB=str(second))
        built = subprocess.run([sys.executable, "-m", "snakemake", "--cores", "1", str(target)],
                               cwd=self.work, env=env, capture_output=True, text=True, timeout=60)
        self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
        self.assertEqual((self.work / target).read_text(), "first first")
        self.assertEqual(alias.resolve(), second)

        # Both a retained analysis with finished jobs and a new workdir now select the second database.
        plan(self.work, second)
        fresh = self.work / "fresh"
        shutil.copytree(self.work / "data/raw", fresh / "data/raw")
        shutil.copyfile(self.work / "Snakefile", fresh / "Snakefile")
        (fresh / alias.name).symlink_to(second, target_is_directory=True)
        plan(fresh, second)

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
        self.assertEqual(args[args.index("--db") + 1], str(self.work / DATABASE))
        self.assertNotIn("--paired", args)

    def test_missing_project_inputs_name_the_path(self):
        for name, target in (("taxo.k2d", "kraken_classify_all"), ("database150mers.kmer_distrib", "multiqc_report")):
            with self.subTest(name=name):
                path = self.work / DATABASE / name
                path.unlink()
                result = self.dry_run(target)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(f"Missing required database file in project ref_dir: {path}", result.stdout)
                path.touch()


if __name__ == "__main__":
    unittest.main()
