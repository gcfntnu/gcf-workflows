"""Exercise report provenance with real repositories and the config-generation CLI.

Run: python -m unittest discover -s misc/multiqc/.tests -p 'test_*.py' -v
Requires Git, pandas, PyYAML, matplotlib and peppy. No sequencing data or containers.
"""

import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import yaml


MULTIQC = Path(__file__).resolve().parents[1]
SCRIPT = MULTIQC / "create_mqc_config.py"
SPEC = importlib.util.spec_from_file_location("create_mqc_config", SCRIPT)
MQC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MQC)
BRANCH = "issues/184-optional-singlecell-notebooks"


class GitRevisionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.repo = self.work / "workflow checkout"
        self.repo.mkdir()
        self.git("init", "-q")
        self.git("symbolic-ref", "HEAD", "refs/heads/bfq-dev")
        self.git("-c", "user.name=Test fixture", "-c", "user.email=fixture@example.invalid",
                 "-c", "commit.gpgsign=false", "commit", "--allow-empty", "-qm", "Fixture commit")
        self.commit = self.git("rev-parse", "HEAD")

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.repo), *args],
                              check=True, capture_output=True, text=True).stdout.strip()

    def expected(self, ref):
        return f"Analysis pipeline: github.com/gcfntnu/gcf-workflows/tree/{ref} commit {self.commit}"

    def assert_revision(self, ref, repo=None):
        self.assertEqual(MQC.get_software_versions(SimpleNamespace(repo_dir=repo or self.repo)), self.expected(ref))

    def test_ordinary_branch(self):
        self.assert_revision("bfq-dev")

    def test_slash_branches(self):
        for branch in (BRANCH, "issues/reporting/nested-branch"):
            with self.subTest(branch=branch):
                self.git("checkout", "-qb", branch)
                self.assert_revision(branch)

    def test_packed_branch_ref(self):
        for branch in ("bfq-dev", BRANCH):
            with self.subTest(branch=branch):
                if branch != "bfq-dev":
                    self.git("checkout", "-qb", branch)
                self.git("pack-refs", "--all", "--prune")
                self.assertFalse((self.repo / ".git/refs/heads" / branch).exists())
                self.assert_revision(branch)

    def test_detached_head(self):
        self.git("checkout", "--detach", "-q", self.commit)
        self.assert_revision(self.commit)

    def test_linked_worktree(self):
        worktree = self.work / "linked workflow"
        self.git("worktree", "add", "-qb", BRANCH, str(worktree))
        self.assertTrue((worktree / ".git").is_file())
        self.assert_revision(BRANCH, worktree)

    def test_branch_name_is_not_changed_by_matching_tag(self):
        self.git("tag", "bfq-dev")
        self.assert_revision("bfq-dev")

    def test_invalid_repository(self):
        for repo in (self.work, self.work / "missing"):
            with self.subTest(repo=repo), self.assertRaisesRegex(RuntimeError, "Cannot resolve Git revision") as caught:
                MQC.get_software_versions(SimpleNamespace(repo_dir=repo))
            self.assertIn(str(repo), str(caught.exception))

    def test_unborn_branch(self):
        self.git("checkout", "--orphan", "no-commits", "-q")
        with self.assertRaisesRegex(RuntimeError, "Cannot resolve Git revision"):
            MQC.get_software_versions(SimpleNamespace(repo_dir=self.repo))

    def test_missing_git(self):
        with patch.dict(os.environ, {"PATH": str(self.work / "no-executables")}):
            with self.assertRaisesRegex(RuntimeError, "Cannot run Git") as caught:
                MQC.get_software_versions(SimpleNamespace(repo_dir=self.repo))
        self.assertIn(str(self.repo), str(caught.exception))

    def test_config_cli(self):
        self.git("checkout", "-qb", BRANCH)
        sample_info = self.work / "sample_info.tsv"
        sample_info.write_text("Sample_ID\tExternal_ID\nS1\tExample\n")
        config = self.work / "config.yaml"
        config.write_text(yaml.safe_dump({"samples": {"S1": {"External_ID": "Example"}}}))
        pep = self.work / "pep_config.yaml"
        pep.write_text("pep_version: 2.0.0\nsample_table: samples.csv\n")
        (self.work / "samples.csv").write_text("sample_name,External_ID\nS1,Example\n")
        output = self.work / "multiqc.yaml"
        for detached in (False, True):
            if detached:
                self.git("checkout", "--detach", "-q", self.commit)
            ref = self.commit if detached else BRANCH
            for project in (config, pep):
                with self.subTest(detached=detached, project=project.name):
                    result = subprocess.run(
                        [sys.executable, str(SCRIPT), "-p", "GCF-TEST", "-S", str(sample_info),
                         "--repo-dir", str(self.repo), "--header-template", str(MULTIQC / "mqc_header.txt"),
                         "--config-template", str(MULTIQC / "multiqc_config-default.yaml"),
                         "--pep", str(project), "--read-geometry", "75,75", "-o", str(output)],
                        cwd=self.work, capture_output=True, text=True, timeout=30,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    report = yaml.safe_load(output.read_text())
                    self.assertEqual(report["title"], "GCF-TEST")
                    self.assertIn(self.expected(ref), report["intro_text"])
                    self.assertEqual(report["custom_data"]["general_statistics"]["data"]["S1"]["External_ID"], "Example")


class SampleMetadataTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        # Deliberately stale metadata: only the TSV may supply sample rows/values.
        self.config = self.work / "active.yaml"
        self.config.write_text(yaml.safe_dump({"samples": {"stale": {"External_ID": "outdated"}}}))
        self.pep = self.work / "pep_config.yaml"
        self.pep.write_text("pep_version: 2.0.0\nsample_table: samples.csv\n")
        (self.work / "samples.csv").write_text("sample_name,External_ID\nstale,outdated\n")
        self.projects = (self.config, self.pep)
        version = patch.object(MQC, "get_software_versions", return_value="Analysis pipeline: synthetic")
        version.start()
        self.addCleanup(version.stop)

    def generate(self, project, contents, descriptors=None):
        config = yaml.safe_load(project.read_text())
        config["descriptors"] = descriptors or {}
        project.write_text(yaml.safe_dump(config))
        info = io.StringIO(contents)
        info.name = "sample_info.tsv"
        with project.open() as pep:
            args = SimpleNamespace(
                pep=pep, sample_info=info, config_template=io.StringIO("extra_fn_clean_exts: []\n"),
                header_template=io.StringIO("Project {pname}"), project_id="GCF-TEST", read_geometry=[75, 75],
                machine="Synthetic", organism="N/A", libkit="default", workflow="default", repo_dir=self.work,
            )
            # Verify the serialized configuration consumed by MultiQC, including key types.
            return yaml.safe_load(yaml.safe_dump(MQC.create_mqc_config(args)))

    def test_tsv_is_authoritative(self):
        for project in self.projects:
            with self.subTest(project=project.name):
                report = self.generate(project, "Sample_ID\tExternal_ID\nS1\tcurrent\nS2\tadditional\n")
                self.assertEqual(report["custom_data"]["general_statistics"]["data"], {
                    "S1": {"External_ID": "current"}, "S2": {"External_ID": "additional"},
                })

    def test_identifiers_are_literal_text(self):
        for project in self.projects:
            for ids in (("001", "002"), ("001", "1"), ("1e3", "1000"), ("NA", "NULL", "NaN", "None", "N/A")):
                with self.subTest(project=project.name, ids=ids):
                    contents = "Sample_ID\tExternal_ID\n" + "".join(f"{sid}\t{sid}\n" for sid in ids)
                    report = self.generate(project, contents)
                    self.assertEqual(report["custom_data"]["general_statistics"]["data"], {
                        sid: {"External_ID": sid} for sid in ids
                    })

    def test_text_metadata_is_not_inferred_as_numbers(self):
        for project in self.projects:
            with self.subTest(project=project.name):
                report = self.generate(
                    project, "Sample_ID\tExternal_ID\tSample_Group\tSubmitted_Code\tCustomer_Comments\nS1\t0007\t01\t1e3\tNA\n",
                    {"External_ID": {"dtype": "numerical"}},
                )
                self.assertEqual(report["custom_data"]["general_statistics"]["data"]["S1"], {
                    "External_ID": "0007", "Sample_Group": "01", "Submitted_Code": "1e3", "Customer_Comments": "NA",
                })

    def test_numeric_descriptors_keep_qc_formatting_and_annotations(self):
        for project in self.projects:
            with self.subTest(project=project.name):
                report = self.generate(
                    project, "Sample_ID\tExternal_ID\tRIN\tConcentration\tEmpty\n001\t0007\t8.126\t12.3456\t\n"
                    "1\t0008\tNA\t<1\t\n", {"RIN": {"dtype": "numerical"}, "Concentration": {"dtype": "numerical"}},
                )
                self.assertEqual(report["custom_data"]["general_statistics"]["data"], {
                    "001": {"External_ID": "0007", "RIN": 8.13, "Concentration": 12.35},
                    "1": {"External_ID": "0008", "RIN": "", "Concentration": "<1"},
                })

    def test_exact_duplicates_are_rejected(self):
        for project in self.projects:
            for sid in ("S1", "001", "NA"):
                with self.subTest(project=project.name, sid=sid), self.assertRaisesRegex(ValueError, "duplicate Sample_ID"):
                    self.generate(project, f"Sample_ID\tExternal_ID\n{sid}\tone\n{sid}\ttwo\n")

    def test_missing_or_blank_sample_ids_are_rejected(self):
        cases = (
            ("External_ID\nExample\n", "missing required column 'Sample_ID'"),
            ("Sample_ID\tExternal_ID\n\tExample\n", "blank Sample_ID"),
            ("Sample_ID\tExternal_ID\n   \tExample\n", "blank Sample_ID"),
        )
        for project in self.projects:
            for contents, message in cases:
                with self.subTest(project=project.name, contents=contents), self.assertRaisesRegex(ValueError, message):
                    self.generate(project, contents)


if __name__ == "__main__":
    unittest.main()
