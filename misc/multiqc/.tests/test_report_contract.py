"""Dry-run the real shared reporting rules with synthetic inputs and no containers."""

import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
PROFILES = ("default", "rnaseq", "microbiome", "singlecell", "smallrna", "wgs", "metagenome")


class ReportContractTests(unittest.TestCase):
    def dry_run(self, profile="singlecell", search_dirs=(), overrides=None, prepared_config=True):
        temporary = tempfile.TemporaryDirectory(prefix="multiqc-contract-")
        self.addCleanup(temporary.cleanup)
        work = Path(temporary.name)
        interim = Path("data/tmp")
        report_dir = interim / profile / "bfq"
        (work / report_dir).mkdir(parents=True)
        (work / interim / "sample_info.tsv").write_text("Sample_ID\tExternal_ID\n001\t0007\n")
        if prepared_config:
            (work / report_dir / ".multiqc_config.yaml").write_text("title: Synthetic\n")
        for directory in search_dirs:
            (work / report_dir / directory).mkdir()
        config = {
            "project_id": ["SYNTHETIC"], "workflow": profile, "organism": "homo_sapiens",
            "libprepkit": "Custom", "read_geometry": [151, 151], "machine": "Synthetic", "skip_peppy": True,
        }
        config.update(overrides or {})
        (work / "config.yaml").write_text(yaml.safe_dump(config))
        # Load the same defaults and rules as a full workflow, without adding scientific jobs.
        (work / "Snakefile").write_text(
            'configfile: "config.yaml"\n'
            f'include: {str(ROOT / "utils.py")!r}\n'
            f'include: {str(ROOT / "postprocess.smk")!r}\n'
            'import json\n'
            'from pathlib import Path\n'
            'report_rule = workflow.get_rule("multiqc_report")\n'
            'config_rule = workflow.get_rule("multiqc_config")\n'
            'Path("contract.json").write_text(json.dumps({\n'
            '    "report_image": report_rule.container_img, "report_output": list(report_rule.output),\n'
            '    "config_image": config_rule.container_img, "sample_info": config_rule.input.sample_info,\n'
            '}))\n'
        )
        env = {key: value for key, value in os.environ.items() if not key.startswith("GCF_")}
        for key, directory in (("TMPDIR", "tmp"), ("XDG_CACHE_HOME", "cache"), ("MPLCONFIGDIR", "matplotlib")):
            (work / directory).mkdir()
            env[key] = str(work / directory)
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--dry-run", "--printshellcmds", "--cores", "1", "multiqc_report"],
            cwd=work, env=env, text=True, capture_output=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        commands = [shlex.split(line.removeprefix("Shell command: ")) for line in result.stdout.splitlines()
                    if line.startswith("Shell command: multiqc ")]
        self.assertEqual(len(commands), 1, result.stdout)
        self.assertFalse((work / report_dir / "multiqc_SYNTHETIC.html").exists(), "A dry-run must not execute MultiQC")
        return commands[0], json.loads((work / "contract.json").read_text()), result.stdout

    def test_profiles_keep_their_reporting_contract(self):
        main_config = yaml.safe_load((ROOT / "main.config").read_text())
        docker_config = yaml.safe_load((ROOT / "docker.config").read_text())["docker"]
        self.assertRegex(docker_config["multiqc"], r"^gcfntnu/gcf-multiqc@sha256:[0-9a-f]{64}$")
        for profile in PROFILES:
            with self.subTest(profile=profile):
                command, contract, _ = self.dry_run(profile, search_dirs=("logs", "figs"))
                modules = main_config["multiqc"][profile]["modules"].split(",")
                module_args = [arg for module in modules for arg in ("-m", module)]
                report_dir = f"data/tmp/{profile}/bfq"
                report = f"{report_dir}/multiqc_SYNTHETIC.html"
                self.assertEqual(command, [
                    "multiqc", "--config", f"{report_dir}/.multiqc_config.yaml", *module_args,
                    "--filename", report, "-f", "--interactive", f"{report_dir}/logs", f"{report_dir}/figs",
                ])
                self.assertEqual(contract["report_image"], "docker://" + docker_config["multiqc"])
                self.assertEqual(contract["config_image"], "docker://" + docker_config["default"])
                self.assertEqual(contract["report_output"], [report])
                self.assertEqual(contract["sample_info"], "data/tmp/sample_info.tsv")
        self.assertIn("parse", main_config["multiqc"]["singlecell"]["modules"].split(","))

    def test_search_path_uses_available_directories_or_bfq_fallback(self):
        report_dir = "data/tmp/singlecell/bfq"
        for directories in ((), ("logs",), ("figs",)):
            with self.subTest(directories=directories):
                command, _, _ = self.dry_run(search_dirs=directories)
                paths = command[command.index("--interactive") + 1:]
                self.assertEqual(paths, [f"{report_dir}/{name}" for name in directories] or [report_dir])

    def test_explicit_project_image_and_module_overrides_remain_authoritative(self):
        image = "example.invalid/multiqc:reviewed-override"
        command, contract, _ = self.dry_run(overrides={
            "docker": {"multiqc": image}, "multiqc": {"singlecell": {"modules": "parse,custom_content"}},
        })
        self.assertEqual(contract["report_image"], "docker://" + image)
        self.assertEqual([command[index + 1] for index, arg in enumerate(command) if arg == "-m"],
                         ["parse", "custom_content"])
        self.assertNotEqual(contract["config_image"], contract["report_image"])

    def test_report_requires_config_generated_from_curated_metadata(self):
        _, _, output = self.dry_run(prepared_config=False)
        commands = [shlex.split(line.removeprefix("Shell command: ")) for line in output.splitlines()
                    if line.startswith("Shell command: python ")]
        self.assertEqual(len(commands), 1, output)
        command = commands[0]
        self.assertEqual(command[1], str(ROOT / "misc/multiqc/create_mqc_config.py"))
        self.assertEqual(command[command.index("-S") + 1], "data/tmp/sample_info.tsv")
        self.assertEqual(command[command.index("--pep") + 1], "config.yaml")
        self.assertEqual(command[command.index("--config-template") + 1],
                         str(ROOT / "misc/multiqc/multiqc_config-singlecell.yaml"))
        self.assertEqual(command[command.index("-o") + 1], "data/tmp/singlecell/bfq/.multiqc_config.yaml")


if __name__ == "__main__":
    unittest.main()
