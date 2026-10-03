"""Check real singlecell DAGs; never execute jobs or download scientific data.

Run with: python -m unittest discover -s singlecell/.tests -p 'test_*.py'
Requires the same Python dependencies as the singlecell dry-run CI job.
"""

import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "singlecell/.tests"
METHODS = ("splitpipe", "cellranger")


class BfqNotebookDefaultsTests(unittest.TestCase):
    def prepare(self, work, method):
        for name in ("data", "pep"):
            shutil.copytree(FIXTURE / name, work / name)
        config = yaml.safe_load((FIXTURE / "config.yaml").read_text())
        config.update(interim_dir="data/tmp", ext_dir="data/ext", fastq_dir="data/raw/fastq")
        if method == "splitpipe":
            config["libprepkit"] = "Parse Biosciences Evercode WT v3"
            config["read_geometry"] = [150, 150]
            config["quant"]["aggregate"]["method"] = "default"
            config["wells"] = {
                "sample_A": {"Sample_ID": "sample_A", "Wells": "A1-A2"},
                "sample_B": {"Sample_ID": "sample_B", "Wells": "A3-A4"},
            }
        (work / "config.yaml").write_text(yaml.safe_dump(config))
        (work / "Snakefile").write_text(
            "pepfile: 'pep/pep_config.yaml'\n"
            "configfile: 'config.yaml'\n"
            f"include: {str(ROOT / 'singlecell/singlecell.smk')!r}\n"
        )

    def snakemake(self, work, target):
        env = {k: v for k, v in os.environ.items() if not k.startswith("GCF_")}
        env["XDG_CACHE_HOME"] = str(work / "cache")
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--cores", "1", "--scheduler", "greedy",
             "--dry-run", "--printshellcmds", "--rerun-incomplete",
             target],
            cwd=work, env=env, capture_output=True, text=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def jobs(self, output):
        return set(re.findall(r"^(?:local)?rule (\w+):", output, re.MULTILINE))

    def test_default_targets_keep_scientific_outputs_without_notebooks(self):
        for method in METHODS:
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                self.prepare(work, method)
                for target in ("bfq_all", "multiqc_report"):
                    output = self.snakemake(work, target)
                    jobs = self.jobs(output)
                    required = {
                        target, "scanpy_aggr_finalize", "autoqc_prepare", "autoqc_mad",
                        "dbl_classification_aggr", "bfq_level1_all", "bfq_level2_exprs",
                        "bfq_level2_logs", f"{method}_quant", f"{method}_aggr",
                    }
                    required.add("bfq_level2_figs" if method == "splitpipe" else "bfq_level2_umap_png")
                    if method == "cellranger":
                        required.add("bfq_level2_data")
                    self.assertTrue(required <= jobs, output)
                    self.assertFalse(any("ipynb" in job or "notebook" in job for job in jobs), output)
                    self.assertIn("all_samples_filtered.h5ad", output)
                    self.assertNotIn("_preprocessed.h5ad", output)

    def test_notebooks_remain_explicit_targets(self):
        for method in METHODS:
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                self.prepare(work, method)
                producer = f"{method}_scanpy_pp_ipynb"
                jobs = self.jobs(self.snakemake(work, "bfq_level2_notebooks"))
                self.assertTrue({producer, producer + "_html", "bfq_level2_notebooks"} <= jobs)
                target = f"data/tmp/singlecell/quant/aggregate/{method}/scanpy/all_samples_preprocessed.h5ad"
                jobs = self.jobs(self.snakemake(work, target))
                self.assertIn(producer, jobs)
                self.assertNotIn(producer + "_html", jobs)
                self.assertNotIn("bfq_level2_notebooks", jobs)

    def test_reused_workdir_only_needs_multiqc(self):
        for method in METHODS:
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                self.prepare(work, method)
                output = self.snakemake(work, "multiqc_report")
                retained = []
                # Synthetic placeholders model completed upstream outputs. Their
                # contents are never read by scientific tools: every call is a dry run.
                # Read the dry-run job outputs rather than --summary, which is
                # broken in the CI-pinned Snakemake 8.1.1 (async generator error).
                rule = None
                for line in output.splitlines():
                    match = re.match(r"^(?:local)?rule (\w+):", line)
                    if match:
                        rule = match.group(1)
                    if not line.startswith("    output: ") or rule == "multiqc_report":
                        continue
                    for filename in line.removeprefix("    output: ").split(", "):
                        path = work / filename
                        if not path.resolve().is_relative_to(work):
                            # Do not create Kraken's temporary /dev/shm copies.
                            self.assertEqual(rule, "langmead_shmem")
                            continue
                        self.assertNotIn("_preprocessed.h5ad", str(path))
                        self.assertNotIn("notebooks", path.parts)
                        path.parent.mkdir(parents=True, exist_ok=True)
                        self.assertFalse(path.exists(), path)
                        path.write_text("Synthetic DAG fixture; not scientific data.\n")
                        retained.append(path)
                now = time.time()
                for path in retained:
                    os.utime(path, (now, now))
                failed_log = work / f"data/tmp/singlecell/quant/aggregate/{method}/scanpy/notebooks/all_samples_pp.ipynb"
                failed_log.parent.mkdir(parents=True, exist_ok=True)
                failed_log.write_text("Synthetic failed notebook log.\n")
                retained.append(failed_log)
                before = {p: (hashlib.sha256(p.read_bytes()).digest(), p.stat().st_mtime_ns) for p in retained}
                self.assertEqual(self.jobs(self.snakemake(work, "multiqc_report")), {"multiqc_report"})
                after = {p: (hashlib.sha256(p.read_bytes()).digest(), p.stat().st_mtime_ns) for p in retained}
                self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()
