"""Execute the real BFQ QC delivery rule with tiny synthetic Snakemake inputs.

Run with: python -m unittest discover -s wgs/.tests -p 'test_bfq_qc.py' -v
Requires Snakemake 9.7.1; no scientific tools, containers or sequencing data.
"""

import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SAMPLES = ("earlier", "later")
METRICS = ("alignment", "wgs", "duplicates")


class BfqQcTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="wgs-qc-")
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.source = self.work / "source"
        for index, sample in enumerate(SAMPLES):
            source = self.source / sample
            qualimap = source / "qualimap"
            (qualimap / "raw_data").mkdir(parents=True)
            for name in ("genome_results.txt", "qualimapReport.html", "raw_data/coverage.txt"):
                (qualimap / name).write_text(f"{sample}: {name}\n")
            (qualimap / ".snakemake_timestamp").touch()
            for metric in METRICS:
                (source / f"{metric}.metrics").write_text(f"{sample}: {metric}\n")
            for suffix in ("bam", "bai"):
                (source / f"aligned.{suffix}").write_text("synthetic alignment placeholder\n")
            # The old aggregate job linked the earlier directory but consumed
            # later inputs too, triggering Snakemake 9's timestamp validation.
            timestamp = time.time() - 120 + index * 30
            for path in sorted(source.rglob("*"), reverse=True):
                os.utime(path, (timestamp, timestamp))
            os.utime(source, (timestamp, timestamp))
        (self.work / "Snakefile").write_text(
            "from os.path import join\n"
            f"SAMPLES = {list(SAMPLES)!r}\n"
            "BFQ_INTERIM = 'delivery'\nBFQ_ALL = []\n"
            "rule qualimap_bamqc:\n"
            "    output: odir=directory('source/{sample}/qualimap')\n"
            "    shell: 'echo synthetic-input-missing >&2; exit 99'\n"
            "rule picard_alignment_summary_metrics:\n"
            "    output: 'source/{sample}/alignment.metrics'\n"
            "    shell: 'echo synthetic-input-missing >&2; exit 99'\n"
            "rule picard_wgs_metrics:\n"
            "    output: 'source/{sample}/wgs.metrics'\n"
            "    shell: 'echo synthetic-input-missing >&2; exit 99'\n"
            "rule picard_mark_duplicates:\n"
            "    output: metrics='source/{sample}/duplicates.metrics', "
            "bam='source/{sample}/aligned.bam', bai='source/{sample}/aligned.bai'\n"
            "    shell: 'echo synthetic-input-missing >&2; exit 99'\n"
            f"include: {str(ROOT / 'wgs/rules/bfq/bfq_level2.smk')!r}\n"
            "rule qc_test_all:\n"
            "    input: expand(rules.bfq_level2_bam_qc.output, sample=SAMPLES)\n"
        )

    def run_workflow(self, *extra_args):
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", "Snakefile", "--cores", "2",
             "--scheduler", "greedy", "--printshellcmds", "qc_test_all", *extra_args],
            cwd=self.work, capture_output=True, text=True, timeout=60,
        )
        return result.returncode, result.stdout + result.stderr

    def snapshot(self, root):
        return {str(path.relative_to(root)): (path.stat().st_mtime_ns, path.read_bytes() if path.is_file() else None)
                for path in root.rglob("*")}

    def test_distinct_sample_timestamps_deliver_without_mutating_inputs(self):
        original = self.snapshot(self.source)
        code, output = self.run_workflow()
        self.assertEqual(code, 0, output)
        self.assertEqual(self.snapshot(self.source), original)
        for sample in SAMPLES:
            source = self.source / sample
            delivered = self.work / "delivery/logs" / sample
            qualimap = delivered / f"{sample}_qualimap"
            self.assertTrue(qualimap.is_dir())
            self.assertFalse(qualimap.is_symlink())
            for name in ("genome_results.txt", "qualimapReport.html", "raw_data/coverage.txt"):
                self.assertEqual((qualimap / name).read_bytes(), (source / "qualimap" / name).read_bytes())
            for metric, name in zip(METRICS, ("picard_alignment_summary.metrics", "picard_wgs.metrics", "picard_rmdup.metrics")):
                self.assertFalse((delivered / name).is_symlink())
                self.assertEqual((delivered / name).read_bytes(), (source / f"{metric}.metrics").read_bytes())
        delivered_snapshot = self.snapshot(self.work / "delivery")
        code, output = self.run_workflow()
        self.assertEqual(code, 0, output)
        self.assertIn("Nothing to be done", output)
        self.assertEqual(self.snapshot(self.work / "delivery"), delivered_snapshot)
        self.assertEqual(self.snapshot(self.source), original)

    def test_forced_rerun_replaces_legacy_symlinks_without_mutating_inputs(self):
        for sample in SAMPLES:
            source = self.source / sample
            delivered = self.work / "delivery/logs" / sample
            delivered.mkdir(parents=True)
            (delivered / f"{sample}_qualimap").symlink_to(source / "qualimap", target_is_directory=True)
            for metric, name in zip(METRICS, ("picard_alignment_summary.metrics", "picard_wgs.metrics", "picard_rmdup.metrics")):
                (delivered / name).symlink_to(source / f"{metric}.metrics")
        original = self.snapshot(self.source)
        code, output = self.run_workflow("--forcerun", "bfq_level2_bam_qc")
        self.assertEqual(code, 0, output)
        self.assertEqual(self.snapshot(self.source), original)
        for sample in SAMPLES:
            delivered = self.work / "delivery/logs" / sample
            for path in delivered.iterdir():
                self.assertFalse(path.is_symlink(), path)
            self.assertEqual((delivered / f"{sample}_qualimap/genome_results.txt").read_bytes(),
                             (self.source / sample / "qualimap/genome_results.txt").read_bytes())

    def test_incomplete_qualimap_report_fails(self):
        (self.source / "earlier/qualimap/qualimapReport.html").unlink()
        code, output = self.run_workflow()
        self.assertNotEqual(code, 0, output)
        self.assertIn("Error in rule bfq_level2_bam_qc", output)

    def test_empty_metrics_fail(self):
        (self.source / "earlier/wgs.metrics").write_text("")
        code, output = self.run_workflow()
        self.assertNotEqual(code, 0, output)
        self.assertIn("Error in rule bfq_level2_bam_qc", output)

    def test_missing_metrics_fail(self):
        (self.source / "earlier/wgs.metrics").unlink()
        code, output = self.run_workflow()
        self.assertNotEqual(code, 0, output)
        self.assertIn("synthetic-input-missing", output)


if __name__ == "__main__":
    unittest.main()
