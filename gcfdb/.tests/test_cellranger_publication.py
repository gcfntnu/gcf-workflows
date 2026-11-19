"""Exercise Cell Ranger index publication across workdirs using only synthetic mkref output."""

import gzip
import json
import os
import shlex
import subprocess
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import build_cellranger_reference as publisher


class CellRangerPublicationTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.source = self.root / "panel"
        self.cache = self.source / "cellranger/builds"
        self.fasta = self.source / "fasta/genome.fa"
        self.gtf = self.source / "anno/genes.gtf"
        self.fasta.parent.mkdir(parents=True)
        self.gtf.parent.mkdir(parents=True)
        self.fasta.write_text(">syn\nACGT\n")
        self.gtf.write_text("# synthetic annotation\n")

    def command(self, genome):
        return ["cellranger", "mkref", "--fasta", str(self.fasta), "--genes", str(self.gtf), "--genome", genome]

    def publish(self, workdir, genome="panel", image="synthetic:v1"):
        receipt = self.root / workdir / "ready"
        publisher.build_reference(self.source, self.cache, receipt, None, image, self.command(genome))
        return receipt

    def stub_mkref(self, command, *, check):
        self.assertEqual(command[:2], ["cellranger", "mkref"])
        self.assertTrue(check)
        directory = Path(command[command.index("--output-dir") + 1])
        for child in ("fasta", "genes", "star"):
            (directory / child).mkdir(parents=True)
        (directory / "reference.json").write_text(json.dumps({"genome": command[command.index("--genome") + 1]}))
        (directory / "fasta/genome.fa").write_bytes(self.fasta.read_bytes())
        (directory / "genes/genes.gtf.gz").write_bytes(gzip.compress(self.gtf.read_bytes()))
        # Distinguish complete snapshots built concurrently with identical inputs.
        for member in ("Genome", "SA", "SAindex"):
            (directory / "star" / member).write_text(directory.parent.name)

    def test_two_publishers_pin_one_complete_generation_through_forced_reruns(self):
        both_staged = threading.Barrier(2)

        def concurrent_mkref(command, *, check):
            self.stub_mkref(command, check=check)
            # Both jobs have already observed a cache miss before either can publish.
            both_staged.wait(timeout=20)

        with patch.object(publisher.subprocess, "run", side_effect=concurrent_mkref) as mkref:
            with ThreadPoolExecutor(max_workers=2) as pool:
                receipts = list(pool.map(self.publish, ("work-a", "work-b")))
        self.assertEqual(mkref.call_count, 2)
        manifest = receipts[0].resolve(strict=True)
        self.assertEqual(receipts[1].resolve(strict=True), manifest)
        star = manifest.parent / "star"
        self.assertEqual(len({(star / member).read_text() for member in ("Genome", "SA", "SAindex")}), 1)
        before = manifest.stat()
        self.assertEqual(list(self.cache.rglob(publisher.MANIFEST)), [manifest])

        # Exercise Snakemake's real pre-job output cleanup. Warm reuse must never run Cell Ranger.
        binary = self.root / "bin/cellranger"
        binary.parent.mkdir()
        binary.write_text("#!/bin/sh\nexit 93\n")
        binary.chmod(0o755)
        env = dict(os.environ, PATH=str(binary.parent) + os.pathsep + os.environ["PATH"])
        for receipt in receipts:
            args = [sys.executable, str(SCRIPTS / "build_cellranger_reference.py"), "--source-dir", str(self.source),
                    "--cache", str(self.cache), "--receipt", "ready", "--image", "synthetic:v1", "--", *self.command("panel")]
            (receipt.parent / "Snakefile").write_text(f'rule pin:\n    output: "ready"\n    shell: {shlex.join(args)!r}\n')
            result = subprocess.run([sys.executable, "-m", "snakemake", "--cores", "1", "--forceall", "pin"],
                                    cwd=receipt.parent, env=env, capture_output=True, text=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(all(r.resolve(strict=True) == manifest for r in receipts))
            current = manifest.stat()
            self.assertEqual((current.st_ino, current.st_mtime_ns), (before.st_ino, before.st_mtime_ns))

    def test_failed_and_distinct_builds_preserve_a_readers_view(self):
        with patch.object(publisher.subprocess, "run", side_effect=self.stub_mkref):
            pinned = self.publish("active")
        index = pinned.parent / "index"
        original = index.resolve(strict=True)
        with (index / "star/SA").open() as reader:
            marker = reader.read()

            def incomplete_mkref(command, *, check):
                self.stub_mkref(command, check=check)
                Path(command[command.index("--output-dir") + 1], "star/SA").unlink()

            with patch.object(publisher.subprocess, "run", side_effect=incomplete_mkref):
                with self.assertRaisesRegex(ValueError, "star/SA"):
                    self.publish("candidate", genome="panel_augmented")
            self.assertFalse((self.root / "candidate/ready").is_symlink())
            self.assertEqual(list(self.cache.rglob(publisher.MANIFEST)), [pinned.resolve(strict=True)])

            with patch.object(publisher.subprocess, "run", side_effect=self.stub_mkref):
                renamed = self.publish("candidate", genome="panel_augmented").resolve(strict=True)
                upgraded = self.publish("upgraded", image="synthetic:v2").resolve(strict=True)
                self.fasta.write_text(">syn\nACGTACGT\n")
                augmented = self.publish("new-sequence").resolve(strict=True)
            self.assertEqual(len({original, renamed.parent, upgraded.parent, augmented.parent}), 4)
            self.assertEqual((original / "fasta/genome.fa").read_text(), ">syn\nACGT\n")
            self.assertEqual((augmented.parent / "fasta/genome.fa").read_text(), self.fasta.read_text())
            self.assertEqual(json.loads(renamed.read_text())["genome"], "panel_augmented")
            self.assertEqual(json.loads(upgraded.read_text())["image"], "synthetic:v2")
            self.assertEqual(index.resolve(strict=True), original)
            self.assertEqual((index / "star/SA").read_text(), marker)
            reader.seek(0)
            self.assertEqual(reader.read(), marker)

    def test_damaged_cache_is_rejected_before_pinning(self):
        with patch.object(publisher.subprocess, "run", side_effect=self.stub_mkref):
            original = self.publish("original").resolve(strict=True)
        (original.parent / "star/SAindex").unlink()
        with patch.object(publisher.subprocess, "run", side_effect=AssertionError("damaged shared index rebuilt in place")):
            with self.assertRaisesRegex(ValueError, "star/SAindex"):
                self.publish("fresh")
        self.assertFalse((self.root / "fresh/ready").is_symlink())


if __name__ == "__main__":
    unittest.main()
