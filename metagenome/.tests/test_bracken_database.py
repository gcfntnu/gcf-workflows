"""Dry-run real database/Bracken rules without downloads or sequencing data.

Run with: python -m unittest discover -s metagenome/.tests -p 'test_*.py'
Requires Snakemake and pandas.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class BrackenDatabaseTests(unittest.TestCase):
    def check_database(self, reference, organism, missing_distribution=False):
        assembly = "k2_pluspf" if reference == "langmead" else "16S"
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            db = work / "ext" / reference / "release-20231009" / "metagenome" / assembly
            db.mkdir(parents=True)
            for name in ("hash.k2d", "opts.k2d", "taxo.k2d", "seqid2taxid.map"):
                (db / name).touch()
            distribution = db / "database150mers.kmer_distrib"
            if not missing_distribution:
                distribution.touch()
            report = work / "quant" / "kraken2" / assembly / "13" / "13.kraken.kreport"
            report.parent.mkdir(parents=True)
            report.write_text("100\t10\t10\tS\t123\tExample species\n")
            settings = {
                "db": {"reference_db": reference,
                       reference: {"release": "20231009", "assembly": assembly}},
                "docker": {name: "example/" + name + ":test" for name in
                           ("kraken2", "bracken", "krona", "kraken-biom", "phyloseq")},
                "winecellar": {"url": "https://example.invalid"},
            }
            snakefile = work / "Snakefile"
            snakefile.write_text(
                "from os.path import join, dirname\n"
                f"config.update({settings!r})\n"
                f"EXT_DIR = {str(work / 'ext')!r}\n"
                f"GCFDB_DIR = {str(ROOT / 'gcfdb')!r}\n"
                f"QUANT_INTERIM = {str(work / 'quant')!r}\n"
                f"INTERIM_DIR = {str(work / 'interim')!r}\n"
                f"ORG = {organism!r}\n"
                "PE = True\nSAMPLES = ['13']\nread_geometry = [151, 151]\n"
                "def get_filtered_fastq(wildcards):\n"
                "    return {'R1': 'unused_R1.fastq.gz', 'R2': 'unused_R2.fastq.gz'}\n"
                f"def src_gcf(path):\n    return join({str(ROOT)!r}, path)\n"
                f"include: {str(ROOT / 'gcfdb' / (reference + '.db'))!r}\n"
                f"include: {str(ROOT / 'metagenome/rules/quant/kraken2.smk')!r}\n"
            )
            result = subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile),
                 "--directory", str(work), "--dry-run", "--printshellcmds", "--cores", "1",
                 str(report.with_name("13.bracken_out"))],
                cwd=work, capture_output=True, text=True, timeout=60,
            )
            output = result.stdout + result.stderr
            self.assertEqual(result.returncode, 0, output)
            self.assertIn(str(distribution), output)
            self.assertIn("-d " + str(db), output)
            if missing_distribution:
                producer = "langmead_kraken_prebuild" if reference == "langmead" else "ncbi_16s_bracken"
                self.assertIn("rule " + producer + ":", output)
            else:
                self.assertNotIn("rule langmead_kraken_prebuild:", output)
                self.assertNotIn("rule ncbi_16s_bracken:", output)

    def test_existing_database_independent_of_organism(self):
        for reference in ("langmead", "ncbi_16s"):
            for organism in ("n/a", "metagenome", "homo_sapiens"):
                with self.subTest(reference=reference, organism=organism):
                    self.check_database(reference, organism)

    def test_missing_distribution_resolves_to_producer(self):
        for reference in ("langmead", "ncbi_16s"):
            with self.subTest(reference=reference):
                self.check_database(reference, "n/a", missing_distribution=True)


if __name__ == "__main__":
    unittest.main()
