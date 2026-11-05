"""miRBase reference rules (#278): release-specific URLs, load-time checks and species FASTA extraction; no downloads."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/mirbase_species_fasta.py'

# Wrapped like miRBase hairpin.fa: sequences continue over several lines.
FASTA = (
    '>cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop\nUACACUGUGG\nAUCCGG\n'
    '>hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop\nUGGGAUGAGG\nUAGUAGGUUG\nUAUAGUU\n'
    '>hsa2-mir-1 MI9999999 Synthetic prefix lookalike\nUUUU\n'
    '>hsa-mir-U1 MI0000061 Homo sapiens synthetic U in name\r\nACGU\r\n'
)


def extract(text, species):
    return subprocess.run([sys.executable, str(SCRIPT), '--species', species, '-'], input=text,
                          capture_output=True, text=True, timeout=30)


class SpeciesFastaTests(unittest.TestCase):
    def test_keeps_whole_wrapped_records_of_the_species_only(self):
        result = extract(FASTA, 'hsa')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, '>hsa-let-7a-1 MI0000060 Homo sapiens let-7a-1 stem-loop\nTGGGATGAGGTAGTAGGTTGTATAGTT\n'
                                        '>hsa-mir-U1 MI0000061 Homo sapiens synthetic U in name\nACGT\n')

    def test_species_without_records_fails(self):
        result = extract(FASTA, 'rno')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('no rno records', result.stderr)


def dry_run(release, organism='homo_sapiens'):
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        settings = {'read_geometry': [50], 'organism': organism, 'db': {'mirbase': {'release': release}},
                    'docker': {'default': 'example/default:1'}}
        (work / 'Snakefile').write_text(f'config.update({settings!r})\n'
                                        f'include: {str(ROOT / "utils.py")!r}\n'
                                        f'include: {str(ROOT / "gcfdb/mirbase.db")!r}\n')
        env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
        result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', 'mirbase_db_all',
                                 'data/ext/mirbase/{}/release-{}/anno/genes.tsv'.format(organism, release)],
                                cwd=work, env=env, capture_output=True, text=True, timeout=60)
        return result.returncode, result.stdout + result.stderr


class MirbaseRuleTests(unittest.TestCase):
    def test_every_served_release_uses_release_specific_files(self):
        for release in ('20', '21', '22', '23.1'):
            with self.subTest(release=release):
                returncode, output = dry_run(release)
                self.assertEqual(returncode, 0, output)
                for url in (f'https://www.mirbase.org/download_version_files/{release}/hairpin.fa -O-',
                            f'https://www.mirbase.org/download_version_files/{release}/mature.fa -O-',
                            f'https://www.mirbase.org/download_version_genome_files/{release}/hsa.gff3'):
                    self.assertIn(url, output)
                self.assertNotIn('gunzip', output)
                self.assertNotIn('miRNA.str', output)
                self.assertIn(f"miRBase v{release}$'", output)
                self.assertEqual(output.count('mirbase_species_fasta.py --species hsa -'), 2)

    def test_release_that_miRBase_does_not_serve_fails_at_load(self):
        for release in ('19', '23', 'CURRENT'):
            with self.subTest(release=release):
                returncode, output = dry_run(release)
                self.assertNotEqual(returncode, 0)
                self.assertIn(f"db.mirbase.release: miRBase serves releases 20, 21, 22, 23.1, got '{release}'", output)

    def test_organism_without_miRBase_code_fails_at_load(self):
        for organism in ('danio_rerio', 'human'):
            with self.subTest(organism=organism):
                returncode, output = dry_run('23.1', organism)
                self.assertNotEqual(returncode, 0)
                self.assertIn(f"miRBase references support organism homo_sapiens, mus_musculus, rattus_norvegicus, got '{organism}'", output)


if __name__ == '__main__':
    unittest.main()
