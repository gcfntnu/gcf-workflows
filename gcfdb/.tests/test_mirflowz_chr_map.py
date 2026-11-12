"""mirflowz_chr_map.py: miRNA GFF chromosome names mapped to genome sequence names; synthetic files."""
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/mirflowz_chr_map.py'


def gff(build, chromosomes):
    header = '##gff-version 3\n# microRNAs:               miRBase v23.1\n'
    header += '# genome-build-id:         {}\n'.format(build) if build else ''
    return header + ''.join('{}\t.\tmiRNA\t10\t31\t.\t+\t.\tID=MIMAT{:07d};Name=x\n'.format(c, i) for i, c in enumerate(chromosomes))


class ChrMapTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.output = self.tmp / 'chr_map.tsv'

    def run_script(self, genome_names, chromosomes, build='GRCh38', extra=(), gzipped=False):
        fasta = ''.join('>{} dna:chromosome\nACGT\n'.format(name) for name in genome_names)
        genome = self.tmp / ('genome.fa.gz' if gzipped else 'genome.fa')
        genome.write_bytes(gzip.compress(fasta.encode()) if gzipped else fasta.encode())
        (self.tmp / 'mirna.gff3').write_text(gff(build, chromosomes))
        result = subprocess.run([sys.executable, str(SCRIPT), '--genome', str(genome), '--mirna-gff', str(self.tmp / 'mirna.gff3'),
                                 '--output', str(self.output), *extra], capture_output=True, text=True, timeout=60)
        rows = [line.split('\t') for line in self.output.read_text().splitlines()] if self.output.exists() else None
        return result, rows

    def test_ucsc_names_map_to_ensembl_names(self):
        result, rows = self.run_script(['1', 'X', 'MT', 'KI270706.1', 'JH584304.1'],
                                       ['chr1', 'chrX', 'chrM', 'chr1_KI270706v1_random', 'chrUn_JH584304'],
                                       extra=('--assembly', 'GRCh38'), gzipped=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [['chr1', '1'], ['chrX', 'X'], ['chrM', 'MT'], ['chr1_KI270706v1_random', 'KI270706.1'],
                                ['chrUn_JH584304', 'JH584304.1']])

    def test_names_already_in_the_genome_map_to_themselves(self):
        result, rows = self.run_script(['1', '10', 'X'], ['1', '10', 'X'], build='GRCm39')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [['1', '1'], ['10', '10'], ['X', 'X']])

    def test_unmapped_chromosome_fails_naming_it(self):
        result, rows = self.run_script(['1'], ['chr1', 'chr2', 'chrUn_GL000220v1'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('chr2, chrUn_GL000220v1', result.stderr)
        self.assertIsNone(rows)

    def test_listed_chromosome_is_removed_and_others_still_checked(self):
        result, rows = self.run_script(['1'], ['chr1', '901866519'], build='Rnor_6.0', extra=('--remove', '901866519'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [['chr1', '1'], ['901866519']])
        result, _ = self.run_script(['1'], ['chr1', 'chr2', '901866519'], build='Rnor_6.0', extra=('--remove', '901866519'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('chr2', result.stderr)

    def test_removal_of_a_chromosome_absent_from_the_gff_fails(self):
        result, _ = self.run_script(['1'], ['chr1'], extra=('--remove', '901866519'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--remove names chromosome(s) not in', result.stderr)

    def test_gff_on_another_assembly_fails(self):
        for build in ('GRCm38', None):
            result, rows = self.run_script(['1'], ['chr1'], build=build, extra=('--assembly', 'GRCm39'))
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('the mirflowz genome is on GRCm39', result.stderr)
            self.assertIsNone(rows)


if __name__ == '__main__':
    unittest.main()
