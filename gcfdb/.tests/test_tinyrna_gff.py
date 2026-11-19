"""tinyrna_gff.py (#373): the miRBase GFF3 renamed to the genome's sequence names for tiny-count."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/tinyrna_gff.py'

GFF = ('##gff-version 3\n# genome-build-id: GRCr8\n'
       'chr1\t.\tmiRNA_primary_transcript\t10\t80\t.\t+\t.\tID=MI0000001;Alias=MI0000001;Name=rno-mir-1\n'
       'chr1\t.\tmiRNA\t30\t51\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=rno-miR-1-5p;Derives_from=MI0000001\n'
       'chrUn_JH584304\t.\tmiRNA\t5\t26\t.\t-\t.\tID=MIMAT0000001_1;Alias=MIMAT0000001;Name=rno-miR-1-5p;Derives_from=MI0000002\n'
       '901866519\t.\tmiRNA\t5\t26\t.\t+\t.\tID=MIMAT0000002;Alias=MIMAT0000002;Name=rno-miR-16-5p;Derives_from=MI0000003\n'
       'chr1\t.\tgene\t1\t100\t.\t+\t.\tID=other\n')
CHR_MAP = 'chr1\t1\nchrUn_JH584304\tJH584304.1\n901866519\n'


class TinyrnaGffTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, gff=GFF, chr_map=CHR_MAP):
        (self.dir / 'mirna.gff3').write_text(gff)
        (self.dir / 'chr_map.tsv').write_text(chr_map)
        return subprocess.run([sys.executable, str(SCRIPT), '--mirna-gff', str(self.dir / 'mirna.gff3'), '--chr-map',
                               str(self.dir / 'chr_map.tsv'), '--output', str(self.dir / 'out.gff3')], capture_output=True, text=True)

    def test_rows_are_renamed_filtered_and_recorded(self):
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = (self.dir / 'out.gff3').read_text().splitlines()
        rows = [line.split('\t') for line in lines if not line.startswith('#')]
        self.assertEqual([(r[0], r[2], r[8].split(';')[0]) for r in rows], [
            ('1', 'miRNA_primary_transcript', 'ID=MI0000001'), ('1', 'miRNA', 'ID=MIMAT0000001'),
            ('JH584304.1', 'miRNA', 'ID=MIMAT0000001_1')])
        self.assertEqual(rows[2][3:7], ['5', '26', '.', '-'])  # coordinates, strand and attributes unchanged
        self.assertEqual(lines[0], '##gff-version 3')
        self.assertIn('# rows miRNA 2', lines)
        self.assertIn('# dropped (chromosome removed by the map) 1', lines)
        self.assertTrue(any(line.startswith('# source chr_map.tsv sha256 ') for line in lines))
        self.assertEqual(sorted(p.name for p in self.dir.iterdir()), ['chr_map.tsv', 'mirna.gff3', 'out.gff3'])  # no .part left
        umask = os.umask(0)
        os.umask(umask)
        self.assertEqual((self.dir / 'out.gff3').stat().st_mode & 0o777, 0o666 & ~umask)  # not mkstemp's 0600

    def test_repeated_mature_loci_keep_alias_and_avoid_native_suffixes(self):
        mature = next(line for line in GFF.splitlines(keepends=True) if '\tmiRNA\t' in line)
        repeated = mature.replace('30\t51', '75\t96')
        # The following native _1 must not be displaced by an earlier repeated ID.
        result = self.run_script(GFF.replace(mature, mature + repeated))
        self.assertEqual(result.returncode, 0, result.stderr)
        features = [row.split('\t')[8] for row in (self.dir / 'out.gff3').read_text().splitlines() if '\tmiRNA\t' in row]
        self.assertEqual([attrs.split(';')[:2] for attrs in features],
                         [[f'ID=MIMAT0000001{suffix}', 'Alias=MIMAT0000001'] for suffix in ('', '_2', '_1')])

    def test_concurrent_builds_of_the_same_reference_do_not_collide(self):
        (self.dir / 'mirna.gff3').write_text(GFF)
        (self.dir / 'chr_map.tsv').write_text(CHR_MAP)
        cmd = [sys.executable, str(SCRIPT), '--mirna-gff', str(self.dir / 'mirna.gff3'), '--chr-map', str(self.dir / 'chr_map.tsv'),
               '--output', str(self.dir / 'out.gff3')]
        runs = [subprocess.Popen(cmd, stderr=subprocess.PIPE, text=True) for _ in range(8)]
        for run in runs:
            self.assertEqual(run.wait(), 0, run.stderr.read())
            run.stderr.close()
        self.assertIn('# rows miRNA 2', (self.dir / 'out.gff3').read_text().splitlines())
        self.assertFalse(list(self.dir.glob('*.part')))

    def test_a_chromosome_missing_from_the_map_stops(self):
        result = self.run_script(chr_map='chr1\t1\n901866519\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('chromosome chrUn_JH584304 is not in the chromosome map', result.stderr)
        self.assertFalse((self.dir / 'out.gff3').exists())

    def test_a_gff_without_mirnas_stops(self):
        result = self.run_script(gff='##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('no miRNA rows', result.stderr)


if __name__ == '__main__':
    unittest.main()
