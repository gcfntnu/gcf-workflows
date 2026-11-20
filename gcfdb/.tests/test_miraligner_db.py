"""miraligner_db.py (#426): the miraligner database, hairpins and the mature miRNAs' positions in them (miRNA.str)."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/miraligner_db.py'

# hsa-mir-1: + strand, 60 nt; hsa-mir-2: - strand, two loci (MI0000002 and MI0000002_1; the second locus's mature rows
# derive from the bare accession), so its positions count from the hairpin's 3' end on the genome.
HAIRPIN = ('>hsa-mir-1 MI0000001 Homo sapiens miR-1 stem-loop\n'
           'AAAAACCCCCTGGAATGTAAAGAAGTATGTATTTTTTGGGGGTACATACTTCTTTACATTCCAGGG\n'
           '>hsa-mir-2 MI0000002 Homo sapiens miR-2 stem-loop\n'
           'GGGGGTAGCAGCACGTAAATATTGGCGTTTTTAAAAACCAATATTTACGTGCTGCTATTTTTT\n'
           '>hsa-mir-3 MI0000003 Homo sapiens miR-3 stem-loop\nACGTACGTACGT\n'
           '>mmu-mir-1 MI0000099 Mus musculus miR-1 stem-loop\nACGT\n')
MATURE = ('>hsa-miR-1-5p MIMAT0000001 Homo sapiens miR-1-5p\nUGGAAUGUAAAGAAGUAUGUAU\n'
          '>hsa-miR-2-5p MIMAT0000002 Homo sapiens miR-2-5p\nUAGCAGCACGUAAAUAUUGGCG\n'
          '>mmu-miR-1-5p MIMAT0000099 Mus musculus miR-1-5p\nACGU\n')


def gff_rows(rows):
    return '##gff-version 3\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows)


# hsa-mir-1 on chr1 1001..1066 (+): miR-1-5p at 1011..1032 is 11-32. hsa-mir-2 on chr2 2001..2063 (-) and chr3 (-):
# miR-2-5p at 2053..2058 would be 6-11 counted from 2063; its sequence is at 6-27, i.e. 2037..2058.
GFF = gff_rows([
    ('chr1', '.', 'miRNA_primary_transcript', 1001, 1066, '.', '+', '.', 'ID=MI0000001;Alias=MI0000001;Name=hsa-mir-1'),
    ('chr1', '.', 'miRNA', 1011, 1032, '.', '+', '.', 'ID=MIMAT0000001;Alias=MIMAT0000001;Name=hsa-miR-1-5p;Derives_from=MI0000001'),
    ('chr2', '.', 'miRNA_primary_transcript', 2001, 2063, '.', '-', '.', 'ID=MI0000002;Alias=MI0000002;Name=hsa-mir-2'),
    ('chr2', '.', 'miRNA', 2037, 2058, '.', '-', '.', 'ID=MIMAT0000002;Alias=MIMAT0000002;Name=hsa-miR-2-5p;Derives_from=MI0000002'),
    ('chr3', '.', 'miRNA_primary_transcript', 3001, 3063, '.', '-', '.', 'ID=MI0000002_1;Alias=MI0000002;Name=hsa-mir-2'),
    ('chr3', '.', 'miRNA', 3037, 3058, '.', '-', '.', 'ID=MIMAT0000002_1;Alias=MIMAT0000002;Name=hsa-miR-2-5p;Derives_from=MI0000002'),
])


class MiralignerDbTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, gff=GFF, hairpin=HAIRPIN, mature=MATURE):
        for name, text in (('hairpin.fa', hairpin), ('mature.fa', mature), ('hsa.gff3', gff)):
            (self.dir / name).write_text(text)
        return subprocess.run([sys.executable, str(SCRIPT), '--species', 'hsa', '--hairpin', str(self.dir / 'hairpin.fa'),
                               '--mature', str(self.dir / 'mature.fa'), '--gff', str(self.dir / 'hsa.gff3'),
                               '--output', str(self.dir / 'db')], capture_output=True, text=True)

    def record(self):
        return dict(line.split('\t', 1) for line in (self.dir / 'db/miraligner_db.tsv').read_text().splitlines()[1:])

    def test_positions_on_both_strands_and_further_loci(self):
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.dir / 'db/miRNA.str').read_text(),
                         '>hsa-mir-1 (0.00)   [hsa-miR-1-5p:11-32]\n>hsa-mir-2 (0.00)   [hsa-miR-2-5p:6-27]\n')
        hairpins = (self.dir / 'db/hairpin.fa').read_text()
        self.assertIn('>hsa-mir-3 MI0000003 Homo sapiens miR-3 stem-loop\nACGTACGTACGT\n', hairpins)
        self.assertNotIn('mmu-', hairpins)
        record = self.record()
        self.assertEqual((record['hairpins'], record['hairpins_with_positions'], record['mature_positions']), ('3', '2', '2'))
        self.assertEqual(record['hairpins_without_gff_rows'], 'hsa-mir-3')
        self.assertEqual((record['positions_from_sequence'], record['matures_left_out']), ('none', 'none'))
        self.assertEqual(len(record['gff_sha256']), 64)

    def test_a_gff_offset_that_misses_the_sequence_takes_the_sequence_position(self):
        """miRBase 23.1: mmu-miR-26a-5p one nt short at mmu-mir-26a-2 in the GFF."""
        result = self.run_script(GFF.replace('1011\t1032', '1011\t1031'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('[hsa-miR-1-5p:11-32]', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['positions_from_sequence'], 'hsa-miR-1-5p@hsa-mir-1:11-32')

    def test_a_locus_whose_sequence_differs_keeps_the_gff_span_or_leaves_the_mature_out(self):
        """miRBase 23.1: mmu-miR-19b-3p, 24 nt, is 23 nt at mmu-mir-19b-1."""
        longer = MATURE.replace('UGGAAUGUAAAGAAGUAUGUAU', 'UGGAAUGUAAAGAAGUAUGUAUC')  # the hairpin has ...UGUAUU
        result = self.run_script(mature=longer)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('[hsa-miR-1-5p:11-32]', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['positions_from_gff_sequence_differs'], 'hsa-miR-1-5p@hsa-mir-1:11-32')
        other = MATURE.replace('UGGAAUGUAAAGAAGUAUGUAU', 'CCCCCCCCCCCCCCCCCCCCCC')
        result = self.run_script(mature=other)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('hsa-mir-1 ', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['matures_left_out'], 'hsa-miR-1-5p@hsa-mir-1:GFF 1011..1032')

    def test_a_gff_span_beyond_the_hairpin_sequence_is_not_kept(self):
        """A GFF locus longer than the hairpin sequence: its span is cut short there, so it cannot stand as a position."""
        gff = GFF.replace('1001\t1066', '1001\t1080').replace('1011\t1032', '1050\t1071')
        # not in the hairpin, but starts with its last 17 nt (from 50), so the cut-short span would pass as a prefix
        short = MATURE.replace('UGGAAUGUAAAGAAGUAUGUAU', 'UUCUUUACAUUCCAGGGAA')
        result = self.run_script(gff, mature=short)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('hsa-mir-1 ', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['matures_left_out'], 'hsa-miR-1-5p@hsa-mir-1:GFF 1050..1071')

    def test_one_note_per_hairpin_and_mature(self):
        """A mature miRNA placed from its GFF offset at one locus is not listed as left out because a row elsewhere could
        not be placed (no locus of its hairpin there, and its sequence twice in the hairpin)."""
        mature = 'UAGCAGCACGUAAAUAUUGGCG'.replace('U', 'T')
        hairpin = HAIRPIN + '>hsa-mir-5 MI0000005 Homo sapiens miR-5 stem-loop\n' + mature + 'AAAAA' + mature + '\n'
        mature_fa = MATURE + '>hsa-miR-5-5p MIMAT0000005 Homo sapiens miR-5-5p\n' + mature + '\n'
        gff = GFF + gff_rows([
            ('chr5', '.', 'miRNA_primary_transcript', 1, 49, '.', '+', '.', 'ID=MI0000005;Alias=MI0000005;Name=hsa-mir-5'),
            ('chr5', '.', 'miRNA', 1, 22, '.', '+', '.', 'ID=MIMAT0000005;Alias=MIMAT0000005;Name=hsa-miR-5-5p;Derives_from=MI0000005'),
            ('chr6', '.', 'miRNA', 1, 22, '.', '+', '.', 'ID=MIMAT0000005_1;Alias=MIMAT0000005;Name=hsa-miR-5-5p;Derives_from=MI0000005'),
        ]).replace('##gff-version 3\n', '')
        result = self.run_script(gff, hairpin=hairpin, mature=mature_fa)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('>hsa-mir-5 (0.00)   [hsa-miR-5-5p:1-22]', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['matures_left_out'], 'none')

    def test_rows_without_a_strand_or_locus_are_placed_by_their_sequence(self):
        """miRBase 23.1: hsa-mir-13168 has no strand and a mature row outside it."""
        gff = GFF.replace('1001\t1066\t.\t+', '1001\t1066\t.\t.').replace('1011\t1032\t.\t+', '1070\t1091\t.\t.')
        result = self.run_script(gff)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('[hsa-miR-1-5p:11-32]', (self.dir / 'db/miRNA.str').read_text())
        self.assertEqual(self.record()['positions_from_sequence'], 'hsa-miR-1-5p@hsa-mir-1:11-32')

    def test_missing_records_stop(self):
        for kwargs, message in ((dict(mature=MATURE.replace('hsa-miR-1-5p MIMAT0000001', 'hsa-miR-1-5p MIMAT0000777')),
                                 'hsa-miR-1-5p (MIMAT0000001) is not in the mature FASTA'),
                                (dict(hairpin=HAIRPIN.replace('>hsa-mir-1 ', '>hsa-mir-9 ')), 'hairpin hsa-mir-1 of mature miRNA'),
                                (dict(gff=GFF.replace('Derives_from=MI0000001', 'Derives_from=MI0000777')), 'derives from MI0000777'),
                                (dict(hairpin=HAIRPIN + '>hsa-mir-3 MI0000003 again\nACGT\n'), 'hsa-mir-3 twice')):
            with self.subTest(message=message):
                result = self.run_script(**kwargs)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)


if __name__ == '__main__':
    unittest.main()
