"""tesmall_db.py: TEsmall database classes from miRBase, Ensembl, GtRNAdb, RepeatMasker and piRNA clusters; synthetic files."""
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'gcfdb/scripts'))
import tesmall_db  # noqa: E402

NAMES = {'1': 1000, '2': 1000, 'MT': 100}


class TesmallDbTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))

    def write(self, name, text):
        path = self.tmp / name
        path.write_text(text)
        return path

    def test_mirbase_features_use_the_chromosome_map(self):
        gff = self.write('mirna.gff3', '##gff-version 3\n'
                         'chr1\t.\tmiRNA_primary_transcript\t11\t80\t.\t+\t.\tID=MI0000001;Name=hsa-mir-1\n'
                         'chr1\t.\tmiRNA\t20\t41\t.\t+\t.\tID=MIMAT0000001;Name=hsa-miR-1-5p;Derives_from=MI0000001\n'
                         'chrUn_X\t.\tmiRNA\t5\t26\t.\t-\t.\tID=MIMAT0000002;Name=hsa-miR-2\n')
        mature, hairpin = tesmall_db.mirna_beds(gff, {'chr1': '1', 'chrUn_X': None})
        self.assertEqual(mature, [('1', 19, 41, 'hsa-miR-1-5p', '+')])
        self.assertEqual(hairpin, [('1', 10, 80, 'hsa-mir-1', '+')])
        with self.assertRaises(SystemExit):
            tesmall_db.mirna_beds(gff, {'chr1': '1'})

    def test_gtf_classes(self):
        line = '{}\tensembl\t{}\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}"; gene_name "{}"; gene_biotype "{}";\n'
        gtf = self.write('genes.gtf', ''.join([
            line.format('1', 'gene', 101, 200, '+', 'G1', 'T1', 'U6', 'snRNA'),
            line.format('1', 'exon', 101, 200, '+', 'G1', 'T1', 'U6', 'snRNA'),
            line.format('MT', 'gene', 1, 50, '+', 'G2', 'T2', 'MT-RNR1', 'Mt_rRNA'),
            line.format('1', 'gene', 301, 700, '-', 'G3', 'T3', 'ABC', 'protein_coding'),
            line.format('1', 'exon', 301, 350, '-', 'G3', 'T3', 'ABC', 'protein_coding'),
            line.format('1', 'exon', 601, 700, '-', 'G3', 'T3', 'ABC', 'protein_coding'),
            line.format('1', 'exon', 401, 420, '+', 'G4', 'T4', 'MIR1', 'miRNA'),
            line.format('KI1', 'exon', 1, 10, '+', 'G5', 'T5', 'X', 'lncRNA')]))
        structural, exons, introns, rdna, dropped = tesmall_db.gtf_beds(gtf, NAMES)
        self.assertEqual(structural, [('1', 100, 200, 'snRNA:U6', '+'), ('MT', 0, 50, 'Mt_rRNA:MT-RNR1', '+')])
        self.assertEqual(exons, [('1', 300, 350, 'ABC', '-'), ('1', 600, 700, 'ABC', '-')])
        self.assertEqual(introns, [('1', 350, 600, 'ABC', '-')])
        self.assertEqual(rdna, [('MT', 0, 50, 'Mt_rRNA:MT-RNR1', '+')])
        self.assertEqual(dropped, 1)

    def test_repeatmasker_keeps_transposon_classes_on_genome_sequences(self):
        row = '0\t1\t0\t0\t0\t{}\t{}\t{}\t-0\t{}\t{}\t{}\t{}\t1\t1\t0\t1\n'
        rmsk = self.write('rmsk.txt', row.format('chr1', 10, 20, '+', 'L1HS', 'LINE', 'L1')
                          + row.format('chr2', 30, 40, '-', 'AluY', 'SINE?', 'Alu')
                          + row.format('chr1', 50, 60, '+', '(CA)n', 'Simple_repeat', 'Simple_repeat')
                          + row.format('chr9_alt', 5, 9, '+', 'L2', 'LINE', 'L2'))
        rows, dropped = tesmall_db.rmsk_rows(rmsk, NAMES)
        self.assertEqual(rows, [('1', 10, 20, 'L1HS:L1:LINE', '+'), ('2', 30, 40, 'AluY:Alu:SINE?', '-')])
        self.assertEqual(dropped, 1)

    def test_gtrnadb_archive_high_confidence_set(self):
        bed = 'chr1\t9\t21\ttRNA-Gly-GCC-1-1\t1000\t+\t9\t21\t0\t1\t12,\t0,\nchr1\t50\t60\ttRX-Arg-NNN-1-1\t0\t+\t50\t60\t0\t1\t10,\t0,\n'
        fa = '>Homo_sapiens_tRNA-Gly-GCC-1-1 (tRNAscan-SE ID: chr1.trna1) Gly (GCC) 12 bp Sc: 70 chr1:10-21 (+)\nacgtacgtacgt\n'
        archive = self.tmp / 'hg38-tRNAs.tar.gz'
        with tarfile.open(archive, 'w:gz') as tar:
            for name, text in (('hg38-tRNAs.bed', bed), ('hg38-tRNAs.fa', fa), ('hg38-mature-tRNAs.fa', '>x\nA\n')):
                info = tarfile.TarInfo(name)
                info.size = len(text.encode())
                tar.addfile(info, io.BytesIO(text.encode()))
        self.assertEqual(tesmall_db.trna_records(archive), [('chr1', 9, 21, 'tRNA-Gly-GCC-1-1', '+', 'ACGTACGTACGT')])

    def test_pirna_cluster_headers(self):
        fa = self.write('piRclDB', '>Hsap_0001_4:101..400\nACGT\n>Rnor_0002_chr1:11..20\nAC\n')
        clusters, declared = tesmall_db.pirna_clusters(fa)
        self.assertEqual([c for c, _ in clusters], ['piRNAclusterDB:Hsap_0001', 'piRNAclusterDB:Rnor_0002'])
        self.assertEqual(declared['piRNAclusterDB:Hsap_0001'], ('4', 100, 400))
        self.assertEqual(declared['piRNAclusterDB:Rnor_0002'], ('chr1', 10, 20))

    def test_clusters_are_placed_from_unique_anchors(self):
        seq = 'A' * 400
        clusters = [('c+', seq), ('c-', seq), ('inner', seq), ('split', seq), ('short', seq), ('none', seq)]
        A = tesmall_db.ANCHOR
        hits = {'c+|s|0': ('1', 100, '+'), 'c+|e|0': ('1', 470, '+'),          # forward: 100..500
                'c-|s|0': ('2', 470, '-'), 'c-|e|0': ('2', 100, '-'),          # reverse: 100..500
                'inner|s|50': ('1', 650, '+'), 'inner|e|100': ('1', 1000 - 100 - A, '+'),          # 600..1000
                'split|s|0': ('1', 10, '+'), 'split|e|0': ('2', 380, '+'),
                'short|s|0': ('1', 10, '+'), 'short|e|0': ('1', 100, '+'),
                'none|s|0': ('1', 10, '+')}
        rows, reasons = tesmall_db.place_clusters(clusters, hits)
        self.assertEqual(rows, [('1', 100, 500, 'c+', '.'), ('2', 100, 500, 'c-', '.'), ('1', 600, 1000, 'inner', '.')])
        self.assertEqual(sum(reasons.values()), 3)

    def test_anchors_from_both_ends(self):
        anchors = dict(tesmall_db.cluster_anchors([('c', 'G' * 30 + 'A' * 100 + 'T' * 30)]))
        self.assertEqual(anchors['c|s|0'], 'G' * 30)
        self.assertEqual(anchors['c|e|0'], 'T' * 30)
        self.assertEqual(anchors['c|s|50'], 'A' * 30)
        self.assertNotIn('c|s|100', anchors)  # 160 nt leave no room for anchors 100 nt in from both ends

    @unittest.skipUnless(shutil.which('samtools') and shutil.which('bowtie') and shutil.which('bowtie-build'),
                         'samtools and bowtie 1 not on PATH')
    def test_build_on_a_synthetic_genome(self):
        import random
        random.seed(3)
        rnd = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
        chr1, chr2 = rnd(3000), rnd(3000)
        genome = self.write('genome.fa', '>1\n{}\n>2\n{}\n'.format(chr1, chr2))
        subprocess.run(['samtools', 'faidx', str(genome)], check=True)
        subprocess.run(['bowtie-build', '-q', str(genome), str(self.tmp / 'genome')], check=True)
        gff = self.write('mirna.gff3', 'chr1\t.\tmiRNA\t101\t122\t.\t+\t.\tID=MIMAT1;Name=syn-miR-1\n'
                         'chr1\t.\tmiRNA_primary_transcript\t91\t160\t.\t+\t.\tID=MI1;Name=syn-mir-1\n')
        chr_map = self.write('chr_map.tsv', 'chr1\t1\n')
        gtf = self.write('genes.gtf', '2\te\tgene\t201\t300\t.\t-\t.\tgene_id "G"; gene_name "RNA5S1"; gene_biotype "rRNA";\n')
        trna = chr2[999:1071]
        archive = self.tmp / 'syn-tRNAs.tar.gz'
        rev = trna[::-1].translate(str.maketrans('ACGT', 'TGCA'))
        files = {'syn-tRNAs.bed': 'chr2\t999\t1071\ttRNA-Gly-GCC-1-1\t1000\t-\t999\t1071\t0\t1\t72,\t0,\n',
                 'syn-tRNAs.fa': '>Synthetic_organism_tRNA-Gly-GCC-1-1 (x) Gly (GCC) 72 bp Sc: 70 chr2:1000-1071 (-)\n{}\n'.format(rev)}
        with tarfile.open(archive, 'w:gz') as tar:
            for name, text in files.items():
                info = tarfile.TarInfo(name)
                info.size = len(text.encode())
                tar.addfile(info, io.BytesIO(text.encode()))
        rmsk = self.write('rmsk.txt', '0\t1\t0\t0\t0\tchr2\t1500\t1800\t-0\t+\tL1HS\tLINE\tL1\t1\t1\t0\t1\n')
        rdna = self.write('rdna.fa', '>U00000.1 synthetic rDNA\nACGTACGT\n')
        clusters = self.write('piRclDB', '>Syn_0001_1:1501..2300\n{}\n'.format(chr1[1500:2300]))
        args = ['--genome', genome, '--genome-index', self.tmp / 'genome', '--gtf', gtf, '--mirna-gff', gff, '--chr-map', chr_map,
                '--trna-archive', archive, '--rmsk', rmsk, '--rdna', rdna, '--pirna-clusters', clusters,
                '--annotation', self.tmp / 'annotation', '--sequence', self.tmp / 'sequence', '--record', self.tmp / 'record.tsv']
        result = subprocess.run([sys.executable, str(ROOT / 'gcfdb/scripts/tesmall_db.py'), *map(str, args)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        annotation = lambda name: (self.tmp / 'annotation' / name).read_text()
        self.assertEqual(annotation('miRNA.bed'), '1\t100\t122\tsyn-miR-1\t0\t+\n')
        self.assertEqual(annotation('piRNA_cluster.bed'), '1\t1500\t2300\tpiRNAclusterDB:Syn_0001\t0\t.\n')
        self.assertIn('2\t999\t1071\ttRNA:tRNA-Gly-GCC-1-1\t0\t-\n', annotation('structural_RNA.bed'))
        self.assertEqual(annotation('TE.bed'), '2\t1500\t1800\tL1HS:L1:LINE\t0\t+\n')
        self.assertEqual(annotation('hairpin.bed'), '1\t90\t160\tsyn-mir-1\t0\t+\n')
        self.assertEqual(annotation('exon.bed') + annotation('intron.bed'), '')
        tdna = (self.tmp / 'sequence' / 'tDNA.fa').read_text().split()
        self.assertEqual(tdna[0], '>tRNA-Gly-GCC-1-1:2:999-1071:-')
        self.assertEqual(''.join(tdna[1:]), rev)
        self.assertEqual((self.tmp / 'sequence' / 'rDNA.fa').read_text().splitlines()[2], '>rRNA:RNA5S1:2:200-300')
        record = (self.tmp / 'record.tsv').read_text()
        self.assertIn('piRNA clusters\tplaced at their published coordinates\t1', record)

        bad = files['syn-tRNAs.fa'].replace(rev, rev[::-1])
        with tarfile.open(archive, 'w:gz') as tar:
            for name, text in (('syn-tRNAs.bed', files['syn-tRNAs.bed']), ('syn-tRNAs.fa', bad)):
                info = tarfile.TarInfo(name)
                info.size = len(text.encode())
                tar.addfile(info, io.BytesIO(text.encode()))
        result = subprocess.run([sys.executable, str(ROOT / 'gcfdb/scripts/tesmall_db.py'), *map(str, args)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('differ from the genome', result.stderr)


if __name__ == '__main__':
    unittest.main()
