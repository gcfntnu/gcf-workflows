"""tesmall_rnacentral.py (#378): RNAcentral's snoRNA, scaRNA, Y RNA, vault RNA and 7SK entries added to a TEsmall
database's structural RNA class, in a database derived from the base one; synthetic files."""
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/tesmall_rnacentral.py'
sys.path.insert(0, str(ROOT / 'gcfdb/scripts'))
import tesmall_rnacentral  # noqa: E402


def exon(chrom, start, end, strand, urs, rna_type, description='an RNA', databases='RFAM'):
    """A GFF3 noncoding_exon line of RNAcentral (1-based start)."""
    return ('{}\tRNAcentral\tnoncoding_exon\t{}\t{}\t.\t{}\t.\tName={};description={};type={};databases={};ID={}.0:ncRNA_exon1;'
            'Parent={}.0\n').format(chrom, start, end, strand, urs, description.replace(',', '%2C'), rna_type, databases, urs, urs)


GFF = ''.join([
    '##gff-version 3\n',
    '1\tRNAcentral\ttranscript\t101\t180\t.\t+\t.\tName=URS0000000001_10090;type=snoRNA;databases=RFAM;ID=URS0000000001_10090.0\n',
    exon('1', 101, 180, '+', 'URS0000000001_10090', 'snoRNA'),                  # new snoRNA
    exon('1', 101, 180, '+', 'URS0000000002_10090', 'snoRNA'),                  # the same exon under another id: one row
    exon('1', 301, 330, '-', 'URS0000000003_10090', 'scaRNA'),                  # overlaps the base snoRNA (other strand)
    exon('chr2', 501, 600, '+', 'URS0000000004_10090', 'Y_RNA'),                # UCSC name renamed
    exon('1', 701, 790, '+', 'URS0000000005_10090', 'vault_RNA'),
    exon('2', 801, 1000, '-', 'URS0000000006_10090', 'snRNA', '7SK RNA'),       # 7SK by description
    exon('2', 201, 300, '+', 'URS0000000007_10090', 'lncRNA', '7SK RNA'),       # a long transcript: not 7SK here
    exon('2', 1, 60, '+', 'URS0000000008_10090', 'snRNA', 'U1 spliceosomal RNA'),   # another snRNA: out of scope
    exon('1', 901, 980, '+', 'URS0000000009_10090', 'snoRNA', 'a snoRNA', 'ENA,miRBase'),  # miRBase lists it: left out
    exon('GL456210.1', 1, 50, '+', 'URS0000000010_10090', 'Y_RNA'),            # not in the genome
    exon('1', 401, 470, '+', 'URS0000000011_10090', 'piRNA'),
])
BASE_STRUCTURAL = '1\t310\t340\tsnoRNA:Snord1\t0\t+\n1\t5\t60\ttRNA:tRNA-Ala-1\t0\t+\n'
MIRNA = ('1\t120\t142\tmmu-miR-1-5p\t0\t+\n'      # inside the new snoRNA
         '1\t170\t192\tmmu-miR-2-3p\t0\t-\n'      # only half inside: not counted
         '2\t880\t902\tmmu-miR-3-5p\t0\t+\n'      # inside the 7SK row on the other strand: counted
         '1\t315\t337\tmmu-miR-4-5p\t0\t+\n')     # inside the base snoRNA only


class TesmallRnacentralTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.base = self.tmp / 'base/genomes/gcf'
        (self.base / 'annotation').mkdir(parents=True)
        (self.base / 'sequence/bowtie_index').mkdir(parents=True)
        (self.base / 'annotation/structural_RNA.bed').write_text(BASE_STRUCTURAL)
        (self.base / 'annotation/miRNA.bed').write_text(MIRNA)
        (self.base / 'annotation/TE.bed').write_text('1\t1\t10\tL1:L1:LINE\t0\t+\n')
        shared = self.tmp / 'genome/genome.fa'
        shared.parent.mkdir()
        shared.write_text('>1\nACGT\n')
        (self.base / 'sequence/genome.fa').symlink_to('../../../../genome/genome.fa')   # a link, as tesmall_genome_links makes
        (self.base / 'sequence/bowtie_index/rDNA.1.ebwt').write_text('index')
        self.fai = self.tmp / 'genome.fa.fai'
        self.fai.write_text('1\t2000\t3\t60\t61\n2\t2000\t2040\t60\t61\n')
        self.gff = self.tmp / 'mus_musculus.GRCm39.gff3.gz'
        with gzip.open(self.gff, 'wt') as fh:
            fh.write(GFF)
        self.record = self.tmp / 'base/tesmall_reference.tsv'
        self.record.write_text('section\tkey\tvalue\nfeatures\tstructural_RNA\t2\nfeatures\tmiRNA\t4\nsource\tname\tgcf\n')

    def run_script(self):
        out = self.tmp / 'derived/genomes/gcf'
        record = self.tmp / 'derived/tesmall_reference.tsv'
        record.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(SCRIPT), '--base', str(self.base), '--base-record', str(self.record), '--rnacentral', str(self.gff),
                        '--genome-fai', str(self.fai), '--out', str(out), '--record', str(record), '--source', 'rnacentral=28.0'], check=True)
        rows = dict(line.rsplit('\t', 1) for line in record.read_text().splitlines()[1:])
        return out, rows

    def test_families_selection_and_naming(self):
        rows, counts = tesmall_rnacentral.read_rnacentral(self.gff, {'1', '2'})
        self.assertEqual(rows, {('1', 100, 180, 'snoRNA', '+'): 'URS0000000001_10090', ('1', 300, 330, 'scaRNA', '-'): 'URS0000000003_10090',
                                ('2', 500, 600, 'Y_RNA', '+'): 'URS0000000004_10090', ('1', 700, 790, 'vault_RNA', '+'): 'URS0000000005_10090',
                                ('2', 800, 1000, '7SK', '-'): 'URS0000000006_10090'})
        self.assertEqual(counts[('snoRNA', 'exons')], 3)
        self.assertEqual(counts[('snoRNA', 'left out: miRBase')], 1)
        self.assertEqual(counts[('Y_RNA', 'dropped: not in the genome')], 1)

    def test_rows_covered_by_the_base_class_are_not_added(self):
        rows, _ = tesmall_rnacentral.read_rnacentral(self.gff, {'1', '2'})
        added, counts = tesmall_rnacentral.added_rows(rows, tesmall_rnacentral.read_bed(self.base / 'annotation/structural_RNA.bed'))
        self.assertNotIn('scaRNA:URS0000000003_10090', [name for _, _, _, name, _ in added])   # overlaps Snord1 on the other strand
        self.assertEqual(counts[('scaRNA', 'covered by the base structural RNA')], 1)
        self.assertEqual(len(added), 4)

    def test_mirnas_inside_added_rows_on_either_strand(self):
        added = [('1', 100, 180, 'snoRNA:URS0000000001_10090', '+'), ('2', 800, 1000, '7SK:URS0000000006_10090', '-')]
        self.assertEqual(tesmall_rnacentral.mirnas_inside(tesmall_rnacentral.read_bed(self.base / 'annotation/miRNA.bed'), added),
                         ['mmu-miR-1-5p', 'mmu-miR-3-5p'])

    def test_derived_database_links_the_base_and_records_the_additions(self):
        out, record = self.run_script()
        structural = (out / 'annotation/structural_RNA.bed').read_text().splitlines()
        self.assertEqual(len(structural), 6)   # two base rows, four added
        self.assertIn('2\t800\t1000\t7SK:URS0000000006_10090\t0\t-', structural)
        self.assertFalse((out / 'annotation/structural_RNA.bed').is_symlink())
        for name in ('annotation/miRNA.bed', 'annotation/TE.bed', 'sequence/genome.fa', 'sequence/bowtie_index/rDNA.1.ebwt'):
            self.assertTrue((out / name).is_symlink(), name)
            self.assertFalse(Path((out / name).readlink()).is_absolute(), name)
            self.assertEqual((out / name).read_text(), (self.base / name).read_text(), name)
        self.assertEqual((self.base / 'annotation/structural_RNA.bed').read_text(), BASE_STRUCTURAL)   # the base is untouched
        self.assertEqual(record['features\tstructural_RNA'], '6')
        self.assertEqual(record['features\tmiRNA'], '4')
        self.assertEqual(record['RNAcentral 7SK\tadded'], '1')
        self.assertEqual(record['RNAcentral snoRNA\tleft out: miRBase'], '1')
        self.assertEqual(record['RNAcentral\tmature miRNAs at least 90% inside an added row'], '2')
        self.assertEqual(record['RNAcentral\tthose miRNAs (the first 20)'], 'mmu-miR-1-5p,mmu-miR-3-5p')
        self.assertEqual(len(record['sha256\tmus_musculus.GRCm39.gff3.gz']), 64)
        self.assertEqual(record['source\trnacentral'], '28.0')
        self.run_script()   # again over its own links: replaced, not stacked
        self.assertEqual((out / 'annotation/TE.bed').read_text(), (self.base / 'annotation/TE.bed').read_text())


if __name__ == '__main__':
    unittest.main()
