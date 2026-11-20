"""miRBase GFF3 correction (#429): reviewed table rows applied, then every other row checked against the sequences."""
import csv
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/mirbase_gff_correct.py'
TABLE = ROOT / 'gcfdb/tables/mirbase_gff_corrections.tsv'
COLUMNS = ['release', 'organism', 'id', 'name', 'old', 'new', 'check', 'reason']

rng = random.Random(429)
H1, H2 = (''.join(rng.choice('ACGT') for _ in range(40)) for _ in range(2))
HAIRPINS = '>hsa-mir-1 MI0000001 synthetic\n{}\n>hsa-mir-2 MI0000002 synthetic\n{}\n'.format(H1, H2)
# 5p of hsa-mir-1 at offset 3, 3p at offset 25 (plus strand, locus chr1:101-140); 5p of hsa-mir-2 at offset 5 (minus
# strand, locus chr2:501-540, so 540 - 5 = 535 is its 5' end).
MATURES = '>hsa-miR-1-5p MIMAT0000001 synthetic\n{}\n>hsa-miR-1-3p MIMAT0000002 synthetic\n{}\n>hsa-miR-2-5p MIMAT0000003 synthetic\n{}\n'.format(
    H1[3:13], H1[25:35], H2[5:15])
HEADER = ['##gff-version 3', '# microRNAs:               miRBase v23.1', '# genome-build-id:         GRCh38']
ROWS = {
    'h1': 'chr1\t.\tmiRNA_primary_transcript\t101\t140\t.\t+\t.\tID=MI0000001;Alias=MI0000001;Name=hsa-mir-1',
    'm1': 'chr1\t.\tmiRNA\t104\t113\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=hsa-miR-1-5p;Derives_from=MI0000001',
    'm2': 'chr1\t.\tmiRNA\t126\t135\t.\t+\t.\tID=MIMAT0000002;Alias=MIMAT0000002;Name=hsa-miR-1-3p;Derives_from=MI0000001',
    'h2': 'chr2\t.\tmiRNA_primary_transcript\t501\t540\t.\t-\t.\tID=MI0000002;Alias=MI0000002;Name=hsa-mir-2',
    'm3': 'chr2\t.\tmiRNA\t526\t535\t.\t-\t.\tID=MIMAT0000003;Alias=MIMAT0000003;Name=hsa-miR-2-5p;Derives_from=MI0000002',
}
NOTE = '# gcf-workflows: '


def gff(*rows, **replace):
    """The GFF3 lines: the header and the named rows, with replace={'m1': row} swapping rows."""
    return HEADER + [replace.get(name, ROWS[name]) for name in rows]


def table_row(row_id, name, old, new, check='sequence', release='23.1', organism='homo_sapiens', reason='test'):
    return [release, organism, row_id, name, old, new, check, reason]


def run(lines, table_rows=(), release='23.1'):
    with tempfile.TemporaryDirectory() as directory:
        d = Path(directory)
        (d / 'in.gff3').write_text('\n'.join(lines) + '\n')
        (d / 'hairpin.fa').write_text(HAIRPINS)
        (d / 'mature.fa').write_text(MATURES)
        with open(d / 'table.tsv', 'w', newline='') as fh:
            fh.write('# synthetic corrections\n')
            writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
            writer.writerow(COLUMNS)
            writer.writerows(table_rows)
        result = subprocess.run([sys.executable, str(SCRIPT), '--gff', str(d / 'in.gff3'), '--hairpins', str(d / 'hairpin.fa'),
                                 '--matures', str(d / 'mature.fa'), '--release', release, '--organism', 'homo_sapiens',
                                 '--table', str(d / 'table.tsv'), '--output', str(d / 'out.gff3'), '--record', str(d / 'record.tsv')],
                                capture_output=True, text=True, timeout=60)
        out = (d / 'out.gff3').read_text().splitlines() if (d / 'out.gff3').exists() else None
        record = (d / 'record.tsv').read_text().splitlines() if (d / 'record.tsv').exists() else None
        return result, out, record


class CorrectionTests(unittest.TestCase):
    def setUp(self):
        for seq, sub in ((H1, H1[3:13]), (H1, H1[25:35]), (H2, H2[5:15])):
            self.assertEqual(seq.count(sub), 1)

    def test_consistent_gff_is_written_unchanged_with_a_note(self):
        lines = gff('h1', 'm1', 'm2', 'h2', 'm3')
        result, out, record = run(lines)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(out[:3], HEADER)
        self.assertTrue(out[3].startswith(NOTE + '0 rows corrected and 0 exempt'), out[3])
        self.assertEqual(out[4:], lines[3:])
        self.assertEqual(record, ['id\tname\told\tnew\tcheck\treason'])

    def test_listed_correction_moves_the_row_and_is_recorded(self):
        wrong = ROWS['m1'].replace('\t104\t113\t', '\t105\t114\t')
        result, out, record = run(gff('h1', 'm1', 'm2', m1=wrong),
                                  [table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:105-114:+', 'chr1:104-113:+', reason='offset')])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(ROWS['m1'], out)
        self.assertTrue(out[3].startswith(NOTE + '1 rows corrected'), out[3])
        self.assertEqual(record[1], 'MIMAT0000001\thsa-miR-1-5p\tchr1:105-114:+\tchr1:104-113:+\tsequence\toffset')

    def test_unlisted_disagreement_stops_naming_the_row(self):
        wrong = ROWS['m1'].replace('\t104\t113\t', '\t105\t114\t')
        result, out, _ = run(gff('h1', 'm1', 'm2', m1=wrong))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('MIMAT0000001 hsa-miR-1-5p chr1:105-114:+: not where its sequence is', result.stderr)
        self.assertIn('sequence at hairpin offsets [3]', result.stderr)
        self.assertIsNone(out)

    def test_minus_strand_offsets_count_from_the_5_prime_end(self):
        wrong = ROWS['m3'].replace('\t526\t535\t', '\t506\t515\t')
        result, _, _ = run(gff('h2', 'm3', m3=wrong))
        self.assertIn('MIMAT0000003 hsa-miR-2-5p chr2:506-515:-: not where its sequence is', result.stderr)
        result, _, _ = run(gff('h2', 'm3'))
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_table_row_at_another_position_stops(self):
        wrong = ROWS['m1'].replace('\t104\t113\t', '\t105\t114\t')
        result, _, _ = run(gff('h1', 'm1', m1=wrong),
                           [table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:106-115:+', 'chr1:104-113:+')])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('no GFF row with this ID and position; the table is stale or the download changed', result.stderr)
        result, _, _ = run(gff('h1', 'm1', m1=wrong),
                           [table_row('MIMAT0000001', 'hsa-miR-1-3p', 'chr1:105-114:+', 'chr1:104-113:+')])
        self.assertIn('the GFF row is named hsa-miR-1-5p, not hsa-miR-1-3p; the table is stale', result.stderr)

    def test_exempt_row_is_left_unchecked(self):
        wrong = ROWS['m1'].replace('\t104\t113\t', '\t105\t114\t')
        result, out, record = run(gff('h1', 'm1', m1=wrong),
                                  [table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:105-114:+', 'chr1:105-114:+', 'exempt', reason='paralog')])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(wrong, out)
        self.assertTrue(out[3].startswith(NOTE + '0 rows corrected and 1 exempt'), out[3])

    def test_hairpin_without_strand_or_with_another_span_stops(self):
        result, _, _ = run(gff('h1', 'h2', h1=ROWS['h1'].replace('\t+\t', '\t.\t'), h2=ROWS['h2'].replace('\t501\t', '\t500\t')))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('MI0000001 hsa-mir-1 chr1:101-140:.: no strand', result.stderr)
        self.assertIn('MI0000002 hsa-mir-2 chr2:500-540:-: spans 41 nt, its sequence has 40', result.stderr)

    def test_mature_outside_its_hairpin_stops(self):
        result, _, _ = run(gff('h1', 'm2', m2=ROWS['m2'].replace('\t126\t135\t', '\t136\t145\t')))
        self.assertIn('MIMAT0000002 hsa-miR-1-3p chr1:136-145:+: outside every locus of its hairpin MI0000001', result.stderr)

    def test_mixed_chromosome_names_stop(self):
        result, _, _ = run(gff('h1', 'm1', 'h2', 'm3', h2=ROWS['h2'].replace('chr2\t', '2\t', 1), m3=ROWS['m3'].replace('chr2\t', '2\t', 1)))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('chromosome names mixed: 2 rows with chr, 2 without', result.stderr)

    def test_mature_at_a_further_locus_of_its_hairpin(self):
        # miRBase names a hairpin's further loci MI0000001_2, ...; their matures derive from the bare accession.
        locus2 = 'chr3\t.\tmiRNA_primary_transcript\t1001\t1040\t.\t+\t.\tID=MI0000001_2;Alias=MI0000001;Name=hsa-mir-1'
        mature2 = 'chr3\t.\tmiRNA\t1004\t1013\t.\t+\t.\tID=MIMAT0000001_1;Alias=MIMAT0000001;Name=hsa-miR-1-5p;Derives_from=MI0000001'
        result, _, _ = run(gff('h1', 'm1') + [locus2, mature2])
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rows_sharing_an_id_are_told_apart_by_position(self):
        # 23.1 gives a mature the same ID at every hairpin it comes from.
        locus2 = 'chr3\t.\tmiRNA_primary_transcript\t1001\t1040\t.\t+\t.\tID=MI0000001_2;Alias=MI0000001;Name=hsa-mir-1'
        wrong2 = 'chr3\t.\tmiRNA\t1005\t1014\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=hsa-miR-1-5p;Derives_from=MI0000001'
        result, out, _ = run(gff('h1', 'm1') + [locus2, wrong2],
                             [table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr3:1005-1014:+', 'chr3:1004-1013:+')])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(ROWS['m1'], out)
        self.assertIn(wrong2.replace('\t1005\t1014\t', '\t1004\t1013\t'), out)

    def test_rows_of_another_release_or_organism_are_not_applied(self):
        wrong = ROWS['m1'].replace('\t104\t113\t', '\t105\t114\t')
        for other in ({'release': '22'}, {'organism': 'mus_musculus'}):
            with self.subTest(**other):
                result, _, _ = run(gff('h1', 'm1', m1=wrong),
                                   [table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:105-114:+', 'chr1:104-113:+', **other)])
                self.assertIn('not where its sequence is', result.stderr)

    def test_table_without_reason_or_with_another_check_stops(self):
        for row in (table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:104-113:+', 'chr1:104-113:+', reason=''),
                    table_row('MIMAT0000001', 'hsa-miR-1-5p', 'chr1:104-113:+', 'chr1:104-113:+', check='skip')):
            result, _, _ = run(gff('h1', 'm1'), [row])
            self.assertIn('check must be sequence or exempt, with a reason', result.stderr)


class ShippedTableTests(unittest.TestCase):
    """The reviewed table: well formed, one row per GFF row, and a sequence row always moves its GFF row."""

    def test_table_is_well_formed(self):
        with open(TABLE, newline='') as fh:
            rows = list(csv.DictReader((line for line in fh if not line.startswith('#')), delimiter='\t'))
        self.assertTrue(rows)
        self.assertEqual(list(rows[0]), COLUMNS)
        keys = set()
        for row in rows:
            with self.subTest(row=row['id'], old=row['old']):
                self.assertIn(row['release'], ('20', '21', '22', '23.1'))
                self.assertIn(row['organism'], ('homo_sapiens', 'mus_musculus', 'rattus_norvegicus'))
                self.assertIn(row['check'], ('sequence', 'exempt'))
                self.assertTrue(row['reason'])
                for position in (row['old'], row['new']):
                    chrom, span, strand = position.rsplit(':', 2)
                    start, end = (int(x) for x in span.split('-'))
                    self.assertTrue(chrom and 0 < start <= end and strand in '+-.', position)
                if row['check'] == 'sequence':
                    self.assertNotEqual(row['old'], row['new'])
                    self.assertIn(row['new'].rsplit(':', 1)[1], '+-')
                key = (row['release'], row['organism'], row['id'], row['old'])
                self.assertNotIn(key, keys)
                keys.add(key)


if __name__ == '__main__':
    unittest.main()
