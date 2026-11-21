"""miRBase reference rules (#278): release-specific URLs, load-time checks and species FASTA extraction; no downloads."""
import hashlib
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


GENES = ROOT / 'gcfdb/scripts/mirbase_gff2genes.py'
GFF = ('##gff-version 3\n# microRNAs:               miRBase v23.1\n'
       '1\t.\tmiRNA_primary_transcript\t100\t180\t.\t+\t.\tID=MI0000001;Alias=MI0000001;Name=mmu-mir-1a\n'
       '1\t.\tmiRNA\t110\t131\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=mmu-miR-1a-5p;Derives_from=MI0000001\n'
       '2\t.\tmiRNA_primary_transcript\t500\t580\t.\t-\t.\tID=MI0000002;Alias=MI0000002;Name=mmu-mir-1a-2\n'
       '2\t.\tmiRNA\t510\t531\t.\t-\t.\tID=MIMAT0000002;Alias=MIMAT0000002;Name=mmu-miR-1a-5p;Derives_from=MI0000002\n')


class MirbaseGenesTests(unittest.TestCase):
    """mirbase_gff2genes.py (#283): runs on current pandas and strips the species prefix, as UNITAS names miRNAs."""

    def test_mature_mirnas_without_species_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            gff = Path(directory) / 'mmu.gff3'
            gff.write_text(GFF)
            result = subprocess.run([sys.executable, '-W', 'error::FutureWarning', str(GENES), str(gff)],
                                    capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0].split('\t'), ['mirna_id', 'seqid', 'source', 'feature', 'start', 'end', 'score', 'strand',
                                                'phase', 'id', 'alias', 'derives_from'])
        rows = [line.split('\t') for line in lines[1:]]
        self.assertEqual([(row[0], row[9]) for row in rows], [('miR-1a-5p', 'MIMAT0000001'), ('miR-1a-1-5p', 'MIMAT0000002')])

def dry_run(release, organism='homo_sapiens', work=None):
    """A dry run of mirbase_db_all and mirbase_genes in a new project: `work`, or an empty temporary folder."""
    if work is None:
        with tempfile.TemporaryDirectory() as directory:
            return dry_run(release, organism, Path(directory))
    settings = {'read_geometry': [50], 'organism': organism, 'db': {'mirbase': {'release': release}},
                'docker': {'default': 'example/default:1'}}
    (work / 'Snakefile').write_text(f'config.update({settings!r})\n'
                                    f'include: {str(ROOT / "utils.py")!r}\n'
                                    f'include: {str(ROOT / "gcfdb/mirbase.db")!r}\n')
    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', 'mirbase_db_all', 'mirbase_genes'],
                            cwd=work, env=env, capture_output=True, text=True, timeout=60)
    return result.returncode, result.stdout + result.stderr


def corrected_id():
    """MIRBASE_GFF_ID (gcfdb/mirbase.db), which names the corrected GFF3 and its record after the table and the script."""
    digest = hashlib.sha256()
    for name in ('gcfdb/tables/mirbase_gff_corrections.tsv', 'gcfdb/scripts/mirbase_gff_correct.py'):
        digest.update(hashlib.sha256((ROOT / name).read_bytes()).digest())
    return digest.hexdigest()[:12]


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
                # The GFF3 consumers read is miRBase's download corrected and checked (#429).
                self.assertIn(f'release-{release}/sources/hsa.gff3', output)
                self.assertIn('scripts/mirbase_gff_correct.py --gff ', output)
                self.assertIn(f'--release {release} --organism homo_sapiens --table ', output)
                self.assertIn('gcfdb/tables/mirbase_gff_corrections.tsv', output)

    def test_a_newer_correction_table_does_not_rebuild_the_corrected_gff3(self):
        # A checkout gives the table the checkout time, later than a corrected GFF3 built before in a shared EXT_DIR;
        # a new project has no metadata to compare the table's content with (#438).
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            ext = work / 'data/ext/mirbase/homo_sapiens/release-23.1'
            gff = ext / 'anno' / f'genes.{corrected_id()}.gff3'
            built = [ext / 'fasta/hairpin.fa', ext / 'fasta/mature.fa', ext / 'sources/hsa.gff3', gff,
                     ext / f'logs/gff_corrections.{corrected_id()}.tsv', ext / f'anno/genes.{corrected_id()}.tsv']
            before = (ROOT / 'gcfdb/tables/mirbase_gff_corrections.tsv').stat().st_mtime - 86400
            for n, path in enumerate(built):  # in build order, all before the table's file time
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('x\n')
                os.utime(path, (before + n, before + n))
            returncode, output = dry_run('23.1', work=work)
            self.assertEqual(returncode, 0, output)
            self.assertNotIn('rule mirbase_gff:', output)
            self.assertIn('Nothing to be done', output)
            # A store with other corrections, or none (genes.gff3 before #429), has no file of this name: the GFF3 is
            # corrected, and its consumers rebuild after it.
            gff.rename(ext / 'anno/genes.gff3')
            returncode, output = dry_run('23.1', work=work)
            self.assertEqual(returncode, 0, output)
            self.assertIn('rule mirbase_gff:', output)
            self.assertIn('rule mirbase_genes:', output)

    def test_the_corrected_gff3_and_its_record_are_named_after_the_table_and_script(self):
        # Any change to either gives new files, so no build with other corrections is taken for this one (#438).
        returncode, output = dry_run('23.1')
        self.assertEqual(returncode, 0, output)
        self.assertIn(f'release-23.1/anno/genes.{corrected_id()}.gff3', output)
        self.assertIn(f'release-23.1/logs/gff_corrections.{corrected_id()}.tsv', output)
        self.assertIn(f'release-23.1/anno/genes.{corrected_id()}.tsv', output)  # and UNITAS's feature table built from it (#440)
        self.assertNotIn('anno/genes.gff3', output)
        self.assertNotIn('anno/genes.tsv', output)

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
