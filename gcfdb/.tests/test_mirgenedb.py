"""UNITAS dumps with MirGeneDB in place of miRBase: dump selection, feature info and db log; no downloads."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
GENES = ROOT / 'gcfdb/scripts/mirgenedb_gff2genes.py'
DBLOG = ROOT / 'gcfdb/scripts/format_unitas_dblog.py'

GFF = ('##gff-version 3\n# microRNAs:    MirGeneDB v3.0\n# genome-build-id:  hg38\n'
       'chr1\t.\tpre_miRNA\t1167124\t1167182\t.\t+\t.\tID=Hsa-Mir-8-P1a_pre;Alias=MI0000342\n'
       'chr1\t.\tmiRNA\t1167124\t1167145\t.\t+\t.\tID=Hsa-Mir-8-P1a_5p*;Alias=MIMAT0004571\n'
       'chr1\t.\tmiRNA\t1167160\t1167182\t.\t+\t.\tID=Hsa-Mir-8-P1a_3p;Alias=MIMAT0000318\n'
       'chr1\t.\tmiRNA\t1169020\t1169040\t.\t-\t.\tID=Hsa-Mir-8-P3a_5p*\n')


def genes(text):
    with tempfile.TemporaryDirectory() as directory:
        gff = Path(directory) / 'genes.gff3'
        gff.write_text(text)
        return subprocess.run([sys.executable, str(GENES), str(gff)], capture_output=True, text=True, timeout=60)


class MirgenedbGenesTests(unittest.TestCase):
    def test_mature_and_star_mirnas_named_as_unitas_reports_them(self):
        result = genes(GFF)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = [line.split('\t') for line in result.stdout.splitlines()]
        self.assertEqual(lines[0], ['mirna_id', 'seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase',
                                    'id', 'alias', 'precursor', 'arm', 'star'])
        self.assertEqual([(row[0], row[9], row[10], row[11], row[12], row[13]) for row in lines[1:]],
                         [('Mir-8-P1a-5p', 'Hsa-Mir-8-P1a_5p*', 'MIMAT0004571', 'Hsa-Mir-8-P1a_pre', '5p', 'True'),
                          ('Mir-8-P1a-3p', 'Hsa-Mir-8-P1a_3p', 'MIMAT0000318', 'Hsa-Mir-8-P1a_pre', '3p', 'False'),
                          ('Mir-8-P3a-5p', 'Hsa-Mir-8-P3a_5p*', '', 'Hsa-Mir-8-P3a_pre', '5p', 'True')])

    def test_rejects_other_gff(self):
        for text, message in [(GFF.replace('MirGeneDB v3.0', 'miRBase v23.1'), 'not a MirGeneDB GFF'),
                              (GFF.replace('##gff-version 3\n', ''), 'not GFF3'),
                              (GFF + GFF.splitlines(True)[4], 'duplicate miRNA Mir-8-P1a-5p'),
                              (GFF.replace('ID=Hsa-Mir-8-P1a_3p', 'ID=Hsa-Mir-8-P1a'), 'not a MirGeneDB miRNA ID'),
                              ('\n'.join(GFF.splitlines()[:4]) + '\n', 'no miRNA features')]:
            with self.subTest(message=message):
                result = genes(text)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)


class DbLogTests(unittest.TestCase):
    def test_mirgenedb_dump_logs_mirgenedb_release(self):
        with tempfile.TemporaryDirectory() as directory:
            info = Path(directory) / 'db_versions.info'
            info.write_text('Ensembl version/date: ..................... Release 116\n'
                            'miRBase database version/date: ............ MirGeneDB 3.0\n'
                            'tRF-1 sequence data version/date: ......... 09.04.2019 (dd.mm.yyyy)\n')
            result = subprocess.run([sys.executable, str(DBLOG), str(info)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ['UNITAS Ensembl,116,NA,NA', 'UNITAS MirGeneDB,3.0,NA,NA',
                                                      'UNITAS tRF-1 sequence data,NA,NA,09-04-2019'])


def dry_run(unitas, target):
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        settings = {'read_geometry': [50], 'organism': 'homo_sapiens', 'db': {'unitas': unitas}, 'winecellar': {'url': 'https://wine.example/gcf'},
                    'docker': {'default': 'example/default:1'}}
        (work / 'Snakefile').write_text(f'config.update({settings!r})\n'
                                        f'include: {str(ROOT / "utils.py")!r}\n'
                                        f'include: {str(ROOT / "gcfdb/unitas.db")!r}\n')
        env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
        result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', target],
                                cwd=work, env=env, capture_output=True, text=True, timeout=60)
        return result.returncode, result.stdout + result.stderr


class DumpSelectionTests(unittest.TestCase):
    RELEASES = {'release': '05-11-2026', 'mirgenedb_release': '06-11-2026-mirgenedb'}

    def test_mirbase_is_the_default(self):
        returncode, output = dry_run(self.RELEASES, 'unitas_db_all')
        self.assertEqual(returncode, 0, output)
        self.assertIn('https://wine.example/gcf/unitas/05-11-2026/UNITAS_refdump_homo_sapiens.tar.gz', output)
        self.assertNotIn('gcf_mirna_db.tsv', output)

    def test_mirgenedb_uses_its_release_checks_the_dump_and_annotates_from_its_gff(self):
        target = 'data/ext/unitas/06-11-2026-mirgenedb/anno/mirgenedb_genes_homo_sapiens.tsv'
        returncode, output = dry_run(dict(self.RELEASES, mirna_db='mirgenedb'), target)
        self.assertEqual(returncode, 0, output)
        self.assertIn('https://wine.example/gcf/unitas/06-11-2026-mirgenedb/UNITAS_refdump_homo_sapiens.tar.gz', output)
        self.assertIn('data/ext/unitas/06-11-2026-mirgenedb/UNITAS_refdump_homo_sapiens/gcf_mirna_db.tsv', output)
        self.assertIn('rule mirgenedb_genes:', output)
        self.assertIn('data/ext/unitas/06-11-2026-mirgenedb/UNITAS_refdump_homo_sapiens/gcf_mirgenedb.gff3', output)

    def test_bad_settings_stop_at_load(self):
        for unitas, message in [(dict(self.RELEASES, mirna_db='mirbas'), "db.unitas.mirna_db: mirbase or mirgenedb, got 'mirbas'"),
                                ({'release': '05-11-2026', 'mirna_db': 'mirgenedb', 'mirgenedb_release': None},
                                 'db.unitas.mirna_db is mirgenedb, but db.unitas.mirgenedb_release is not set')]:
            with self.subTest(message=message):
                returncode, output = dry_run(unitas, 'unitas_db_all')
                self.assertNotEqual(returncode, 0)
                self.assertIn(message, output)


if __name__ == '__main__':
    unittest.main()
