"""mirge3_library.py: download, checksum, file checks and installation of a miRge3 library; synthetic file:// store."""
import hashlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/mirge3_library.py'
spec = importlib.util.spec_from_file_location('mirge3_library', SCRIPT)
mirge3_library = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mirge3_library)

SPECIES, ORG, RELEASE = 'mus_musculus', 'mouse', '09-11-2026'
ARCHIVE = 'miRge3_lib_{}.tar.gz'.format(SPECIES)
CONTENT = {'gcf_mirna_db.tsv': 'database\trelease\nmirbase\t23.1\nmirgenedb\t3.0\n',
           'gcf_manifest.tsv': 'key\tvalue\nspecies\tmus_musculus\nensembl_release\t116\n',
           'gcf_sources.tsv': 'file\turl\tsha256\tretrieved\nmature.fa\thttps://example.org/mature.fa\t{}\tNA\n'.format('a' * 64)}
SPIKEIN = ['index.Libs/mouse_spike-in.1.ebwt', 'index.Libs/mouse_spike-in.rev.2.ebwt']


def library(drop=(), extra=()):
    files = {name: CONTENT.get(name, 'synthetic ' + name) for name in mirge3_library.library_files(ORG)}
    files.update({name: 'synthetic' for name in SPIKEIN + list(extra)})
    return {name: text for name, text in files.items() if name not in drop}


def archive_bytes(files, top=ORG):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode='w:gz') as tar:
        for name, text in sorted(files.items()):
            data = text.encode()
            info = tarfile.TarInfo('{}/{}'.format(top, name))
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


class LibraryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.store = self.tmp / 'store/mirge3' / RELEASE
        self.store.mkdir(parents=True)
        self.root = self.tmp / 'ext/mirge3' / RELEASE
        self.output = self.root / 'verify' / 'miRge3_lib_{}.tsv'.format(SPECIES)
        self.publish(library())

    def publish(self, files, top=ORG, checksum=None, name=ARCHIVE):
        data = archive_bytes(files, top)
        (self.store / ARCHIVE).write_bytes(data)
        (self.store / 'SHA256SUMS').write_text('{}  {}\n'.format(checksum or hashlib.sha256(data).hexdigest(), name))
        return data

    def run_script(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--release-url', self.store.as_uri(), '--release', RELEASE,
                                 '--species', SPECIES, '--organism-name', ORG, '--libraries-path', str(self.root),
                                 '--output', str(self.output)], capture_output=True, text=True, timeout=60)
        record = dict(line.split('\t') for line in self.output.read_text().splitlines()) if self.output.exists() else {}
        return result, record

    def assertNotInstalled(self, record):
        self.assertEqual(record, {})
        self.assertFalse((self.root / ORG).exists())
        if self.root.exists():
            self.assertEqual([p.name for p in self.root.iterdir() if p.name.startswith('.mirge3-')], [])

    def test_published_library_is_installed_without_its_spikein_index(self):
        data = self.publish(library())
        result, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record['status'], 'verified')
        self.assertEqual(record['published_sha256'], hashlib.sha256(data).hexdigest())
        self.assertEqual((record['mirbase_release'], record['mirgenedb_release'], record['ensembl_release']), ('23.1', '3.0', '116'))
        self.assertEqual(record['removed_spikein_index'], 'mouse_spike-in.1.ebwt,mouse_spike-in.rev.2.ebwt')
        lib = self.root / ORG
        self.assertTrue((lib / 'index.Libs/mouse_mirna_miRBase.1.ebwt').is_file())
        self.assertTrue((lib / 'gcf_manifest.tsv').is_file())
        self.assertEqual(list((lib / 'index.Libs').glob('mouse_spike-in.*')), [])
        self.assertEqual(sorted(p.name for p in self.root.iterdir()), [ORG, 'verify'])

    def test_human_library_needs_the_trna_fragment_annotation(self):
        self.assertIn('annotation.Libs/human_tRF_infor.csv', mirge3_library.library_files('human'))
        self.assertNotIn('annotation.Libs/mouse_tRF_infor.csv', mirge3_library.library_files(ORG))

    def test_archive_differing_from_the_published_checksum_is_not_unpacked(self):
        self.publish(library(), checksum='0' * 64)
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('publishes {}'.format('0' * 64), result.stderr)
        self.assertNotInstalled(record)

    def test_release_without_checksums_fails(self):
        (self.store / 'SHA256SUMS').unlink()
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('publishes no SHA256SUMS', result.stderr)
        self.assertNotInstalled(record)

    def test_archive_not_listed_in_the_checksums_fails(self):
        self.publish(library(), name='miRge3_lib_homo_sapiens.tar.gz')
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not list {}'.format(ARCHIVE), result.stderr)
        self.assertNotInstalled(record)

    def test_missing_library_or_provenance_files_fail_naming_them(self):
        self.publish(library(drop=('index.Libs/mouse_mirna_MirGeneDB.1.ebwt', 'gcf_manifest.tsv')))
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('index.Libs/mouse_mirna_MirGeneDB.1.ebwt', result.stderr)
        self.assertIn('gcf_manifest.tsv', result.stderr)
        self.assertNotInstalled(record)

    def test_archive_for_another_organism_name_fails(self):
        self.publish(library(), top='human')
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not only the mouse/ directory', result.stderr)
        self.assertNotInstalled(record)

    def test_existing_library_without_record_is_not_replaced(self):
        (self.root / ORG).mkdir(parents=True)
        (self.root / ORG / 'local.txt').write_text('kept')
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exists but has no verification record', result.stderr)
        self.assertEqual(record, {})
        self.assertEqual((self.root / ORG / 'local.txt').read_text(), 'kept')

    def test_unset_store_url_fails_clearly(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--release-url', 'mirge3/' + RELEASE, '--release', RELEASE,
                                 '--species', SPECIES, '--organism-name', ORG, '--libraries-path', str(self.root),
                                 '--output', str(self.output)], capture_output=True, text=True, timeout=60)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('is winecellar.url set?', result.stderr)


if __name__ == '__main__':
    unittest.main()
