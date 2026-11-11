"""unitas_verify_dump.py: a cached UNITAS dump against its published archive checksum; synthetic file:// store."""
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/unitas_verify_dump.py'
DUMP = 'UNITAS_refdump_homo_sapiens'
FILES = {'db_versions.info': 'miRBase database version/date: ............ Release 23.1\n',
         'miRBase_mature': '>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p\nUGAGGUAGUAGGUUGUAUAGUU\n'}


def archive_bytes(files):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode='w:gz') as tar:
        for name, text in sorted(files.items()):
            data = text.encode()
            info = tarfile.TarInfo('{}/{}'.format(DUMP, name))
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


class VerifyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.store = self.tmp / 'store/unitas/05-11-2026'
        self.store.mkdir(parents=True)
        self.archive = archive_bytes(FILES)
        (self.store / (DUMP + '.tar.gz')).write_bytes(self.archive)
        self.publish(hashlib.sha256(self.archive).hexdigest())
        self.dump = self.tmp / 'ext/unitas/05-11-2026' / DUMP
        self.dump.mkdir(parents=True)
        for name, text in FILES.items():
            (self.dump / name).write_text(text)
        self.output = self.tmp / 'ext/unitas/05-11-2026/verify' / (DUMP + '.tsv')

    def publish(self, checksum, name=DUMP + '.tar.gz'):
        (self.store / 'SHA256SUMS').write_text('{}  {}\n'.format(checksum, name))

    def run_script(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--dump', str(self.dump), '--release-url', self.store.as_uri(),
                                 '--output', str(self.output)], capture_output=True, text=True, timeout=60)
        record = dict(line.split('\t') for line in self.output.read_text().splitlines()) if self.output.exists() else {}
        return result, record

    def test_cached_dump_matching_the_published_archive_is_verified(self):
        result, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record['status'], 'verified')
        self.assertEqual(record['published_sha256'], hashlib.sha256(self.archive).hexdigest())
        self.assertEqual(record['files_compared'], '2')
        self.assertEqual(sorted(p.name for p in self.output.parent.iterdir()), [self.output.name])

    def test_snakemake_directory_marker_is_not_a_difference(self):
        (self.dump / '.snakemake_timestamp').write_text('')
        result, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record['files_compared'], '2')

    def test_archive_differing_from_the_published_checksum_fails(self):
        self.publish('0' * 64)
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('but {} publishes {}'.format((self.store / 'SHA256SUMS').as_uri(), '0' * 64), result.stderr)
        self.assertEqual(record, {})

    def test_cache_differing_from_the_archive_fails_naming_the_files(self):
        (self.dump / 'miRBase_mature').write_text('>edited\nACGT\n')
        (self.dump / 'extra.fas').write_text('>x\nA\n')
        (self.dump / 'db_versions.info').unlink()
        result, record = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        for line in ['missing from the cache: db_versions.info', 'not in the archive: extra.fas',
                     'differs from the archive: miRBase_mature']:
            self.assertIn(line, result.stderr)
        self.assertEqual(record, {})
        self.assertEqual((self.dump / 'miRBase_mature').read_text(), '>edited\nACGT\n')

    def test_release_without_published_checksums_is_recorded_as_not_verified(self):
        (self.store / 'SHA256SUMS').unlink()
        result, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record['status'], 'not verified: the release publishes no SHA256SUMS')
        self.assertEqual(record['published_sha256'], 'NA')

    def test_archive_not_listed_is_recorded_as_not_verified(self):
        self.publish('0' * 64, name='UNITAS_refdump_mus_musculus.tar.gz')
        result, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(record['status'], 'not verified: SHA256SUMS does not list {}.tar.gz'.format(DUMP))

    def test_listed_archive_that_cannot_be_downloaded_fails(self):
        (self.store / (DUMP + '.tar.gz')).unlink()
        result, _ = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('is listed in', result.stderr)

    def test_unset_store_url_fails_clearly(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--dump', str(self.dump), '--release-url', 'unitas/05-11-2026',
                                 '--output', str(self.output)], capture_output=True, text=True, timeout=60)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('is winecellar.url set?', result.stderr)


if __name__ == '__main__':
    unittest.main()
