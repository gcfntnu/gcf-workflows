"""ERCC92 rules (#336): maintained download URLs and SHA-256 checks; a stand-in wget serves synthetic files, no downloads."""
import hashlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[2]
CHECKSUMS = ROOT / 'gcfdb/checksums/ERCC92.sha256'
FASTA = '>ERCC-00002\nTCCAGATTACTTCCATTTCCGCCCAAGCTGCTCACAGTATACGGGCGTCGGCATCCAGACCGTCGGCTGATCGTGGTTTTACTAGGCTAGACTAGCGTACGAGCACTATGGTCA\n'
GTF = 'ERCC-00002\tERCC\texon\t1\t1061\t0.000000\t+\t.\tgene_id "ERCC-00002"; transcript_id "DQ459430";\n'
ANNO = 'Re-sort ID\tERCC ID\tsubgroup\tconcentration in Mix 1 (attomoles/ul)\n1\tERCC-00130\tA\t30000\n'

# Serves `-O <target> <url>` from $SERVE by the URL's file name and records the URL.
WGET = '''#!/bin/sh
while [ $# -gt 0 ]; do case "$1" in -O) target="$2"; shift 2;; -*) shift;; *) url="$1"; shift;; esac; done
echo "$url" >> "$SERVE/urls"
cp "$SERVE/$(basename "$url")" "$target"
'''

SNAKEFILE = '''from os.path import join
config = {{'proxy': {{}}}}
EXT_DIR = 'ext'
GCFDB_DIR = {gcfdb!r}
def src_gcf(path): return join(GCFDB_DIR, path)
include: join(GCFDB_DIR, 'ercc.db')
'''


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


class ErccChecksumFileTests(unittest.TestCase):
    def test_lists_the_published_files(self):
        rows = [line.split() for line in CHECKSUMS.read_text().splitlines()]
        self.assertEqual([name for _, name in rows], ['ERCC92.fa', 'ERCC92.gtf', 'cms_095046.txt'])
        self.assertTrue(all(len(digest) == 64 for digest, _ in rows))


class ErccRuleTests(unittest.TestCase):
    """The rules run with a copy of gcfdb whose checksum file matches synthetic files, or deliberately does not."""

    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.serve = self.tmp / 'serve'
        self.serve.mkdir()
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w') as archive:
            archive.writestr('ERCC92.fa', FASTA)
            archive.writestr('ERCC92.gtf', GTF)
        (self.serve / 'ERCC92.zip').write_bytes(buffer.getvalue())
        (self.serve / 'cms_095046.txt').write_text(ANNO)
        bin_dir = self.tmp / 'bin'
        bin_dir.mkdir()
        (bin_dir / 'wget').write_text(WGET)
        (bin_dir / 'wget').chmod(0o755)
        self.gcfdb = self.tmp / 'gcfdb'
        (self.gcfdb / 'checksums').mkdir(parents=True)
        (self.gcfdb / 'scripts').mkdir()
        (self.gcfdb / 'ercc.db').write_text((ROOT / 'gcfdb/ercc.db').read_text())
        (self.gcfdb / 'scripts/fix_ercc_gtf.py').write_text((ROOT / 'gcfdb/scripts/fix_ercc_gtf.py').read_text())
        self.project = self.tmp / 'project'
        self.project.mkdir()
        (self.project / 'Snakefile').write_text(SNAKEFILE.format(gcfdb=str(self.gcfdb)))
        self.env = dict(os.environ, PATH='{}:{}'.format(bin_dir, os.environ['PATH']), SERVE=str(self.serve))

    def publish(self, fasta=FASTA):
        (self.gcfdb / 'checksums/ERCC92.sha256').write_text(
            '{}  ERCC92.fa\n{}  ERCC92.gtf\n{}  cms_095046.txt\n'.format(sha256(fasta), sha256(GTF), sha256(ANNO)))

    def run_rules(self, *targets):
        return subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', *targets], cwd=self.project, env=self.env,
                              capture_output=True, text=True, timeout=300)

    def leftovers(self):
        return sorted(str(p.relative_to(self.project / 'ext')) for p in (self.project / 'ext').rglob('*')
                      if p.is_file() and p.parent.name != 'logs')

    def test_matching_files_are_put_in_place_from_the_maintained_urls(self):
        self.publish()
        result = self.run_rules('ext/spikein/fasta/ERCC92.fa', 'ext/spikein/anno/ERCC92.anno', 'ext/spikein/anno/ERCC92.gtf')
        self.assertEqual(result.returncode, 0, result.stderr[-3000:])
        self.assertEqual((self.project / 'ext/spikein/fasta/ERCC92.fa').read_text(), FASTA)
        self.assertEqual((self.project / 'ext/spikein/anno/ERCC92.anno').read_text(), ANNO)
        self.assertIn('gene_biotype "spikein"', (self.project / 'ext/spikein/anno/ERCC92.gtf').read_text())
        self.assertEqual(self.leftovers(), ['spikein/anno/ERCC92.anno', 'spikein/anno/ERCC92.gtf', 'spikein/fasta/ERCC92.fa'])
        self.assertEqual(sorted((self.serve / 'urls').read_text().split()),
                         ['https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip',
                          'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'])

    def test_a_checksum_mismatch_fails_and_leaves_nothing_in_place(self):
        self.publish(fasta=FASTA + 'A\n')
        result = self.run_rules('ext/spikein/fasta/ERCC92.fa')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('ERCC92.fa: FAILED', result.stdout + result.stderr)
        self.assertEqual(self.leftovers(), [])


if __name__ == '__main__':
    unittest.main()
