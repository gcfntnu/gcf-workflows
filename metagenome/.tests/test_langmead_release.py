"""Real DAG checks for reference identity; no downloads or classification."""
import hashlib
import os
import shlex
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
MERS = (50, 75, 100, 150, 200, 250, 300)


class LangmeadReleaseTests(unittest.TestCase):
    def test_acquisition_checksums_reject_corrupt_or_missing_files(self):
        for defect in ('none', 'corrupt', 'missing'):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                source = work / 'archive-files'
                source.mkdir()
                names = [line.split()[1] for line in (ROOT / 'gcfdb/checksums/k2_pluspf_20260626.md5').read_text().splitlines()]
                manifest = work / 'gcfdb/checksums/k2_pluspf_20260626.md5'
                manifest.parent.mkdir(parents=True)
                for name in names:
                    (source / name).write_text('synthetic ' + name)
                manifest.write_text(''.join(f'{hashlib.md5((source / name).read_bytes()).hexdigest()}  {name}\n' for name in names))
                if defect == 'corrupt':
                    (source / 'hash.k2d').write_text('damaged')
                elif defect == 'missing':
                    (source / 'database150mers.kmer_distrib').unlink()
                archive = work / 'synthetic.tar.gz'
                with tarfile.open(archive, 'w:gz') as tar:
                    for path in source.iterdir():
                        tar.add(path, arcname=path.name)
                binary = work / 'bin/wget'
                binary.parent.mkdir()
                binary.write_text('#!/bin/sh\ncat ' + shlex.quote(str(archive)) + '\n')
                binary.chmod(0o755)
                old = work / 'ext/langmead/release-20231009/metagenome/k2_pluspf/hash.k2d'
                old.parent.mkdir(parents=True)
                old.write_text('retained old database')
                (work / 'Snakefile').write_text(
                    'from os.path import join\n'
                    "config.update({'db': {'langmead': {'release': '20260626', 'assembly': 'k2_pluspf'}}})\n"
                    f'EXT_DIR = {str(work / "ext")!r}\nGCFDB_DIR = {str(work / "gcfdb")!r}\n'
                    f'include: {str(ROOT / "gcfdb/langmead.db")!r}\n')
                result = subprocess.run([sys.executable, '-m', 'snakemake', '--cores', '1', 'langmead_kraken_prebuild'],
                                        cwd=work, capture_output=True, text=True, timeout=60,
                                        env=dict(os.environ, PATH=str(binary.parent) + os.pathsep + os.environ['PATH']))
                output = result.stdout + result.stderr
                self.assertEqual(result.returncode == 0, defect == 'none', output)
                self.assertIn('hash.k2d: ' + ('FAILED' if defect == 'corrupt' else 'OK'), output)
                self.assertEqual(old.read_text(), 'retained old database')

    def test_full_default_and_published_manifest_contract(self):
        settings = yaml.safe_load((ROOT / 'main.config').read_text())
        self.assertEqual(settings['db']['langmead'], {'release': '20260626', 'assembly': 'k2_pluspf'})
        required = {'hash.k2d', 'opts.k2d', 'taxo.k2d', 'seqid2taxid.map'}
        required.update(f'database{n}mers.kmer_distrib' for n in MERS)
        for assembly in ('k2_pluspf', 'k2_pluspf_16_GB'):
            manifest = ROOT / 'gcfdb/checksums' / f'{assembly}_20260626.md5'
            entries = dict(line.split()[::-1] for line in manifest.read_text().splitlines())
            self.assertTrue(required <= entries.keys())
            self.assertEqual(len(entries), 17)
            self.assertTrue(all(len(v) == 32 and int(v, 16) >= 0 for v in entries.values()))

    def test_screening_and_metagenome_use_release_scoped_ram(self):
        for consumer in ('screening', 'metagenome'):
            for release, assembly in (('20231009', 'k2_pluspf'), ('20260626', 'k2_pluspf'), ('20260626', 'k2_pluspf_16_GB')):
                with self.subTest(consumer=consumer, release=release, assembly=assembly), tempfile.TemporaryDirectory() as directory:
                    work = Path(directory)
                    db = work / 'ext/langmead' / f'release-{release}' / 'metagenome' / assembly
                    db.mkdir(parents=True)
                    for name in ['hash.k2d', 'opts.k2d', 'taxo.k2d', 'seqid2taxid.map', *[f'database{n}mers.kmer_distrib' for n in MERS]]:
                        (db / name).touch()
                    for name in ('R1.fastq.gz', 'R2.fastq.gz'):
                        (work / name).touch()
                    settings = {'db': {'reference_db': 'langmead', 'langmead': {'release': release, 'assembly': assembly}},
                                'libprepkit': 'Synthetic paired library', 'winecellar': {'url': 'https://example.invalid'},
                                'docker': {name: 'example/' + name + ':test' for name in ('kraken2', 'bracken', 'krona', 'kraken-biom', 'phyloseq', 'seqkit')}}
                    text = ('import os\nfrom os.path import join, dirname\n'
                            f'config.update({settings!r})\nEXT_DIR = {str(work / "ext")!r}\n'
                            f'GCFDB_DIR = {str(ROOT / "gcfdb")!r}\n'
                            "FILTER_INTERIM = 'filter'\nQC_INTERIM = 'qc'\nQUANT_INTERIM = 'quant'\nINTERIM_DIR = 'interim'\n"
                            "PE = True\nSAMPLES = ['synthetic']\nread_geometry = [151, 151]\n"
                            "def get_filtered_fastq(wildcards):\n    return {'R1': 'R1.fastq.gz', 'R2': 'R2.fastq.gz'}\n"
                            f'def src_gcf(path):\n    return join({str(ROOT)!r}, path)\n')
                    if consumer == 'screening':
                        text += f'include: {str(ROOT / "common/kraken2.smk")!r}\n'
                        target = 'qc/kraken2/synthetic.kraken.kreport'
                    else:
                        text += f'include: {str(ROOT / "gcfdb/langmead.db")!r}\n'
                        text += f'include: {str(ROOT / "metagenome/rules/quant/kraken2.smk")!r}\n'
                        target = f'quant/kraken2/{assembly}/synthetic/synthetic.kraken.kreport'
                    (work / 'Snakefile').write_text(text)
                    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', target],
                                            cwd=work, capture_output=True, text=True, timeout=60)
                    output = result.stdout + result.stderr
                    self.assertEqual(result.returncode, 0, output)
                    self.assertIn(f'--db /dev/shm/langmead/release-{release}/metagenome/{assembly}', output)
                    self.assertNotIn(f'--db /dev/shm/{assembly} ', output)
                    self.assertIn('rule langmead_shmem:', output)
                    self.assertNotIn('rule langmead_kraken_prebuild:', output)


if __name__ == '__main__':
    unittest.main()
