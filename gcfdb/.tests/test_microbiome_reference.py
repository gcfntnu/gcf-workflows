"""Dry-run the actual classifiers and QC references selected by the microbiome workflow."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'microbiome/.tests'
KIT = 'QIAseq 16S ITS Region Panels PE'
REF = Path('refs/panel v1')
EXPORT = REF / 'qiime2/classifiers/export'


class MicrobiomeReferenceTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        shutil.copytree(FIXTURE / 'data', self.work / 'data')
        self.config = yaml.safe_load((FIXTURE / 'config.yaml').read_text())
        self.config.update(ref_dir=str(REF), skip_peppy=True)
        self.primers = yaml.safe_load((ROOT / 'libprep.config').read_text())[KIT]['db']['primers']
        (self.work / 'Snakefile').write_text(
            f'configfile: "config.yaml"\ninclude: {str(ROOT / "microbiome/microbiome.smk")!r}\n'
        )

    def configure(self, provider='silva', level='99', **db):
        self.config['db'] = dict(reference_db=provider, classifier_level=level, **db)
        (self.work / 'config.yaml').write_text(yaml.safe_dump(self.config))

    def prepare_archives(self, regions, level='99'):
        archives = []
        for region in regions:
            path = EXPORT / f'{level}_{self.primers[region]}' / 'sklearn_pipeline.tar'
            (self.work / path).parent.mkdir(parents=True, exist_ok=True)
            (self.work / path).write_bytes(b'synthetic classifier')
            archives.append(path)
        return archives

    def dry_run(self, target='qiime2_run_regions', success=True):
        env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
        command = [sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', target]
        result = subprocess.run(command, cwd=self.work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60)
        self.assertEqual(result.returncode == 0, success, result.stdout)
        return result.stdout

    def test_unite_requires_selected_kit_archives_and_level(self):
        # The runtime reads primer pairs from the kit file passed to --libprep-config.
        self.configure('unite', level='dynamic', primers={'ITS1': 'IGNORED-OVERRIDE'})
        archives = self.prepare_archives(['ITS1', 'ITS2'], level='dynamic')
        (self.work / archives[1]).unlink()
        missing = self.dry_run(success=False)
        self.assertIn('Missing project classifier archive: ' + str(archives[1]), missing)

        (self.work / archives[1]).touch()
        planned = self.dry_run()
        for path in archives:
            self.assertIn(str(path), planned)
        commands = [line for line in planned.splitlines() if 'run_qiime2.py ' in line]
        self.assertEqual(len(commands), 1, planned)
        args = shlex.split(commands[0])
        self.assertEqual(args[args.index('--taxonomy-db') + 1], 'unite')
        self.assertEqual(args[args.index('--classifier-dir') + 1], str(EXPORT))
        self.assertEqual(args[args.index('--classifier-level') + 1], 'dynamic')

    def test_newer_classifier_replans_completed_taxonomy(self):
        self.configure()
        archives = self.prepare_archives([r for r in self.primers if r.startswith('V')])
        # Simulate completed inputs and a retained taxonomy from an earlier run.
        interim = Path('data/tmp/microbiome')
        quant = interim / 'quant/qiime2/silva'
        for sample in self.config['samples']:
            for read in ('R1', 'R2'):
                path = self.work / interim / 'filter/fastp' / f'{sample}_{read}.fastq'
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
        (self.work / quant).mkdir(parents=True)
        (self.work / quant / 'q2_sample_info.tsv').touch()
        target = quant / 'taxonomy.qza'
        (self.work / target).touch()
        timestamp = time.time() - 120
        for path in self.work.rglob('*'):
            os.utime(path, (timestamp, timestamp))
        os.utime(self.work / target, (timestamp + 1, timestamp + 1))
        self.assertIn('Nothing to be done', self.dry_run(str(target)))

        changed = self.work / archives[0]
        changed.write_bytes(b'updated synthetic classifier')
        os.utime(changed, (timestamp + 2, timestamp + 2))
        planned = self.dry_run(str(target))
        self.assertIn('checkpoint qiime2_run_regions:', planned)
        self.assertIn('Updated input files: ' + str(archives[0]), planned)

    def test_independent_ensembl_qc_reference_in_reporting(self):
        self.configure(ensembl={'release': None})
        self.prepare_archives([r for r in self.primers if r.startswith('V')])
        planned = self.dry_run('multiqc_report')
        release = yaml.safe_load((ROOT / 'main.config').read_text())['db']['ensembl']['release']
        self.assertIn('rule multiqc_report:', planned)
        self.assertIn(f'data/ext/ensembl/release-{release}/homo_sapiens/GRCh38/index/genome/bowtie2/genome.1.bt2', planned)
        self.assertIn(f'ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz', planned)
        self.assertNotIn('release-None', planned)
