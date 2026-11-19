"""Check the project reference passed to the microbiome QIIME2 checkpoint."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'microbiome/.tests'


class MicrobiomeReferenceTests(unittest.TestCase):
    def test_explicit_classifier_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            shutil.copytree(FIXTURE / 'data', work / 'data')
            config = yaml.safe_load((FIXTURE / 'config.yaml').read_text())
            config.update(ref_dir='refs/panel v1', skip_peppy=True)
            config.setdefault('db', {})['reference_db'] = 'silva'
            (work / 'config.yaml').write_text(yaml.safe_dump(config))
            (work / 'Snakefile').write_text(
                f'configfile: "config.yaml"\ninclude: {str(ROOT / "microbiome/microbiome.smk")!r}\n'
            )
            classifiers = Path(config['ref_dir']) / 'qiime2/classifiers/export'
            command = [sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', 'qiime2_run_regions']
            env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
            run_options = dict(cwd=work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60)

            missing = subprocess.run(command, **run_options)
            self.assertNotEqual(missing.returncode, 0, missing.stdout)
            self.assertIn('MissingInputException', missing.stdout)
            self.assertIn(str(classifiers), missing.stdout)

            (work / classifiers).mkdir(parents=True)
            planned = subprocess.run(command, **run_options)
            self.assertEqual(planned.returncode, 0, planned.stdout)
            commands = [line for line in planned.stdout.splitlines() if 'run_qiime2.py ' in line]
            self.assertEqual(len(commands), 1, planned.stdout)
            args = shlex.split(commands[0])
            self.assertEqual(args[args.index('--taxonomy-db') + 1], 'silva')
            self.assertEqual(args[args.index('--classifier-dir') + 1], str(classifiers))
