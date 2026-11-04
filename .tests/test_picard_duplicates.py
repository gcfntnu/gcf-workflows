"""Instrument-aware Picard MarkDuplicates options shared by all workflows (#275)."""
import ast
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_helper():
    """Execute only the helper and its constants from utils.py, without Snakemake."""
    tree = ast.parse((ROOT / 'utils.py').read_text())
    names = {'PATTERNED_FLOWCELL_MACHINES', 'NONPATTERNED_FLOWCELL_MACHINES', 'picard_duplicate_args'}
    nodes = [node for node in tree.body
             if (isinstance(node, ast.FunctionDef) and node.name in names)
             or (isinstance(node, ast.Assign) and {t.id for t in node.targets if isinstance(t, ast.Name)} & names)]
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'utils.py', 'exec'), namespace)
    return namespace['picard_duplicate_args']


picard_duplicate_args = load_helper()


class PicardDuplicateArgsTests(unittest.TestCase):
    def test_optical_distance_follows_flowcell_type(self):
        cases = {'NovaSeq 6000': 2500, 'HiSeq 4000': 2500, 'HiSeq X': 2500,
                 'NextSeq 500': 100, 'MiSeq NTNU': 100, 'MiniSeq NTNU': 100, 'HiSeq 2500': 100, '': 100, None: 100}
        for machine, distance in cases.items():
            with self.subTest(machine=machine):
                args = picard_duplicate_args(machine).split()
                self.assertIn('TAGGING_POLICY=All', args)
                self.assertIn(f'OPTICAL_DUPLICATE_PIXEL_DISTANCE={distance}', args)
                self.assertNotIn('REMOVE_SEQUENCING_DUPLICATES=TRUE', args)

    def test_removal_applies_to_every_instrument(self):
        for machine in ('NovaSeq 6000', 'MiSeq NTNU', ''):
            with self.subTest(machine=machine):
                self.assertIn('REMOVE_SEQUENCING_DUPLICATES=TRUE', picard_duplicate_args(machine, True).split())

    def test_unknown_instrument_fails(self):
        with self.assertRaisesRegex(ValueError, "no optical-duplicate distance for 'Unknown Sequencer'"):
            picard_duplicate_args('Unknown Sequencer')


class PicardCommandTests(unittest.TestCase):
    """The delivered BAM rules receive the options; a NovaSeq run gets the patterned-flowcell distance."""

    def plan(self, workflow, machine):
        with tempfile.TemporaryDirectory(prefix='picard-') as temp:
            work = Path(temp) / workflow
            shutil.copytree(ROOT / workflow / '.tests', work, ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
            settings = yaml.safe_load((work / 'config.yaml').read_text())
            settings['machine'] = machine
            (work / 'config.yaml').write_text(yaml.safe_dump(settings))
            snakefile = (work / 'Snakefile').read_text()
            (work / 'Snakefile').write_text(snakefile.replace(f"'../{workflow}.smk'", repr(str(ROOT / workflow / f'{workflow}.smk'))))
            env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
            env.update(TMPDIR=temp, XDG_CACHE_HOME=str(Path(temp) / '.cache'))
            result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '1', 'multiqc_report'],
                                    cwd=work, env=env, capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return [line for line in result.stdout.splitlines() if 'picard MarkDuplicates' in line]

    def test_novaseq_bams_use_patterned_flowcell_distance(self):
        for workflow in ('smallrna', 'rnaseq', 'wgs'):
            with self.subTest(workflow=workflow):
                commands = self.plan(workflow, 'NovaSeq 6000')
                self.assertTrue(commands, 'no MarkDuplicates job planned')
                for command in commands:
                    self.assertIn('TAGGING_POLICY=All OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500', command)


if __name__ == '__main__':
    unittest.main()
