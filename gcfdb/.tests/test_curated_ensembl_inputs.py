"""Curated source files may occupy the same paths as the Ensembl getters."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[2]
DB = yaml.safe_load((ROOT / 'main.config').read_text())['db']
MAIN = DB['ensembl']['release']
DIVISION = DB['ensemblgenomes']['release']


class CuratedEnsemblTests(unittest.TestCase):
    def test_overlapping_source_files_have_no_getter(self):
        for organism, assembly, release, cache in [
            ('salmo_salar', 'Ssal_v3.1', MAIN, f'ensembl/release-{MAIN}'),
            ('homo_sapiens', 'GRCh38', MAIN, f'ensembl/release-{MAIN}'),
            ('arabidopsis_thaliana', 'TAIR10', DIVISION, f'ensemblgenomes/release-{DIVISION}/plants'),
            ('arabidopsis_thaliana', 'TAIR10', MAIN, f'ensembl/release-{MAIN}'),
        ]:
            with self.subTest(organism=organism, cache=cache), tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                ref = Path('data/ext') / cache / organism / assembly
                fasta, gtf = ref / 'fasta/genome.fa', ref / 'anno/genes.gtf'
                (work / fasta.parent).mkdir(parents=True)
                (work / gtf.parent).mkdir(parents=True)
                (work / fasta).write_text('>synthetic\nACGT\n')
                (work / gtf).write_text('# curated panel\n')
                os.utime(work / gtf, (100, 100))
                (work / gtf.parent / f'{organism}.{assembly}.{release}.gtf').write_text('# newer stock annotation\n')
                (work / 'panel').symlink_to(ref, target_is_directory=True)

                settings = {'organism': organism, 'read_geometry': [75], 'ref_dir': 'panel'}
                (work / 'Snakefile').write_text(
                    f'config.update({settings!r})\n'
                    f'include: {str(ROOT / "utils.py")!r}\n'
                    f'include: {str(ROOT / "gcfdb/ensembl.db")!r}\n'
                    'rule consume:\n'
                    f'    input: {str(fasta)!r}, {str(gtf)!r}\n'
                    '    output: "result.txt"\n'
                    '    shell: "cat {input:q} > {output:q}"\n'
                )
                env = {key: value for key, value in os.environ.items() if not key.startswith('GCF_')}
                # A forced downstream rebuild must still treat both files as source inputs.
                result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--forceall', '--cores', '1', 'consume'],
                                        cwd=work, env=env, capture_output=True, text=True, timeout=60)
                output = result.stdout + result.stderr
                self.assertEqual(result.returncode, 0, output)
                self.assertIn(f'cat {fasta} {gtf} > result.txt', output)
                for getter in ('ensembl_symlink_genes', 'ensembl_genome', 'human_ensembl_genome', 'ensembl_gtf'):
                    self.assertNotIn(f'rule {getter}:', output)
