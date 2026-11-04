"""Dry-run each custom-reference reader; no scientific tools or downloads."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class CustomReferenceTests(unittest.TestCase):
    def test_project_then_environment_then_generated_reference(self):
        readers = (
            ('gcfdb/reference_db.db', 'ensembl', 'data/ext/ensembl/release-116/homo_sapiens/'),
            ('metagenome/rules/reference_db.smk', 'langmead', 'data/ext/langmead/20260626'),
            ('microbiome/rules/reference_db.smk', 'silva', 'data/ext/silva/132'),
        )
        cases = ((None, None), ('project-ref', None), (None, 'env-ref'), ('project-ref', 'env-ref'),
                 ('', 'env-ref'), ('', ''), (None, ''))
        for reader, database, default in readers:
            for project, environment in cases:
                with self.subTest(reader=reader, project=project, environment=environment), tempfile.TemporaryDirectory() as directory:
                    work = Path(directory)
                    settings = {'read_geometry': [75], 'organism': 'homo_sapiens', 'db': {'reference_db': database},
                                'ref_dir': project}
                    expected = project or environment or default
                    (work / 'Snakefile').write_text(
                        f'config.update({settings!r})\n'
                        f'include: {str(ROOT / "utils.py")!r}\n'
                        f'include: {str(ROOT / reader)!r}\n'
                        f'assert REF_DIR == {expected!r}, REF_DIR\n'
                        'rule reference_probe:\n    output: "probe.txt"\n    params: reference=REF_DIR\n'
                        '    shell: "printf %s {params.reference:q} > {output:q}"\n'
                    )
                    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
                    if environment is not None:
                        env['GCF_REFDIR'] = environment
                    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', 'reference_probe'],
                                            cwd=work, env=env, capture_output=True, text=True, timeout=60)
                    output = result.stdout + result.stderr
                    self.assertEqual(result.returncode, 0, output)
                    self.assertIn('rule reference_probe:', output)
                    self.assertIn(expected, output)
                    self.assertFalse((work / 'probe.txt').exists())


if __name__ == '__main__':
    unittest.main()
