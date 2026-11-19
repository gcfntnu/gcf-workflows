"""Dry-run each custom-reference reader; no scientific tools or downloads."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class CustomReferenceTests(unittest.TestCase):
    def test_project_directory_or_configured_provider(self):
        readers = (
            ('gcfdb/reference_db.db', 'ensembl', 'data/ext/ensembl/release-116/homo_sapiens/GRCh38'),
            ('metagenome/rules/reference_db.smk', 'langmead', 'data/ext/langmead/release-20260626/metagenome/k2_pluspf'),
            ('microbiome/rules/reference_db.smk', 'silva', 'data/ext/silva/132'),
        )
        for reader, database, default in readers:
            for project in ({}, {'ref_dir': None}, {'ref_dir': ''}, {'ref_dir': 'refs/panel-v1'}):
                with self.subTest(reader=reader, project=project), tempfile.TemporaryDirectory() as directory:
                    work = Path(directory)
                    settings = {'read_geometry': [75], 'organism': 'homo_sapiens', 'db': {'reference_db': database}, **project}
                    expected = project.get('ref_dir') or default
                    if project.get('ref_dir'):
                        (work / expected).mkdir(parents=True)
                        if database == 'langmead':
                            expected = str((work / expected).resolve())
                    # The shared reader follows the Ensembl lookup in the genomic workflows.
                    ensembl_include = f'include: {str(ROOT / "gcfdb/ensembl.db")!r}\n' if database == 'ensembl' else ''
                    (work / 'Snakefile').write_text(
                        f'config.update({settings!r})\n'
                        f'include: {str(ROOT / "utils.py")!r}\n'
                        f'{ensembl_include}'
                        f'include: {str(ROOT / reader)!r}\n'
                        f'assert REF_DIR == {expected!r}, REF_DIR\n'
                        f'assert config["db"]["reference_db"] == {database!r}\n'
                        'assert DB_CONF == config["db"][config["db"]["reference_db"]]\n'
                        'rule reference_probe:\n    output: "probe.txt"\n    params: reference=REF_DIR\n'
                        '    shell: "printf %s {params.reference:q} > {output:q}"\n'
                    )
                    env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
                    result = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', 'reference_probe'],
                                            cwd=work, env=env, capture_output=True, text=True, timeout=60)
                    output = result.stdout + result.stderr
                    self.assertEqual(result.returncode, 0, output)
                    self.assertIn('rule reference_probe:', output)
                    self.assertIn(expected, output)
                    self.assertFalse((work / 'probe.txt').exists())


if __name__ == '__main__':
    unittest.main()
