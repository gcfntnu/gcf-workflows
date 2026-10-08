"""Exercise configmaker's generated entry point with private inputs and dry runs only."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml
from configmaker.configmaker import add_workflow

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'singlecell/.tests'
PARSE_METHODS = ('splitpipe', 'parsebio_starsolo', 'parsebio_starsolo_rt')
METHODS = ('cellranger', '10x_starsolo', *PARSE_METHODS)


class ConfigurationLoadingTests(unittest.TestCase):
    def prepare(self, method='cellranger', overrides=None, defaults_edit=None):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        work = Path(temporary.name)
        workflow = work / 'src/gcf-workflows'
        shutil.copytree(ROOT, workflow, ignore=shutil.ignore_patterns('.git', '.dev', '.tests', '.snakemake', '*.swp'))
        for name in ('data', 'pep'):
            shutil.copytree(FIXTURE / name, work / name)
        config = yaml.safe_load((FIXTURE / 'config.yaml').read_text())
        config.update(interim_dir='data/tmp', ext_dir='data/ext', fastq_dir='data/raw/fastq')
        config['quant']['method'] = method
        config['quant']['aggregate']['method'] = 'default'
        if method in PARSE_METHODS:
            config['libprepkit'] = 'Parse Biosciences Evercode WT v3'
            config['read_geometry'] = [150, 150]
            config['wells'] = {sample: {'Sample_ID': sample, 'Wells': wells}
                               for sample, wells in zip(('sample_A', 'sample_B'), ('A1-A2', 'A3-A4'))}
        if overrides:
            for key, value in overrides.items():
                config['quant'][key] = value
        original_cwd = Path.cwd()
        try:
            os.chdir(work)
            # This is the real producer used by BFQ; the workflow tree is pre-staged,
            # so add_workflow neither clones nor contacts a service.
            add_workflow(config)
        finally:
            os.chdir(original_cwd)
        (work / 'config.yaml').write_text(yaml.safe_dump(config))
        with (work / 'Snakefile').open('a') as handle:
            handle.write("\nimport json\nwith open('effective-config.json', 'w') as handle:\n    json.dump(config, handle)\n")
        defaults = workflow / 'singlecell/singlecell.config'
        if defaults_edit:
            defaults_edit(defaults)
        return work

    def dry_run(self, work, success=True, target='multiqc_report'):
        env = {key: value for key, value in os.environ.items() if not key.startswith('GCF_')}
        env['XDG_CACHE_HOME'] = str(work / 'cache')
        result = subprocess.run([sys.executable, '-m', 'snakemake', '--dry-run', '--printshellcmds',
                                 '--cores', '1', '--scheduler', 'greedy', target],
                                cwd=work, env=env, capture_output=True, text=True, timeout=60)
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode == 0, success, output)
        if success:
            self.assertIn(f'rule {target}:', output)
        return output

    def edit_defaults(self, path, edit):
        defaults = yaml.safe_load(path.read_text())
        edit(defaults)
        path.write_text(yaml.safe_dump(defaults))

    def test_generated_entry_point_loads_all_supported_methods(self):
        for method in METHODS:
            with self.subTest(method=method):
                work = self.prepare(method)
                output = self.dry_run(work)
                if 'starsolo' in method:
                    output += self.dry_run(work, target='quant_all')
                effective = json.loads((work / 'effective-config.json').read_text())
                self.assertEqual(effective['quant']['method'], method)
                starsolo = effective['quant']['starsolo']
                self.assertEqual(starsolo['10x_starsolo']['umi_dedup'], '1MM_CR')
                self.assertEqual(starsolo['parsebio_starsolo']['umi_dedup'], '1MM_All')
                self.assertEqual(starsolo['feature_count'], 'GeneFull_Ex50pAS')
                if 'starsolo' in method:
                    section = '10x_starsolo' if method == '10x_starsolo' else 'parsebio_starsolo'
                    self.assertIn('--soloUMIdedup ' + starsolo[section]['umi_dedup'], output)
                    self.assertIn('--soloUMIfiltering ' + starsolo[section]['umi_filtering'], output)
                self.assertTrue(effective['preprocessing']['enabled'])
                self.assertIn('_preprocessed.h5ad', output)

    def test_missing_defaults_reproduce_original_failure_with_actionable_error(self):
        work = self.prepare(defaults_edit=lambda path: path.unlink())
        output = self.dry_run(work, success=False)
        self.assertIn('cannot load required defaults', output)
        self.assertIn('singlecell/singlecell.config', output)
        self.assertIn("quant.method='cellranger'", output)
        self.assertNotIn('KeyError', output)

    def test_invalid_defaults_have_file_and_method_diagnostics(self):
        for content in ('[not: valid', '[]', '', '{}'):
            with self.subTest(content=content):
                work = self.prepare(defaults_edit=lambda path: path.write_text(content))
                output = self.dry_run(work, success=False)
                self.assertIn('Defaults file:', output)
                self.assertIn('singlecell/singlecell.config', output)
                self.assertIn("quant.method='cellranger'", output)
                self.assertNotIn('KeyError', output)

    def test_inactive_starsolo_settings_are_not_required(self):
        for method in ('cellranger', 'splitpipe'):
            with self.subTest(method=method):
                work = self.prepare(method, defaults_edit=lambda path: self.edit_defaults(
                    path, lambda defaults: defaults['quant'].pop('starsolo')))
                self.dry_run(work)

    def test_inactive_method_overrides_are_not_validated(self):
        for method, inactive in (('cellranger', '10x_starsolo'), ('splitpipe', 'parsebio_starsolo'),
                                 ('10x_starsolo', 'parsebio_starsolo'), ('parsebio_starsolo', '10x_starsolo')):
            with self.subTest(method=method):
                work = self.prepare(method, {'starsolo': {inactive: None}})
                self.dry_run(work)
                effective = json.loads((work / 'effective-config.json').read_text())
                self.assertIsNone(effective['quant']['starsolo'][inactive])

    def test_explicit_null_starsolo_is_preserved_for_inactive_methods(self):
        for method in ('cellranger', 'splitpipe'):
            with self.subTest(method=method):
                work = self.prepare(method, {'starsolo': None})
                self.dry_run(work)
                effective = json.loads((work / 'effective-config.json').read_text())
                self.assertIsNone(effective['quant']['starsolo'])

    def test_explicit_missing_active_settings_are_diagnosed(self):
        for starsolo, key in ((None, 'quant.starsolo'), ({'feature_count': None}, 'quant.starsolo.feature_count')):
            with self.subTest(key=key):
                work = self.prepare('10x_starsolo', {'starsolo': starsolo})
                output = self.dry_run(work, success=False)
                self.assertIn(key, output)
                self.assertIn("quant.method='10x_starsolo'", output)
                self.assertNotIn('KeyError', output)
                self.assertNotIn('AttributeError', output)

    def test_required_active_settings_name_key_and_method(self):
        cases = (('10x_starsolo', '10x_starsolo', 'umi_dedup'),
                 ('parsebio_starsolo', 'parsebio_starsolo', 'umi_filtering'))
        for method, section, key in cases:
            with self.subTest(method=method):
                work = self.prepare(method, defaults_edit=lambda path: self.edit_defaults(
                    path, lambda defaults: defaults['quant']['starsolo'][section].pop(key)))
                output = self.dry_run(work, success=False)
                self.assertIn(f'quant.starsolo.{section}.{key}', output)
                self.assertIn(f"quant.method='{method}'", output)
                self.assertIn('singlecell/singlecell.config', output)
                self.assertNotIn('KeyError', output)

    def test_starsolo_cellbender_keeps_count_filter_rule(self):
        work = self.prepare('10x_starsolo', {'cellbender': {'enabled': True}})
        output = self.dry_run(work, target='quant_all')
        self.assertIn('rule cellbender_filter_starsolo_counts:', output)

    def test_explicit_starsolo_overrides_preserve_parameters(self):
        for method, section in (('10x_starsolo', '10x_starsolo'), ('parsebio_starsolo', 'parsebio_starsolo')):
            with self.subTest(method=method):
                starsolo = {'feature_count': 'GeneFull', section: {'umi_dedup': 'Exact', 'umi_filtering': '-'}}
                work = self.prepare(method, {'starsolo': starsolo})
                output = self.dry_run(work)
                if 'starsolo' in method:
                    output += self.dry_run(work, target='quant_all')
                effective = json.loads((work / 'effective-config.json').read_text())['quant']['starsolo']
                self.assertEqual(effective['feature_count'], 'GeneFull')
                self.assertEqual(effective[section], starsolo[section])
                self.assertIn('--soloUMIdedup Exact', output)
                self.assertIn('--soloUMIfiltering -', output)

    def test_inconsistent_active_umi_options_have_diagnostic(self):
        work = self.prepare('10x_starsolo', {'starsolo': {'10x_starsolo': {'umi_dedup': 'Exact'}}})
        output = self.dry_run(work, success=False)
        self.assertIn('umi_filtering=MultiGeneUMI_CR requires umi_dedup=1MM_CR', output)
        self.assertIn("quant.method='10x_starsolo'", output)

    def test_cellranger_cellbender_does_not_require_starsolo(self):
        work = self.prepare(overrides={'cellbender': {'enabled': True}}, defaults_edit=lambda path: self.edit_defaults(
            path, lambda defaults: defaults['quant'].pop('starsolo')))
        self.dry_run(work)


if __name__ == '__main__':
    unittest.main()
