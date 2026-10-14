"""Configuration contract tests; no sequencing data or server access required."""
import copy
import importlib.util
from pathlib import Path
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('immune_config', ROOT / 'singlecell/rules/immune_config.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
validate = MODULE.validate_singlecell_assay


def tcr_config():
    preset = yaml.safe_load((ROOT / 'libprep.config').read_text())['Parse Biosciences Evercode TCR v4 PE']
    preset.update(organism='homo_sapiens', read_geometry=[242, 58], samples={'sublib1': {}, 'sublib2': {}})
    preset['immune']['wt_parent'] = {
        'root': '/external/wt', 'sublibraries': {'sublib1': 'splitpipe/sublib1', 'sublib2': 'splitpipe/sublib2'},
        'combined': 'aggregate/splitpipe',
    }
    preset['immune']['sublibrary_order'] = ['sublib1', 'sublib2']
    return preset


class ImmuneConfigTests(unittest.TestCase):
    def test_existing_wt_default_is_unchanged(self):
        for config in ({}, {'quant': {'method': 'cellranger'}}, {'quant': {'method': 'splitpipe', 'assay': 'wt'}}):
            original = copy.deepcopy(config)
            self.assertEqual(validate(config), 'wt')
            self.assertEqual(config, original)

    def test_valid_parent_linked_tcr_and_imgt(self):
        config = tcr_config()
        original = copy.deepcopy(config)
        self.assertEqual(validate(config), 'tcr')
        self.assertEqual(config, original)
        config['immune']['reference'] = 'imgt'
        self.assertEqual(validate(config), 'tcr')

    def test_preset_requires_project_parent(self):
        config = tcr_config()
        config['immune']['wt_parent']['root'] = None
        with self.assertRaisesRegex(ValueError, 'absolute path'):
            validate(config)

    def test_no_silent_wt_routing(self):
        for assay in ('TCR', 'bcr', None):
            with self.subTest(assay=assay), self.assertRaisesRegex(ValueError, 'quant.assay'):
                validate({'quant': {'assay': assay}})
        with self.assertRaisesRegex(ValueError, 'immune configuration requires'):
            validate({'immune': {'species': 'human'}})

    def test_mapping_must_cover_samples_and_be_unique(self):
        for mapping in ({'sublib1': 'one'}, {'sublib1': 'one', 'sublib2': 'one'}, {'sublib1': 'one', 'wrong': 'two'}):
            config = tcr_config()
            config['immune']['wt_parent']['sublibraries'] = mapping
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                validate(config)

    def test_paths_stay_under_parent_root(self):
        for path in ('/absolute', '../outside', '.', '', None):
            config = tcr_config()
            config['immune']['wt_parent']['sublibraries']['sublib1'] = path
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, 'WT parent path'):
                validate(config)

    def test_combine_order_is_complete_and_unique(self):
        for order in ([], ['sublib1', 'sublib1'], ['sublib2'], 'sublib1,sublib2'):
            config = tcr_config()
            config['immune']['sublibrary_order'] = order
            with self.subTest(order=order), self.assertRaisesRegex(ValueError, 'sublibrary_order'):
                validate(config)

    def test_wt_processing_options_rejected(self):
        for key, value in [('use_velo', True), ('cellbender', {'enabled': True}), ('aggregate', {'skip': False})]:
            config = tcr_config()
            config['quant'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate(config)

    def test_unsupported_assay_inputs_rejected(self):
        for section, key, value in [('quant', 'method', 'cellranger'), ('quant', 'kit', 'TCR'),
                                    ('quant', 'chemistry', 'v3'), ('immune', 'species', 'mouse'),
                                    ('immune', 'reference', 'unknown')]:
            config = tcr_config()
            config[section][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate(config)
        for key, value in [('read_geometry', [242]), ('read_geometry', [242, 0]), ('organism', 'mus_musculus')]:
            config = tcr_config()
            config[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate(config)


if __name__ == '__main__':
    unittest.main()
