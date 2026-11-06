"""Real configuration-loader and reporting DAG checks, using synthetic fixtures only."""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml
from configmaker.configmaker import add_workflow
from configmaker.libprep import LibprepConfig

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ('default', 'rnaseq', 'smallrna', 'singlecell', 'metagenome', 'microbiome', 'wgs')
MRNA_KIT = 'Illumina Stranded mRNA Prep'
QIASEQ_KIT = 'QIAseq 16S ITS Region Panels'


def merge(target, source):
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            merge(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def leaves(settings, prefix=()):
    for key, value in settings.items():
        path = (*prefix, key)
        if isinstance(value, dict):
            yield from leaves(value, path)
        else:
            yield path


class ConfigurationTests(unittest.TestCase):
    def load_fixture(self, workflow, overrides=None, initialize=False, secrets=None, success=True):
        temporary = tempfile.TemporaryDirectory(prefix='workflow-config-')
        self.addCleanup(temporary.cleanup)
        work = Path(temporary.name)
        shutil.copytree(ROOT / workflow / '.tests', work, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('.snakemake', '__pycache__'))
        settings = yaml.safe_load((work / 'config.yaml').read_text())
        merge(settings, overrides or {})
        if initialize:
            # Exercise configmaker's actual producer API before the workflow consumes its output.
            source = work / 'src/gcf-workflows'
            source.mkdir(parents=True)
            previous = Path.cwd()
            try:
                os.chdir(work)
                add_workflow(settings, src_dir=work / 'src', libprep_config=ROOT / 'libprep.config')
            finally:
                os.chdir(previous)
        if len(settings['read_geometry']) == 2:
            for sample in settings['samples'].values():
                if not sample.get('R2'):
                    sample['R2'] = sample['R1'].replace('_R1.', '_R2.')
                for read in ('R1', 'R2'):
                    for filename in sample[read].split(','):
                        path = work / 'data/raw/fastq' / filename
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.touch()
        (work / 'config.yaml').write_text(yaml.safe_dump(settings))
        (work / 'Snakefile').write_text(
            'pepfile: "pep/pep_config.yaml"\nconfigfile: "config.yaml"\n'
            f'include: {str(ROOT / workflow / (workflow + ".smk"))!r}\n'
            'import json\nfrom pathlib import Path\n'
            'Path("effective.json").write_text(json.dumps(config, default=str))\n'
        )
        env = {key: value for key, value in os.environ.items() if not key.startswith('GCF_')}
        env.update(TMPDIR=str(work), XDG_CACHE_HOME=str(work / '.cache'))
        if secrets:
            (work / 'secrets.yaml').write_text(yaml.safe_dump(secrets))
            env['GCF_SECRET'] = str(work / 'secrets.yaml')
        # Use the installed command: python -m snakemake does not expose sibling-script import shadowing.
        result = subprocess.run([str(Path(sys.executable).with_name('snakemake')), '-n', '-p', '--cores', '1', 'multiqc_report'],
                                cwd=work, env=env, capture_output=True, text=True, timeout=120)
        output = result.stdout + result.stderr
        self.assertEqual(result.returncode == 0, success, output)
        if not success:
            self.assertNotIn('rule multiqc_report:', output)
            return output
        self.assertIn('rule multiqc_report:', output)
        return json.loads((work / 'effective.json').read_text()), output

    def test_configuration_files_have_no_duplicate_keys(self):
        # PyYAML silently keeps the last duplicate, so a repeated key hides which value is in effect.
        duplicates = []

        class Loader(yaml.SafeLoader):
            pass

        def construct_mapping(loader, node, deep=False):
            seen = set()
            for key_node, _ in node.value:
                key = loader.construct_object(key_node, deep=deep)
                if key in seen:
                    duplicates.append(f'{Path(loader.name).name}:{key_node.start_mark.line + 1}: {key}')
                seen.add(key)
            return yaml.SafeLoader.construct_mapping(loader, node, deep)

        Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)
        files = [ROOT / name for name in ('main.config', 'libprep.config', 'docker.config', 'demultiplexing.config')]
        files += sorted(ROOT.glob('*/*.config'))
        for path in files:
            with path.open() as stream:
                yaml.load(stream, Loader=Loader)
        self.assertEqual(duplicates, [])

    def test_shared_and_workflow_defaults_have_no_overlapping_leaves(self):
        shared = set(leaves(yaml.safe_load((ROOT / 'main.config').read_text())))
        for workflow in WORKFLOWS:
            path = ROOT / workflow / (workflow + '.config')
            specific = set(leaves(yaml.safe_load(path.read_text()) or {})) if path.exists() else set()
            with self.subTest(workflow=workflow):
                conflicts = {(a, b) for a in shared for b in specific if a[:len(b)] == b or b[:len(a)] == a}
                self.assertEqual(conflicts, set())

    def test_rnaseq_defaults_and_explicit_precedence(self):
        for project, secret, expected in (({}, {}, 'fastq'), ({}, {'assembly': {'transcriptome': {'merge': {'type': 'bam'}}}}, 'bam'),
                                           ({'assembly': {'transcriptome': {'merge': {'type': 'fastq'}}}},
                                            {'assembly': {'transcriptome': {'merge': {'type': 'bam'}}}}, 'fastq')):
            with self.subTest(expected=expected, project=project, secret=secret):
                config, _ = self.load_fixture('rnaseq', project, secrets=secret)
                self.assertEqual(config['assembly']['transcriptome']['merge']['type'], expected)
                self.assertEqual(config['assembly']['transcriptome']['assembler'], 'skip')
                self.assertEqual(config['assembly']['genome']['assembler'], 'skip')
                self.assertEqual(config['align']['genome']['aligner'], 'star')
                self.assertEqual(config['qc']['bam_subset'], 'housekeeping_genes')
                self.assertNotIn('base_reference_db', config['assembly'])
                self.assertNotIn('base_ref_dir', config['assembly'])
                self.assertNotIn('method', config['align'])
                self.assertNotIn('transcriptome', config['align'])

    def test_smallrna_defaults_keep_existing_kit_trimming(self):
        config, _ = self.load_fixture('smallrna')
        self.assertEqual(config['align']['genome']['aligner'], 'star')
        self.assertIsNone(config['filter']['spikein']['ref'])
        self.assertEqual(config['filter']['ribosomal']['ref'], 'ncbi')
        self.assertEqual(config['filter']['contaminants']['ref'], 'univec_core')
        before = LibprepConfig.load(ROOT / 'libprep.config').select(config['libprepkit'], [75]).parameters
        self.assertEqual(config['filter']['trim']['fastp']['args'], before['filter']['trim']['fastp']['args'])
        self.assertNotIn('trim_args', config['filter']['trim']['fastp'])

    def test_singlecell_sample_id_is_a_live_qc_default(self):
        config, output = self.load_fixture('singlecell')
        self.assertEqual(config['sample_id'], 'Sample_ID')
        self.assertIn('--qc-sample Sample_ID', output)
        self.assertNotIn('bam_subset', config['qc'])

    def test_configmaker_accepted_kit_spellings_plan_rnaseq_jobs(self):
        for name in (MRNA_KIT, MRNA_KIT + ' PE', MRNA_KIT.lower()):
            with self.subTest(name=name):
                config, output = self.load_fixture('rnaseq', {'libprepkit': name, 'read_geometry': [75, 150]}, initialize=True)
                self.assertEqual(config['libprepkit'], MRNA_KIT)
                self.assertEqual(config['libprep_selection']['entry'], MRNA_KIT + ' PE')
                self.assertEqual(config['libprep_selection']['kit'], name)
                # postprocess.smk restores reporting geometry after utils applies the kit delta.
                self.assertEqual(config['read_geometry'], [75, 150])
                self.assertIn('rule fastp:', output)

    def test_unknown_or_geometry_mismatched_kit_fails(self):
        for name, message in (('Synthetic unknown kit', 'Unsupported library kit'), (MRNA_KIT + ' SE', 'conflicts with PE')):
            with self.subTest(name=name):
                output = self.load_fixture('rnaseq', {'libprepkit': name, 'read_geometry': [75, 150]}, success=False)
                self.assertIn(message, output)

    def test_recorded_selection_mismatch_fails(self):
        selection = LibprepConfig.load(ROOT / 'libprep.config').select(MRNA_KIT, [75, 150]).diagnostics()
        for field, value in (('entry', 'Custom PE'), ('read_geometry', [50, 50]), ('workflow', 'default')):
            with self.subTest(field=field):
                recorded = dict(selection, **{field: value})
                output = self.load_fixture('rnaseq', {'libprepkit': MRNA_KIT, 'read_geometry': [75, 150],
                                                     'libprep_selection': recorded}, success=False)
                self.assertIn('libprep_selection.' + field, output)

    def test_qiaseq_retains_dormant_settings_and_pluspf_demultiplexing(self):
        config, output = self.load_fixture('metagenome', {'libprepkit': QIASEQ_KIT.lower() + ' PE'})
        self.assertEqual(config['db']['reference_db'], 'langmead')
        self.assertEqual(config['db']['langmead']['assembly'], 'k2_pluspf')
        self.assertEqual(config['quant']['method'], 'qiime2')
        self.assertIn('dada2', config['quant']['qiime2'])
        self.assertEqual(config['db']['classifier_level'], 99)
        self.assertIn('rule cutadapt_demultiplex:', output)
        self.assertIn("--libkit 'QIAseq 16S ITS Region Panels PE'", output)
        self.assertNotIn('rule qiime2_', output)


if __name__ == '__main__':
    unittest.main()
