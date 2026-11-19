"""Project FastQ Screen selects and tracks a complete Bowtie2 index of the required size."""
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
FIXTURE = ROOT / 'rnaseq/.tests'
PARTS = ('1', '2', '3', '4', 'rev.1', 'rev.2')
SMALL_MAX_BYTES = 4 * 1024**3 - 200


class ProjectBowtie2Tests(unittest.TestCase):
    def check_plan(self, size, cached, expected_extension, needs_build):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            ref = Path('refs/panel one')
            prefix = ref / 'index/genome/bowtie2/genome'
            (work / prefix.parent).mkdir(parents=True)
            fasta = ref / 'fasta/genome.fa'
            (work / fasta.parent).mkdir()
            (work / fasta).write_text('>synthetic\nACGT\n')
            for extension, parts in cached.items():
                for part in parts:
                    (work / f'{prefix}.{part}.{extension}').touch()
            config = yaml.safe_load((FIXTURE / 'config.yaml').read_text())
            config.update(ref_dir=str(ref), organism='synthetic_organism', skip_peppy=True)
            (work / 'config.yaml').write_text(yaml.safe_dump(config))
            shutil.copytree(FIXTURE / 'data', work / 'data')
            # Simulate only the FASTA size reported to the planner. The fixture stays four bases long.
            (work / 'Snakefile').write_text(
                'import os\n'
                'real_size = os.path.getsize\n'
                f'os.path.getsize = lambda path: {size} if os.fspath(path) == {str(fasta)!r} else real_size(path)\n'
                'configfile: "config.yaml"\n'
                f'include: {str(ROOT / "rnaseq/rnaseq.smk")!r}\n'
            )
            sample = next(iter(config['samples']))
            target = f'data/tmp/rnaseq/qc/fastq_screen/{sample}_screen.txt'
            env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
            proc = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', target],
                                  cwd=work, env=env, text=True, capture_output=True, timeout=60)
            log = proc.stdout + proc.stderr
            self.assertEqual(proc.returncode, 0, log)
            jobs = log.splitlines()
            config_job = next(i for i, line in enumerate(jobs) if line == 'rule fastq_screen_config:')
            dependencies = jobs[config_job + 1].strip().removeprefix('input: ').split(', ')
            expected_parts = [f'{prefix}.{part}.{expected_extension}' for part in PARTS]
            self.assertEqual([path for path in dependencies if path.startswith(str(prefix))], expected_parts)
            builds = [shlex.split(line) for line in jobs if 'bowtie2-build ' in line and str(fasta) in line]
            self.assertEqual(len(builds), int(needs_build), log)
            if needs_build:
                self.assertEqual('--large-index' in builds[0], expected_extension == 'bt2l')
                self.assertEqual(builds[0][-2:], [str(fasta), str(prefix)])
                self.assertIn('output: ' + ', '.join(expected_parts), log)
            self.assertIn('DATABASE Synthetic_Organism project_bowtie2_index/genome', log)
            screen = next(shlex.split(line) for line in jobs if '&& fastq_screen ' in line)
            self.assertIn(f'{work / prefix}.*.{expected_extension}', screen)

    def test_fresh_build_at_and_above_small_index_limit(self):
        self.check_plan(SMALL_MAX_BYTES, {}, 'bt2', True)
        self.check_plan(SMALL_MAX_BYTES + 1, {}, 'bt2l', True)

    def test_complete_forced_large_index_and_incomplete_small_index(self):
        self.check_plan(20, {'bt2l': PARTS}, 'bt2l', False)
        self.check_plan(20, {'bt2': PARTS[:1]}, 'bt2', True)

    def test_large_fasta_with_both_index_formats(self):
        self.check_plan(SMALL_MAX_BYTES + 1, {'bt2': PARTS, 'bt2l': PARTS}, 'bt2l', False)
