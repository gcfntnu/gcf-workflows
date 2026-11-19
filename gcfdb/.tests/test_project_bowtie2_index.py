"""Use the selected project FASTA bytes for FastQ Screen's nonhuman Bowtie2 index."""
import hashlib
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
BUILDER = ROOT / 'gcfdb/scripts/build_project_screen_index.sh'
PARTS = ('1', '2', '3', '4', 'rev.1', 'rev.2')
SMALL_MAX_BYTES = 4 * 1024**3 - 200
SCREEN = Path('data/tmp/rnaseq/qc/fastq_screen')


class ProjectBowtie2Tests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.ref = Path('refs/panel one')
        self.fasta = self.ref / 'fasta/genome.fa'

    def project(self, size=None, human_cache=None):
        if human_cache:
            (self.work / human_cache).mkdir(parents=True)
            (self.work / self.ref.parent).mkdir(parents=True)
            (self.work / self.ref).symlink_to(self.work / human_cache, target_is_directory=True)
        (self.work / self.fasta.parent).mkdir(parents=True)
        (self.work / self.fasta).write_text('>synthetic\nACGT\n')
        config = yaml.safe_load((FIXTURE / 'config.yaml').read_text())
        config.update(ref_dir=str(self.ref), organism='synthetic_organism', skip_peppy=True)
        (self.work / 'config.yaml').write_text(yaml.safe_dump(config))
        shutil.copytree(FIXTURE / 'data', self.work / 'data', dirs_exist_ok=True)
        # Change only the FASTA size seen by the planner, using a four-base fixture.
        mock_size = '' if size is None else (
            'import os\n'
            'real_size = os.path.getsize\n'
            f"fasta_inode = os.stat({str(self.fasta)!r}).st_ino\n"
            f'os.path.getsize = lambda p: {size} if isinstance(p, int) and os.fstat(p).st_ino == fasta_inode else real_size(p)\n'
        )
        (self.work / 'Snakefile').write_text(mock_size + 'configfile: "config.yaml"\n'
                                             + f'include: {str(ROOT / "rnaseq/rnaseq.smk")!r}\n')
        self.target = SCREEN / f'{next(iter(config["samples"]))}_screen.txt'

    def index_dir(self):
        return SCREEN / 'indexes' / hashlib.sha256((self.work / self.fasta).read_bytes()).hexdigest()

    def seed_index(self, directory, extension, parts=PARTS):
        (self.work / directory).mkdir(parents=True, exist_ok=True)
        for part in parts:
            (self.work / directory / f'genome.{part}.{extension}').touch()

    def check_plan(self, extension, needs_build):
        env = {k: v for k, v in os.environ.items() if not k.startswith('GCF_')}
        proc = subprocess.run([sys.executable, '-m', 'snakemake', '-n', '-p', '--cores', '1', self.target],
                              cwd=self.work, env=env, text=True, capture_output=True, timeout=60)
        log = proc.stdout + proc.stderr
        self.assertEqual(proc.returncode, 0, log)
        jobs = log.splitlines()
        config_job = jobs.index('rule fastq_screen_config:')
        dependencies = jobs[config_job + 1].strip().removeprefix('input: ').split(', ')
        directory = self.index_dir()
        expected = [str(directory / f'genome.{part}.{extension}') for part in PARTS]
        self.assertEqual([path for path in dependencies if path.startswith(str(SCREEN / 'indexes'))], expected)
        builds = [shlex.split(line.removeprefix('Shell command: ')) for line in jobs if str(BUILDER) in line]
        self.assertEqual(len(builds), int(needs_build), log)
        if needs_build:
            self.assertEqual(builds[0], ['bash', str(BUILDER), str(self.fasta), directory.name,
                                         extension, '1', str(directory)])
            self.assertIn('output: ' + ', '.join(expected), log)
        self.assertIn('DATABASE Synthetic_Organism project_bowtie2_index/genome', log)
        screen = next(shlex.split(line) for line in jobs if '&& fastq_screen ' in line)
        self.assertIn(f'{self.work / directory}/genome.*.{extension}', screen)
        return log

    def test_human_control_uses_stock_sequence_when_project_occupies_its_cache(self):
        release = yaml.safe_load((ROOT / 'main.config').read_text())['db']['ensembl']['release']
        cache_path = Path('ensembl') / f'release-{release}' / 'homo_sapiens/GRCh38'
        self.project(human_cache=Path('data/ext') / cache_path)
        log = self.check_plan('bt2', True)
        control = Path('data/ext/fastq_screen') / cache_path
        fasta, index = control / 'fasta/genome.fa', control / 'index/genome/bowtie2/genome'
        url = f'https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
        self.assertIn(f'{url} | gunzip -c > {fasta}', log)
        build = next(shlex.split(line) for line in log.splitlines() if 'bowtie2-build ' in line and str(fasta) in line)
        self.assertEqual(build[-2:], [str(fasta), str(index)])
        records = next(shlex.split(line) for line in log.splitlines() if 'DATABASE Human ' in line)
        human = next(record for record in records if record.startswith('DATABASE Human '))
        self.assertTrue(human.endswith('/' + str(index)), human)

    def test_partial_small_index_at_size_boundary(self):
        self.project(SMALL_MAX_BYTES)
        self.seed_index(self.index_dir(), 'bt2', PARTS[:1])
        self.check_plan('bt2', True)

    def test_large_index_above_size_boundary(self):
        self.project(SMALL_MAX_BYTES + 1)
        self.check_plan('bt2l', True)
        self.seed_index(self.index_dir(), 'bt2l')
        self.check_plan('bt2l', False)

    def test_changed_fasta_with_preserved_mtime_bypasses_old_complete_indexes(self):
        self.project()
        self.seed_index(self.ref / 'index/genome/bowtie2', 'bt2')
        self.check_plan('bt2', True)
        old = self.index_dir()
        self.seed_index(old, 'bt2')
        self.check_plan('bt2', False)
        fasta = self.work / self.fasta
        stat = fasta.stat()
        fasta.write_text('>synthetic\nAAAA\n')
        os.utime(fasta, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertNotEqual(self.index_dir(), old)
        self.check_plan('bt2', True)
        self.assertTrue((self.work / old / 'genome.1.bt2').is_file())

    def test_builder_uses_verified_private_copy(self):
        self.project()
        directory = self.work / self.index_dir()
        directory.mkdir(parents=True)
        executable = self.work / 'bowtie2-build'
        # Modify the curated FASTA after verification, then record the actual builder input.
        executable.write_text(
            '#!/usr/bin/env bash\nset -eu\n'
            '[[ "$1" == --large-index ]]\n'
            'echo changed > "$CURATED"\n'
            'for part in 1 2 3 4 rev.1 rev.2; do cp -- "${@: -2:1}" "${@: -1}.$part.bt2l"; done\n'
        )
        executable.chmod(0o755)
        fasta = self.work / self.fasta
        expected = fasta.read_bytes()
        env = dict(os.environ, PATH=str(self.work) + os.pathsep + os.environ['PATH'], CURATED=str(fasta))
        args = ['bash', BUILDER, fasta, directory.name, 'bt2l', '1', directory]
        subprocess.run(args, env=env, check=True)
        for part in PARTS:
            self.assertEqual((directory / f'genome.{part}.bt2l').read_bytes(), expected)
        changed = subprocess.run(args, env=env, capture_output=True, text=True)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn(f'FASTA changed after planning: {fasta}', changed.stderr)
        self.assertFalse(list(directory.glob('.build-*')))
