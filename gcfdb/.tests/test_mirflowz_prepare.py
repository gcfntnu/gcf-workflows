"""mirflowz_prepare.py: mirflowz's preparation rules run without its top-level Snakefile; snakemake is mocked."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'gcfdb/scripts/mirflowz_prepare.py'
GFF = ('##gff-version 3\n'
       'chr1\t.\tmiRNA_primary_transcript\t1\t60\t.\t+\t.\tID=MI0014047;Alias=MI0014047;Name=mmu-mir-3084-1\n'
       'chr1\t.\tmiRNA\t10\t31\t.\t+\t.\tID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0014047\n'
       'chr2\t.\tmiRNA_primary_transcript\t1\t60\t.\t+\t.\tID=MI0020958;Alias=MI0020958;Name=mmu-mir-3084-2\n'
       'chr2\t.\tmiRNA\t10\t31\t.\t+\t.\tID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020958\n'
       'chr3\t.\tmiRNA\t10\t31\t.\t+\t.\tID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020959\n')
GENOME = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
          'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv')
ANNOTATION = ('extended_mirna_annotation_6_nt.gff3', 'extended_primir_annotation_6_nt.gff3')

# Records its arguments and configuration, checks the inputs as mirflowz reads them, writes the prepared files.
FAKE_SNAKEMAKE = '''#!/usr/bin/env python3
import gzip, json, os, sys
from pathlib import Path
args = sys.argv[1:]
config = dict(line.split(': ', 1) for line in open(args[args.index('--configfile') + 1]).read().splitlines())
config = {k: json.loads(v) for k, v in config.items()}
work = Path(args[args.index('--directory') + 1])
out = Path(config['intermediates_dir'])
read = lambda key, opener: opener(config[key], 'rt').read() if config[key] != 'NA' else None
(work / 'call.json').write_text(json.dumps({'args': args, 'config': config,
    'snakefile': Path(args[args.index('-s') + 1]).read_text(),
    'genome': read('genome_file', gzip.open), 'gtf': read('gtf_file', gzip.open), 'gff': read('mirna_file', open),
    'linked': sorted(p.name for p in out.iterdir() if p.is_symlink()) if out.exists() else []}))
import shutil; shutil.copy(work / 'call.json', os.environ['CALL_COPY'])
if os.environ.get('FAIL'):
    Path(config['local_log']).mkdir(parents=True); (Path(config['local_log']) / 'genome_process.log').write_text('segemehl broke')
    sys.exit(1)
out.mkdir(parents=True, exist_ok=True)
for name in os.environ['WRITE'].split(','):
    (out / name).write_text('prepared ' + name)
'''


class PrepareTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.mirflowz = self.tmp / 'opt/mirflowz'
        (self.mirflowz / 'workflow/rules').mkdir(parents=True)
        (self.mirflowz / 'workflow/rules/prepare.smk').write_text('# mirflowz prepare rules\n')
        (self.tmp / 'bin').mkdir()
        (self.tmp / 'bin/snakemake').write_text(FAKE_SNAKEMAKE)
        (self.tmp / 'bin/snakemake').chmod(0o755)
        (self.tmp / 'genome.fa').write_text('>1\nACGT\n')
        (self.tmp / 'genes.gtf.gz').write_bytes(gzip.compress(b'1\tensembl\texon\t1\t4\n'))
        (self.tmp / 'mirna.gff3.gz').write_bytes(gzip.compress(GFF.encode()))
        (self.tmp / 'chr_map.tsv').write_text('chr1\t1\n')
        (self.tmp / 'chr_size.tsv').write_text('1\t4\n')
        self.ref = self.tmp / 'ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus'

    def run_script(self, part, write=None, fail=False, mirflowz_dir=None):
        env = dict(os.environ, PATH='{}:{}'.format(self.tmp / 'bin', os.environ['PATH']), MIRFLOWZ_DIR=str(mirflowz_dir or self.mirflowz),
                   WRITE=','.join(write or (GENOME if part == 'genome' else ANNOTATION)), CALL_COPY=str(self.tmp / 'call.json'),
                   FAIL='1' if fail else '')
        if part == 'genome':
            output, record_path = self.ref / 'prepared', self.ref / 'mirflowz_reference.tsv'
            options = ['--genome', str(self.tmp / 'genome.fa'), '--gtf', str(self.tmp / 'genes.gtf.gz'), '--ensembl-release', '116',
                       '--assembly', 'GRCm39', '--threads', '4', '--mem-mb', '30720']
        else:
            output, record_path = self.ref / 'gff-abc/annotation', self.ref / 'gff-abc/mirflowz_annotation.tsv'
            options = ['--mirna-gff', str(self.tmp / 'mirna.gff3.gz'), '--chr-map', str(self.tmp / 'chr_map.tsv'),
                       '--chr-size', str(self.tmp / 'chr_size.tsv')]
        result = subprocess.run([sys.executable, str(SCRIPT), part, *options, '--mirbase-release', '23.1', '--output', str(output),
                                 '--record', str(record_path)], capture_output=True, text=True, env=env, timeout=60)
        call = json.loads((self.tmp / 'call.json').read_text()) if (self.tmp / 'call.json').exists() else None
        record = {r['field']: r['value'] for r in csv.DictReader(record_path.read_text().splitlines(), delimiter='\t')} if record_path.exists() else {}
        return result, call, record, output

    def targets(self, call):
        args = call['args']
        return [Path(a).name for a in args[args.index('none') + 1:args.index('--drop-metadata')]]

    def test_genome_part_builds_the_genome_files_without_the_gff(self):
        result, call, record, output = self.run_script('genome')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(call['snakefile'], "include: '{}'\n".format(self.mirflowz / 'workflow/rules/prepare.smk'))
        self.assertEqual(self.targets(call), list(GENOME))
        args = call['args']
        self.assertNotIn('--allowed-rules', args)
        self.assertEqual(args[args.index('--cores') + 1], '4')
        self.assertIn('generate_segemehl_index_genome=4', args)
        self.assertIn('generate_segemehl_index_genome:mem=30', args)
        self.assertIn('generate_segemehl_index_transcriptome:mem=10', args)
        self.assertEqual(args[args.index('--resources') + 1:args.index('--resources') + 3], ['threads=4', 'mem=30'])
        self.assertEqual(call['config']['scripts_dir'], str(self.mirflowz / 'scripts') + '/')
        self.assertEqual((call['config']['mirna_file'], call['config']['map_chr_file']), ('NA', 'NA'))  # the GFF is the other part's
        # mirflowz zcats the genome and GTF: inputs are converted where needed.
        self.assertEqual((call['genome'], call['gtf']), ('>1\nACGT\n', '1\tensembl\texon\t1\t4\n'))
        self.assertEqual(sorted(p.name for p in output.iterdir()), sorted(GENOME))
        self.assertEqual([p.name for p in self.ref.iterdir() if p.name.startswith('.')], [])
        self.assertEqual((record['ensembl_release'], record['assembly'], record['prepared_files']), ('116', 'GRCm39', ','.join(GENOME)))
        self.assertEqual(len(record['genome_sha256']), 64)
        self.assertNotIn('mirna_gff', record)

    def test_annotation_part_runs_only_the_gff_rules_on_the_genome_parts_sizes(self):
        result, call, record, output = self.run_script('annotation')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.targets(call), list(ANNOTATION))
        args = call['args']
        self.assertEqual(args[args.index('--allowed-rules') + 1:], ['map_chr_names', 'extend_mirs_annotations'])
        self.assertEqual(call['linked'], ['chr_size.tsv'])  # the genome part's, so no genome rule could run
        self.assertEqual((call['config']['genome_file'], call['config']['gtf_file']), ('NA', 'NA'))
        self.assertEqual(call['config']['map_chr_file'], str(self.tmp / 'chr_map.tsv'))
        self.assertEqual(call['config']['extension'], 6)
        # Repeated mature IDs get miRBase's replica suffix after their first occurrence; nothing else changes.
        expected = GFF.replace('ID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020958',
                               'ID=MIMAT0014876_1;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020958')
        expected = expected.replace('ID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020959',
                                    'ID=MIMAT0014876_2;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020959')
        self.assertEqual(call['gff'], expected)
        self.assertEqual(sorted(p.name for p in output.iterdir()), sorted(ANNOTATION))  # not the linked sizes
        self.assertEqual((record['gff_ids_suffixed'], record['extension']), ('2', '6'))
        self.assertEqual([len(record[k]) for k in ('mirna_gff_sha256', 'chr_map_sha256', 'chr_size_sha256')], [64] * 3)
        self.assertEqual([p.name for p in output.parent.iterdir() if p.name.startswith('.')], [])

    def test_missing_prepared_file_fails_naming_it(self):
        result, _, record, output = self.run_script('genome', write=[n for n in GENOME if n != 'segemehl_genome_index.idx'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('wrote no segemehl_genome_index.idx', result.stderr)
        self.assertEqual(record, {})
        self.assertFalse(output.exists())
        result, _, record, output = self.run_script('annotation', write=ANNOTATION[:1])
        self.assertIn('wrote no ' + ANNOTATION[1], result.stderr)
        self.assertFalse(output.exists())

    def test_failed_preparation_shows_the_tool_log(self):
        result, _, record, output = self.run_script('genome', fail=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('mirflowz preparation failed', result.stderr)
        self.assertIn('segemehl broke', result.stderr)
        self.assertFalse(output.exists())

    def test_outside_the_mirflowz_image_fails_clearly(self):
        for part in ('genome', 'annotation'):
            result, _, _, _ = self.run_script(part, mirflowz_dir=self.tmp / 'nowhere')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('run in the mirflowz image', result.stderr)
