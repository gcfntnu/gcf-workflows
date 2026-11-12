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
PREPARED = ('genome_processed.fa', 'segemehl_genome_index.idx', 'transcriptome_trimmed_id.fa',
            'segemehl_transcriptome_index.idx', 'exons.bed', 'genome_header.sam', 'chr_size.tsv',
            'extended_mirna_annotation_6_nt.gff3', 'extended_primir_annotation_6_nt.gff3')

# Records its arguments and configuration, checks the inputs as mirflowz reads them, writes the prepared files.
FAKE_SNAKEMAKE = '''#!/usr/bin/env python3
import gzip, json, os, sys
from pathlib import Path
args = sys.argv[1:]
config = dict(line.split(': ', 1) for line in open(args[args.index('--configfile') + 1]).read().splitlines())
config = {k: json.loads(v) for k, v in config.items()}
work = Path(args[args.index('--directory') + 1])
(work / 'call.json').write_text(json.dumps({'args': args, 'config': config,
    'snakefile': Path(args[args.index('-s') + 1]).read_text(),
    'genome': gzip.open(config['genome_file'], 'rt').read(), 'gtf': gzip.open(config['gtf_file'], 'rt').read(),
    'gff': open(config['mirna_file']).read()}))
import shutil; shutil.copy(work / 'call.json', os.environ['CALL_COPY'])
if os.environ.get('FAIL'):
    Path(config['local_log']).mkdir(parents=True); (Path(config['local_log']) / 'genome_process.log').write_text('segemehl broke')
    sys.exit(1)
out = Path(config['intermediates_dir']); out.mkdir(parents=True)
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
        self.ref = self.tmp / 'ext/mirflowz/mirbase-23.1_ensembl-116/mus_musculus'

    def run_script(self, write=PREPARED, fail=False):
        env = dict(os.environ, PATH='{}:{}'.format(self.tmp / 'bin', os.environ['PATH']), MIRFLOWZ_DIR=str(self.mirflowz),
                   WRITE=','.join(write), CALL_COPY=str(self.tmp / 'call.json'), FAIL='1' if fail else '')
        result = subprocess.run([sys.executable, str(SCRIPT), '--genome', str(self.tmp / 'genome.fa'), '--gtf', str(self.tmp / 'genes.gtf.gz'),
                                 '--mirna-gff', str(self.tmp / 'mirna.gff3.gz'), '--chr-map', str(self.tmp / 'chr_map.tsv'),
                                 '--mirbase-release', '23.1', '--ensembl-release', '116', '--assembly', 'GRCm39', '--threads', '4', '--mem-gb', '30',
                                 '--output', str(self.ref / 'prepared'), '--record', str(self.ref / 'mirflowz_reference.tsv')],
                                capture_output=True, text=True, env=env, timeout=60)
        call = json.loads((self.tmp / 'call.json').read_text()) if (self.tmp / 'call.json').exists() else None
        record_path = self.ref / 'mirflowz_reference.tsv'
        record = {r['field']: r['value'] for r in csv.DictReader(record_path.read_text().splitlines(), delimiter='\t')} if record_path.exists() else {}
        return result, call, record

    def test_preparation_runs_only_the_prepare_rules_and_keeps_its_output(self):
        result, call, record = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(call['snakefile'], "include: '{}'\n".format(self.mirflowz / 'workflow/rules/prepare.smk'))
        self.assertEqual(call['args'][call['args'].index('--nolock') + 3], 'finish_prepare')
        args = call['args']
        self.assertEqual(args[args.index('--cores') + 1], '4')
        self.assertIn('generate_segemehl_index_genome=4', args)
        self.assertIn('generate_segemehl_index_genome:mem=30', args)
        self.assertIn('generate_segemehl_index_transcriptome:mem=10', args)
        self.assertEqual(args[args.index('--resources') + 1:args.index('--resources') + 3], ['threads=4', 'mem=30'])
        self.assertEqual(call['config']['scripts_dir'], str(self.mirflowz / 'scripts') + '/')
        self.assertEqual(call['config']['map_chr_file'], str(self.tmp / 'chr_map.tsv'))
        self.assertEqual(call['config']['extension'], 6)
        # mirflowz zcats the genome and GTF and reads the GFF as text: inputs are converted where needed.
        self.assertEqual((call['genome'], call['gtf']), ('>1\nACGT\n', '1\tensembl\texon\t1\t4\n'))
        # Repeated mature IDs get miRBase's replica suffix after their first occurrence; nothing else changes.
        expected = GFF.replace('ID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020958',
                               'ID=MIMAT0014876_1;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020958')
        expected = expected.replace('ID=MIMAT0014876;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020959',
                                    'ID=MIMAT0014876_2;Alias=MIMAT0014876;Name=mmu-miR-3084-5p;Derives_from=MI0020959')
        self.assertEqual(call['gff'], expected)
        self.assertEqual(record['gff_ids_suffixed'], '2')
        self.assertEqual(sorted(p.name for p in (self.ref / 'prepared').iterdir()), sorted(PREPARED))
        self.assertEqual([p.name for p in self.ref.iterdir() if p.name.startswith('.')], [])
        self.assertEqual((record['ensembl_release'], record['assembly'], record['extension']), ('116', 'GRCm39', '6'))
        self.assertEqual(len(record['genome_sha256']), 64)

    def test_missing_prepared_file_fails_naming_it(self):
        result, _, record = self.run_script(write=[n for n in PREPARED if n != 'segemehl_genome_index.idx'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('wrote no segemehl_genome_index.idx', result.stderr)
        self.assertEqual(record, {})
        self.assertFalse((self.ref / 'prepared').exists())

    def test_failed_preparation_shows_the_tool_log(self):
        result, _, record = self.run_script(fail=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('mirflowz preparation failed', result.stderr)
        self.assertIn('segemehl broke', result.stderr)
        self.assertFalse((self.ref / 'prepared').exists())

    def test_outside_the_mirflowz_image_fails_clearly(self):
        env = dict(os.environ, MIRFLOWZ_DIR=str(self.tmp / 'nowhere'))
        result = subprocess.run([sys.executable, str(SCRIPT), '--genome', 'g', '--gtf', 'g', '--mirna-gff', 'g', '--chr-map', 'c',
                                 '--mirbase-release', '23.1', '--ensembl-release', '116', '--assembly', 'GRCm39', '--mem-gb', '30',
                                 '--output', str(self.tmp / 'p'), '--record', str(self.tmp / 'r.tsv')],
                                capture_output=True, text=True, env=env, timeout=60)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('run in the mirflowz image', result.stderr)


if __name__ == '__main__':
    unittest.main()
