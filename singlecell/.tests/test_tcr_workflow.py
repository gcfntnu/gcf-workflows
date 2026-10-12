import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'rules/immune/scripts/tcr_workflow.py'
SPEC = importlib.util.spec_from_file_location('tcr_workflow', SCRIPT)
w = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(w)


class TcrWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def parents(self):
        config = {'immune': {'sublibrary_order': ['sublib1', 'sublib2'], 'wt_parent': {
            'root': str(self.root), 'sublibraries': {'sublib1': 'sub1', 'sublib2': 'sub2'}, 'combined': 'combined'}}}
        for name in ('sub1', 'sub2', 'combined'):
            p = self.root / name
            cells = 'bc_wells,sample\n01_01_01,S1\n' if name != 'combined' else 'bc_wells,sample\n01_01_01__s1,S1\n01_01_01__s2,S1\n'
            for sample in ('all-sample', 'S1'):
                d = p / sample / 'DGE_filtered'; d.mkdir(parents=True)
                (d / 'cell_metadata.csv').write_text(cells)
            proc = {'header': {'ver_number': '1.8.2'}, 'kit': {'kit_name': 'WT'}, 'chemistry': {'name': 4},
                    'run': {'use_case': 'normal'}, 'sample': {'samples': [
                        {'name': 'S1', 'wspec': 'A1-A2', 'meta': False}, {'name': 'all-sample', 'wspec': 'A1-A2', 'meta': True}]},
                    'combine': {'sublib_number': 2, 'sublib_s1_name': 'sub1', 'sublib_s2_name': 'sub2'}}
            w.save_json(p / 'process/run_proc_def.json', proc)
        return config

    def test_parent_union_and_order(self):
        config = self.parents()
        w.validate_parents(config, self.root / 'valid.json')
        self.assertEqual(json.loads((self.root / 'valid.json').read_text())['combined']['cells'], 2)
        config['immune']['sublibrary_order'].reverse()
        with self.assertRaisesRegex(ValueError, 'order mismatch'):
            w.validate_parents(config, self.root / 'bad.json')

    def test_parent_sample_assignment_mismatch(self):
        config = self.parents()
        p = self.root / 'combined/all-sample/DGE_filtered/cell_metadata.csv'
        p.write_text(p.read_text().replace('S1', 'S2'))
        with self.assertRaisesRegex(ValueError, 'sample assignments'):
            w.validate_parents(config, self.root / 'bad.json')

    def test_link_retains_wt_cells_and_counts_pairing(self):
        meta = self.root / 'cells.csv'
        meta.write_text('bc_wells,sample\n01__s1,S1\n02__s1,S1\n')
        airr = self.root / 'airr.tsv.gz'
        with gzip.open(airr, 'wt') as f:
            f.write('cell_id\tsequence_id\tlocus\tproductive\n01__s1\tc1\tTRA\tT\n01__s1\tc2\tTRB\tT\n')
        w.link_wt(airr, meta, self.root / 'linked.tsv.gz', self.root / 'qc.json')
        qc = json.loads((self.root / 'qc.json').read_text())
        self.assertEqual(qc['tcr_cells'], 1)
        self.assertEqual(qc['wt_cells_without_tcr'], 1)
        with gzip.open(self.root / 'linked.tsv.gz', 'rt') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['paired_TRA_TRB'], '1')
        self.assertEqual(rows[1]['tcr_contigs'], '0')
        with gzip.open(airr, 'wt') as f:
            f.write('cell_id\tsequence_id\tlocus\tproductive\n99__s1\tc1\tTRA\tT\n')
        with self.assertRaisesRegex(ValueError, 'outside the WT parent'):
            w.link_wt(airr, meta, self.root / 'linked2.tsv.gz', self.root / 'qc2.json')
        self.assertEqual(json.loads((self.root / 'qc2.json').read_text())['unmatched_tcr_cells'], ['99__s1'])

    @unittest.skipUnless(importlib.util.find_spec('pysam'), 'requires pinned immune image')
    def test_full_fastq_pairing_and_truncation(self):
        paths = [str(self.root / 'R1.gz'), str(self.root / 'R2.gz')]
        for i, p in enumerate(paths, 1):
            with gzip.open(p, 'wt') as f:
                f.write(f'@read1/{i}\nACGT\n+\nFFFF\n')
        digests = [hashlib.md5(Path(p).read_bytes()).hexdigest() for p in paths]
        w.check_fastqs(paths, digests, self.root / 'fq.json')
        self.assertEqual(json.loads((self.root / 'fq.json').read_text())['read_pairs'], 1)
        with gzip.open(paths[1], 'wt') as f:
            f.write('@different/2\nACGT\n+\nFFFF\n')
        digests = [hashlib.md5(Path(p).read_bytes()).hexdigest() for p in paths]
        with self.assertRaisesRegex(ValueError, 'IDs differ'):
            w.check_fastqs(paths, digests, self.root / 'bad.json')
        with self.assertRaisesRegex(ValueError, 'MD5 mismatch'):
            w.check_fastqs(paths, ['0' * 32] * 2, self.root / 'bad.json')
        with gzip.open(paths[1], 'wt') as f:
            f.write('@read1/2\nACGT\n+\nFF\n')
        digests = [hashlib.md5(Path(p).read_bytes()).hexdigest() for p in paths]
        with self.assertRaises(ValueError):
            w.check_fastqs(paths, digests, self.root / 'bad.json')

    def test_run_refuses_false_success_and_existing_results(self):
        output = self.root / 'run'
        output.mkdir()  # Snakemake creates output parents before invoking the command.
        args = argparse.Namespace(output_dir=str(output), parent='/parent', threads=1, sublibraries=None,
                                  name='sub1', fastqs=['a', 'b'], reference='default', image='image@digest',
                                  workflow_commit='commit', dryrun=True, provenance=self.root / 'provenance.json',
                                  success=output / 'complete.json')
        def fake_run(*a, **kw):
            (output / 'split-pipe_v1_8_2.log').write_text('Problems encountered')
            return argparse.Namespace(returncode=0)
        with patch.object(w.subprocess, 'run', side_effect=fake_run), self.assertRaisesRegex(RuntimeError, 'lacks successful'):
            w.run_splitpipe(args)
        self.assertEqual(json.loads(args.provenance.read_text())['exit_status'], 1)
        self.assertFalse(args.success.exists())
        with self.assertRaisesRegex(ValueError, 'Refusing to reuse'):
            w.run_splitpipe(args)


if __name__ == '__main__':
    unittest.main()
