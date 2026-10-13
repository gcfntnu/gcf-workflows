#!/usr/bin/env python3
"""Parent-linked Split-Pipe orchestration checks; all paths supplied explicitly."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def file_hash(path, algorithm):
    digest = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def sha256(path):
    return file_hash(path, 'sha256')


def read_cells(path):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not {'bc_wells', 'sample'} <= rows[0].keys():
        raise ValueError(f'Missing WT cells or columns: {path}')
    if len({row['bc_wells'] for row in rows}) != len(rows):
        raise ValueError(f'Duplicate WT barcodes: {path}')
    return rows


def validate_parents(settings, output):
    immune = settings['immune']
    root = Path(immune['wt_parent']['root'])
    order = immune['sublibrary_order']
    manifest = {'settings': settings, 'parents': {}, 'sample_wells': None}
    all_cells = {}
    for index, name in enumerate(order, 1):
        parent = root / immune['wt_parent']['sublibraries'][name]
        definition = parent / 'process/run_proc_def.json'
        proc = json.loads(definition.read_text())
        if proc['header']['ver_number'] != '1.8.2' or proc['kit']['kit_name'] != 'WT' or proc['chemistry']['name'] != 4:
            raise ValueError(f'Unsupported parent version/kit/chemistry: {parent}')
        if proc['run']['use_case'] != 'normal':
            raise ValueError(f'Expected a WT parent: {parent}')
        sample_wells = {s['name']: s['wspec'] for s in proc['sample']['samples'] if not s['meta']}
        if manifest['sample_wells'] is None:
            manifest['sample_wells'] = sample_wells
        if sample_wells != manifest['sample_wells']:
            raise ValueError(f'Inconsistent WT sample/well definitions: {parent}')
        files = [definition]
        for sample in proc['sample']['samples']:
            path = parent / sample['name'] / 'DGE_filtered/cell_metadata.csv'
            if not path.is_file():
                raise ValueError(f'Missing WT filtered metadata: {path}')
            files.append(path)
        cells = read_cells(parent / 'all-sample/DGE_filtered/cell_metadata.csv')
        for cell in cells:
            all_cells[f"{cell['bc_wells']}__s{index}"] = cell['sample']
        manifest['parents'][name] = {'path': str(parent), 'cells': len(cells),
                                     'files': {str(p): sha256(p) for p in files}}
    combined = root / immune['wt_parent']['combined']
    definition = combined / 'process/run_proc_def.json'
    proc = json.loads(definition.read_text())
    comb = proc['combine']
    recorded = [comb[f'sublib_s{i}_name'] for i in range(1, comb['sublib_number'] + 1)]
    expected = [Path(immune['wt_parent']['sublibraries'][name]).name for name in order]
    if recorded != expected:
        raise ValueError(f'WT combine order mismatch: {recorded} != {expected}')
    path = combined / 'all-sample/DGE_filtered/cell_metadata.csv'
    combined_cells = {r['bc_wells']: r['sample'] for r in read_cells(path)}
    if combined_cells != all_cells:
        raise ValueError('WT combined barcodes/sample assignments do not match the ordered parents')
    manifest['combined'] = {'cells': len(combined_cells), 'files': {str(p): sha256(p) for p in (definition, path)}}
    save_json(output, manifest)


def run_splitpipe(args):
    output = Path(args.output_dir).resolve()
    if output.exists() and any(not p.is_dir() or p.is_symlink() for p in output.rglob('*')):
        raise ValueError(f'Refusing to reuse output directory: {output}; preserve/rename the failed attempt first')
    output.parent.mkdir(parents=True, exist_ok=True)
    base = ['split-pipe', '--mode', 'comb' if args.sublibraries else 'all', '--output_dir', str(output),
            '--parent_dir', args.parent, '--nthreads', str(args.threads), '--rseed', '1234']
    if args.sublibraries:
        base += ['--sublibraries', *args.sublibraries]
    else:
        base += ['--tcr_analysis', '--immune_genome', 'human', '--run_name', args.name,
                 '--fq1', args.fastqs[0], '--fq2', args.fastqs[1]]
        if args.reference == 'imgt':
            base += ['--use_imgt_db']
    record = {'commands': [], 'image': args.image, 'workflow_commit': args.workflow_commit,
              'started': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    code = 1
    try:
        # Combine setup requires completed sublibrary outputs: check it immediately before combine.
        commands = [base + ['--dryrun'], base] if args.sublibraries and not args.dryrun else [base + (['--dryrun'] if args.dryrun else [])]
        for command in commands:
            print(json.dumps(command), flush=True)
            record['commands'].append(command)
            result = subprocess.run(command, check=False)
            if result.returncode:
                code = result.returncode
                raise RuntimeError(f'Split-Pipe exited {result.returncode}')
            log = output / 'split-pipe_v1_8_2.log'
            text = log.read_text()
            if 'Finished successfully (0 problems)' not in text:
                raise RuntimeError('Split-Pipe log lacks successful completion')
            if 'Transcriptome file provided is incorrect' in text:
                raise RuntimeError('Split-Pipe fell back from parent filtering')
        if not args.dryrun:
            for fraction in ('filtered', 'unfiltered'):
                folder = output / 'all-sample' / f'TCR_{fraction}'
                airr = folder / 'tcr_annotation_airr.tsv.gz'
                if not airr.is_file():
                    raise RuntimeError(f'Missing expected compressed AIRR output: {airr}')
                with gzip.open(airr, 'rt') as handle:
                    fields = next(csv.reader(handle, delimiter='\t'))
                if not {'cell_id', 'sequence_id', 'locus', 'productive'} <= set(fields):
                    raise RuntimeError(f'Missing AIRR columns: {airr}')
                for filename in ('barcode_report.tsv', 'clonotype_frequency.tsv', 'tcr_contigs.fa.gz'):
                    if not (folder / filename).is_file():
                        raise RuntimeError(f'Missing immune output: {folder / filename}')
            if not (output / 'all-sample_analysis_summary.html').is_file():
                raise RuntimeError('Missing all-sample report')
        code = 0
    finally:
        record['exit_status'] = code
        record['finished'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        save_json(args.provenance, record)
    if code == 0:
        save_json(args.success, {'output_dir': str(output), 'dryrun': args.dryrun})


def link_wt(airr, metadata, output, summary):
    cells = read_cells(metadata)
    wt = {row['bc_wells']: row for row in cells}
    counts = collections.defaultdict(lambda: collections.Counter())
    unmatched = set()
    seen = set()
    with gzip.open(airr, 'rt') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'cell_id', 'sequence_id', 'locus', 'productive'} <= set(reader.fieldnames or []):
            raise ValueError('Missing required AIRR columns')
        for row in reader:
            if row['sequence_id'] in seen:
                raise ValueError('Duplicate AIRR sequence_id')
            seen.add(row['sequence_id'])
            barcode = row['cell_id']
            if barcode not in wt:
                unmatched.add(barcode)
            counts[barcode]['tcr_contigs'] += 1
            if row['productive'].lower() in ('true', 't', '1'):
                counts[barcode]['productive_contigs'] += 1
                if row['locus'] in ('TRA', 'TRB'):
                    counts[barcode][f"productive_{row['locus']}"] += 1
    detected = set(counts)
    extras = ['tcr_contigs', 'productive_contigs', 'productive_TRA', 'productive_TRB', 'paired_TRA_TRB']
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cells[0]) + extras, delimiter='\t')
        writer.writeheader()
        for row in cells:
            values = counts[row['bc_wells']]
            row.update({key: values[key] for key in extras[:-1]})
            row['paired_TRA_TRB'] = int(values['productive_TRA'] > 0 and values['productive_TRB'] > 0)
            writer.writerow(row)
    save_json(summary, {'wt_cells': len(wt), 'tcr_cells': len(detected), 'matched_tcr_cells': len(detected & set(wt)),
                        'unmatched_tcr_cells': sorted(unmatched), 'wt_cells_without_tcr': len(set(wt) - detected)})
    if unmatched:
        raise ValueError('Filtered TCR has cells outside the WT parent; inspect the saved overlap report')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('validate-parents'); p.add_argument('--settings', required=True); p.add_argument('--output', required=True)
    p = sub.add_parser('run')
    for option in ('output-dir', 'parent', 'provenance', 'success', 'image', 'workflow-commit'):
        p.add_argument('--' + option, required=True)
    p.add_argument('--threads', type=int, default=12); p.add_argument('--name'); p.add_argument('--fastqs', nargs=2)
    p.add_argument('--reference', choices=['default', 'imgt'], default='default'); p.add_argument('--sublibraries', nargs='+')
    p.add_argument('--dryrun', action='store_true')
    p = sub.add_parser('link-wt')
    for option in ('airr', 'metadata', 'output', 'summary'):
        p.add_argument('--' + option, required=True)
    args = parser.parse_args()
    if args.action == 'validate-parents': validate_parents(json.loads(args.settings), args.output)
    elif args.action == 'run': run_splitpipe(args)
    elif args.action == 'link-wt': link_wt(args.airr, args.metadata, args.output, args.summary)


if __name__ == '__main__':
    main()
