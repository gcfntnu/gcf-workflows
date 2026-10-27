#!/usr/bin/env python
"""Shared paired QIAseq region assignment for microbiome and metagenome."""

import argparse
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

from Bio.Seq import Seq
import yaml


def get_16s_its_primers(args):
    with open(args.libprepconfig) as handle:
        config = yaml.safe_load(handle)
    forward, reverse = {}, {}
    for region, primers in config[args.libkit]['db']['primers'].items():
        # Region names become filenames and Cutadapt adapter names. Spaces and
        # shell punctuation are literal; path separators and templates are not.
        if not isinstance(region, str) or not region or any(c in region for c in '/\\{}=\r\n\0'):
            raise ValueError(f'Invalid region name: {region!r}')
        if region == 'unknown':
            raise ValueError("The region name 'unknown' is reserved for unmatched reads")
        forward[region], reverse[region] = primers.split('-')
    if not forward:
        raise ValueError('At least one primer pair is required')
    return {'forward': forward, 'reverse': reverse}


def remove_unknown_marker(source, destination):
    """Stream the previous unmatched pairs without accumulating suffixes."""
    with source.open('rb') as reader, destination.open('wb') as writer:
        for line in reader:
            writer.write(line.replace(b':region=no_adapter', b''))


def concatenate(sources, destination, compress):
    with destination.open('wb') as output:
        if not compress:
            for source in sources:
                with source.open('rb') as handle:
                    shutil.copyfileobj(handle, output)
            return
        # One compressor for the complete stream, preserving pigz's level and
        # thread count. No shell pipeline can mask a failed producer/compressor.
        with subprocess.Popen(['pigz', '-6', '-p', '2'], stdin=subprocess.PIPE, stdout=output) as process:
            try:
                for source in sources:
                    with source.open('rb') as handle:
                        shutil.copyfileobj(handle, process.stdin)
                process.stdin.close()
            except BaseException:
                process.kill()
                raise
            returncode = process.wait()
            if returncode:
                raise subprocess.CalledProcessError(returncode, process.args)


def cutadapt_worker(args, primers):
    sample = Path(args.in1).name.replace('_R1.fastq.gz', '')
    output_dir, log_dir = Path(args.output_dir), Path(args.log_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    reverse_logs = output_dir / 'rev_comp_log'
    reverse_logs.mkdir(exist_ok=True)

    forward, reverse = primers['forward'], primers['reverse']
    log_path = log_dir / f'{sample}_qiaseq_demultiplex.log'
    # Only this invocation owns these intermediates. A retry starts from the
    # supplied reads, regardless of stale files from older implementations.
    with TemporaryDirectory(prefix='.qiaseq-', dir=output_dir) as directory, log_path.open('w') as log:
        work = Path(directory).resolve()

        def pair(region, prefix=''):
            return [work / f'{prefix}region_{region}_R{mate}.fastq' for mate in (1, 2)]

        reads = [Path(args.in1).resolve(), Path(args.in2).resolve()]
        unknown = pair('unknown')
        regions = ['unknown']
        for region, primer in forward.items():
            if len(regions) > 1:
                reads = pair('unknown', prefix='input_')
                for source, destination in zip(unknown, reads):
                    remove_unknown_marker(source, destination)

            command = [
                'cutadapt', '-g', f'{region}={primer};min_overlap={len(primer.replace("N", ""))};max_error_rate=0.25',
                '-G', f'{region}={reverse[region]};min_overlap={len(reverse[region].replace("N", ""))};max_error_rate=0.3',
                '--pair-adapters', '--untrimmed-output', unknown[0].name, '--untrimmed-paired-output', unknown[1].name,
                '--suffix', ':region={name}', '-o', pair('{name}')[0].name, '-p', pair('{name}')[1].name,
                '--minimum-length', '20', *map(str, reads),
            ]
            # Cutadapt expands {name} throughout output templates. Keep sample
            # names and parent directories out of those templates entirely.
            subprocess.run(command, cwd=work, stdout=log, check=True)

            assigned = pair(region)
            if assigned[0].exists():
                trimmed = pair(region, prefix='tmp_')
                command = [
                    'cutadapt', '-a', str(Seq(reverse[region]).reverse_complement()),
                    '-A', str(Seq(primer).reverse_complement()), '--pair-adapters',
                    '-o', trimmed[0].name, '-p', trimmed[1].name, '--minimum-length', '20', *map(str, assigned),
                ]
                with (reverse_logs / f'{sample}_{region}_rev_comp.log').open('w') as reverse_log:
                    subprocess.run(command, cwd=work, stdout=reverse_log, check=True)
                for source, destination in zip(trimmed, assigned):
                    source.replace(destination)
            regions.append(region)

        pairs = [pair(region) for region in regions if pair(region)[0].exists()]
        suffix = '.gz' if args.compress else ''
        staged = work / 'final'
        staged.mkdir()
        combined = [staged / f'{sample}_R{mate}.fastq{suffix}' for mate in (1, 2)]
        for mate, destination in enumerate(combined):
            concatenate([reads[mate] for reads in pairs], destination, args.compress)
        # Publish only after both mates and compression have succeeded.
        for source in combined:
            source.replace(output_dir / source.name)
    print(f'finished sample: {sample}')


def demultiplex_16s_its(args):
    cutadapt_worker(args, get_16s_its_primers(args))


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--libkit', required=True, help='Library preparation kit name.')
    parser.add_argument('--libprepconfig', required=True, help='Path to libprepconfig.')
    parser.add_argument('--in1', required=True, help='Input R1.')
    parser.add_argument('--in2', required=True, help='Input R2.')
    parser.add_argument('--compress', action='store_true', help='Compress output files with pigz.')
    parser.add_argument('--output-dir', required=True, help='Path to output directory.')
    parser.add_argument('--log-dir', required=True, help='Path to log directory.')
    demultiplex_16s_its(parser.parse_args())


if __name__ == '__main__':
    main()
