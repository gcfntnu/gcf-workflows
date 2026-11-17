#!/usr/bin/env python
"""Spike-in reads per sample for MultiQC's General Statistics (#336): reads matching a spike-in and their percentage of
the filtered reads, from the configured spike-in route.

    spikein_mqc.py bowtie --counts SAMPLE.tsv [...] --logs SAMPLE.log [...] > spikein_mqc.yaml
        Bowtie route: the sum of each sample's count table (every spike-in and __ambiguous), of the reads bowtie
        processed (`# reads processed: N` in its log), sample names from the file names.
    spikein_mqc.py unitas SAMPLE/unitas.annotation_summary.txt [...] > spikein_mqc.yaml
        Native route (filter.spikein.quantifier: native), from UNITAS, which delivery uses: the `Spikein` class of each summary, of the sum of its top-level classes (every read UNITAS
        annotated or not), sample names from the summaries' directories.

Writes MultiQC custom content (plot_type: generalstats) as YAML; values are plain numbers, names quoted.
"""
import argparse
import json
from pathlib import Path
import re
import sys

COLUMNS = {
    'bowtie': 'Filtered reads matching a spike-in exactly, forward strand (bowtie)',
    'unitas': 'Filtered reads UNITAS annotated as spike-in (mismatches and 3\' tails allowed)',
}


def bowtie_totals(counts, logs):
    processed = {}
    for path in logs:
        match = re.search(r'^# reads processed: (\d+)$', Path(path).read_text(), re.M)
        if not match:
            raise ValueError('{}: no "# reads processed" line'.format(path))
        processed[Path(path).name[:-len('.log')]] = int(match.group(1))
    totals = {}
    for path in counts:
        sample = Path(path).name[:-len('.tsv')]
        lines = Path(path).read_text().splitlines()
        if not lines or lines[0] != 'spikein\tcount':
            raise ValueError('{}: not a spike-in count table'.format(path))
        if sample not in processed:
            raise ValueError('{}: no bowtie log for sample {}'.format(path, sample))
        totals[sample] = (sum(int(line.split('\t')[1]) for line in lines[1:]), processed[sample])
    return totals


def unitas_totals(summaries):
    totals = {}
    for path in summaries:
        spikein, total = None, 0.0
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line[0].isspace():
                continue
            name, value = line.rsplit('\t', 1)
            total += float(value)
            if name == 'Spikein':
                spikein = float(value)
        totals[Path(path).parent.name] = (spikein or 0.0, total)
    return totals


def mqc(route, totals):
    data = {sample: {'spikein_reads': round(n, 2), 'spikein_percent': round(100.0 * n / total, 2) if total else 0.0}
            for sample, (n, total) in sorted(totals.items())}
    lines = ['id: smallrna_spikein', 'section_name: Spike-ins', 'plot_type: generalstats', 'pconfig:',
             '  - spikein_reads:', '      title: Spike-in reads', '      description: {}'.format(json.dumps(COLUMNS[route])),
             "      format: '{:,.0f}'", '      scale: Purples', '      placement: 116',
             '  - spikein_percent:', "      title: '% Spike-in'",
             '      description: Spike-in reads as a percentage of the filtered reads', '      min: 0', '      max: 100',
             "      suffix: '%'", "      format: '{:,.1f}'", '      scale: Purples', '      placement: 117', 'data:']
    lines += ['  {}: {}'.format(json.dumps(sample), json.dumps(values)) for sample, values in data.items()]
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='route', required=True)
    bowtie = sub.add_parser('bowtie')
    bowtie.add_argument('--counts', nargs='+', required=True)
    bowtie.add_argument('--logs', nargs='+', required=True)
    sub.add_parser('unitas').add_argument('summaries', nargs='+')
    args = parser.parse_args()
    try:
        totals = bowtie_totals(args.counts, args.logs) if args.route == 'bowtie' else unitas_totals(args.summaries)
    except (OSError, ValueError) as error:
        sys.exit('spikein_mqc.py: {}'.format(error))
    sys.stdout.write(mqc(args.route, totals))


if __name__ == '__main__':
    main()
