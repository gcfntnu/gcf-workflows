#!/usr/bin/env python
"""Merge a sample's two cutadapt JSON reports (3' adapter pass, filter pass) into one report in cutadapt's schema.

Read counts: `input` from the adapter pass, `output` from the filter pass, filtered categories summed over both passes,
`read1_with_adapter` (3' adapter found) from the adapter pass; base pair counts likewise; `adapters_read1` (trimmed
length statistics) from the adapter pass only, as the filter pass's 5' adapter matches only flag reads for discarding.
`input.path1` is `<sample>.fastq`, so report tools name the sample after it. The report records that it is a merge
(`tag`: "Cutadapt report, merged from the passes in gcf_passes"; `gcf_passes`: each pass's command line and file),
and the reports of the passes are kept unchanged next to it.

Usage: merge_cutadapt.py --sample NAME ADAPTER_JSON FILTER_JSON > SAMPLE_JSON
"""
import argparse
import json
import sys

COUNTS = ('too_short', 'too_long', 'too_many_n', 'too_many_expected_errors', 'casava_filtered', 'discard_trimmed',
          'discard_untrimmed')


def add(a, b):
    return None if a is None and b is None else (a or 0) + (b or 0)


def merge(adapter, filtered, sample, paths):
    for name, report in (('adapter', adapter), ('filter', filtered)):
        if report.get('tag') != 'Cutadapt report' or 'read_counts' not in report:
            raise ValueError('{} report is not a cutadapt JSON report: {}'.format(name, paths[name]))
    if adapter['read_counts']['output'] != filtered['read_counts']['input']:
        raise ValueError('the filter pass read {} reads, the adapter pass wrote {}: not the same sample run'.format(
            filtered['read_counts']['input'], adapter['read_counts']['output']))
    merged = json.loads(json.dumps(filtered))
    merged['tag'] = 'Cutadapt report, merged from the passes in gcf_passes'
    merged['input'] = dict(adapter['input'], path1='{}.fastq'.format(sample))
    merged['command_line_arguments'] = adapter['command_line_arguments']
    merged['gcf_passes'] = [{'pass': name, 'file': paths[name], 'command_line_arguments': report['command_line_arguments']}
                            for name, report in (('adapter', adapter), ('filter', filtered))]
    counts = merged['read_counts']
    counts['input'] = adapter['read_counts']['input']
    counts['filtered'] = {key: add(adapter['read_counts']['filtered'].get(key), filtered['read_counts']['filtered'].get(key))
                          for key in COUNTS}
    counts['read1_with_adapter'] = adapter['read_counts']['read1_with_adapter']
    bp = merged['basepair_counts']
    bp['input'] = adapter['basepair_counts']['input']
    bp['input_read1'] = adapter['basepair_counts'].get('input_read1')
    for key in ('quality_trimmed', 'quality_trimmed_read1', 'poly_a_trimmed', 'poly_a_trimmed_read1'):
        if key in bp:
            bp[key] = add(adapter['basepair_counts'].get(key), filtered['basepair_counts'].get(key))
    merged['adapters_read1'] = adapter['adapters_read1']
    return merged


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--sample', required=True)
    parser.add_argument('adapter_json')
    parser.add_argument('filter_json')
    args = parser.parse_args()
    paths = {'adapter': args.adapter_json, 'filter': args.filter_json}
    try:
        reports = {name: json.load(open(path)) for name, path in paths.items()}
        merged = merge(reports['adapter'], reports['filter'], args.sample, paths)
    except (OSError, ValueError, KeyError) as error:
        sys.exit('merge_cutadapt.py: {}'.format(error))
    json.dump(merged, sys.stdout, indent=2)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
