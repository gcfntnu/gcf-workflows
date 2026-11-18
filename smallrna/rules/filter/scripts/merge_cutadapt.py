#!/usr/bin/env python
"""Merge a sample's three cutadapt JSON reports (3' adapter pass, filter pass, artefact pass) into one report in
cutadapt's schema.

Read counts: `input` from the adapter pass, `output` from the artefact pass, filtered categories summed over the
adapter and filter passes, and the artefact pass's discarded reads (`discard_trimmed` there) split by the record that
matched (`adapters_read1`, exact, as that pass only flags reads): reads still containing 3' adapter or Illumina primer
sequence (#353) as `gcf_three_prime_artefacts`, reads of a reagent oligo of the kit (records `reagent_oligo_*`, #380) as
`gcf_reagent_oligos` (a read matching both counts for the record cutadapt found best); any other record fails. So
`discard_trimmed` stays the 5' adapter artefacts. `read1_with_adapter` (3' adapter found) comes from the adapter pass,
and base pair counts likewise; `adapters_read1` (trimmed length statistics) from the adapter pass only, as the later
passes only flag reads for discarding.
`input.path1` is `<sample>.fastq`, so report tools name the sample after it. The report records that it is a merge
(`tag`: "Cutadapt report, merged from the passes in gcf_passes"; `gcf_passes`: each pass's command line and file),
and the reports of the passes are kept unchanged next to it.

Usage: merge_cutadapt.py --sample NAME ADAPTER_JSON FILTER_JSON ARTEFACTS_JSON > SAMPLE_JSON
"""
import argparse
import json
import sys

COUNTS = ('too_short', 'too_long', 'too_many_n', 'too_many_expected_errors', 'casava_filtered', 'discard_trimmed',
          'discard_untrimmed')
PASSES = ('adapter', 'filter', 'artefacts')
# record name prefixes of the artefact pass (cutadapt_config.THREE_PRIME_ARTEFACT_RECORD, REAGENT_OLIGO_RECORD)
THREE_PRIME_ARTEFACT_RECORD, REAGENT_OLIGO_RECORD = 'three_prime_artefact_', 'reagent_oligo_'


def add(*values):
    return None if all(v is None for v in values) else sum(v or 0 for v in values)


def merge(adapter, filtered, artefacts, sample, paths):
    reports = dict(zip(PASSES, (adapter, filtered, artefacts)))
    for name, report in reports.items():
        if report.get('tag') != 'Cutadapt report' or 'read_counts' not in report:
            raise ValueError('{} report is not a cutadapt JSON report: {}'.format(name, paths[name]))
    for before, after in zip(PASSES, PASSES[1:]):
        written, read = reports[before]['read_counts']['output'], reports[after]['read_counts']['input']
        if written != read:
            raise ValueError('the {} pass read {} reads, the {} pass wrote {}: not the same sample run'.format(
                after, read, before, written))
    extra = {k: v for k, v in artefacts['read_counts']['filtered'].items() if v and k != 'discard_trimmed'}
    if extra:
        raise ValueError('the artefact pass filtered reads other than by discard_trimmed: {}'.format(extra))
    merged = json.loads(json.dumps(artefacts))
    merged['tag'] = 'Cutadapt report, merged from the passes in gcf_passes'
    merged['input'] = dict(adapter['input'], path1='{}.fastq'.format(sample))
    merged['command_line_arguments'] = adapter['command_line_arguments']
    merged['gcf_passes'] = [{'pass': name, 'file': paths[name], 'command_line_arguments': report['command_line_arguments']}
                            for name, report in reports.items()]
    counts = merged['read_counts']
    counts['input'] = adapter['read_counts']['input']
    counts['filtered'] = {key: add(adapter['read_counts']['filtered'].get(key), filtered['read_counts']['filtered'].get(key))
                          for key in COUNTS}
    discarded = artefacts['read_counts']['filtered'].get('discard_trimmed') or 0
    records = artefacts.get('adapters_read1') or []
    if any(not isinstance(r.get('name'), str) or not isinstance(r.get('total_matches'), int) for r in records):
        raise ValueError('the artefact pass report has records without a name or total_matches: {}'.format(paths['artefacts']))
    unknown = [r['name'] for r in records if not r['name'].startswith((THREE_PRIME_ARTEFACT_RECORD, REAGENT_OLIGO_RECORD))]
    if unknown:
        raise ValueError('the artefact pass has records of unknown kind: {}: {}'.format(', '.join(unknown), paths['artefacts']))
    matched = sum(r['total_matches'] for r in records)
    if matched != discarded:
        raise ValueError('the artefact pass discarded {} reads, but its records matched {}: {}'.format(discarded, matched, paths['artefacts']))
    oligos = sum(r['total_matches'] for r in records if r['name'].startswith(REAGENT_OLIGO_RECORD))
    counts['filtered']['gcf_three_prime_artefacts'] = discarded - oligos
    counts['filtered']['gcf_reagent_oligos'] = oligos
    counts['read1_with_adapter'] = adapter['read_counts']['read1_with_adapter']
    bp = merged['basepair_counts']
    bp['input'] = adapter['basepair_counts']['input']
    bp['input_read1'] = adapter['basepair_counts'].get('input_read1')
    for key in ('quality_trimmed', 'quality_trimmed_read1', 'poly_a_trimmed', 'poly_a_trimmed_read1'):
        if key in bp:
            bp[key] = add(*(report['basepair_counts'].get(key) for report in reports.values()))
    merged['adapters_read1'] = adapter['adapters_read1']
    return merged


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--sample', required=True)
    parser.add_argument('adapter_json')
    parser.add_argument('filter_json')
    parser.add_argument('artefacts_json')
    args = parser.parse_args()
    paths = dict(zip(PASSES, (args.adapter_json, args.filter_json, args.artefacts_json)))
    try:
        reports = {name: json.load(open(path)) for name, path in paths.items()}
        merged = merge(reports['adapter'], reports['filter'], reports['artefacts'], args.sample, paths)
    except (OSError, ValueError, KeyError) as error:
        sys.exit('merge_cutadapt.py: {}'.format(error))
    json.dump(merged, sys.stdout, indent=2)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
