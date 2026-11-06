#!/usr/bin/env python
"""Reconcile sequential SE fastp reports without losing stage provenance."""

import argparse
import copy
import html
import json
from pathlib import Path
import sys


def merge_reports(adapter, final):
    stages = [('adapter', adapter), ('filter', final)]
    for name, report in stages:
        counts = report['filtering_result']
        before = report['summary']['before_filtering']['total_reads']
        after = report['summary']['after_filtering']['total_reads']
        if counts['passed_filter_reads'] != after or sum(counts.values()) != before:
            raise ValueError(f'{name} fastp counts do not reconcile with its input/output reads')
    for key in ('total_reads', 'total_bases'):
        if adapter['summary']['after_filtering'][key] != final['summary']['before_filtering'][key]:
            raise ValueError(f'fastp stages do not form a sequential chain: {key}')

    merged = copy.deepcopy(final)
    merged['summary']['before_filtering'] = copy.deepcopy(adapter['summary']['before_filtering'])
    merged['filtering_result'] = {
        key: sum(report['filtering_result'].get(key, 0) for _, report in stages)
        for key in adapter['filtering_result'].keys() | final['filtering_result'].keys()
        if key != 'passed_filter_reads'
    }
    merged['filtering_result']['passed_filter_reads'] = final['filtering_result']['passed_filter_reads']
    for key in ('adapter_cutting', 'read1_before_filtering'):
        merged[key] = copy.deepcopy(adapter[key])
    merged['command'] = '\n'.join(report['command'].strip() for _, report in stages)
    merged['processing_stages'] = [
        dict(stage=name, command=report['command'], summary=copy.deepcopy(report['summary']),
             filtering_result=copy.deepcopy(report['filtering_result'])) for name, report in stages
    ]
    merged['report_type'] = 'smallrna_sequential_fastp'
    # The first pass measured duplication on the original input.
    if 'duplication' in adapter:
        merged['duplication'] = copy.deepcopy(adapter['duplication'])
    else:
        merged.pop('duplication', None)
    return merged


def render_html(report):
    summary = report['summary']
    rows = []
    for stage in report['processing_stages']:
        counts = stage['filtering_result']
        cells = ''.join(f'<tr><td>{html.escape(key)}</td><td>{value}</td></tr>' for key, value in sorted(counts.items()))
        rows.append(f'<h2>{html.escape(stage["stage"])} stage</h2>'
                    f'<pre>{html.escape(stage["command"])}</pre><table>{cells}</table>')
    return ('<!doctype html><html lang="en"><meta charset="utf-8"><title>Small-RNA fastp stages</title>'
            '<style>body{font-family:sans-serif;max-width:1000px;margin:2em auto}pre{white-space:pre-wrap}'
            'td{padding:.3em 1em;border-bottom:1px solid #ddd}</style><h1>Small-RNA fastp stages</h1>'
            f'<p>Original reads: {summary["before_filtering"]["total_reads"]}; '
            f'final retained reads: {summary["after_filtering"]["total_reads"]}.</p>'
            '<p>Adapter removal precedes final trimming and filtering. Each rejection belongs to the stage shown below. '
            'Empty reads can be rejected even when length filtering is disabled.</p>'
            '<p>Native per-sample adapter and filter JSON/HTML reports are delivered under provenance/fastp/. '
            'The combined JSON supplies original-input/final-output metrics to MultiQC.</p>'
            + ''.join(rows) + '</html>\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('adapter_report', type=Path)
    parser.add_argument('filter_report', type=Path)
    parser.add_argument('--html', type=Path)
    args = parser.parse_args()
    report = merge_reports(json.loads(args.adapter_report.read_text()), json.loads(args.filter_report.read_text()))
    if args.html:
        args.html.write_text(render_html(report))
    json.dump(report, sys.stdout, sort_keys=True)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
