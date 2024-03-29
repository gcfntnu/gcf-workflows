#!/usr/bin/env python


import sys
import json
import copy


if __name__ == '__main__':
    with open(sys.argv[1]) as fh:
        trim_json = json.load(fh)
    with open(sys.argv[2]) as fh:
        filter_json = json.load(fh)


    filter_json['summary']['before_filtering'] = copy.copy(trim_json['summary']['before_filtering'])

    for k in ['low_quality_reads', 'too_many_N_reads', 'low_complexity_reads', 'too_short_reads', 'too_long_reads']:
        if k in  filter_json['filtering_result'] and k in trim_json['filtering_result']:
            filter_json['filtering_result'][k] += trim_json['filtering_result'][k]
        elif k in trim_json['filtering_result']:
            filter_json['filtering_result'][k] = trim_json['filtering_result'][k]
        else:
            pass

    filter_json['adapter_cutting'] = copy.copy(trim_json['adapter_cutting'])
    filter_json['read1_before_filtering'] = copy.copy(trim_json['read1_before_filtering'])
    if 'read2_before_filtering' in filter_json.keys():
        filter_json['read2_before_filtering'] = copy.copy(trim_json['read2_before_filtering'])

    json.dump(filter_json, sys.stdout)
