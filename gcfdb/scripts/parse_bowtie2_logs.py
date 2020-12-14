#!/usr/bin/env python
"""
Mostly borrowed code from MultiQC's parser. 
https://github.com/ewels/MultiQC/blob/master/multiqc/modules/bowtie2/bowtie2.py
"""

import sys
import re

regexes = {
    'unpaired': {
        'unpaired_aligned_none': r"(\d+) \([\d\.]+%\) aligned 0 times",
        'unpaired_aligned_one': r"(\d+) \([\d\.]+%\) aligned exactly 1 time",
        'unpaired_aligned_multi': r"(\d+) \([\d\.]+%\) aligned >1 times"
    },
    'paired': {
        'paired_aligned_none': r"(\d+) \([\d\.]+%\) aligned concordantly 0 times",
        'paired_aligned_one': r"(\d+) \([\d\.]+%\) aligned concordantly exactly 1 time",
        'paired_aligned_multi': r"(\d+) \([\d\.]+%\) aligned concordantly >1 times",
        'paired_aligned_discord_one': r"(\d+) \([\d\.]+%\) aligned discordantly 1 time",
        'paired_aligned_discord_multi': r"(\d+) \([\d\.]+%\) aligned discordantly >1 times",
        'paired_aligned_mate_one': r"(\d+) \([\d\.]+%\) aligned exactly 1 time",
        'paired_aligned_mate_multi': r"(\d+) \([\d\.]+%\) aligned >1 times",
        'paired_aligned_mate_none': r"(\d+) \([\d\.]+%\) aligned 0 times"
    }
}


def parse_log(fh):
    parsed_data = {}
    num_se, num_pe = 0.0, 0.0
    for l in fh:
        total = re.search(r"(\d+) reads; of these:", l)
        if total:
            parsed_data['total_reads'] = int(total.group(1))
        unpaired = re.search(r"(\d+) \([\d\.]+%\) were unpaired; of these:", l)

        if unpaired:
            parsed_data['unpaired_total'] = int(unpaired.group(1))
            num_se += 1
            l = fh.readline()
            while l.startswith('    '):
                for k, r in regexes['unpaired'].items():
                    match = re.search(r, l)
                    if match:
                        parsed_data[k] = int(match.group(1))
                l = fh.readline()

        paired = re.search(r"(\d+) \([\d\.]+%\) were paired; of these:", l)
        if paired:
            parsed_data['paired_total'] = int(paired.group(1))
            num_pe += 1

            # Do nested loop whilst we have this level of indentation
            l = fh.readline()
            while l.startswith('    '):
                for k, r in regexes['paired'].items():
                    match = re.search(r, l)
                    if match:
                        parsed_data[k] = int(match.group(1))
                l = fh.readline()

        overall = re.search(r"([\d\.]+)% overall alignment rate", l)
        
        if overall:
            parsed_data['overall_alignment_rate'] = float(overall.group(1))
            m_keys = ['paired_aligned_mate_multi', 'paired_aligned_mate_none', 'paired_aligned_mate_one']
            for k in m_keys:
                if k in parsed_data:
                    parsed_data['{}_halved'.format(k)] = float(parsed_data[k]) / 2.0
    parsed_data['num_se'] = num_se
    parsed_data['num_pe'] = num_pe
    
    return parsed_data

if __name__ == '__main__':
    log_data = parse_log(sys.stdin)
    num_se = log_data.pop('num_se')
    num_pe = log_data.pop('num_pe')
    sep = '\t'
    if num_se > 0:
        header = ['total_reads', 'mapped_unique', 'multimapped', 'not_mapped']
        keys =  ['total_reads', 'unpaired_aligned_one', 'unpaired_aligned_multi', 'unpaired_aligned_none']
        vals = [str(log_data[k]) for k in keys]
    if num_pe > 0:
        header = ['total_reads', 'mapped_unique', 'mapped_discordantly_unique', 'one_mate_mapped_unique', 'multimapped', 'dicordant_multimapped', 'one_mate_mutimapped', 'not_mapped']
        keys =  ['total_reads', 'paired_aligned_one', 'paired_aligned_discord_one', 'paired_aligned_mate_one_halved', 'paired_aligned_multi', 'paired_aligned_discord_multi', 'paired_aligned_mate_multi_halved', 'paired_aligned_mate_none_halved']
        vals = [str(log_data.get(k, 0)) for k in keys]
    sys.stdout.write(sep.join(header) + '\n')
    sys.stdout.write(sep.join(vals) + '\n')
