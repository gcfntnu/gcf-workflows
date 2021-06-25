import sys
import os
import re
import glob
import argparse
import pandas as pd
import yaml
import json
import subprocess

QC_PLACEMENT = {
    'External_ID': 0,
    'Sample_Biosource': 10,
    'Sample_Group': 20,
    'Customer_Comment': 30,
    'Fragment_Length': 35,
    '260/230': 40,
    '260/280': 50,
    'Concentration': 50,
    'RIN': 60
}

def get_software_versions(args):
    versions = {}
    with open(os.path.join(args.repo_dir, ".git", "HEAD"), 'r') as head_fh:
        branch = head_fh.read().split('/')[-1].rstrip()
    with open(os.path.join(args.repo_dir, ".git", "refs", "heads", branch), 'r') as commit_fh:
        commit = commit_fh.read().rstrip()
    versions["Analysis pipeline"] = "github.com/gcfntnu/gcf-workflows/tree/{} commit {}".format(branch, commit).encode()
    #version["Workflow"] = args.workflow
    software = '\n'.join('{}: {}'.format(key,val) for (key,val) in versions.items())
    return software


def str_read_geometry(read_geometry):
    if len(read_geometry) == 1:
        rg_str = 'Single end - read length (R1): {}'.format(read_geometry[0])
    else:
        rg_str = 'Paired end - forward read length (R1): {}, reverse read length (R2): {}'.format(read_geometry[0], read_geometry[1])
    return rg_str


def create_mqc_config(args):
    mqc_conf = yaml.load(args.config_template)
    mqc_conf['title'] = args.project_id
    read_geometry = args.read_geometry.split(",")
    #pipeline = config.get("Options","pipeline")

    header_text = args.header_template.read()
    mqc_conf['intro_text'] = header_text.format(pname=args.project_id)
    software = get_software_versions(args)
    format_software = '<br/>'.join(["<strong>Software versions</strong>"] + software.split("\n"))
    mqc_conf['intro_text'] = '<br/><br/>'.join([mqc_conf['intro_text'],format_software])

    # ommit {'Contact E-mail': contact},
    report_header = [
    {'Sequencing Platform': args.machine},
    {'Read Geometry': str_read_geometry(read_geometry)},
    {'Organism': args.organism},
    {'Lib prep kit': args.libkit},
    ]

    mqc_conf['report_header_info'] = report_header

    if len(read_geometry) == 1:
        mqc_conf['extra_fn_clean_exts'].append('_R1')

    s_df = pd.read_csv(args.sample_info, sep='\t')
    s_df.index = s_df['Sample_ID']

    max_260_230 = float(s_df['260/230'].max()) if '260/230' in s_df.columns else 3
    max_260_280 = float(s_df['260/280'].max()) if '260/280' in s_df.columns else 3

    MAX = {
        'RIN': 10,
        '260/230': max_260_230,
        '260/280': max_260_280
        }

    COL_SCALE = {
        'RIN': 'OrRd',
        '260/230': 'PuBu',
        '260/280': 'BuGn'
    }

    s_df.drop(['Sample_ID'], axis=1,inplace=True)
    s_df.dropna(how='all', axis=1, inplace=True)
    s_df = s_df.round(2)
    s_dict = s_df.to_dict(orient='index')

    pconfig = {}
    for col in list(s_df.columns.values):
        if col not in QC_PLACEMENT.keys():
            continue
        pconfig[col] = {'format': '{}', 'min': 0, 'placement': QC_PLACEMENT[col]}
        pconfig[col]['max'] = MAX.get(col,None)
        pconfig[col]['scale'] = COL_SCALE.get(col,False)

    data = s_dict

    general_statistics = {
        'plot_type': 'generalstats',
        'pconfig': [pconfig],
        'data': data
    }
    custom_data = {'general_statistics': general_statistics}

    mqc_conf['custom_data'] = custom_data

    return mqc_conf


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--project-id", help="GCF project ID", required=True)
    parser.add_argument("-o", "--output", default=".multiqc_config.yaml", help="Output config file", type=argparse.FileType('w'), required=True)
    parser.add_argument("-S", "--sample-info", type=argparse.FileType('r'), help="Sample info in tsv format", required=True)
    parser.add_argument("--organism",  help="Organism (if applicable to all samples). Overrides value from samplesheet.", required=True)
    parser.add_argument("--libkit",  help="Library preparation kit name. (if applicable for all samples). Overrides value from samplesheet.", required=True)
    parser.add_argument("--machine",  help="Sequencer model.", required=True)
    parser.add_argument("--read-geometry",  help="Read geometry.", required=True)
    parser.add_argument("--repo-dir",  help="Path to git repo of workflow.", required=True)
    parser.add_argument("--header-template",  help="Path to multiqc header template.", type=argparse.FileType('r'), required=True)
    parser.add_argument("--config-template",  help="Path to multiqc config template.", type=argparse.FileType('r'), required=True)

    args = parser.parse_args()

    mqc_conf = create_mqc_config(args)
    yaml.dump(mqc_conf, args.output)


