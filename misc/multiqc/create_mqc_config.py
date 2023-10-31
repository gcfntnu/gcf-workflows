import sys
import os
import re
import glob
import argparse
import pandas as pd
import yaml
import json
import subprocess

import peppy
from matplotlib import cm, colors


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
    pep = peppy.Project(args.pep.name)
    mqc_conf = yaml.load(args.config_template, Loader=yaml.Loader)
    title = ','.join(pep.config.get('Project_ID', [args.project_id]))
    mqc_conf['title'] = title

    header_text = args.header_template.read()
    mqc_conf['intro_text'] = header_text.format(pname=title)
    software = get_software_versions(args)
    format_software = '<br/>'.join(["<strong>Software versions</strong>"] + software.split("\n"))
    mqc_conf['intro_text'] = '<br/><br/>'.join([mqc_conf['intro_text'],format_software])

    # ommit {'Contact E-mail': contact},
    report_header = [
        {'Sequencing Platform': pep.config.get('machine', args.machine)},
        {'Read Geometry': str_read_geometry(pep.config.read_geometry)},
        {'Organism': pep.config.get('organism', args.organism).replace('_', ' ').title()},
        {'Lib prep kit': pep.config.get('libprepkit', args.libkit)},
        {'Workflow': pep.config.get('workflow', args.workflow)}
    ]

    mqc_conf['report_header_info'] = report_header

    if len(pep.config.read_geometry) == 1:
        mqc_conf['extra_fn_clean_exts'].append('_R1')

 
    s_df = pep.sample_table
    
    s_df = s_df.rename(columns={'sample_name': 'Sample_ID'}).set_index('Sample_ID')
    drop = list(set(['Flowcell_Name', 'Project_ID', 'R1', 'R2', 'subsample_name', 'sample_name', 'Flowcell_ID', 'Lane', 'lane', 'run_number', 'I1', 'I1_md5sum', 'R1_md5sum', 'R2_md5sum', 'I2', 'I2_md5sum']).intersection(s_df.columns))
    if 'Organism' in s_df.columns and len(set(s_df['Organism'])) == 1:
        drop.append('Organism')
    s_df = s_df.drop(drop, axis=1)

    na_vals = ['nan', 'NAN', 'na', 'NA', 'n/a', 'N/A', 'None', 'none', '<na>', '<NA>']
    s_df = s_df.replace(na_vals, pd.NA)
    s_df.dropna(how='all', axis=1, inplace=True)
    s_df = s_df.round(2)
    s_df = s_df.fillna('')
    
    COL_SCALE = {
        'RIN': 'RdYlGn',
        '260/230': 'BuGn',
        '260/280': 'BuGn',
        'Concentration': 'BuGn'
    }
    
    def _get_colors(df, col_name, scale='pairs'):
        if not col_name in df.columns:
            return None
        col = df[col_name]
        levels = col.astype('category').cat.categories
        if scale == 'pairs':
            cols = list(map(colors.to_hex, cm.tab20.colors))[1:15:2]
        else:
            if len(levels) <= 10:
                cols = list(map(colors.to_hex, cm.tab10.colors))
            else:
                cols = list(map(colors.to_hex, cm.tab20.colors))
                
        
        levels = col.astype('category').cat.categories
        return {k:cols[i] for i,k in enumerate(levels)}

    BGCOLS = {}
    for col_name in s_df.columns:
        if col_name in ['Sample_Group']:
            BGCOLS[col_name] = _get_colors(s_df, col_name, scale='mqc')
        elif col_name in ['Sample_Biosource']:
            BGCOLS[col_name] = _get_colors(s_df, col_name, scale='pairs')
        else:
            pass
    

    s_dict = s_df.to_dict(orient='index')

    pconfig = {}
    #pconfig['title'] = 'GCF'
    desc = pep.config.get('descriptors', {})
    for col in list(s_df.columns.values):
        pconfig[col] = {'format': '{}', 'namespace': 'gcf'}
        if col not in desc.keys():
            continue
        if 'max' in desc[col]:
            pconfig[col]['max'] = desc[col]['max']
        if 'min' in desc[col]:
            pconfig[col]['min'] = desc[col]['min']
        if 'placement' in desc[col]:
            pconfig[col]['placement'] = desc[col]['placement']
        if 'display_name' in desc[col]:
            pconfig[col]['title'] = desc[col]['display_name']
        if 'description' in desc[col]:
            pconfig[col]['description'] = desc[col]['description']
        if 'suffix' in desc[col]:
            pconfig[col]['suffix'] = ' ' + desc[col]['suffix']
        if col in COL_SCALE:
            pconfig[col]['scale'] = COL_SCALE[col]
        if col in BGCOLS:
            pconfig[col]['bgcols'] = BGCOLS[col]

    general_statistics = {
        'plot_type': 'generalstats',
        'pconfig': [pconfig],
        'data': s_df.to_dict(orient='index')
    }
    custom_data = {'general_statistics': general_statistics}

    mqc_conf['custom_data'] = custom_data

    return mqc_conf


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--project-id", help="GCF project ID", required=True)
    parser.add_argument("-o", "--output", default=".multiqc_config.yaml", help="Output config file", type=argparse.FileType('w'), required=True)
    parser.add_argument("-S", "--sample-info", type=argparse.FileType('r'), help="Sample info in tsv format", required=True)
    parser.add_argument("--organism",  help="Organism (if applicable to all samples). Overrides value from samplesheet.", default='N/A')
    parser.add_argument("--libkit",  help="Library preparation kit name. (if applicable for all samples). Overrides value from samplesheet.", default='default')
    parser.add_argument("--workflow",  help="Snakemake workflow.", default='default')    
    parser.add_argument("--machine",  help="Sequencer model.", default='')
    parser.add_argument("--read-geometry",  help="Read geometry.", default=[75])
    parser.add_argument("--repo-dir",  help="Path to git repo of workflow.", required=True)
    parser.add_argument("--header-template",  help="Path to multiqc header template.", type=argparse.FileType('r'), required=True)
    parser.add_argument("--config-template",  help="Path to multiqc config template.", type=argparse.FileType('r'), required=True)
    parser.add_argument("--pep",  help="Path to peppy project", type=argparse.FileType('r'), required=True)

    args = parser.parse_args()

    mqc_conf = create_mqc_config(args)
    yaml.dump(mqc_conf, args.output)


