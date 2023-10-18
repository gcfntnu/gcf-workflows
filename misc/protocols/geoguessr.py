import sys
import os
import glob
import string
import re
import argparse
import copy

import pandas as pd
import openpyxl
from openpyxl import load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import peppy


METADATA_SHEET = '2. Metadata Template'
MD5_SHEET = '3. MD5 Checksums'
STUDY_START = 11
SAMPLES_START = 34
PROTOCOLS_START = 45
PAIRED_START = 65
MD5_START = 9

def _dedup_by_counter(x):
    """
    Make list unique by adding underscore + counter on duplicated values
    """
    index = pd.Index(x)
    if index.is_unique:
        return index
    from collections import Counter
    values = index.values.copy()
    indices_dup = index.duplicated(keep=False)
    values_dup = values[indices_dup]
    values_set = set(values)
    counter = Counter()
    for i, v in enumerate(values_dup):
        while True:
            counter[v] += 1
            tentative_new_name = v + '_' + str(counter[v])
            if tentative_new_name not in values_set:
                values_set.add(tentative_new_name)
                values_dup[i] = tentative_new_name
                break
    values[indices_dup] = values_dup
    index = pd.Index(values, name=index.name)
    return index.values

    
def insert_samples(pep,  workbook, md5, supplementary_file=None):
    sheet = workbook[METADATA_SHEET]
    header = [c.value for c in sheet[SAMPLES_START-1]]
    col_map = dict(zip(header, list(string.ascii_uppercase)))
    
    df = pep.sample_table[['sample_name', 'Sample_Group', 'R1']].copy()
    for n, row in df.iterrows():
        sample_id = row['sample_name']
        _md5 = md5.loc[md5.sample_id==sample_id]

        
    df.R1 = [os.path.basename(i) for i in df['R1']]
    df.columns = ['*library name', 'treatment', '*raw file']
    if 'Organism' in df.columns:
        organism = df['Organism']
    else:
        organism = df['*library name'].copy()
        organism.values[:] = pep.config.get('organism')
    df['*organism'] = organism
    title = pep.sample_table['Sample_Biosource'] + '_' + pep.sample_table['Sample_Group']
    title = title.str.replace(' ', '_')
    title = _dedup_by_counter(title)
    df['*title'] = title
    df['*instrument model'] = pep.config['machine']
    df['*molecule'] = pep.config['molecule']
    df['*single or paired-end'] = 'paired' if len(pep.config['read_geometry']) > 1 else 'single'
    if supplementary_file is not None:
        supplementary_file = [os.path.basename(fn) for fn in supplementary_file]
        if len(supplementary_file) == 1:
            df['*processed data file '] = supplementary_file[0]
        elif len(supplementary_file) == 2:
            df['*processed data file '] = supplementary_file[0]
            df['processed data file '] = supplementary_file[1]
        else:
            raise ValueError
    if len(pep.config['read_geometry']) > 1 :
        df['raw file'] = [os.path.basename(i) for i in df['R2']]

    # create room for samples
    n_rows = df.shape[0]
    if n_rows > 8: #available rows in template
        add_rows = n_rows - 8
        sheet.insert_rows(SAMPLES_START, amount=add_rows)

    for col_name in header:
        if col_name in df.columns:
            for i, val in enumerate(df[col_name].values):
                row_num = SAMPLES_START + i
                col_letter = col_map[col_name]
                sheet[f'{col_letter}{row_num}'].value = val

    return workbook

def insert_study(pep, workbook, supplementary_file=None, experimental_design=None):
    supplementary_file = copy.deepcopy(supplementary_file)
    sheet = workbook[METADATA_SHEET]
    project_id =  pep.config.get('project_id', ['NA'])[0]
    pi_filled = False
    for row in sheet.iter_rows(min_row=STUDY_START, max_col=2, max_row=STUDY_START+20, values_only=False):
        name, val = row
        #print("|{}|".format(name.value))
        if not name.value:
            if len(supplementary_file) > 0:
                val.value = os.path.basename(supplementary_file.pop(0))
                name.value = 'supplementary file'
            else:
                break
        if name.value == '*title':
            experiment_title = pep.config.get('experiment_title')
            if isinstance(experiment_title, (list, tuple)):
                experiment_title = experiment_title[0]
            experiment_title = experiment_title or project_id
            val.value = experiment_title
        elif name.value == '*summary (abstract)':
            experiment_summary = pep.config.get('experiment_summary')
            if experiment_summary in [None, '', 'NA']:
                experiment_summary = '***'
            val.value = experiment_summary
        elif name.value == 'contributor' and not pi_filled:
            PI = pep.config.get('experiment_principal_inverstigator')
            if PI:
                if PI == 'NA':
                    PI = '***'
                val.value = PI
            pi_filled = True
        elif name.value == '*experimental design':
            experimental_design = experimental_design or '***'
            val.value = experimental_design
        elif name.value == 'supplementary file':
            if len(supplementary_file) > 0:
                val.value = os.path.basename(supplementary_file.pop(0))
    
    return workbook
    

    
def insert_protocols(pep, workbook, demultiplex_protocol=None, workflow_protocols=None, processed_data_format=None, assembly=None):
    sheet = workbook[METADATA_SHEET]
    protocol_rownum = {
        'extract_protocol': 2,
        'library_construction_protocol': 3,
        'library_strategy': 4,
        'demultiplex_protocol': 6,
        'genome_build': 11,
        'processed_data_content': 12
    }
    extraction_protocol = pep.config.get('protocols', {}).get('extraction')
    if extraction_protocol:
        row_num = PROTOCOLS_START + protocol_rownum['extract_protocol']
        sheet[f'B{row_num}'].value = extraction_protocol

    library_construction = pep.config.get('protocols', {}).get('extraction')
    if library_construction:
        row_num = PROTOCOLS_START + protocol_rownum['library_construction_protocol']
        sheet[f'B{row_num}'].value = library_construction
    
    library_strategy = pep.config.get('library_strategy')
    if library_strategy:
        row_num = PROTOCOLS_START + protocol_rownum['library_strategy']
        sheet[f'B{row_num}'].value = library_strategy
    
    if demultiplex_protocol:
        row_num = PROTOCOLS_START + protocol_rownum['demultiplex_protocol']
        sheet[f'B{row_num}'].value = demultiplex_protocol

    if workflow_protocols:
        #fixme: insert rows if number of workflow protocols > 5
        for i, (k, v) in enumerate(workflow_protocols.items()):
            row_num = PROTOCOLS_START + protocol_rownum['demultiplex_protocol'] + i
            sheet[f'B{row_num}'].value = v
    if assembly:
        row_num = PROTOCOLS_START + protocol_rownum['genome_build']
        sheet[f'B{row_num}'].value = assembly
    
    if processed_data_format:
        row_num = PROTOCOLS_START + protocol_rownum['processed_data_content']
        sheet[f'B{row_num}'].value = processed_data_format
    
    return workbook

def insert_matched_pe_fastq(pep, workbook):
    sheet = workbook[METADATA_SHEET]
    header = [c.value for c in sheet[PAIRED_START-1]]
    col_map = dict(zip(header, list(string.ascii_uppercase)))
    df = pep.sample_table[['R1', 'R2']].copy()
    df['R1'] = [os.path.basename(i) for i in df['R1']]
    df['R2'] = [os.path.basename(i) for i in df['R2']]
    df.columns = header
    for col_name in header:
        if col_name in df.columns:
            for i, val in enumerate(df[col_name].values):
                row_num = PAIRED_START + i
                col_letter = col_map[col_name]
                sheet[f'{col_letter}{row_num}'].value = val
    return workbook


def read_md5(fn):
    md5 = pd.read_table(fn, sep="\s+", engine='python', names=['file checksum', 'path'])
    md5['file name'] = [os.path.basename(i) for i in md5.path]
    md5 = md5.drop('path', axis=1)
    sample_id = []
    for f in md5['file name']:
        m = re.match('(.*)_([RI][1-2]).fastq.gz$', f)
        if m:
            sample_id.append(m.groups()[0])
        else:
            msg = 'failed to identify sample_id from path `{}` '.format(f)
            raise ValueError(msg)
    md5['sample_id'] = sample_id
    
    return md5
    

def insert_md5sums(pep, workbook, md5, supplementary_file=None):
    sheet = workbook[MD5_SHEET]
    df = md5.loc[md5['sample_id'].isin(pep.sample_table.index)]
    raw_processed_index = {c.value:i for i, c in enumerate(sheet[MD5_START-2]) if c.value}
    header = ['file name', 'file checksum']
    for col_index, col_name in enumerate(header):
        for row_index, val in enumerate(df[col_name].values):
            row_index = MD5_START + row_index
            col_index = col_index + raw_processed_index['RAW FILES']
            col_letter = string.ascii_uppercase[col_index]
            sheet[f'{col_letter}{row_index}'].value = val

    if supplementary_file:
        import hashlib
        for i, fn in enumerate(supplementary_file):
            if os.path.exists(fn):
                with open(fn, "rb") as fh:
                    bytes = fh.read()
                    checksum = hashlib.md5(bytes).hexdigest()
                    row = [os.path.basename(fn), checksum]
                for col_index, col_name in enumerate(header):
                    val = row[col_index]
                    col_index = col_index + raw_processed_index['PROCESSED DATA FILES']
                    col_letter = string.ascii_uppercase[col_index]
                    row_index = MD5_START + i
                    sheet[f'{col_letter}{row_index}'].value = val
    return workbook

def get_workflow_protocols(repo_dir, workflow):
    protocols = {}
    with open(os.path.join(repo_dir, ".git", "HEAD"), 'r') as head_fh:
        branch = head_fh.read().split('/')[-1].rstrip()
    with open(os.path.join(repo_dir, ".git", "refs", "heads", branch), 'r') as commit_fh:
        commit = commit_fh.read().rstrip()
    workflow_version = "github.com/gcfntnu/gcf-workflows/tree/{} commit {}".format(branch, commit)
    protocols['gcf-workflow'] = f'FASTQ files were processsed by the {workflow} pipline of gcf-tools ({workflow_version}) ' 
    return protocols

def create_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-S", "--sample-info", required=True, 
                        help="Path to sample meta data")
    
    parser.add_argument("-T", "--template", required=True, 
                        help="Path to GEO seq_template.xlsx")

    parser.add_argument("-P", "--pep", type=argparse.FileType('r'), required=True,
                        help="Path to peppy project")

    parser.add_argument("-A", "--assembly",
                        help="", default=None)

    parser.add_argument("--read-geometry",
                        help="", required=True)

    parser.add_argument("--md5sum-file",
                        help="", required=True)

    parser.add_argument("--processed-data", default=None, nargs='+',
                        help="")

    parser.add_argument("--demultiplex-protocol", default=None,
                        help="")

    parser.add_argument("--repo-dir", default=None,
                        help="")

    parser.add_argument("-o", "--output", required=True,
                        help="output excel file")
    return parser

if __name__ == '__main__':
    args = create_parser().parse_args()
    pep = peppy.Project(args.pep.name)
    md5 = read_md5(args.md5sum_file)

    experimental_design = pep.config.get('experimental_design', '***')
    supplementary_file = args.processed_data
    if supplementary_file is not None:
        description_file = 'tab separated file with header and row names in first column' 
        protocol_file = ''
    demultiplex_protocol = args.demultiplex_protocol
    if demultiplex_protocol is not None:
        with open(demultiplex_protocol) as fh:
            demultiplex_protocol = fh.read().splitlines()
    else:
        demultiplex_protocol = "Sequencing basecalls were demultiplexed and converted into FASTQ using bcl2fastq v2.20.0.422"
    workflow_protocols = get_workflow_protocols(args.repo_dir, pep.config['workflow'])
    processed_data_format = None
    if processed_data_format is None:
        processed_data_format = "Tab seprated files with header and rownames in first column"
    if args.assembly == 'NA':
        assembly = None
    else:
        assembly = args.assembly
    
    wb = load_workbook(args.template)
    wb = insert_study(pep, wb, experimental_design=experimental_design, supplementary_file=supplementary_file)
    wb = insert_samples(pep, wb, md5, supplementary_file=supplementary_file)
    wb = insert_protocols(pep, wb, demultiplex_protocol, workflow_protocols, processed_data_format, assembly)
    if len(pep.config['read_geometry']) > 1:
        wb = insert_matched_pe_fastq(pep, workbook)
    wb = insert_md5sums(pep, wb, md5, supplementary_file=supplementary_file)

    wb.save(args.output)
    
    
