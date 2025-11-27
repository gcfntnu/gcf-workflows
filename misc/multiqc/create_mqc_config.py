import sys
import os
import re
import glob
import argparse
import pandas as pd
import yaml
import json
import subprocess

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
    versions["Analysis pipeline"] = "github.com/gcfntnu/gcf-workflows/tree/{} commit {}".format(branch, commit)
    software = '\n'.join(f"{key}: {val}" for key, val in versions.items())
    return software


def str_read_geometry(read_geometry):
    if len(read_geometry) == 1:
        rg_str = 'Single end - read length (R1): {}'.format(read_geometry[0])
    else:
        rg_str = 'Paired end - forward read length (R1): {}, reverse read length (R2): {}'.format(read_geometry[0], read_geometry[1])
    return rg_str


class DummyPep:
    def __init__(self):
        self.config = {}
        self.samples = {}
        self.sample_table = None

        # flags
        self.subsamples = False
        self.multiple_flowcells = False
        self.multiple_projects = False
        self.single_cell = False

    def _empty_col(self, col):
        """Return True if a pandas Series is entirely empty / NA / ''."""
        # col == '' will be False for non-string entries
        return ((col == '') | col.isna()).all()

    def set_sample_table(self):
        samples = self.config.get("samples", {}) or {}
        self.samples = samples

        subsamples = False
        multiple_flowcells = False
        multiple_projects = False
        single_cell = False

        for sample_id, sample_conf in samples.items():
            # Treat presence of I1 as single-cell (or adjust logic as needed)
            if sample_conf.get("I1"):
                single_cell = True

            r1 = sample_conf.get("R1", "")
            r2 = sample_conf.get("R2", "")
            proj = sample_conf.get("Project_ID", "")
            flow = sample_conf.get("Flowcell_ID", "")

            # Any comma in R1/R2 means subsamples
            if "," in r1 or "," in r2:
                subsamples = True

            # Multiple projects if project list has >1 unique entries
            if "," in proj:
                if len(set(proj.split(","))) > 1:
                    multiple_projects = True

            # Multiple flowcells if flowcell list has >1 unique entries
            if "," in flow:
                if len(set(flow.split(","))) > 1:
                    multiple_flowcells = True

        self.subsamples = subsamples
        self.multiple_flowcells = multiple_flowcells
        self.multiple_projects = multiple_projects
        self.single_cell = single_cell

    def from_configfile(self, fn):
        with open(fn) as fh:
            self.config = yaml.safe_load(fh)

        self.set_sample_table()

        # samples is a dict mapping Sample_ID -> dict of attributes
        df = pd.DataFrame.from_dict(self.samples, orient="index")

        # If the index is actually Sample_ID, keep that knowledge
        if df.index.name is None:
            df.index.name = "Sample_ID"

        if self.subsamples:
            drop_cols = [c for c in df.columns if c.endswith("_md5sum")]

            # Flowcell columns
            if "Flowcell_Name" in df.columns or "Flowcell_ID" in df.columns:
                if self.multiple_flowcells:
                    # Can't safely collapse; just drop these
                    drop_cols.extend(
                        [c for c in ("Flowcell_Name", "Flowcell_ID") if c in df.columns]
                    )
                else:
                    # Collapse comma-separated flowcell entries to first element
                    if "Flowcell_Name" in df.columns:
                        df["Flowcell_Name"] = df["Flowcell_Name"].astype(str).str.split(",").str[0]
                    if "Flowcell_ID" in df.columns:
                        df["Flowcell_ID"] = df["Flowcell_ID"].astype(str).str.split(",").str[0]

            # Project column
            if "Project_ID" in df.columns:
                if self.multiple_projects:
                    drop_cols.append("Project_ID")
                else:
                    df["Project_ID"] = df["Project_ID"].astype(str).str.split(",").str[0]

            df = df.drop(columns=drop_cols, errors="ignore")

        # Normalize to have a 'sample_name' column
        if "Sample_ID" in df.columns:
            # Sample_ID was in columns (unlikely with orient="index", but safe)
            df = df.set_index("Sample_ID")
        else:
            # Use index as Sample_ID
            df = df.copy()
            df.index.name = "Sample_ID"

        df = df.reset_index().rename(columns={"Sample_ID": "sample_name"})
        self.sample_table = df
        return self

def create_mqc_config(args):
    pep_path = args.pep.name
    if os.path.basename(pep_path) == "pep_config.yaml":
        import peppy
        pep = peppy.Project(pep_path)
    else:
        pep = DummyPep().from_configfile(pep_path)
    
    mqc_conf = yaml.load(args.config_template, Loader=yaml.Loader)
    project_id = pep.config.get('Project_ID', args.project_id)
    if isinstance(project_id, str):
        title = project_id
    else:
        title = ",".join(map(str, project_id))
    mqc_conf['title'] = title
    rg = pep.config.get("read_geometry", args.read_geometry)
    # normalize to list of ints/strings
    if isinstance(rg, str):
        # e.g. "151,151" → ["151","151"]
        rg = [x.strip() for x in rg.split(",")]
    elif not isinstance(rg, (list, tuple)):
        rg = [rg]
    
    header_text = args.header_template.read()
    mqc_conf['intro_text'] = header_text.format(pname=title)
    software = get_software_versions(args)
    format_software = '<br/>'.join(["<strong>Software versions</strong>"] + software.split("\n"))
    mqc_conf['intro_text'] = '<br/><br/>'.join([mqc_conf['intro_text'],format_software])

    # ommit {'Contact E-mail': contact},
    report_header = [
        {'Sequencing Platform': pep.config.get('machine', args.machine)},
        {'Read Geometry': str_read_geometry(rg)},
        {'Organism': pep.config.get('organism', args.organism).replace('_', ' ').title()},
        {'Lib prep kit': pep.config.get('libprepkit', args.libkit)},
        {'Workflow': pep.config.get('workflow', args.workflow)}
    ]

    mqc_conf['report_header_info'] = report_header

    if len(rg) == 1:
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
        if col_name not in df.columns:
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

        return {k: cols[i] for i, k in enumerate(levels)}


    BGCOLS = {}
    for col_name in s_df.columns:
        if col_name in ['Sample_Group']:
            BGCOLS[col_name] = _get_colors(s_df, col_name, scale='mqc')
        elif col_name in ['Sample_Biosource']:
            BGCOLS[col_name] = _get_colors(s_df, col_name, scale='pairs')
        else:
            pass

    if s_df.shape[-1] > 0: # only if any extra cols exists
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
    parser.add_argument("--pep",  help="Path to peppy project (can also be config.yaml)", type=argparse.FileType('r'), required=True)

    args = parser.parse_args()

    mqc_conf = create_mqc_config(args)
    yaml.dump(mqc_conf, args.output)


