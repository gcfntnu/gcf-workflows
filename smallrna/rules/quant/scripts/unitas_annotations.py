#!/usr/bin env python
import sys
import os
import glob
import argparse
import io
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import pandas as pd

import unitas_misc_rna


def move_y_rna(df, moved, source):
    """The sample's class table with the Y RNA that Ensembl files under another class moved to misc_RNA (#423): `moved`
    is {class: reads} from its hits_per_target table (unitas_misc_rna.y_rna_by_class). The classes still add up; a
    class left without reads (rat's Y_RNA) is dropped."""
    moved = {name: reads for name, reads in moved.items() if reads > 0}
    if not moved:
        return df
    df = df.copy()
    df['Count'] = df['Count'].astype(float)
    top = (df.Level2 == '') & (df.Level3 == '')
    for name, reads in moved.items():
        rows = df.index[top & (df.Level1 == name)]
        if len(rows) != 1 or df.loc[rows[0], 'Count'] < reads * (1 - unitas_misc_rna.TOTAL_REL_TOL):
            raise ValueError('{}: its {} targets that are Y RNA have {!r} reads, but the summary has {} for the class'.format(
                source, name, reads, 'no single row' if len(rows) != 1 else repr(float(df.loc[rows[0], 'Count']))))
        left = df.loc[rows[0], 'Count'] - reads
        if left <= reads * unitas_misc_rna.TOTAL_REL_TOL:
            df = df[df.Level1 != name]
            top = (df.Level2 == '') & (df.Level3 == '')
        else:
            df.loc[rows[0], 'Count'] = left
    total = sum(moved.values())
    rows = df.index[top & (df.Level1 == 'misc_RNA')]
    if len(rows):
        df.loc[rows[0], 'Count'] += total
    else:
        df = pd.concat([df, pd.DataFrame([{'Level1': 'misc_RNA', 'Level2': '', 'Level3': '', 'Count': total}])], ignore_index=True)
    return df


def samplesheet_ids(fn, sep='\t'):
    sample_ids = []
    with open(fn) as fh:
        txt = fh.read().splitlines()
        header = txt.pop(0).split(sep)
        if not 'Sample_ID' in header:
            raise ValueError('`Sample_ID` column not found in samplesheet')
        for line in txt:
            sample_ids.append(line.split('\t')[0])
        return sample_ids

def fill_col(col):
    P = None
    for i, e in enumerate(col):
        if pd.isnull(e):
            col.iloc[i] = P
        else:
            P = e
    return col

def fill_row(row):
    P = None
    for i, e in enumerate(row):
        if pd.isnull(e):
            row.iloc[i] = row.iloc[i-1]
        else:
            P = e
    return row

def file2pandas(fn):
    with open(fn) as fh:
        txt = fh.read()
    txt = txt.replace('   ', ';')
    header = '\t'.join(['Level1', 'Level2', 'Level3', 'Count'])
    new_lines = [header]
    parents = ['', '', '']
    for line in txt.split('\n'):
        if not line:
            break
        desc, count = line.split('\t')
        level = desc.count(';')
        # Each row names its parent levels, so the table does not depend on row order.
        parents[level] = desc.replace(';', '')
        parents[level + 1:] = [''] * (2 - level)
        new_lines.append('\t'.join(parents) + '\t' + count)
    txt = '\n'.join(new_lines)
    df = pd.read_csv(io.StringIO(txt), sep='\t', keep_default_na=False, dtype={'Level1': str, 'Level2': str, 'Level3': str})
    #df.Level1 = fill_col(df.Level1)
    #df = df.apply(fill_row, axis=1)
    return df

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output prefix')
    parser.add_argument('--hits-per-target', nargs='+', default=[], dest='hits',
                        help="each sample's unitas.hits_per_target.txt: its Y RNA outside the misc_RNA class is moved there (#423)")
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    
    import pandas as pd
    df_list = []
    

    
    hits = {os.path.dirname(fn).split(os.path.sep)[-1]: fn for fn in args.hits}
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = file2pandas(fn)
        if sample_id in hits:
            try:
                df = move_y_rna(df, unitas_misc_rna.y_rna_by_class(hits[sample_id]), hits[sample_id])
            except (OSError, ValueError) as error:
                sys.exit('unitas_annotations.py: {}'.format(error))
        df['Sample_ID'] = [sample_id] * df.shape[0]
        df_list.append(df)
    DF = pd.concat(df_list, axis=0, join='outer', sort=False)
    DF.sort_values(['Sample_ID', 'Level1', 'Level2', 'Level3'], kind='mergesort', inplace=True)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = os.path.join(args.output, 'annotations.tsv')
                                
    DF.to_csv(out_fn, sep='\t', index=False)
    
