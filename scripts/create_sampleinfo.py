import sys
import yaml
import peppy
import pandas as pd
import peppy

def col_is_empty(col):
    if all(col==''):
        return True
    if col.isna().all():
        return True
    if col.isnull().all():
        return True
    return False


def sampleinfo_from_yaml(fn):
    with open(fn) as fh:
        c = yaml.safe_load(fh)
        assert('samples' in c)
    df = pd.DataFrame.from_dict(c['samples'], orient='index')
    df.index.name = 'Sample_ID'
    if not 'Sample_ID' in df.columns:
        df = df.reset_index()
    cols = list(df.columns.copy())
    cols.pop(cols.index('Sample_ID'))
    cols.insert(0, 'Sample_ID')
    df = df.loc[:,cols]
    empty = df.apply(col_is_empty, axis=0) # rm empty cols
    df = df.loc[:,empty==False]   
    return df

def sampleinfo_from_peppy(fn):
    pep = peppy.Project(fn)
    df = pep.sample_table
    df = df.loc[:,df.convert_dtypes().dtypes != 'object']
    if 'sample_name' in df.columns and 'Sample_ID' not in df.columns:
        df = df.rename(columns={'sample_name':'Sample_ID'})
    return df

    
if __name__ == '__main__':
    fn = sys.argv[1]
    if len(sys.argv) > 2:
        out = sys.argv[2]
    else:
        out = sys.stdout

    if fn.endswith('pep_config.yaml'):
        df = sampleinfo_from_peppy(fn)
    elif fn.endswith('config.yaml'):
        df = sampleinfo_from_yaml(fn)
    else:
        raise ValueError
    
    df.to_csv(out, sep='\t', index=False)
    
        
