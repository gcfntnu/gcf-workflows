import sys
import yaml

import pandas as pd


def col_is_empty(col):
    if all(col==''):
        return True
    if col.isna().all():
        return True
    if col.isnull().all():
        return True
    return False

if __name__ == '__main__':
    with open(sys.argv[1]) as fh:
        c = yaml.load(fh)
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
    df.to_csv(sys.stdout, sep='\t', index=False)
    
        
