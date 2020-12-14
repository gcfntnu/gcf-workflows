import sys
import yaml

import pandas as pd


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
    df.to_csv(sys.stdout, sep='\t', index=False)
    
        
