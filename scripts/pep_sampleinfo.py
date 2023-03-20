import sys

import peppy
import pandas as pd


if __name__ == '__main__':
    pep = peppy.Project(sys.argv[1])
    df = pep.sample_table
    df = df.loc[:,df.convert_dtypes().dtypes != 'object']
    if 'sample_name' in df.columns:
        df = df.rename(columns={'sample_name':'Sample_ID'}).set_index('Sample_ID')
    df.to_csv(sys.argv[2], sep='\t')
