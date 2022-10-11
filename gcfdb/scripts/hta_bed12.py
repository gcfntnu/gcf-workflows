import sys

import pandas as pd

if __name__ == '__main__':
    bed = sys.argv[1]
    hta = sys.argv[2]
    output_fn = sys.argv[3]


    bed = pd.read_csv(bed, sep='\t', header=None)
    hk = pd.read_csv(hta, sep=";").filter(regex="Ensembl.*")
    hk = set(hk.iloc[:,0])
    # keep original sorting of bed file
    bed.index = bed.iloc[:,3]
    common = hk.intersection(bed.iloc[:,3])
    keep = [tx for tx in bed.index if tx in common]
    bed_hk = bed.loc[keep,:]
    print('found {} out of {} houskeeping genes in bed file'.format(len(keep), len(hk)))
    bed_hk.to_csv(output_fn, sep='\t', index=False, header=False)
    
