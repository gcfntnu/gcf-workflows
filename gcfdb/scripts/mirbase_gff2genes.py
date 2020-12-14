import sys
from collections import defaultdict

import pandas as pd

if __name__ == '__main__':
    
    df = pd.read_csv(sys.argv[1], sep='\t', header=None, comment='#')
    df.columns = ['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attrs']

    data = {}
    delete_keys = []
    i = 0
    hei = 0
    for name, row in df.iterrows():
        if row.feature == 'miRNA':
            attrs = row[8].split(';')
            for a in attrs:
                key, value = a.split('=')
                key = key.lower()
                if key == 'name':
                    key = 'mirna_id'
                row[key] = value
                if key == 'mirna_id':
                    if value in data:
                        for i in range(1, 11):
                            if value.endswith('p'):
                                base = '-'.join(value.split('-')[:-1])
                                ext = '-' + value.split('-')[-1]
                            else:
                                base = value
                                ext = ''
                            new_value = base + '-' + str(i) + ext
                            if hei:
                                pass
                                #print(new_value)
                            if new_value not in data:
                                row[key] = new_value
                                if hei:
                                    print("new name: " + new_value)
                                break
                    else:
                        row[key] = value
                    
                    i = 0
            data[row.mirna_id] = row
            if hei:
                pass
        del row['attrs']
    
    genes = pd.DataFrame.from_dict(data, orient='index')
    genes.set_index('mirna_id', inplace=True, drop=True)
    genes.index = genes.index.str.replace('^[a-z]{3}\-', '')
    genes.to_csv(sys.stdout, sep='\t', index=True)
