import sys
from gtfparse import read_gtf

ADD_TXVER = False

df = read_gtf(sys.argv[1])
df_tx = df[df['feature'] == 'transcript'].copy()
if ADD_TXVER:
    tx = df_tx['transcript_id']
    tx_ver = df_tx['transcript_version']
    df_tx.loc[:,'transcript_id'] = ['{}.{}'.format(i, j) for i,j in zip(tx, tx_ver)]
cols = df_tx.columns.copy()
cols = cols.insert(0, 'transcript_id')
df_tx = df_tx[cols]
df_tx.to_csv(sys.stdout, sep='\t', index=False)
