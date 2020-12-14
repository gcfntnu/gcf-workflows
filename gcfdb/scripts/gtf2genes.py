import sys
from gtfparse import read_gtf

df = read_gtf(sys.argv[1])
df_genes = df[df['feature'] == 'gene']
cols = df_genes.columns.copy()
cols = cols.insert(0, 'gene_id')
df_genes = df_genes[cols]
df_genes.to_csv(sys.stdout, sep='\t', index=False)
