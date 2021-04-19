import sys
from gtfparse import read_gtf

df = read_gtf(sys.argv[1])
df_gene = df[df['feature'] == 'gene']
cols = df_gene.columns.copy()
cols = cols.insert(0, 'gene_id')
df_gene = df_gene[cols]

# add mitochondrial_protein_coding as a biotype category
if 'gene_biotype' in df_gene.columns and 'seqname' in df_gene.columns:
    biotypes = df_gene.gene_biotype.copy()
    biotypes[(df_gene.seqname=="MT") & (df_gene.gene_biotype=="protein_coding")] = "Mt_protein_coding"
    df_gene.gene_biotype = biotypes

df_gene.to_csv(sys.stdout, sep='\t', index=False)
